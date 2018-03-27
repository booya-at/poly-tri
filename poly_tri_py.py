import numpy as np
import copy

def make_key(i1, i2):
    """
    Make a tuple key such at i1 < i2
    """
    if i1 < i2:
        return (i1, i2)
    return (i2, i1)


class PolyTri(object):

    small = 1e-15

    def __init__(self, pts, boundaries=None, delaunay=True, holes=True, border=[]):

        # data structures
        self.pts = pts[:]  # copy
        self.tris = []  # cells
        self.edge2tris = {}  # edge to triangle(s) map
        self.pnt2tris = {}
        self.boundary_edges = set()
        self.delaunay = delaunay
        self.boundaries = boundaries
        self.border = border

        # compute center of gravity

        cg = sum(pts, np.zeros(2, np.float64)) / len(pts)
        i = 0
        while True:
            self.pts = pts[:]
            dist = self.pts - cg
            square_dist = dist.T[0]**2 + dist.T[1]**2
            self.order = np.argsort(square_dist)
            self.pts = self.pts[self.order]
    
            # create first triangle, make sure we're getting a non-zero area
            # otherwise drop the pts
    
            tri = [0, 1, 2]
            self.make_counter_clockwise(tri)
            if self.get_area(*tri) > self.small:
                self.tris.append(tri)
                break
            else:
                cg = self.pts[i]
                i += 1
        self.unorder = np.argsort(self.order)
            
        # boundary edges
        e01 = (tri[0], tri[1])
        e12 = (tri[1], tri[2])
        e20 = (tri[2], tri[0])

        self.boundary_edges.add(e01)
        self.boundary_edges.add(e12)
        self.boundary_edges.add(e20)
        
        self.edge2tris[make_key(*e01)] = [0]
        self.edge2tris[make_key(*e12)] = [0]
        self.edge2tris[make_key(*e20)] = [0]

        for i in tri:
            self.pnt2tris[i] = set([0])

        # add additional pts
        for i in range(3, len(self.pts)):
            self.add_point(i)

        if self.boundaries:
            self.constraint_boundaries()
            if holes:
                self.remove_empty()
                self.update_mapping()
                self.remove_holes()
    
    def get_tris(self):
        return [self.order[tri] for tri in self.tris]

    def get_area(self, ip0, ip1, ip2):
        """
        Compute the parallelipiped area
        @param ip0 index of first vertex
        @param ip1 index of second vertex
        @param ip2 index of third vertex
        """
        d1 = self.pts[ip1] - self.pts[ip0]
        d2 = self.pts[ip2] - self.pts[ip0]
        return (d1[0]*d2[1] - d1[1]*d2[0])

    def is_edge_visible(self, ip, edge):
        """
        Return true iff the point lies to its right when the edge pts down
        @param ip point index
        @param edge (2 point indices with orientation)
        @return True if visible
        """
        area = self.get_area(ip, edge[0], edge[1])
        if area < -self.small:
            return True
        return False

    def make_counter_clockwise(self, ips):
        """
        Re-order nodes to ensure positive area (in-place operation)
        """
        area = self.get_area(ips[0], ips[1], ips[2])
        if area < -self.small:
            ip1, ip2 = ips[1], ips[2]
            # swap
            ips[1], ips[2] = ip2, ip1

    def constraint_edge(self, cb):
        # start with first point in edge:
        if cb in self.edge2tris.keys():
            return
        pt0, pt1 = cb

        # find first intesecting edge
        for tri in self.pnt2tris[pt1]:
            edge = copy.copy(self.tris[tri])
            edge.remove(pt1)
            edge = make_key(*edge)
            if self.is_intersecting(edge, cb):
                break
        else:
            # no intesection
            # something went wrong
            return
        
        # now we know the first intersecting edge
        # and flip this edge
        edges = self.flipOneEdge(edge, delaunay=False, check_self_intersection=True)
        # 4 edges should be returned (theoretically only 2 are possebly intersecting)
        while True:
            for e in edges:
                if self.is_intersecting(e, cb):
                    edge = e
                    break
            else:
                # no intersection means we are done with constraining
                break
            edges = self.flipOneEdge(edge, delaunay=False, check_self_intersection=True)
            if not edges:
                if cb in self.edge2tris.keys():
                    break
                else:
                    self.constraint_edge(make_key(edge[0], pt0))
                    self.constraint_edge(make_key(pt0, pt1))


    def flipOneEdge(self, edge, delaunay=True, check_self_intersection=False):
        """
        Flip one edge then update the data structures
        @return set of edges to interate over at next iteration
        """

        # start with empty set
        res = set()

        # assume edge is sorted
        _tris = self.edge2tris.get(edge, [])
        if len(_tris) < 2:
                # nothing to do, just return
                return res

        iTri1, iTri2 = _tris
        tri1 = self.tris[iTri1]
        tri2 = self.tris[iTri2]

        # find the opposite vertices, not part of the edge
        iOpposite1 = -1
        iOpposite2 = -1
        for i in range(3):
            if not tri1[i] in edge:
                iOpposite1 = tri1[i]
            if not tri2[i] in edge:
                iOpposite2 = tri2[i]

        if check_self_intersection:
            diagonal2 = make_key(iOpposite1, iOpposite2)
            if not self.is_intersecting(edge, diagonal2):
                return set()

        if delaunay:
            # compute the 2 angles at the opposite vertices
            da1 = self.pts[edge[0]] - self.pts[iOpposite1]
            db1 = self.pts[edge[1]] - self.pts[iOpposite1]
            da2 = self.pts[edge[0]] - self.pts[iOpposite2]
            db2 = self.pts[edge[1]] - self.pts[iOpposite2]
            crossProd1 = self.get_area(iOpposite1, edge[0], edge[1])
            crossProd2 = self.get_area(iOpposite2, edge[1], edge[0])
            dotProd1 = np.dot(da1, db1)
            dotProd2 = np.dot(da2, db2)
            angle1 = abs(np.arctan2(crossProd1, dotProd1))
            angle2 = abs(np.arctan2(crossProd2, dotProd2))

            # Delaunay's test
            if not (angle1 + angle2 > np.pi*(1.0 + self.small)):
                return res

        # flip the tris
        #                         / ^ \                                        / b \
        # iOpposite1 + a|b + iOpposite2    =>     + - > +
        #                         \     /                                        \ a /

        newTri1 = [iOpposite1, edge[0], iOpposite2]  # triangle a
        newTri2 = [iOpposite1, iOpposite2, edge[1]]  # triangle b

        # update the triangle data structure
        self.tris[iTri1] = newTri1
        self.tris[iTri2] = newTri2

        # now handle the topolgy of the edges

        # remove this edge
        del self.edge2tris[edge]

        # add new edge
        e = make_key(iOpposite1, iOpposite2)
        self.edge2tris[e] = [iTri1, iTri2]
        self.pnt2tris[e[0]] |= set([iTri1, iTri2])
        self.pnt2tris[e[1]] |= set([iTri1, iTri2])

        # modify two edge entries which now connect to
        # a different triangle
        e = make_key(iOpposite1, edge[1])
        v = self.edge2tris[e]
        for i in range(len(v)):
            if v[i] == iTri1:
                v[i] = iTri2
        res.add(e)



        e = make_key(iOpposite2, edge[0])
        v = self.edge2tris[e]
        for i in range(len(v)):
            if v[i] == iTri2:
                v[i] = iTri1
        res.add(e)

        # updating the remining points
        # assume iTri2 is not connected to edge U newTri1
        for i in newTri1:
            if i in edge:
                tr = list(self.pnt2tris[i])
                for j in range(len(tr)):
                    if tr[j] == iTri2:
                        tr[j] = iTri1
                self.pnt2tris[i] = set(tr)

        # assume iTri1 is not connected to edge U newTri2
        for i in newTri2:
            if i in edge:
                tr = list(self.pnt2tris[i])
                for j in range(len(tr)):
                    if tr[j] == iTri1:
                        tr[j] = iTri2
                self.pnt2tris[i] = set(tr)


        # these two edges might need to be flipped at the
        # next iteration
        res.add(make_key(iOpposite1, edge[0]))
        res.add(make_key(iOpposite2, edge[1]))

        return res

    def flip_edges(self):
        edgeSet = set(self.edge2tris.keys())

        continueFlipping = True

        while continueFlipping:
            newEdgeSet = set()
            for edge in edgeSet:
                edges2proceed = self.flipOneEdge(edge)
                for add_edge in edges2proceed:
                    newEdgeSet.add(add_edge)
            edgeSet = copy.copy(newEdgeSet)
            continueFlipping = (len(edgeSet) > 0)

    def add_point(self, ip):
        """
        Add point
        @param ip point index
        """

        # collection for later updates
        boundary_edges2remove = set()
        boundary_edges2add = set()

        for edge in self.boundary_edges:
            if self.is_edge_visible(ip, edge):

                # create new triangle
                newTri = [edge[0], edge[1], ip]
                newTri.sort()
                self.make_counter_clockwise(newTri)
                self.tris.append(newTri)
                iTri = len(self.tris) - 1

                # add the two boundary edges
                e0 = make_key(*edge)
                e1 = make_key(ip, edge[0])
                e2 = make_key(edge[1], ip)
                for e in (e0, e1, e2):
                    v = self.edge2tris.get(e, [])
                    v.append(iTri)
                    self.edge2tris[e] = v

                # add point to triangle information
                for i in newTri:
                    p2t = self.pnt2tris.get(i, set())
                    p2t.add(iTri)
                    self.pnt2tris[i] = p2t

                # keep track of the boundary edges to update
                boundary_edges2remove.add(edge)
                boundary_edges2add.add((edge[0], ip))
                boundary_edges2add.add((ip, edge[1]))

        # update the boundary edges
        for bedge in boundary_edges2remove:
            self.boundary_edges.remove(bedge)
        
        for bedge in boundary_edges2add:
            bEdgeSorted = make_key(*bedge)
            if len(self.edge2tris[bEdgeSorted]) == 1:
                # only add boundary edge if it does not appear
                # twice in different order
                self.boundary_edges.add(bedge)

        if self.delaunay:  # recursively flip edges
            self.flip_edges()

    
    def create_boundary_list(self, border=None, create_key=True):
        constrained_boundary = []
        for k, boundary in enumerate(self.boundaries):
            if border and k not in border:
                continue
            b = self.unorder[boundary]
            for i, j in zip(b[:-1], b[1:]):
                item = (i, j)
                if create_key:
                    item = make_key(*item)
                constrained_boundary.append(item)
        return constrained_boundary

    def constraint_boundaries(self):
        boundary = self.create_boundary_list()
        for cb in boundary:
            self.constraint_edge(cb)

    def update_mapping(self):
        self.edge2tris = {}
        self.pnt2tris = {}
        for i, tri in enumerate(self.tris):
            for edge in self.tri2edges(tri):
                e2t = self.edge2tris.get(edge, [])
                self.edge2tris[edge] = e2t + [i]
            for point in tri:
                p2t = self.pnt2tris.get(point, set())
                p2t.add(i)
                self.pnt2tris[point] = p2t
                

    def remove_empty(self):
        tris2remove = []
        for i, tri in enumerate(self.tris):
            self.make_counter_clockwise(tri)
            area = self.get_area(*tri)
            if area < 1e-10:
                tris2remove.append(i)

        tris2remove.sort()
        tris2remove.reverse()
        for i in tris2remove:
            self.tris.pop(i)

    def tri2edges(self, tri, create_key=True):
        _tri = tri + [tri[0]]
        if create_key:
            return [make_key(*edge) for edge in zip(_tri[:-1], _tri[1:])]
        else:
            return [tuple(edge) for edge in zip(_tri[:-1], _tri[1:])]

    def is_intersecting(self, edge1, edge2):
        if edge1[0] in edge2 or edge1[1] in edge2:
            return False
        p11 = self.pts[edge1[0]]
        p12 = self.pts[edge1[1]]
        p21 = self.pts[edge2[0]]
        p22 = self.pts[edge2[1]]
        t = p12 - p11
        s = p22 - p21
        r = p21 - p11
        try:
            c1, c2 = np.linalg.inv(np.array([t, -s]).T) @ r
        except np.linalg.linalg.LinAlgError:
            return False
        return (0 < c1 < 1) and (0 < c2 < 1)

    def remove_holes(self):
        bs = self.create_boundary_list(self.border)
        o_bs = self.create_boundary_list(self.border, create_key=False)
        remove_edges = set()
        for b, o_b in zip(bs, o_bs):
            _tris = self.edge2tris[b]
            for tri in _tris:
                if o_b in self.tri2edges(self.tris[tri], create_key=False):
                    edges = self.tri2edges(self.tris[tri])
                    for edge in edges:
                        remove_edges.add(edge)
        for b in bs:
            if b in remove_edges:
                remove_edges.remove(b)
        tris2remove = set()
        for edge in remove_edges:
            for tri in self.edge2tris[edge]:
                tris2remove.add(tri)
        num_tris2remove = len(tris2remove)
        while True:
            for tri in tris2remove:
                for _edge in self.tri2edges(self.tris[tri], create_key=True):
                    remove_edges.add(_edge)
            for b in bs:
                if b in remove_edges:
                    remove_edges.remove(b)
            for edge in remove_edges:
                for tri in self.edge2tris[edge]:
                    tris2remove.add(tri)
            if num_tris2remove == len(tris2remove):
                break
            else:
                num_tris2remove = len(tris2remove)
        tris2remove = list(tris2remove)
        tris2remove.sort()
        tris2remove.reverse()
        for i in tris2remove:
            self.tris.pop(i)
