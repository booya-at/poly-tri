import numpy as np
import copy
import numba

small = 1e-15

@numba.jit
def make_key(i1, i2):
    """
    Make a tuple key such at i1 < i2
    """
    if i1 < i2:
        return (i1, i2)
    return (i2, i1)

@numba.jit
def get_area(pts, ip0, ip1, ip2):
    """
    Compute the parallelipiped area
    """
    d1 = pts[ip1] - pts[ip0]
    d2 = pts[ip2] - pts[ip0]
    return (d1[0]*d2[1] - d1[1]*d2[0])

@numba.jit
def is_edge_visible(pts, ip, edge):
    """
    Return true if the point lies to its right when the edge pts down
    """
    area = get_area(pts, ip, edge[0], edge[1])
    if area < -small:
        return True
    return False

@numba.jit
def make_counter_clockwise(pts, ips):
    """
    Re-order nodes to ensure positive area (in-place operation)
    """
    area = get_area(pts, ips[0], ips[1], ips[2])
    if area < -small:
        ip1, ip2 = ips[1], ips[2]
        # swap
        ips[1], ips[2] = ip2, ip1

@numba.jit
def is_intersecting(pts, edge1, edge2):
    """
    checks if two edges are intersecting
    """
    if edge1[0] in edge2 or edge1[1] in edge2:
        return False
    p11 = pts[edge1[0]]
    p12 = pts[edge1[1]]
    p21 = pts[edge2[0]]
    p22 = pts[edge2[1]]
    t = p12 - p11
    s = p22 - p21
    r = p21 - p11
    a = t[0]
    b = -s[0]
    c = t[1]
    d = -s[1]
    e = r[0]
    f = r[1]
    g = (a*d - b*c)
    if abs(g) < small:
        return False  # parallel
    # computing the inverse
    c1 = (-b * f + d * e) / g
    c2 = (a * f - c * e) / g
    return (0 <= c1 <= 1) and (0 <= c2 <= 1)

@numba.jit
def tri2edges(tri, create_key=True):
    """
    creates edges from a triangle
    """
    _tri = tri + [tri[0]]
    edges = []
    for edge in zip(_tri[:-1], _tri[1:]):
        if create_key:
            edges.append(make_key(*edge))
        else:
            edges.append((edge[0], edge[1]))
    return edges

@numba.jit
def delaunay_check(pts, edge, iOpposite1, iOpposite2):
    # compute the 2 angles at the opposite vertices
    da1 = pts[edge[0]] - pts[iOpposite1]
    db1 = pts[edge[1]] - pts[iOpposite1]
    da2 = pts[edge[0]] - pts[iOpposite2]
    db2 = pts[edge[1]] - pts[iOpposite2]
    crossProd1 = get_area(pts, iOpposite1, edge[0], edge[1])
    crossProd2 = get_area(pts, iOpposite2, edge[1], edge[0])
    dotProd1 = np.dot(da1, db1)
    dotProd2 = np.dot(da2, db2)
    angle1 = abs(np.arctan2(crossProd1, dotProd1))
    angle2 = abs(np.arctan2(crossProd2, dotProd2))
    return (angle1 + angle2 > np.pi*(1.0 + small))


def create_boundary_list(boundaries, unorder, border, create_key=True):
    constrained_boundary = []
    for k, boundary in enumerate(boundaries):
        if k in border:
            b = unorder[boundary]
            for i, j in zip(b[:-1], b[1:]):
                item = (i, j)
                if create_key:
                    item = make_key(*item)
                constrained_boundary.append(item)
    return constrained_boundary


def poly_tri(pts, boundaries=[], delaunay=True, holes=True, border=[]):

    # data structures
    _pts = pts[:]  # copy
    tris = []  # cells
    edge2tris = {}  # edge to triangle(s) map
    pnt2tris = {}
    boundary_edges = set()
    border = border or list(range(len(boundaries)))

    # compute center of gravity

    cg = sum(pts, np.zeros(2, np.float64)) / len(pts)
    i = 0
    while True:
        pts = _pts[:]
        dist = pts - cg
        square_dist = dist.T[0]**2 + dist.T[1]**2
        order = np.argsort(square_dist)
        pts = pts[order]

        # create first triangle, make sure we're getting a non-zero area
        # otherwise drop the pts

        tri = [0, 1, 2]
        make_counter_clockwise(pts, tri)
        if get_area(pts, *tri) > small:
            tris.append(tri)
            break
        else:
            cg = pts[i]
            i += 1
    unorder = np.argsort(order)
        
    # boundary edges
    e01 = (tri[0], tri[1])
    e12 = (tri[1], tri[2])
    e20 = (tri[2], tri[0])

    boundary_edges.add(e01)
    boundary_edges.add(e12)
    boundary_edges.add(e20)
    
    edge2tris[make_key(*e01)] = [0]
    edge2tris[make_key(*e12)] = [0]
    edge2tris[make_key(*e20)] = [0]

    for i in tri:
        pnt2tris[i] = set([0])

    # add additional pts
    for i in range(3, len(pts)):
        add_point(pts, tris, edge2tris, pnt2tris, i, boundary_edges, delaunay)

    if boundaries:
        boundary = create_boundary_list(boundaries, unorder,
                                        border=range(len(boundaries)),
                                        create_key=True)
        for cb in boundary:
            constraint_edge(pts, tris, edge2tris, pnt2tris, cb)
        if holes:
            remove_empty(pts, tris, edge2tris, pnt2tris)
            update_mapping(pts, tris, edge2tris, pnt2tris)
            remove_holes(tris, edge2tris, pnt2tris,  boundaries, border, unorder)
    
    output = []
    for tri in tris:
        output.append(order[tri])
    return output

def constraint_edge(pts, tris, edge2tris, pnt2tris, cb):
    # start with first point in edge:
    if cb in edge2tris.keys():
        return
    pt0, pt1 = cb

    # find first intesecting edge
    for tri in pnt2tris[pt1]:
        edge = copy.copy(tris[tri])
        edge.remove(pt1)
        edge = make_key(*edge)
        if is_intersecting(pts, edge, cb):
            break
    else:
        # no intesection
        # something went wrong
        return

    # now we know the first intersecting edge
    # and flip this edge
    edges = flip_one_edge(pts, tris, edge2tris, pnt2tris, edge, delaunay=False,
                          check_self_intersection=True)

    # 4 edges should be returned (theoretically only 2 are possebly intersecting)
    while True:
        for e in edges:
            if is_intersecting(pts, e, cb):
                edge = e
                break
        else:
            # no intersection means we are done with constraining
            break
        edges = flip_one_edge(pts, tris, edge2tris, pnt2tris, edge,
                              delaunay=False, check_self_intersection=True)
        if not edges:
            if cb in edge2tris.keys():
                break
            else:
                constraint_edge(pts, tris, edge2tris, pnt2tris,
                                make_key(edge[0], pt0))
                constraint_edge(pts, tris, edge2tris, pnt2tris,
                                make_key(pt0, pt1))

def flip_one_edge(pts, tris, edge2tris, pnt2tris, edge, delaunay=True,
                check_self_intersection=False):
    """
    Flip one edge then update the data structures
    """
    # start with empty set
    res = set()
    proceed = True

    # assume edge is sorted
    _tris = edge2tris.get(edge, [])
    if len(_tris) < 2:
        proceed = False

    if proceed:
        iTri1, iTri2 = _tris
        tri1 = tris[iTri1]
        tri2 = tris[iTri2]

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
            if not is_intersecting(pts, edge, diagonal2):
                proceed = False
    if proceed:
        if delaunay and not delaunay_check(pts, edge, iOpposite1, iOpposite2):
            proceed = False

    if proceed:
            # flip the tris
        #                         / ^ \                                        / b \
        # iOpposite1 + a|b + iOpposite2    =>     + - > +
        #                         \     /                                        \ a /

        newTri1 = [iOpposite1, edge[0], iOpposite2]  # triangle a
        newTri2 = [iOpposite1, iOpposite2, edge[1]]  # triangle b

        # update the triangle data structure
        tris[iTri1] = newTri1
        tris[iTri2] = newTri2

        # now handle the topolgy of the edges

        # remove this edge
        del edge2tris[edge]

        # add new edge
        e = make_key(iOpposite1, iOpposite2)
        edge2tris[e] = [iTri1, iTri2]
        pnt2tris[e[0]] |= set([iTri1, iTri2])
        pnt2tris[e[1]] |= set([iTri1, iTri2])

        # modify two edge entries which now connect to
        # a different triangle
        e = make_key(iOpposite1, edge[1])
        v = edge2tris[e]
        for i in range(len(v)):
            if v[i] == iTri1:
                v[i] = iTri2
        res.add(e)



        e = make_key(iOpposite2, edge[0])
        v = edge2tris[e]
        for i in range(len(v)):
            if v[i] == iTri2:
                v[i] = iTri1
        res.add(e)

        # updating the remining points
        # assume iTri2 is not connected to edge U newTri1
        for i in newTri1:
            if i in edge:
                tr = list(pnt2tris[i])
                for j in range(len(tr)):
                    if tr[j] == iTri2:
                        tr[j] = iTri1
                pnt2tris[i] = set(tr)

        # assume iTri1 is not connected to edge U newTri2
        for i in newTri2:
            if i in edge:
                tr = list(pnt2tris[i])
                for j in range(len(tr)):
                    if tr[j] == iTri1:
                        tr[j] = iTri2
                pnt2tris[i] = set(tr)


        # these two edges might need to be flipped at the
        # next iteration
        res.add(make_key(iOpposite1, edge[0]))
        res.add(make_key(iOpposite2, edge[1]))

    return res

def flip_edges(pts, tris, edge2tris, pnt2tris):
    edgeSet = set(edge2tris.keys())

    continueFlipping = True

    while continueFlipping:
        newEdgeSet = set()
        for edge in edgeSet:
            edges2proceed = flip_one_edge(pts, tris, edge2tris, pnt2tris, edge)
            for add_edge in edges2proceed:
                newEdgeSet.add(add_edge)
        edgeSet = copy.copy(newEdgeSet)
        continueFlipping = (len(edgeSet) > 0)

def add_point(pts, tris, edge2tris, pnt2tris, ip, boundary_edges, delaunay):
    """
    Add point
    """

    # collection for later updates
    boundary_edges2remove = set()
    boundary_edges2add = set()

    for edge in boundary_edges:
        if is_edge_visible(pts, ip, edge):

            # create new triangle
            newTri = [edge[0], edge[1], ip]
            newTri.sort()
            make_counter_clockwise(pts, newTri)
            tris.append(newTri)
            iTri = len(tris) - 1

            # add the two boundary edges
            e0 = make_key(*edge)
            e1 = make_key(ip, edge[0])
            e2 = make_key(edge[1], ip)
            for e in (e0, e1, e2):
                v = edge2tris.get(e, [])
                v.append(iTri)
                edge2tris[e] = v

            # add point to triangle information
            for i in newTri:
                p = pnt2tris.get(i, set())
                p.add(iTri)
                pnt2tris[i] = p

            # keep track of the boundary edges to update
            boundary_edges2remove.add(edge)
            boundary_edges2add.add((edge[0], ip))
            boundary_edges2add.add((ip, edge[1]))

    # update the boundary edges
    for bedge in boundary_edges2remove:
        boundary_edges.remove(bedge)
    
    for bedge in boundary_edges2add:
        bEdgeSorted = make_key(*bedge)
        if len(edge2tris[bEdgeSorted]) == 1:
            # only add boundary edge if it does not appear
            # twice in different order
            boundary_edges.add(bedge)

    if delaunay:  # recursively flip edges
        flip_edges(pts, tris, edge2tris, pnt2tris)

def update_mapping(pts, tris, edge2tris, pnt2tris):
    edge2tris = {}
    pnt2tris = {}
    for i, tri in enumerate(tris):
        for edge in tri2edges(tri):
            e2t = edge2tris.get(edge, [])
            edge2tris[edge] = e2t + [i]
        for point in tri:
            p2t = pnt2tris.get(point, set())
            p2t.add(i)
            pnt2tris[point] = p2t

def remove_empty(pts, tris, edge2tris, pnt2tris):
    tris2remove = []
    for i, tri in enumerate(tris):
        make_counter_clockwise(pts, tri)
        area = get_area(pts, *tri)
        if area < 1e-10:
            tris2remove.append(i)

    tris2remove.sort()
    tris2remove.reverse()
    for i in tris2remove:
        tris.pop(i)

def remove_holes(tris, edge2tris, pnt2tris, boundaries, border, unorder):
    bs = create_boundary_list(boundaries, unorder, border)
    o_bs = create_boundary_list(boundaries, unorder, border, create_key=False)
    remove_edges = set()
    for b, o_b in zip(bs, o_bs):
        _tris = edge2tris[b]
        for tri in _tris:
            if o_b in tri2edges(tris[tri], create_key=False):
                edges = tri2edges(tris[tri])
                for edge in edges:
                    remove_edges.add(edge)
    for b in bs:
        if b in remove_edges:
            remove_edges.remove(b)
    tris2remove = set()
    for edge in remove_edges:
        for tri in edge2tris[edge]:
            tris2remove.add(tri)
    num_tris2remove = len(tris2remove)
    while True:
        for tri in tris2remove:
            for _edge in tri2edges(tris[tri], create_key=True):
                remove_edges.add(_edge)
        for b in bs:
            if b in remove_edges:
                remove_edges.remove(b)
        for edge in remove_edges:
            for tri in edge2tris[edge]:
                tris2remove.add(tri)
        if num_tris2remove == len(tris2remove):
            break
        else:
            num_tris2remove = len(tris2remove)
    tris2remove = list(tris2remove)
    tris2remove.sort()
    tris2remove.reverse()
    for i in tris2remove:
        tris.pop(i)


class PolyTri(object):
    def __init__(self, pts, boundaries=None, delaunay=True, holes=True, border=[]):
        self.tris = poly_tri(pts, boundaries, delaunay, holes, border)
    def get_tris(self):
        return self.tris
        