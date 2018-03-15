import numpy as np
import copy
import numba

small = 1e-10

@numba.jit
def make_key(i1, i2):
    """
    Make a tuple key such at i1 < i2
    """
    if i1 < i2:
        return (i1, i2)
    return (i2, i1)


def poly_tri(pts_in, boundaries=None, delaunay=True, holes=True, boarder=[]):

    # data structures
    triangles = []  # cells
    edge2tris = {}  # edge to triangle(s) map
    point2triangles = {}
    boundary_edges = set()

    # compute center of gravity

    cg = sum(pts_in, np.zeros(2, np.float64)) / len(pts_in)
    i = 0
    while True:
        pts = pts_in[:]
        dist = pts - cg
        square_dist = dist.T[0]**2 + dist.T[1]**2
        order = np.argsort(square_dist)
        pts = pts[order]

        # create first triangle, make sure we're getting a non-zero area
        # otherwise drop the pts

        index = 0

        tri = [index, index + 1, index + 2]
        makeCounterClockwise(pts, tri)
        if getArea(pts, *tri) > small:
            triangles.append(tri)
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

    # add additional pts
    for i in range(3, len(pts)):
        addPoint(i, pts, edge2tris, triangles, boundary_edges,
                 point2triangles, delaunay)

    update_edge2tri(triangles, edge2tris, point2triangles)

    if boundaries:
        constraintBoundaries(pts, edge2tris, point2triangles, triangles, 
                             boundaries, unorder)
        if holes:
            remove_empty(pts, triangles)
            update_edge2tri(triangles, edge2tris, point2triangles)
            removeHoles(edge2tris, triangles, boundaries, boarder, unorder)
    
    return [order[tri] for tri in triangles]

@numba.jit
def getArea(pts, ip0, ip1, ip2):
    """
    Compute the parallelipiped area
    @param ip0 index of first vertex
    @param ip1 index of second vertex
    @param ip2 index of third vertex
    """
    d1 = pts[ip1] - pts[ip0]
    d2 = pts[ip2] - pts[ip0]
    return (d1[0]*d2[1] - d1[1]*d2[0])

@numba.jit
def edgeArea(pts, ordered_edge):
    d1, d2 = ordered_edge
    d1 = pts[d1]
    d2 = pts[d2]
    return (d1[0]*d2[1] - d1[1]*d2[0])

@numba.jit
def isEdgeVisible(pts, ip, edge):
    """
    Return true iff the point lies to its right when the edge pts down
    @param ip point index
    @param edge (2 point indices with orientation)
    @return True if visible
    """
    area = getArea(pts, ip, edge[0], edge[1])
    if area < small:
        return True
    return False

@numba.jit
def makeCounterClockwise(pts, ips):
    """
    Re-order nodes to ensure positive area (in-place operation)
    """
    area = getArea(pts, ips[0], ips[1], ips[2])
    if area < -small:
        ip1, ip2 = ips[1], ips[2]
        # swap
        ips[1], ips[2] = ip2, ip1

@numba.jit
def flipOneEdge(edge, pts, edge2tris, triangles):
    """
    Flip one edge then update the data structures
    @return set of edges to interate over at next iteration
    """

    # start with empty set
    res = set()

    # assume edge is sorted
    tris = edge2tris.get(edge, [])
    if len(tris) < 2:
            # nothing to do, just return
            return res

    iTri1, iTri2 = tris
    tri1 = triangles[iTri1]
    tri2 = triangles[iTri2]

    # find the opposite vertices, not part of the edge
    iOpposite1 = -1
    iOpposite2 = -1
    for i in range(3):
        if not tri1[i] in edge:
            iOpposite1 = tri1[i]
        if not tri2[i] in edge:
            iOpposite2 = tri2[i]

    # compute the 2 angles at the opposite vertices
    da1 = pts[edge[0]] - pts[iOpposite1]
    db1 = pts[edge[1]] - pts[iOpposite1]
    da2 = pts[edge[0]] - pts[iOpposite2]
    db2 = pts[edge[1]] - pts[iOpposite2]
    crossProd1 = getArea(pts, iOpposite1, edge[0], edge[1])
    crossProd2 = getArea(pts, iOpposite2, edge[1], edge[0])
    dotProd1 = np.dot(da1, db1)
    dotProd2 = np.dot(da2, db2)
    angle1 = abs(np.arctan2(crossProd1, dotProd1))
    angle2 = abs(np.arctan2(crossProd2, dotProd2))

    # Delaunay's test
    if angle1 + angle2 > np.pi*(1.0 + small):

        # flip the triangles
        #                         / ^ \                                        / b \
        # iOpposite1 + a|b + iOpposite2    =>     + - > +
        #                         \     /                                        \ a /

        newTri1 = [iOpposite1, edge[0], iOpposite2]  # triangle a
        newTri2 = [iOpposite1, iOpposite2, edge[1]]  # triangle b

        # update the triangle data structure
        triangles[iTri1] = newTri1
        triangles[iTri2] = newTri2

        # now handle the topolgy of the edges

        # remove this edge
        del edge2tris[edge]

        # add new edge
        e = make_key(iOpposite1, iOpposite2)
        edge2tris[e] = [iTri1, iTri2]

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

        # these two edges might need to be flipped at the
        # next iteration
        res.add(make_key(iOpposite1, edge[0]))
        res.add(make_key(iOpposite2, edge[1]))

    return res

@numba.jit
def flip_edges(pts, edge2tris, triangles):
    edgeSet = set(edge2tris.keys())

    continueFlipping = True

    while continueFlipping:
        newEdgeSet = set()
        for edge in edgeSet:
            newEdgeSet |= flipOneEdge(edge, pts, edge2tris, triangles)
        edgeSet = copy.copy(newEdgeSet)
        continueFlipping = (len(edgeSet) > 0)

@numba.jit
def addPoint(ip, pts, edge2tris, triangles, boundary_edges, point2triangles, delaunay):
    """
    Add point
    @param ip point index
    """

    # collection for later updates
    boundary_edges2remove = set()
    boundary_edges2add = set()

    for edge in boundary_edges:
        if isEdgeVisible(pts, ip, edge):

            # create new triangle
            newTri = [edge[0], edge[1], ip]
            newTri.sort()
            makeCounterClockwise(pts, newTri)
            triangles.append(newTri)
            iTri = len(triangles) - 1

            # add the two boundary edges
            e0 = make_key(*edge)
            e1 = make_key(ip, edge[0])
            e2 = make_key(edge[1], ip)
            for e in (e0, e1, e2):
                v = edge2tris.get(e, [])
                v.append(iTri)
                edge2tris[e] = v

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
        flip_edges(pts, edge2tris, triangles)

def create_boundary_list(boundaries, unorder, boarder=None, create_key=True):
    constrained_boundary = []
    for k, boundary in enumerate(boundaries):
        if boarder and k not in boarder:
            continue
        b = unorder[boundary]
        for i, j in zip(b[:-1], b[1:]):
            item = [i, j]
            if create_key:
                item = make_key(*item)
            constrained_boundary.append(item)
    return constrained_boundary

def constraintBoundaries(pts, edge2tris, point2triangles, triangles, 
                         boundaries, unorder):
    boundary = create_boundary_list(boundaries, unorder)
    tris2remove = set()  # nr
    tris2add = []
    for cb in boundary:
        removed_edges = set()
        if cb not in edge2tris.keys():
            pt1, pt2 = cb
            for tri in point2triangles[pt1]:
                edge = copy.copy(triangles[tri])
                edge.remove(pt1)
                edge = make_key(*edge)
                if isIntersecting(pts, edge, cb):
                    tris2remove.add(tri)
                    break
            edges = tri2edges(triangles[tri], create_key=True)
            edges.remove(edge)
            removed_edges.add(edges[0])
            removed_edges.add(edges[1])
            
            if pt2 in triangles[tri]:
                break
            for tri in edge2tris[edge]:
                if tri not in tris2remove:
                    break
            while True:
                tris2remove.add(tri)

                edges = tri2edges(triangles[tri])
                edges.remove(edge)
                for edge in edges:
                    if not pt2 in triangles[tri] and isIntersecting(pts, edge, cb):
                        tris2remove.add(tri)
                        edge_to_proceed = edge
                    else:
                        removed_edges.add(edge)
                if pt2 in triangles[tri]:
                    break
                edge = edge_to_proceed
                for tri in edge2tris[edge]:
                    if tri not in tris2remove:
                       break

            
            ul, ll = create_loop(pts, removed_edges, cb[0], cb[1])
            if len(ul) == 3:
                tris2add.append(ul)
            else:
                ul_cb = list(range(len(ul))) + [0]
                ul, ll = np.array(ul), np.array(ll)
                tris = poly_tri(pts[ul], [ul_cb], holes=True, delaunay=False)
                tris = [ul[tri] for tri in tris]
                tris2add += tris
            if len(ll) == 3:
                tris2add.append(ll)
            else:
                ll_cb = list(range(len(ll))) + [0]
                tris = poly_tri(pts[ll], [ll_cb], holes=True, delaunay=False)
                tris = [ll[tri] for tri in tris]
                tris2add += tris              
            
            # somehow we have to gather the edges!
            # and find a closed loop. Divide them by the edge which is
            # inserted and fill up the hole -> a full PolyTri with edge
            # constraining and boundary removal.
    tris2remove = list(tris2remove)
    tris2remove.sort()
    tris2remove.reverse()
    for tri in tris2remove:
        triangles.pop(tri)
    for tri in tris2add:
        triangles.append(list(tri))
    for tri in triangles:
        makeCounterClockwise(pts, tri)
    update_edge2tri(triangles, edge2tris, point2triangles)

@numba.jit
def update_edge2tri(triangles, edge2tris, point2triangles):
    edge2tris.clear()
    point2triangles.clear()
    for i, tri in enumerate(triangles):
        for edge in tri2edges(tri):
            e2t = edge2tris.get(edge, [])
            edge2tris[edge] = e2t + [i]
        for point in tri:
            p2t = point2triangles.get(point, [])
            point2triangles[point] = p2t + [i]
            
def remove_empty(pts, triangles):
    remove_tris = []
    for i, tri in enumerate(triangles):
        makeCounterClockwise(pts, tri)
        area = getArea(pts, *tri)
        if area < 1e-10:
            remove_tris.append(i)

    remove_tris.sort()
    remove_tris.reverse()
    for i in remove_tris:
        triangles.pop(i)

def tri2edges(tri, create_key=True):
    _tri = tri + [tri[0]]
    if create_key:
        return [make_key(*edge) for edge in zip(_tri[:-1], _tri[1:])]
    else:
        return [[*edge] for edge in zip(_tri[:-1], _tri[1:])]

def isIntersecting(pts, edge1, edge2):
    p11 = pts[edge1[0]]
    p12 = pts[edge1[1]]
    p21 = pts[edge2[0]]
    p22 = pts[edge2[1]]
    t = p12 - p11
    s = p22 - p21
    r = p21 - p11
    try:
        c1, c2 = np.linalg.inv(np.array([t, -s]).T) @ r
    except np.linalg.linalg.LinAlgError:
        return False
    return (0 < c1 < 1) and (0 < c2 < 1)

def removeHoles(edge2tris, triangles, boundaries, boarder, unorder):
    bs = create_boundary_list(boundaries, unorder, boarder)
    o_bs = create_boundary_list(boundaries, unorder, boarder, create_key=False)
    remove_edges = set()
    for b, o_b in zip(bs, o_bs):
        tris = edge2tris[b]
        for tri in tris:
            if o_b in tri2edges(triangles[tri], create_key=False):
                edges = tri2edges(triangles[tri])
                for edge in edges:
                    remove_edges.add(edge)
    for b in bs:
        if b in remove_edges:
            remove_edges.remove(b)
    tris2remove = set()
    for edge in remove_edges:
        for tri in edge2tris[edge]:
            tris2remove.add(tri)
    tris2remove = list(tris2remove)
    tris2remove.sort()
    tris2remove.reverse()
    for i in tris2remove:
        triangles.pop(i)

def create_loop(pts, edges, start, end):
    loop = []
    points = set()
    point2edge = {}
    for edge in edges:
        points.add(edge[0])
        points.add(edge[1])
        p2e = point2edge.get(edge[0], [])
        point2edge[edge[0]] = p2e + [edge]
        p2e = point2edge.get(edge[1], [])
        point2edge[edge[1]] = p2e + [edge]
    e0, e1 = point2edge[start]
    loop.append(start)
    p0 = start
    while True:
        # p1 is the other point of the edge
        p1 = e0[0]
        if p1 == p0:
            p1 = e0[1]
        loop.append(p1)
        e1 = point2edge[p1][0]
        if e1 == e0:
            e1 = point2edge[p1][1]
        e0 = e1
        p0 = p1
        if p0==start:
            break

    # loop.append(start)
    area = sum([edgeArea(pts, e) for e in zip(loop[:-1], loop[1:])])
    if area > 0:
        loop.reverse()
    
    lower_loop = []
    upper_loop = []
    insert_upper = True
    for i in loop:
        if insert_upper:
            upper_loop.append(i)
            if i == end:
                insert_upper=False
                lower_loop.append(i)
        else:
            lower_loop.append(i)
    return upper_loop, lower_loop
