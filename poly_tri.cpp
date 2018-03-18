#include "poly_tri.h"
#include <algorithm>
#include <math.h>

// TODO:
// constraint boundaries
// remove_holes


double small = 1e-10;


Edge kedge2edge(kEdge edge)
{
    Edge e;
    int j = 0;
    for (auto i_pnt: edge)
    {
        e[j] = i_pnt;
        j++;
    }
    return e;
}

kEdge edge2kedge(Edge edge)
{
    return kEdge{edge[0], edge[1]};
}


std::vector<int> range(int l)
{
    std::vector<int> y(l);
    int n(0);
    std::generate(std::begin(y), std::end(y), [&]{ return n++; });
    return y;
}

PolyTri::PolyTri(Vecs pts, Boundaries boundaries, bool remove_holes, bool delaunay, std::list<int> is_border)
{
    pts = pts;
    boundaries = boundaries;
    _remove_holes = remove_holes;
    delaunay = delaunay;
    is_border = is_border;
    
// 1 order points by distance to a point
    Vec cg(0, 0);
    for (auto pnt: pts)
    cg += pnt;
    cg /= pts.size();
    
    Triangle tri0;
    bool found = false;
    for (int k=0; k<=pts.size(); k++)
    {
    std::vector<int> _order = sort_pts(cg);
    int i = 0;
 
// 2 create a map which points back to the input points
    for (auto element: _order)
    {
        order[i] = element;
        unorder[element] = i;
    }
 
// 3 first triangle
    tri0 = {0, 1, 2};
    
// 3.1 check if area > 0
    make_counter_clockwise(tri0);
    double area = get_area(tri0);
    if (area > small)
    {
        found = true;
        tris.push_back(tri0);
        break;
    }
    else
        cg = pts[k];
    }
    if (! found)
    std::cout << "throw exception here" << std::endl;
    return;
    
    
// 4 add_triangle
    // boundary edges
    Edge e01 = {tri0[0], tri0[1]};
    Edge e12 = {tri0[1], tri0[2]};
    Edge e20 = {tri0[2], tri0[0]};
    
    boundary_edges.insert(e01);
    boundary_edges.insert(e12);
    boundary_edges.insert(e20);

    kEdge ke01 = {tri0[0], tri0[1]};
    kEdge ke12 = {tri0[1], tri0[2]};
    kEdge ke20 = {tri0[2], tri0[0]};
 
    edge2tris[ke01] = std::vector<int> {0};
    edge2tris[ke12] = std::vector<int> {0};
    edge2tris[ke20] = std::vector<int> {0};
    
    
// 5 add other points
    for (int i=3; i < pts.size(); i++)
    {
        add_point(i);
    }
    
// 6 constraine boundaries
// 7 remove triangles
}

void PolyTri::add_point(int id_pnt)
{
    
    Edges boundary_edges2remove;
    Edges boundary_edges2add;
    for (Edge edge: boundary_edges)
    {
    if (is_edge_visible(edge, id_pnt))
    {
        // create new triangle
        Triangle new_tri{edge[0], edge[1], id_pnt};
        std::sort(new_tri.begin(), new_tri.end());
        make_counter_clockwise(new_tri);
        tris.push_back(new_tri);
        int iTri = tris.size() - 1;
        
        // update the edge to triangle map
            kEdge e0{edge[0], edge[1]};  // set
        kEdge e1{id_pnt, edge[0]};
        kEdge e2{edge[1], id_pnt};
        append_tri(e0, iTri);
        append_tri(e1, iTri);
        append_tri(e2, iTri);
        
        // keep track of the boundary edges to update
        boundary_edges2remove.insert(edge);
        boundary_edges2add.insert(Edge{edge[0], id_pnt});
        boundary_edges2add.insert(Edge{id_pnt, edge[1]});
    }
    }
    // update the boundary edges
    for (auto bedge: boundary_edges2remove)
    {
    boundary_edges.erase(bedge);
    }
    for (auto bedge: boundary_edges2add)
    {
    if (edge2tris[kEdge{bedge[0], bedge[1]}].size() == 1)
    {
        // only add boundary edge if it does not appear
        // twice in different order
            boundary_edges.insert(bedge);
    }
    }
    if (delaunay)  // recursively flip edges
    {
    flip_edges();
    
        for (int i=0; i < pts.size(); i++)
    {
            pnt2tris[i] = std::vector<int>{};
    }

        for (int j=0; j < tris.size(); j++)
    {
            for (int k: tris[j])
        {
                pnt2tris[k].push_back(j);
        }
    }
    }
}

void PolyTri::append_tri(kEdge edge, int iTri)
{
    if (edge2tris.find(edge) != edge2tris.end())
        edge2tris[edge].push_back(iTri);
    else
        edge2tris[edge] = std::vector<int>{iTri};  // empty triangle
}


kEdges PolyTri::flip_one_edge(kEdge edge)
{
    // start with empty set
    kEdges res;
    std::vector<int> v;
    kEdge _edge;
    std::array<int, 2> v_edge;
    
    int i = 0;
    for (int element: edge)
    {
    v_edge[i++] = element;
    }

    // assume edge is sorted
    if (! edge2tris.count(edge))
    {
    return res;
    }
    std::vector<int> i_tris = edge2tris[edge];
    
    if (i_tris.size() < 2)
    {
    // nothing to do, just return
    return res;
    }
    int iTri1 = i_tris[0];
    int iTri2 = i_tris[1];

    Triangle tri1 = tris[iTri1];
    Triangle tri2 = tris[iTri2];
    
    int iOpposite1 = -1;
    int iOpposite2 = -1;
    for (int i=0; i < 3; i++)
    {
    if (edge.find(tri1[i]) == edge.end())
        iOpposite1 = tri1[i];
    if (edge.find(tri2[i]) == edge.end())
        iOpposite2 = tri2[i];
    }
    // compute the 2 angles at the opposite vertices
    Vec da1 = pts[v_edge[0]] - pts[iOpposite1];
    Vec db1 = pts[v_edge[1]] - pts[iOpposite1];
    Vec da2 = pts[v_edge[0]] - pts[iOpposite2];
    Vec db2 = pts[v_edge[1]] - pts[iOpposite2];
    double crossProd1 = get_edge_area(v_edge, pts[iOpposite1]);
    double crossProd2 = - get_edge_area(v_edge, pts[iOpposite2]);
    double dotProd1 = da1.dot(db1);
    double dotProd2 = da2.dot(db2);
    double angle1 = abs(atan2(crossProd1, dotProd1));
    double angle2 = abs(atan2(crossProd2, dotProd2));
    
    // Delaunay's test
    if (angle1 + angle2 > M_PI*(1.0 + small))
    {

    // flip the triangles
    //                       / ^ \                                        / b \
    // iOpposite1 + a|b + iOpposite2    =>     + - > +
    //                      \     /                                        \ a /

    Triangle newTri1 {iOpposite1, v_edge[0], iOpposite2};  // triangle a
    Triangle newTri2 {iOpposite1, iOpposite2, v_edge[1]};  // triangle b

    // update the triangle data structure
    tris[iTri1] = newTri1;
    tris[iTri2] = newTri2;

    // now handle the topolgy of the edges

    // remove this edge
    edge2tris.erase(edge);

    // add new edge
    _edge = kEdge{iOpposite1, iOpposite2};
    edge2tris[_edge] = std::vector<int> {iTri1, iTri2};

    // modify two edge entries which now connect to
    // a different triangle
    _edge = kEdge{iOpposite1, v_edge[1]};
    v = edge2tris[_edge];
    for (int & i : v)
    {
        if (i == iTri1)
        i = iTri2;
    }
    res.insert(_edge);

    _edge = kEdge{iOpposite2, v_edge[0]};
    v = edge2tris[_edge];
    for (int & i : v)
    {
        if (i == iTri2)
        i = iTri1;
    }
    res.insert(_edge);

    // these two edges might need to be flipped at the
    // next iteration
    res.insert(kEdge{iOpposite1, v_edge[0]});
    res.insert(kEdge{iOpposite2, v_edge[1]});
    }
    return res;
}

void PolyTri::flip_edges()
{
    kEdges edge_set;
    bool continueFlipping = true;
    
    for (auto const& item: edge2tris)
    edge_set.insert(item.first);


    while (continueFlipping)
    {
    kEdges new_edge_set;
    for (auto edge: edge_set)
    {
        kEdges flipped_edges = flip_one_edge(edge);
        new_edge_set.insert(flipped_edges.begin(), flipped_edges.end());
    }

    edge_set = new_edge_set;
    continueFlipping = (edge_set.size() > 0);
    }
}

kEdges PolyTri::create_boundary_list(std::vector<int> border)
{
    kEdges constrained_boundary;
    for (int k; k < boundaries.size(); k++)
    {
        if (std::find(border.begin(), border.end(), k) != border.end())
            continue;
        Boundary boundary = boundaries[k];
        int b0 = unorder[boundary[0]];
        int b1;
        for (int i=1; i<boundaries.size(); i++)
        {
            b1 = unorder[boundary[i]];
            constrained_boundary.insert(kEdge{b0, b1});
            b0 = b1;
        }
    }
    return constrained_boundary;
}

Edges PolyTri::create_ordered_boundary_list(std::vector<int> border)
{
    Edges constrained_boundary;
    for (int k; k < boundaries.size(); k++)
    {
        if (std::find(border.begin(), border.end(), k) != border.end())
            continue;
        Boundary boundary = boundaries[k];
        int b0 = unorder[boundary[0]];
        int b1;
        for (int i=1; i<boundary.size(); i++)
        {
            b1 = unorder[boundary[i]];
            constrained_boundary.insert(Edge{b0, b1});
            b0 = b1;
        }
    }
    return constrained_boundary;
}


bool PolyTri::is_intersecting(Edge edge1, Edge edge2)
{
    Vec p11 = pts[edge1[0]];
    Vec p12 = pts[edge1[1]];
    Vec p21 = pts[edge2[0]];
    Vec p22 = pts[edge2[1]];
    Vec t = p12 - p11;
    Vec s = p22 - p21;
    Vec r = p21 - p11;
    Vec c;
    if ((1 - std::abs(t.dot(s) / t.norm() / s.norm())) < small)
    return false;  // parallel
    else
    {
    Eigen::Matrix<double, 2, 2> mat;
    mat << t, -s;
    mat = mat.transpose();  // necessary?
    c = mat.inverse() * r;
    }
    return (0 < c(0) < 1) and (0 < c(1) < 1);
}

bool PolyTri::is_intersecting(kEdge edge1, kEdge edge2)
{
    Edge e1 = kedge2edge(edge1);
    Edge e2 = kedge2edge(edge2);
    return is_intersecting(e1, e2);
}

double PolyTri::get_area(Triangle& tri)
{
    Vec p1 = pts[tri[0]];
    Vec p2 = pts[tri[1]];
    Vec p3 = pts[tri[2]];
    Vec d1 = p2 - p1;
    Vec d2 = p3 - p1;
    return d1[0]*d2[1] - d1[1]*d2[0];
}

double PolyTri::get_edge_area(Edge edge)
{
    Vec d1 = pts[edge[0]];
    Vec d2 = pts[edge[1]];
    return d1[0]*d2[1] - d1[1]*d2[0];
}

double PolyTri::get_edge_area(Edge edge, Vec pnt)
{
    Vec d1 = pts[edge[0]];
    Vec d2 = pts[edge[1]];
    d1 -= pnt;
    d2 -= pnt;
    return d1[0]*d2[1] - d1[1]*d2[0];
}

bool PolyTri::is_edge_visible(Edge edge, int id_pnt)
{
    Vec pnt = pts[id_pnt];
    double area = get_edge_area(edge, pnt);
        if (area < small)
            return true;
        return false;
}

void PolyTri::make_counter_clockwise(Triangle& tri)
{
    double area = get_area(tri);
    if (area < -small)
    {
    int i1, i2;
    i1 = tri[1];
    i2 = tri[2];
    tri[2] = i1;
    tri[1] = i2;
    }
}

void PolyTri::update_mapping()
{
    std::vector<int> _tris;
    edge2tris.clear();
    pnt2tris.clear();
    int i = 0;
    for (auto tri: tris)
    {
        for (auto edge: tri2edges(tri))
        {
            if (edge2tris.find(edge) != edge2tris.end())
                edge2tris[edge].push_back(i);
            else
                edge2tris[edge] = std::vector<int>{i};
        }
        for (auto point: tri)
        {
            if (pnt2tris.find(point) != pnt2tris.end())
                pnt2tris[point].push_back(i);
            else
                pnt2tris[point] = std::vector<int>{i};
        }
        i++;
    }
}

kEdges PolyTri::tri2edges(Triangle tri)
{
    kEdges edges;
    edges.insert(kEdge{tri[0], tri[1]});
    edges.insert(kEdge{tri[1], tri[2]});
    edges.insert(kEdge{tri[2], tri[0]});
    return edges;
}

void PolyTri::remove_empty()
{
    std::set<int>tris2remove;
    int i = 0;
    double area;
    for (auto tri: tris)
    {
        make_counter_clockwise(tri);
        area = get_area(tri);
        if (area < small)
            tris2remove.insert(i);
    }
    std::vector<int> reverse_vec(tris2remove.rbegin(), tris2remove.rend());
    for (int j : reverse_vec)
        tris.erase(tris.begin() + j);
}


std::vector<int> PolyTri::sort_pts(Vec p0)
{
    int i = 0;
    std::vector<double> diff(pts.size());
    std::generate(diff.begin(), 
          diff.end(), 
          [&]{return (pts[i++] - p0).norm(); } ); //generate distance
    std::vector<int> sort_ind(pts.size());
    std::sort(std::begin(sort_ind), 
          std::end(sort_ind), 
          [&](int i1, int i2) {return diff[i1] < diff[i2]; } );
    Vecs sorted_pts(pts.size());
    i = 0;
    std::generate(sorted_pts.begin(), sorted_pts.end(), [&]{return pts[sort_ind[i++]]; } );
    pts = sorted_pts;
    return sort_ind;    
}


Boundaries PolyTri::create_loop(kEdges edges, int start, int end)
{
        std::vector<int> loop;
        std::map<int, std::vector<std::vector<int>>>point2edge;
        for (auto edge: edges)
        {
            for (auto p: edge)
            {
                if (point2edge.find(p) != point2edge.end())
                    point2edge[p].push_back(std::vector<int>(edge.begin(), edge.end()));
                else
                    point2edge[p] = std::vector<std::vector<int>>{std::vector<int>(edge.begin(), edge.end())};
            }            
        }
        int p0 = start;
        std::vector<int> e0 = point2edge[p0][0];
        loop.push_back(start);
        int p1;
        std::vector<int> e1;
        double area = 0;
        while (true)
        {
            // p1 is the other point of the edge
            p1 = e0[0];
            if (p1 == p0)
                p1 = e0[1];
            loop.push_back(p1);
            e1 = point2edge[p1][0];
            if (e1 == e0)
                e1 = point2edge[p1][1];
            e0 = e1;
            p0 = p1;
            area += get_edge_area(Edge{p0, p1}, Vec{0, 0});
            if (p0 == start)
                break;
        }

        // loop.append(start)
        if (area > 0)
            loop = std::vector<int>(loop.rbegin(), loop.rend());
        std::vector<int> lower_loop, upper_loop;
        bool insert_upper = true;
        for (int i: loop)
        {
            if (insert_upper)
            {
                upper_loop.push_back(i);
                if (i == end)
                {
                    insert_upper = false;
                    lower_loop.push_back(i);
                }
            }
            else
                lower_loop.push_back(i);
        }
        return Boundaries{upper_loop, lower_loop};
}

PolyTri::remove_holes()
{
        kEdges bs = create_boundary_list(is_border);
        Edges o_bs = create_ordered_boundary_list(is_border);
        kEdges remove_edges;
        for (int i=0; i<bs.size(); i++)
        {
            Edge ob = 
            tris = self.edge2tris[b]
            for tri in tris:
                if o_b in self.tri2edges(self.tris[tri], create_key=False):
                    edges = self.tri2edges(self.tris[tri])
                    for edge in edges:
                        remove_edges.add(edge)
        }
        for b in bs:
            if b in remove_edges:
                remove_edges.remove(b)
        tris2remove = set()
        for edge in remove_edges:
            for tri in self.edge2tris[edge]:
                tris2remove.add(tri)
        tris2remove = list(tris2remove)
        tris2remove.sort()
        tris2remove.reverse()
        for i in tris2remove:
            self.tris.pop(i)
}
