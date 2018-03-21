#include "poly_tri.h"
#include <algorithm>
#include <math.h>


double small = 1e-10;


Edge kedge2edge(kEdge edge)
{
    Edge e;
    int j = 0;
    for (auto i_pnt: edge)
        e[j++] = i_pnt;
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

Vec add(Vec vec1, Vec vec2)
{
    Vec result;
    result[0] = vec1[0] + vec2[0];
    result[1] = vec1[1] + vec2[1];
    return result;
}

Vec sub(Vec vec1, Vec vec2)
{
    Vec result;
    result[0] = vec1[0] - vec2[0];
    result[1] = vec1[1] - vec2[1];
    return result;
}


void iadd(Vec& vec1, Vec vec2)
{
    vec1[0] += vec2[0];
    vec1[1] += vec2[1];
}

void isub(Vec& vec1, Vec vec2)
{
    vec1[0] -= vec2[0];
    vec1[1] -= vec2[1];
}


void idiv(Vec& vec1, double scalar)
{
    vec1[0] /= scalar;
    vec1[1] /= scalar;
}

double dot(Vec vec1, Vec vec2)
{
    return vec1[0] * vec2[0] + vec1[1] * vec2[1];
}

double norm(Vec vec)
{
    return std::sqrt(dot(vec, vec));
}

void print_pts(Vecs pts)
{
    for (Vec pnt: pts)
    {
        for (auto p: pnt)
            std::cout << p << ", ";
        std::cout << std::endl;
    }
}

PolyTri::PolyTri(Vecs pts, Boundaries boundaries, bool holes, 
                 bool delaunay, std::vector<int> is_border)
{
    this->boundaries = boundaries;
    this->delaunay = delaunay;
    this->is_border = is_border;
    
// 1 order points by distance to a point
    Vec cg{0, 0};
    for (auto pnt: this->pts)
        iadd(cg, pnt);
    idiv(cg, pts.size());
    
    Triangle tri0;
    bool found = false;
    for (int k=0; k<=pts.size(); k++)
    {
        this->pts = pts;
        std::vector<int> _order = sort_pts(cg);
        int i = 0;
    
// 2 create a map which points back to the input points
        for (auto element: _order)
        {
            std::cout << "order: " << element;
            order[i] = element;
            unorder[element] = i;
            i++;
        }
    
// 3 first triangle
        tri0 = {0, 1, 2};
        
// 3.1 check if area > 0
        make_counter_clockwise(tri0);
        double area = get_area(tri0);
        std::cout << "area = ";
        std::cout << area << std::endl;
        if (area > small)
        {
            found = true;
            tris.push_back(tri0);
            break;
        }
        else
            cg = this->pts[k];
    }
    if (! found)
    {
        std::cout << "throw exception here" << std::endl;
        return;
    }
    
    
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
    
    std::cout << "cg" <<  cg[0] << cg[1] << std::endl;
    print_pts(this->pts);
    std::cout << "tri0: " << tri0[0] << tri0[1] << tri0[2] << std::endl;
 
    edge2tris[ke01] = std::vector<int> {0};
    edge2tris[ke12] = std::vector<int> {0};
    edge2tris[ke20] = std::vector<int> {0};
    
    std::cout << "size of points: ";
    std::cout << pts.size() << std::endl;
    
// 5 add other points
    for (int i=3; i < pts.size(); i++)
    {
        add_point(i);
    }
    remove_empty();
    update_mapping();
// 6 constraine boundaries
    if (boundaries.size() > 0)
    {
        constraint_boundaries();
        if (holes)
        {
            remove_empty();
            update_mapping();
            remove_holes();
        }
    }
// 7 remove triangles
}

void PolyTri::add_point(int id_pnt)
{
    
    Edges boundary_edges2remove;
    Edges boundary_edges2add;
    for (Edge edge: boundary_edges)
    {
//         std::cout << "all: " << id_pnt << edge[0] << edge[1] << std::endl;
        if (is_edge_visible(edge, id_pnt))
        {
//             std::cout << "vis: " << id_pnt << edge[0] << edge[1] << std::endl;
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
        kEdge kbedge(bedge.begin(), bedge.end());
        if (edge2tris[kbedge].size() == 1)
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
        edge2tris[edge] = std::vector<int>{iTri};
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
    Vec da1 = sub(pts[v_edge[0]], pts[iOpposite1]);
    Vec db1 = sub(pts[v_edge[1]], pts[iOpposite1]);
    Vec da2 = sub(pts[v_edge[0]], pts[iOpposite2]);
    Vec db2 = sub(pts[v_edge[1]], pts[iOpposite2]);
    double crossProd1 = get_edge_area(v_edge, pts[iOpposite1]);
    double crossProd2 = - get_edge_area(v_edge, pts[iOpposite2]);
    double dotProd1 = dot(da1, db1);
    double dotProd2 = dot(da2, db2);
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
    for (int k=0; k < boundaries.size(); k++)
    {
        if (std::find(border.begin(), border.end(), k) != border.end())
            continue;
        Boundary boundary = boundaries[k];
        std::cout << "boundary" << std::endl;
        for (auto j: boundary)
            std::cout << j << std::endl;
        int b0 = unorder[boundary[0]];
        int b1;
        for (int i=1; i<boundary.size(); i++)
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
    Vec t = sub(p12, p11);
    Vec s = sub(p22, p21);
    Vec r = sub(p21, p11);
    Vec c;
    if ((1 - std::abs(dot(t, s) / std::sqrt(dot(t, t)) / std::sqrt(dot(s, s)))) < small)
        return false;  // parallel
    else
    {
        double a, b, c, d, e, f, c0, c1;
        a = t[0];
        b = -s[0];
        c = t[1];
        d = -s[1];
        // maybe transpose?
        e = r[0];
        f = r[1];
        c0 = -b*f/(a*d - b*c) + d*e/(a*d - b*c);
        c1 = a*f/(a*d - b*c) - c*e/(a*d - b*c);
        return 0 < c0 && c0 < 1 && 0 < c1 && c1 < 1;
    }
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
    Vec d1 = sub(p2, p1);
    Vec d2 = sub(p3, p1);
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
    isub(d1, pnt);
    isub(d2, pnt);
    return d1[0]*d2[1] - d1[1]*d2[0];
}

bool PolyTri::is_edge_visible(Edge edge, int id_pnt)
{
    Vec pnt = pts[id_pnt];
    double area = get_edge_area(edge, pnt);
    if (area < -small)
        return true;
    else
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
        for (auto edge: tri2ordered_edges(tri))
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

Edges PolyTri::tri2edges(Triangle tri)
{
    Edges edges;
    edges.insert(Edge{tri[0], tri[1]});
    edges.insert(Edge{tri[1], tri[2]});
    edges.insert(Edge{tri[2], tri[0]});
    return edges;
}

kEdges PolyTri::tri2ordered_edges(Triangle tri)
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
    std::generate(diff.begin(), diff.end(), [&]{return norm(sub(pts[i++], p0)); } ); //generate distance
    std::vector<int> sort_ind = range(pts.size());
    std::sort(sort_ind.begin(), sort_ind.end(), [&](int i1, int i2) {return diff[i1] < diff[i2]; } );
    
    Vecs sorted_pts(pts.size());
    i = 0;
    std::generate(sorted_pts.begin(), sorted_pts.end(), [&]{return pts[sort_ind[i++]]; });

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

void PolyTri::remove_holes()
{
        kEdges _o_bs = create_boundary_list(is_border);
        Edges _bs = create_ordered_boundary_list(is_border);
        // create a vector from the set to use [] operator
        std::vector<kEdge> o_bs(_o_bs.begin(), _o_bs.end());
        std::vector<Edge> bs(_bs.begin(), _bs.end());
        kEdges remove_edges;
        for (int i=0; i<bs.size(); i++)
        {
            kEdge ob = o_bs[i];
            Edge b = bs[i];
            std::vector<int> _tris = edge2tris[ob];
            for (int tri: _tris)
            {
                Edges ordered_edges = tri2edges(tris[tri]);
                if (ordered_edges.find(b) != ordered_edges.end())
                {
                    for (auto edge: tri2ordered_edges(tris[tri]))
                        remove_edges.insert(edge);
                }
            }
        }
        
        // we don't want to remove edges on the boundary
        for (auto b: o_bs)
        {
            if (remove_edges.find(b) != remove_edges.end())
                remove_edges.erase(b);
        }
        std::set<int> tris2remove;
        for (auto edge: remove_edges)
        {
            for (auto tri: edge2tris[edge])
                tris2remove.insert(tri);
        }
        std::vector<int> reversed_tris2remove(tris2remove.rbegin(), tris2remove.rend());
        for (auto i: reversed_tris2remove)
            tris.erase(tris.begin() + i);
}

void PolyTri::constraint_boundaries()
{
    kEdges boundary = create_boundary_list(std::vector<int> {});
    std::set<int> tris2remove;
    Triangles tris2add;
    kEdges edges;
    Triangle tr;
    int tri_index;
    kEdge edge, edge2proceed;
    
//     std::cout << "size of pts: " << pnt2tris.size() << std::endl;
//     std::cout << "size of tris: " << tris.size() << std::endl;
//     for (const auto &test: pnt2tris)
//         std::cout << test.first << std::endl;
    for (kEdge cb: boundary)
    {
        kEdges removed_edges;
        // check existence of every boundary
        if (edge2tris.find(cb) == edge2tris.end())   
        {

            // if the edge is not existent, we have to delete intersecting edges
            // and the connected triangles
            // afterwards we create new triangles from the remaining empty loop
            std::vector<int> _cb(cb.begin(), cb.end());
            
            std::cout << "cb: ";
            for (int e: _cb)
                std::cout << e << ", ";
            std::cout << std::endl;
            
            // we start removing triangles at this point
            int pt1 = _cb[0];
     
            // once we reach this point we are done with tri removal
            int pt2 = _cb[1];

            for (auto tri: pnt2tris[pt1])
            {
                // get the edge which is not connected to pt1
                tr = tris[tri];
                edge.clear();
                int index = 0;
                for (auto pt: tr)
                {
                    if (pt != pt1)
                        edge.insert(pt);
                }
                
                // if this edge is intersecting with the cb
                // (cb = not yet applied), then we can proceed
                // otherwise choose another triangle connected to
                // pt1
                if (is_intersecting(edge, cb))
                {
                    tris2remove.insert(tri);
                    tri_index = tri;
                    break;
                }
            }
            // one triangle was added to the removal list
            // now we have to 
            edges = tri2ordered_edges(tr);
            edges.erase(edge);
            removed_edges.insert(edges.begin(), edges.end());
            if (std::find(tr.begin(), tr.end(), pt2) != tr.end()) 
                {
                    // if we find the end point we can stop the triangle removal for this edge
                    break;
                }
            for (auto tri: edge2tris[edge])
                {
                    if (tris2remove.find(tri) == tris2remove.end())
                    {
                        tri_index = tri;
                        break;
                    }
                }
                
            while (true)
            {
//                 std::cout << "currrent edge: ";
//                 for (int e: edge)
//                     std::cout << e << ", ";
//                 std::cout << std::endl;
                tris2remove.insert(tri_index);
                tr = tris[tri_index];
                
                std::cout << "current triangle: ";
                for (int e: tr)
                    std::cout << e << ", ";
                std::cout << std::endl;
                
                edges = tri2ordered_edges(tr);
                edges.erase(edge);

                // removed triangles form a polygon which can be 
                // represented by a closed loop (connected edges)
                // here we gather these edges
                for (kEdge edge: edges)
                {
                    if (is_intersecting(edge, cb))
                        edge2proceed = edge;
                    else
                        removed_edges.insert(edge);
                }
                
                if (std::find(tr.begin(), tr.end(), pt2) != tr.end()) 
                {
                    // if we find the end point we can stop the triangle removal for this edge
                    break;
                }
                else
                {
                    // as we haven't found the end point of the edge cb (constraint boundary)
                    // we have to look at the next triangle which has a intersecting edge
                    edge = edge2proceed;
                }
                
                for (auto tri: edge2tris[edge])
                {
                    if (tris2remove.find(tri) == tris2remove.end())
                    {
                        tri_index = tri;
                        break;
                    }
                }
            }
            Boundaries loops = create_loop(removed_edges, _cb[0], _cb[1]);
            for (auto loop: loops)
            {
                if (loop.size() == 3)
                {
                    Triangle _tr;
                    std::copy_n(loop.begin(), 3, _tr.begin());
                    tris2add.push_back(_tr);
                }
                else
                {
                    std::vector<int> loop_cb = range(loop.size());
                    std::cout << "loop= " << std::endl;
                    loop_cb.push_back(0);
                    for (auto l: loop_cb)
                        std::cout << l << ", " << std::endl;
                    Vecs new_points;
                    for (auto el_i: loop)
                        new_points.push_back(pts[el_i]);
                    for (Triangle tri: PolyTri(new_points, Boundaries{loop_cb}, true, false, std::vector<int>{}).get_tris())
                    {
                        for (int &t: tri)
                            t = loop[t];
                        tris2add.push_back(tri);
                    }
                }
            }
        }
    }
            
    // somehow we have to gather the edges!
    // and find a closed loop. Divide them by the edge which is
    // inserted and fill up the hole -> a full PolyTri with edge
    // constraining and boundary removal.
    std::vector<int> _tris2remove(tris2remove.rbegin(), tris2remove.rend());
    for (auto tri: _tris2remove)
        tris.erase(tris.begin() + tri);
    for (auto tri: tris2add)
        tris.push_back(tri);
    for (auto tri: tris)
        make_counter_clockwise(tri);
    update_mapping();
}

Triangles PolyTri::get_tris()
{
    Triangles reorder;
    for (Triangle tri: tris)
    {
        Triangle reordered_tri;
        int i = 0;
        for (int index: tri)
            reordered_tri[i++] = order[index];
        reorder.push_back(reordered_tri);
    }
    return reorder;
}
