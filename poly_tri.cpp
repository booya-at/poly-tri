#include "poly_tri.h"
#include <algorithm>
#include <math.h>

double small = 1e-10;

std::vector<int> range(int l)
{
    std::vector<int> y(l);
    int n(0);
    std::generate(std::begin(y), std::end(y), [&]{ return n++; });
    return y;
}

PolyTri::PolyTri(Vecs pts, Boundary boundary, bool remove_holes, bool delaunay, std::list<int> isBoarder)
{
    pts = pts;
    boundary = boundary;
    remove_holes = remove_holes;
    delaunay = delaunay;
    isBoarder = isBoarder;
    
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
    
    kEdges boundary_edges2remove;
    kEdges boundary_edges2add;
    for (Edge edge: boundary_edges)
    {
	if (is_edge_visible(edge, id_pnt))
	{
	    // create new triangle
	    Triangle new_tri{edge[0], edge[1], id_pnt};
	    std::sort(new_tri.begin(), new_tri.end());
	    make_counter_clockwise(new_tri);
	    tris.push_back(new_tri);
	    
	    // update the edge to triangle map
            kEdge e{edge[0], edge[1]};  // set
	    int iTri = tris.size() - 1;
	    edge2tris[e].push_back(iTri);
	    
	    // add the two boundary edges
	    kEdge e1{id_pnt, edge[0]};
	    kEdge e2{edge[1], id_pnt};
	    append_tri(e1, iTri);
	    append_tri(e2, iTri);
	    
	    // keep track of the boundary edges to update
	    boundary_edges2remove.insert(kEdge{edge[0], edge[1]});
	    boundary_edges2add.insert(kEdge{edge[0], id_pnt});
	    boundary_edges2add.insert(kEdge{id_pnt, edge[1]});
	}
    }
    // update the boundary edges
    for (auto bedge: boundary_edges2remove)
    {
	boundary_edges.erase(bedge);
    }
    for (auto bedge: boundary_edges2add)
    {
	if (edge2tris[bedge].size() == 1) // does work? we removed the sorting here!
	{
	    // only add boundary edge if it does not appear
	    // twice in different order
            self.boundary_edges.add(bedge)
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

void PolyTri::flipEdges()
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