#include "poly_tri.h"
#include <algorithm>

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
    
    bool found = false;
    for (int k=0; k<pts.size(); k++)
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
	Triangle tri0{0, 1, 2};
    
// 3.1 check if area > 0
	makeCounterClockwise(tri0);
	double area = get_area(tri0);
	if (area > small)
	    found = true;
	    break;
	else
	    cg = pts[k];
    }
    if (! found)
	std::cout << "throw exception here" << std::endl;
	return;
    
    
// 4 add_point
// 5 constraine boundaries
// 6 remove triangles
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

void PolyTri::makeCounterClockwise(Triangle& tri)
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