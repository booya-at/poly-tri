#include <iostream>
#include <map>
#include <set>
#include <list>
#include <vector>
#include <array>
#include <cmath>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Geometry"

typedef std::array<int, 3> Tri;
typedef std::list<Tri> Tris;
typedef Eigen::Matrix<double, 2, 1> Vec;
typedef Eigen::Matrix<double, -1, 2> Vecs;
typedef std::array<int, 2> Edge;  // ordered edges

typedef std::set<int> kEdge;   //key edge
typedef std::set<kEdge> Edges;

typedef std::map<Edge, Tris> Edge2Tris;
typedef std::map<int, Tris> Point2Tris;

double small = 1e-10;


Tris polytris(Vecs pts, std::list<std::list<int>> boundaries, bool remove_holes, bool delaunay, 
	      std::list<int> isBoarder)
{
// 1 order points by distance to a point
// 2 create a map which points back to the input points
// 3 first triangle
// 4 add_point
// 5 constraine boundaries
// 6 remove triangles
}


bool isIntersecting(Vecs pts, Edge edge1, Edge edge2)
{
    Vec p11, p12, p21, p22, c;
    p11 << pts(edge1[0]);
    p12 << pts(edge1[1]);
    p21 << pts(edge2[0]);
    p22 << pts(edge2[1]);
    Vec t = p12 - p11;
    Vec s = p22 - p21;
    Vec r = p21 - p11;
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


