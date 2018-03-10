#include <iostream>
#include <map>
#include <set>
#include <list>
#include <vector>
#include <array>
#include <cmath>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Geometry"

typedef std::array<int, 3> Triangle;
typedef std::list<Triangle> Triangles;
typedef Eigen::Vector2d Vec;
typedef std::vector<Eigen::Vector2d> Vecs;
typedef std::array<int, 2> Edge;  // ordered edges

typedef std::set<int> kEdge;   //key edge
typedef std::set<kEdge> Edges;

typedef std::list<std::list<int>> Boundary;

typedef std::map<Edge, Triangles> Edge2Tris;
typedef std::map<int, Triangles> Point2Tris;


class PolyTri
{
public:
    Triangles tris;
    Vecs pts;
    Boundary boundary;
    bool remove_holes;
    bool delaunay;
    std::list<int> isBoarder;
    std::map<int, int> order;
    std::map<int, int> unorder;

    PolyTri(Vecs pts, Boundary boundaries, bool remove_holes,
	     bool delaunay, std::list<int> isBoarder);
    bool is_intersecting(Edge edge1, Edge edge2);
    bool is_edge_visible(Edge edge, int id_pnt);
    double get_edge_area(Edge edge, Vec pnt);
    double get_edge_area(Edge edge);
    double get_area(Triangle& tri);
    void makeCounterClockwise(Triangle& tri);
    std::vector<int> sort_pts(Vec p0);
    
};