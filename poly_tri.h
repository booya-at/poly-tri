#include <iostream>
#include <map>
#include <set>
#include <list>
#include <vector>
#include <array>
#include <cmath>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Geometry"

typedef Eigen::Vector2d Vec;
typedef std::array<int, 3> Triangle;
typedef std::vector<Triangle> Triangles;
typedef std::vector<Eigen::Vector2d> Vecs;
typedef std::array<int, 2> Edge;  // ordered edges
typedef std::set<int> kEdge;   //key edge
typedef std::set<Edge> Edges;
typedef std::set<kEdge> kEdges;
typedef std::vector<int> Boundary;
typedef std::vector<std::vector<int>> Boundaries;
typedef std::map<kEdge, std::vector<int>> Edge2Tris;
typedef std::map<int, std::vector<int>> Point2Tris;


class PolyTri
{
public:
    Triangles tris;
    Vecs pts;
    Boundaries boundaries;
    bool remove_holes;
    bool delaunay;
    Edges boundary_edges;
    std::vector<int> isBoarder;
    std::map<int, int> order;
    std::map<int, int> unorder;
    Edge2Tris edge2tris;
    Point2Tris point2triangles;

    PolyTri(Vecs pts, Boundaries boundaries, bool remove_holes,
	     bool delaunay, std::list<int> isBoarder);
    void add_point(int id_pnt);
    bool is_intersecting(Edge edge1, Edge edge2);
    bool is_edge_visible(Edge edge, int id_pnt);
    double get_edge_area(Edge edge, Vec pnt);
    double get_edge_area(Edge edge);
    double get_area(Triangle& tri);
    void make_counter_clockwise(Triangle& tri);
    std::vector<int> sort_pts(Vec p0);
    void flip_edges();
    kEdges flip_one_edge(kEdge edge);
    void append_tri(kEdge edge, int iTri);
    kEdges create_boundary_list(std::vector<int> boarder);
    Edges create_ordered_boundary_list(std::vector<int> boarder);
};