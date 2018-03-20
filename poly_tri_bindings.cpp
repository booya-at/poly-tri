#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "poly_tri.h"

PYBIND11_MODULE(poly_tri_cpp, m) {
    pybind11::class_<PolyTri>(m, "PolyTri")
        .def(pybind11::init<Vecs, Boundaries, bool, bool, std::vector<int>>(),
                pybind11::arg("pts"),
                pybind11::arg("boundaries")=Boundaries{},
                pybind11::arg("holes")=false,
                pybind11::arg("delaunay")=false,
                pybind11::arg("borders")=std::vector<int>{} )            
        .def("get_tris", &PolyTri::get_tris);
}