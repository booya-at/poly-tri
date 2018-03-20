cmake_minimum_required(VERSION 2.6)
project(poly_tri)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANGXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANGXX)


find_package(PythonInterp ${py} REQUIRED)
find_package(PythonLibs ${py} REQUIRED)
find_package(pybind11 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})


pybind11_add_module(poly_tri_cpp poly_tri_bindings.cpp poly_tri.cpp poly_tri.h)
add_executable(main main.cpp poly_tri.cpp)


execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
"from distutils.sysconfig import get_python_lib; import os; print(get_python_lib())"
OUTPUT_VARIABLE _PYTHON_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE )

string(REPLACE "\\" "/" PYTHON_INSTALL_DIR ${_PYTHON_INSTALL_DIR})


install(TARGETS poly_tri_cpp DESTINATION ${PYTHON_INSTALL_DIR})