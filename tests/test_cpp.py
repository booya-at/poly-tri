# -*- coding: utf-8 -*-

import unittest
import copy
import numpy as np
import matplotlib.pyplot as plt
from poly_tri_cpp import PolyTri

class TriangleTests(unittest.TestCase):
    def test_easy_2(self):
        pts = [[-1, 0], [1, 0], [0., 0.5], [0., -0.5]]
        edge = [[0, 1]]
        tri = PolyTri(pts, edge, delaunay=True)
        plt.triplot(*np.array(pts).T, tri.get_tris())
        plt.show()

    def test_constraint_edge(self):
        n = 3
        x = np.linspace(0, 1, n)
        y = np.array([0] * n)
        pts = np.array([x, y]).T
        pts_upper = pts[1:-1].copy()
        pts_upper[:,1] += 0.1
        pts[1:-1, 1] -= 0.1
        cb = [[0, n-1]]
        additional_pts = np.array([[-1., 0], [2, 0.]])
        pts = np.array(list(pts) + list(pts_upper) + list(additional_pts))
        tri = PolyTri(pts.tolist(),  delaunay=True)
        plt.triplot(*pts.T, tri.get_tris())
        for i, p in enumerate(pts):
            plt.annotate(str(i), p)
        plt.show()
        
    def test_easy_1(self):
        ul = [[1, 2, 3, 4, 0, 1]]
        pts = [[0.,  0. ],
               [0.2, 0.1],
               [0.5, 0.1],
               [0.8, 0.1],
               [1.,  0. ]]
        tri = PolyTri(pts, ul, delaunay=False)
        plt.triplot(*np.array(pts).T, tri.get_tris())
        print(tri.get_tris())
        for i, p in enumerate(pts):
            plt.annotate(str(i), p)
        plt.show()
        
    def test_easy3(self):
        
        ul = [[0, 1, 2, 3, 4, 5, 6, 7, 0]]
        pts =[
            [ 0.,          0.        ],
            [ 0.6981317,  -0.45721239],
            [ 1.04719755, -0.2339746 ],
            [ 1.3962634,  -0.11519225],
            [ 1.74532925, -0.11519225],
            [ 2.0943951,  -0.2339746 ],
            [ 2.44346095, -0.45721239],
            [ 3.14159265,  0.        ]]
        tri = PolyTri(pts, ul, holes=False, delaunay=False)
        plt.triplot(*np.array(pts).T, tri.get_tris())
        plt.show()
        
    def test_easy4(self):
        ul = [[0, 1]]
        pts = [[0., 0. ],
               [1., 0.],
               [0.2, 0.1],
               [0.4, 0.1],
               [0.5, 0.1],
               [0.8, 0.1],
               [0.2, -0.1],
               [0.4, -0.1],
               [0.5, -0.1],
               [0.8, -0.1]]
        tri = PolyTri(pts, ul, holes=False, delaunay=False)
        print("#############", tri.get_tris())
        plt.triplot(*np.array(pts).T, tri.get_tris())
        plt.show()

    def test_ellipse(self):
        pts =  [[1, 0],
                [0, 1],
                [-1, 0],
                [0, -1],
                [2, 0],
                [0, 2],
                [-2, 0],
                [0, -2]]

        # tri = PolyTri(pts.tolist(), [list(range(len(inner_pts))) + [0]], holes=False, delaunay=False)
        tri = PolyTri(pts, [[0, 1, 2, 3], [4, 5, 6, 7]], holes=True, delaunay=False)
        plt.figure(figsize=(10, 10))
        plt.triplot(*np.array(pts).T, tri.get_tris())
        plt.show()

if __name__ == '__main__':
    unittest.main()