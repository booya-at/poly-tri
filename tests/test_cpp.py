# -*- coding: utf-8 -*-

import unittest
import copy
import numpy as np
import matplotlib.pyplot as plt
from poly_tri_cpp import PolyTri

class TriangleTests(unittest.TestCase):
    
    @unittest.skip
    def test_easy_2(self):
        pts = [[-1, 0], [1, 0], [0., 0.5], [0., -0.5]]
        edge = [[0, 1]]
        tri = PolyTri(pts, edge, delaunay=True)
        plt.triplot(*np.array(pts).T, tri.get_tris())
        plt.show()

    @unittest.skip
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


    @unittest.skip
    def test_easy_1(self):
        ul = [[1, 2, 3, 4, 0, 1]]
        pts = [[0.,  0. ],
               [0.2, 0.1],
               [0.5, 0.1],
               [0.8, 0.1],
               [1.,  0. ]]
        tri = PolyTri(pts, ul, delaunay=True)
        plt.triplot(*np.array(pts).T, tri.get_tris())
        print(tri.get_tris())
        for i, p in enumerate(pts):
            plt.annotate(str(i), p)
        plt.show()

    def test_easy3(self):
        pts =[
            [ 0., 0. ],
            [ 0.2, 0.5],
            [ 0.4, 0.7],
            [ 0.6, 0.7],
            [ 0.8, 0.5],
            [ 1.0, 0.]]
        tri = PolyTri(pts, holes=False, delaunay=False)
        plt.triplot(*np.array(pts).T, tri.get_tris())
        plt.show()


    @unittest.skip
    def test_easy4(self):
        ul = [[0, 1]]
        pts = [[0., 0. ],
               [1., 0.],
               [0.2, 0.1],
               [0.5, 0.1],
               [0.8, 0.1],
               [0.2, -0.1],
               [0.5, -0.1],
               [0.8, -0.1]]
        tri = PolyTri(pts, ul, holes=False, delaunay=True)
        print("#############", tri.get_tris())
        plt.triplot(*np.array(pts).T, tri.get_tris())
        plt.show()


    def test_ellipse(self):
        phi = np.linspace(0, np.pi * 2, 200)[:-1]
        outer_pts = np.array([np.cos(phi), np.sin(phi)]).T
        inner_pts = copy.copy(outer_pts)
        outer_pts *= np.array([1., 1.])
        inner_pts *= 0.5
        pts = np.array(list(inner_pts) + list(outer_pts))
        tri = PolyTri(pts.tolist(),[list(range(len(inner_pts))) + [0]], delaunay=True, holes=True)
        plt.figure(figsize=(10, 10))
        plt.triplot(*pts.T, tri.get_tris())
        plt.show()


    @unittest.skip
    def test_easy_5(self):
        pts = [[0, 0],
               [1, 2],
               [0, 1],
               [-1, 2]]
        tri = PolyTri(pts, boundaries=[[0, 3, 2, 1, 0], [4, 5, 6, 7, 4]], borders=[0], delaunay=True, holes=True)
        plt.triplot(*np.array(pts).T, tri.get_tris())
        for i, p in enumerate(pts):
            plt.annotate(str(i), p)
        plt.show()


    @unittest.skip
    def test_profile(self):
        profile = [0, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10,
                     9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
        hole = [24, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24]
        hole.reverse()
        pts = np.array([[1., 0.], 
                    [0.90196078, 0.01922532], 
                    [0.80392157, 0.0377041 ],
                    [0.70588235, 0.05526607], 
                    [0.60784314, 0.07167488], 
                    [0.50980392, 0.08654832], 
                    [0.41176471, 0.09930623],
                    [0.31372549, 0.10875978], 
                    [0.21568627, 0.11197072], 
                    [0.12254902, 0.10145591], 
                    [0.05446623, 0.07889333], 
                    [0.01361656, 0.04499772], 
                    [0., 0.], 
                    [ 0.01361656, -0.02548616], 
                    [ 0.05446623, -0.04430561], 
                    [ 0.12254902, -0.05642389], 
                    [ 0.21568627, -0.06297434], 
                    [ 0.31372549, -0.0661538 ], 
                    [ 0.41176471, -0.06701299], 
                    [ 0.50980392, -0.06625902], 
                    [ 0.60784314, -0.0629272 ], 
                    [ 0.70588235, -0.05471569], 
                    [ 0.80392157, -0.04154998], 
                    [ 0.90196078, -0.02334161], 
                    [0.22404773, 0.02395152], 
                    [0.21754519, 0.04396428], 
                    [0.20052133, 0.05633284], 
                    [0.17947867, 0.05633284], 
                    [0.16245481, 0.04396428], 
                    [0.15595227, 0.02395152], 
                    [0.16245481, 0.00393877], 
                    [ 0.17947867, -0.00842979], 
                    [ 0.20052133, -0.00842979], 
                    [0.21754519, 0.00393877]])
        tri = PolyTri(pts.tolist(), [profile, hole], delaunay=True, holes=True)
        plt.triplot(*pts.T, tri.get_tris())
        for i, p in enumerate(pts):
            plt.annotate(str(i), p)
        plt.show()

if __name__ == '__main__':
    unittest.main()