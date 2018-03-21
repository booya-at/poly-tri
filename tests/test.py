# -*- coding: utf-8 -*-

import unittest
import copy
import numpy as np
import matplotlib.pyplot as plt
from poly_tri import PolyTri

class TriangleTests(unittest.TestCase):
    @property
    def an_int(self):
        return np.random.randint(10, 30)

    def setUp(self):
        self.phi = np.linspace(0, 2 * np.pi, self.an_int)[:-1]
        self.outer = np.array([np.cos(self.phi), np.sin(self.phi)]).T
        self.inner = (self.outer * 0.5)
        self.points = np.array(list(self.outer) + list(self.inner))
        self.inner_bound = np.array(list(range(len(self.inner))) + [0])
        self.outer_bound = np.array(list(range(len(self.outer))) + [0])
        self.inner_bound +=  max(self.outer_bound) + 1
        self.outer_bound = self.outer_bound[::-1]
    
    def test_triangulization(self):
        tri = PolyTri(self.points, [self.inner_bound, self.outer_bound], delaunay=True, holes=True)
        plt.figure(figsize=(10, 10))
        plt.triplot(*tri.pts.T, tri.tris)
        plt.show()
    
    def test_constraint_edge(self):
        n = self.an_int
        x = np.linspace(0, 1, n)
        y = np.array([0] * n)
        pts = np.array([x, y]).T
        pts_upper = pts[1:-1].copy()
        pts_upper[:,1] += 0.1
        pts[1:-1, 1] -= 0.1
        cb = [[0, n-1]]
        additional_pts = np.array([[-1., 0], [2, 0.]])
        pts = np.array(list(pts) + list(pts_upper) + list(additional_pts))
        tri = PolyTri(pts, cb, holes=False, delaunay=True)
        plt.triplot(*pts.T, tri.get_tris())
        for i, p in enumerate(pts):
            plt.annotate(str(i), p)
        plt.show()

    def test_constraint_edge_1(self):
        n = self.an_int
        x = np.linspace(0, np.pi, n)
        y = abs(np.sin(x)) - 1.1
        pts_inner = list(np.array([x, y]).T)
        pts_inner.reverse()
        
        pts_outer = list(np.array(pts_inner) * np.array([1., -1.]))
        pts_outer.reverse()
        pts = pts_inner + pts_outer
        cb = [list(range(len(pts_inner)))]
        cb += [list(np.array(range(len(pts_outer))) + max(cb[0]) + 1)]
        pts += list(np.array([[0., 0.], [np.pi, 0.]]))
        cb += [[len(pts) - 2, len(pts) - 1]]
        pts = np.array(pts)
        tri = PolyTri(np.array(pts), boundaries=cb, holes=True,
                      border=[0, 1], delaunay=False)
        plt.triplot(*pts.T, tri.get_tris())
        plt.show()
    
    def test_constraint_edge_2(self):
        n = self.an_int
        x = np.linspace(0, np.pi, n)
        y = abs(np.sin(x)) - 1.1
        pts = list(np.array([x, y]).T)
        pts.reverse()
        pts += list(np.array([[0., 0.], [np.pi, 0.]]))
        pts = np.array(pts)
        tri = PolyTri(pts, boundaries=[list(range(len(pts))) + [0]], 
                      holes=True, delaunay=False)
        plt.triplot(*pts.T, tri.get_tris())
        plt.show()
    
    def test_easy(self):
        
        ul = np.array([0, 1, 2, 3, 4, 5, 6, 7, 0])
        pts = np.array([
                        [ 0.,          0.        ],
                        [ 0.6981317,  -0.45721239],
                        [ 1.04719755, -0.2339746 ],
                        [ 1.3962634,  -0.11519225],
                        [ 1.74532925, -0.11519225],
                        [ 2.0943951,  -0.2339746 ],
                        [ 2.44346095, -0.45721239],
                        [ 3.14159265,  0.        ]])
        tri = PolyTri(pts, [ul], holes=True, delaunay=False)
        plt.triplot(*pts.T, tri.get_tris())
        plt.show()
    
    def test_easy_1(self):
        ul = [1, 2, 3, 4, 5, 0, 1]
        pts = np.array([[0.,  0. ],
                        [0.2, 0.1],
                        [0.4, 0.1],
                        [0.6, 0.1],
                        [0.8, 0.1],
                        [1.,  0. ]])
        tri = PolyTri(pts, holes=False)
        plt.triplot(*pts.T, tri.get_tris())
        plt.show()
    
    def test_easy_2(self):
        pts = np.array([[-1, 0], [1, 0], [0., 0.5], [0., -0.5]])
        edge = [np.array([2, 3])]
        tri = PolyTri(pts, edge, holes=False, delaunay=False)
        plt.triplot(*pts.T, tri.get_tris())
        plt.show()
        
    def test_easy_3(self):
        pts = np.array([[-1., 0.], [1., 0.], [0., 0.5], [0., -0.5], [0., 1.]])
        edge = [np.array([0, 1])]
        tri = PolyTri(pts, holes=False, delaunay=False)
        plt.triplot(*pts.T, tri.get_tris())
        plt.show()
    
    def test_ellipse(self):
        outer_pts = np.array([np.cos(self.phi), np.sin(self.phi)]).T
        inner_pts = copy.copy(outer_pts)
        outer_pts *= np.array([2., 1.])
        inner_pts *= 0.5
        pts = np.array(list(inner_pts) + list(outer_pts))
        tri = PolyTri(pts,[list(range(len(inner_pts))) + [0]], delaunay=True, holes=True)
        plt.figure(figsize=(10, 10))
        plt.triplot(*pts.T, tri.get_tris())
        plt.show()
    
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
        tri = PolyTri(pts, [profile, hole], delaunay=False, holes=True)
        plt.triplot(*pts.T, tri.get_tris())
        for i, p in enumerate(tri.pts):
            plt.annotate(str(i), p)
        plt.show()
if __name__ == '__main__':
    unittest.main()