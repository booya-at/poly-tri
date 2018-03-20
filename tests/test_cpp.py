# -*- coding: utf-8 -*-

import unittest
import copy
import numpy as np
#import matplotlib.pyplot as plt
from poly_tri import PolyTri

class TriangleTests(unittest.TestCase):
    def test_easy_2(self):
        pts = [[-1, 0], [1, 0], [0., 0.5], [0., -0.5]]
        edge = [[0, 1]]
        tri = PolyTri(pts)
        print(tri.get_tris())

    def test_constraint_edge(self):
        n = 10
        x = np.linspace(0, 1, n)
        y = np.array([0] * n)
        pts = np.array([x, y]).T
        pts_upper = pts[1:-1].copy()
        pts_upper[:,1] += 0.1
        pts[1:-1, 1] -= 0.1
        cb = [[0, n-1]]
        additional_pts = np.array([[-1., 0], [2, 0.]])
        pts = np.array(list(pts) + list(pts_upper) + list(additional_pts))
        print(np.array(pts))
        tri = PolyTri(pts.tolist())
        print(tri.get_tris())
        #plt.triplot(*pts.T, tri.get_tris())
        #for i, p in enumerate(pts):
            #plt.annotate(str(i), p)
        #plt.show()

if __name__ == '__main__':
    unittest.main()