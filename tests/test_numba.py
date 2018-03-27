# -*- coding: utf-8 -*-

import unittest
import copy
import numpy as np
import matplotlib.pyplot as plt
from poly_tri_numba import poly_tri

class TriangleTests(unittest.TestCase):
    @property
    def an_int(self):
        return np.random.randint(9, 200)

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
        tri = poly_tri(self.points, [self.inner_bound, self.outer_bound], delaunay=False, holes=True)
        plt.figure(figsize=(10, 10))
        plt.triplot(*self.points.T, tri)
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
        tri = poly_tri(pts, cb, holes=False, delaunay=False)
        plt.triplot(*pts.T, tri)
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
        try:
            tri = poly_tri(np.array(pts), cb, border=[0, 1], holes=True, delaunay=False)
        except Exception as e:
            print(pts, cb)
            raise e
            
        plt.triplot(*pts.T, tri)
        plt.show()


    def test_constraint_edge_2(self):
        n = self.an_int
        x = np.linspace(0, np.pi, n)
        y = abs(np.sin(x)) - 1.1
        pts = list(np.array([x, y]).T)
        pts.reverse()
        pts += list(np.array([[0., 0.], [np.pi, 0.]]))
        pts = np.array(pts)
        tri = poly_tri(pts, boundaries=[list(range(len(pts))) + [0]], 
                      holes=True, delaunay=False)
        plt.triplot(*pts.T, tri)
        plt.show()

    def test_easy3(self):
        pts =np.array([
            [ 0., 0. ],
            [ 0.2, 0.5],
            [ 0.4, 0.7],
            [ 0.6, 0.7],
            [ 0.8, 0.5],
            [ 1.0, 0.]])
        tri = poly_tri(pts, holes=False, delaunay=False)
        plt.triplot(*pts.T, tri)
        plt.show()
       

    def test_easy_1(self):
        ul = [1, 2, 3, 4, 0, 1]
        pts = np.array([[0.,  0. ],
                        [0.2, 0.1],
                        [0.5, 0.1],
                        [0.8, 0.1],
                        [1.,  0. ]])
        tri = poly_tri(pts, [ul], holes=False, delaunay=True)
        plt.triplot(*pts.T, tri)
        plt.show()


    def test_easy_2(self):
        pts = np.array([[-1, 0], [1, 0], [0., 0.5], [0., -0.5]])
        edge = [np.array([2, 3])]
        tri = poly_tri(pts, edge, holes=False, delaunay=False)
        plt.triplot(*pts.T, tri)
        plt.show()


    def test_easy_3(self):
        pts = np.array([[-1., 0.], [1., 0.], [0., 0.5], [0., -0.5], [0., 1.]])
        edge = [np.array([0, 1])]
        tri = poly_tri(pts, holes=False, delaunay=False)
        plt.triplot(*pts.T, tri)
        plt.show()

    def test_ellipse(self):
        outer_pts = np.array([np.cos(self.phi), np.sin(self.phi)]).T
        inner_pts = copy.copy(outer_pts)
        outer_pts *= np.array([1., 1.])
        inner_pts *= 0.5
        pts = np.array(list(inner_pts) + list(outer_pts))
        tri = poly_tri(pts,[list(range(len(inner_pts))) + [0]], delaunay=True, holes=True)
        plt.figure(figsize=(10, 10))
        plt.triplot(*pts.T, tri)
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
        tri = poly_tri(pts, [hole], delaunay=True, holes=False)
        plt.triplot(*pts.T, tri)
        for i, p in enumerate(pts):
            plt.annotate(str(i), p)
        plt.show()
    

    def test_profile2(self):
        pts = np.array([
                 [ 1.        ,  0.        ],
                 [ 0.93079585,  0.01364559],
                 [ 0.8615917 ,  0.02693948],
                 [ 0.79238754,  0.03983196],
                 [ 0.72318339,  0.05224763],
                 [ 0.65397924,  0.06411011],
                 [ 0.58477509,  0.07533044],
                 [ 0.51557093,  0.08571854],
                 [ 0.44636678,  0.09508809],
                 [ 0.37716263,  0.10309254],
                 [ 0.30795848,  0.10914777],
                 [ 0.23875433,  0.11209604],
                 [ 0.16961802,  0.10886508],
                 [ 0.10855553,  0.09829529],
                 [ 0.06106249,  0.08212217],
                 [ 0.02713888,  0.06031021],
                 [ 0.00678472,  0.03289393],
                 [ 0.        , 0.        ],
                 [ 0.00678472, -0.01869534],
                 [ 0.02713888, -0.03406389],
                 [ 0.06106249, -0.04609344],
                 [ 0.10855553, -0.05478371],
                 [ 0.16961802, -0.06033743],
                 [ 0.23875433, -0.06399067],
                 [ 0.30795848, -0.06603619],
                 [ 0.37716263, -0.06691995],
                 [ 0.44636678, -0.06690913],
                 [ 0.51557093, -0.06617305],
                 [ 0.58477509, -0.06412989],
                 [ 0.65397924, -0.0596694 ],
                 [ 0.72318339, -0.05275214],
                 [ 0.79238754, -0.04337563],
                 [ 0.8615917 , -0.03145862],
                 [ 0.93079585, -0.0170007 ],
                 [ 0.22424641,  0.02420158],
                 [ 0.21770592,  0.04433111],
                 [ 0.20058272,  0.05677185],
                 [ 0.17941728,  0.05677185],
                 [ 0.16229408,  0.04433111],
                 [ 0.15575359,  0.02420158],
                 [ 0.16229408,  0.00407205],
                 [ 0.17941728, -0.00836869],
                 [ 0.20058272, -0.00836869],
                 [ 0.21770592,  0.00407205],
                 [ 0.43347349,  0.01676735],
                 [ 0.42708062,  0.03644257],
                 [ 0.41034388,  0.04860253],
                 [ 0.38965612,  0.04860253],
                 [ 0.37291938,  0.03644257],
                 [ 0.36652651,  0.01676735],
                 [ 0.37291938, -0.00290788],
                 [ 0.38965612, -0.01506784],
                 [ 0.41034388, -0.01506784],
                 [ 0.42708062, -0.00290788],
                 [ 0.67016282,  0.0026589 ],
                 [ 0.66535715,  0.01744923],
                 [ 0.65277574,  0.02659016],
                 [ 0.63722426,  0.02659016],
                 [ 0.62464285,  0.01744923],
                 [ 0.61983718,  0.0026589 ],
                 [ 0.62464285, -0.01213144],
                 [ 0.63722426, -0.02127236],
                 [ 0.65277574, -0.02127236],
                 [ 0.66535715, -0.01213144]])
        boundaries = [[0, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 
                       15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0], [34, 35, 36, 37, 38,
                       39, 40, 41, 42, 43, 34], [44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 44], [54, 
                       55, 56, 57, 58, 59, 60, 61, 62, 63, 54]]

        tri = poly_tri(pts, boundaries, delaunay=False, holes=True)
        plt.triplot(*pts.T, tri)
        plt.show()


    def test_profile3(self):
        pts = np.array([
            [ 1.        ,  0.        ],
            [ 0.88235294,  0.02298626],
            [ 0.76470588,  0.0448548 ],
            [ 0.64705882,  0.06526509],
            [ 0.52941176,  0.08372709],
            [ 0.41176471,  0.09930623],
            [ 0.29411765,  0.11005194],
            [ 0.17647059,  0.10961876],
            [ 0.07843137,  0.08925774],
            [ 0.01960784,  0.05267707],
            [ 0.        ,  0.        ],
            [ 0.01960784, -0.02975812],
            [ 0.07843137, -0.04997027],
            [ 0.17647059, -0.06078761],
            [ 0.29411765, -0.06573109],
            [ 0.41176471, -0.06701299],
            [ 0.52941176, -0.06593259],
            [ 0.64705882, -0.06023152],
            [ 0.76470588, -0.04741638],
            [ 0.88235294, -0.02739101],
            [ 0.22420494,  0.02415623],
            [ 0.21767237,  0.04426139],
            [ 0.20056991,  0.05668706],
            [ 0.17943009,  0.05668706],
            [ 0.16232763,  0.04426139],
            [ 0.15579506,  0.02415623],
            [ 0.16232763,  0.00405107],
            [ 0.17943009, -0.0083746 ],
            [ 0.20056991, -0.0083746 ],
            [ 0.21767237,  0.00405107],
            [ 0.43345312,  0.016748  ],
            [ 0.42706414,  0.03641125],
            [ 0.41033758,  0.04856381],
            [ 0.38966242,  0.04856381],
            [ 0.37293586,  0.03641125],
            [ 0.36654688,  0.016748  ],
            [ 0.37293586, -0.00291525],
            [ 0.38966242, -0.01506781],
            [ 0.41033758, -0.01506781],
            [ 0.42706414, -0.00291525],
            [ 0.67018389,  0.00262844],
            [ 0.6653742 ,  0.01743116],
            [ 0.65278225,  0.02657975],
            [ 0.63721775,  0.02657975],
            [ 0.6246258 ,  0.01743116],
            [ 0.61981611,  0.00262844],
            [ 0.6246258 , -0.01217428],
            [ 0.63721775, -0.02132286],
            [ 0.65278225, -0.02132286],
            [ 0.6653742 , -0.01217428]
            ])
        boundaries = [[0, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0], [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 20],
                      [30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 30], [40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 40]]
        tri = poly_tri(pts, boundaries, delaunay=True, holes=True)
        plt.triplot(*pts.T, tri)
        plt.show()

if __name__ == '__main__':
    unittest.main()


