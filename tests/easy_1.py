# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from poly_tri_py import PolyTri

ul = [[1, 2, 3, 4, 0, 1]]
pts = np.array([[0.,  0. ],
        [0.2, 0.1],
        [0.5, 0.1],
        [0.8, 0.1],
        [1.,  0. ]])
tri = PolyTri(pts, delaunay=False)
print(tri.get_tris())