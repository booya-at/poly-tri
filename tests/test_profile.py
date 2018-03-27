from poly_tri_py import PolyTri
import matplotlib.pyplot as plt
import numpy as np

pts = np.array([
 [ 1.        ,  0.        ],
 [ 0.94117647,  0.01161274],
 [ 0.88235294,  0.02298626],
 [ 0.82352941,  0.03408675],
 [ 0.76470588,  0.0448548 ],
 [ 0.70588235,  0.05526607],
 [ 0.64705882,  0.06526509],
 [ 0.58823529,  0.0747821 ],
 [ 0.52941176,  0.08372709],
 [ 0.47058824,  0.09194376],
 [ 0.41176471,  0.09930623],
 [ 0.35294118,  0.1054852 ],
 [ 0.29411765,  0.11005194],
 [ 0.23529412,  0.11214142],
 [ 0.17647059,  0.10961876],
 [ 0.12254902,  0.10145591],
 [ 0.07843137,  0.08925774],
 [ 0.04411765,  0.0730155 ],
 [ 0.01960784,  0.05267707],
 [ 0.00490196,  0.02828742],
 [ 0.        ,  0.        ],
 [ 0.00490196, -0.016024  ],
 [ 0.01960784, -0.02975812],
 [ 0.04411765, -0.04108445],
 [ 0.07843137, -0.04997027],
 [ 0.12254902, -0.05642389],
 [ 0.17647059, -0.06078761],
 [ 0.23529412, -0.06385111],
 [ 0.29411765, -0.06573109],
 [ 0.35294118, -0.06672856],
 [ 0.41176471, -0.06701299],
 [ 0.47058824, -0.06672403],
 [ 0.52941176, -0.06593259],
 [ 0.58823529, -0.0639694 ],
 [ 0.64705882, -0.06023152],
 [ 0.70588235, -0.05471569],
 [ 0.76470588, -0.04741638],
 [ 0.82352941, -0.03831459],
 [ 0.88235294, -0.02739101],
 [ 0.94117647, -0.0146261 ],
 [ 0.22433824,  0.02435338],
 [ 0.19      ,  0.05869161],
 [ 0.15566176,  0.02435338],
 [ 0.19      , -0.00998486],
 [ 0.43349963,  0.01679296],
 [ 0.4       ,  0.05029259],
 [ 0.36650037,  0.01679296],
 [ 0.4       , -0.01670667],
 [ 0.67019211,  0.00261792],
 [ 0.645     ,  0.02781002],
 [ 0.61980789,  0.00261792],
 [ 0.645     , -0.02257419]])

boundaries = [[0, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
			  [40, 41, 42, 43, 40],
			  [44, 45, 46, 47, 44],
			  [48, 49, 50, 51, 48]]

tris = PolyTri(pts, boundaries, holes=True, delaunay=False).get_tris()
plt.triplot(*pts.T, tris)