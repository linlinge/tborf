import numpy as np
from sklearn.neighbors.kde import KernelDensity,KDTree
import open3d as o3d
import matplotlib.pyplot as plt
from numba import cuda
from KDEpy import FFTKDE


# estimate pdf using KDE with gaussian kernel
# @cuda.jit
def pts_entropy(ctmp):
    kde = KernelDensity(kernel='gaussian',bandwidth=0.8).fit(ctmp)
    log_p = kde.score_samples(ctmp)  # returns log(p) of data sample
    p = np.exp(log_p)                # estimate p of data sample
    entropy = -np.sum(p*log_p)       # evaluate entropy
    return entropy

def pts_fft_kde(ctmp):
    grid,point= FFTKDE(kernel='gaussian').fit(ctmp).evaluate(ctmp)
    print(grid)

# @cuda.jit('void(int32[:], int32[:])',device=True)
def sigularity(ply_pts,K=30):
    tree=KDTree(ply_pts)
    N=ply_pts.shape[0]
    for i in range(1,N+1):
        dist,idx=tree.query(ply_pts[:i],K)
        ctmp=ply_pts[idx,:][0]
        # print(ctmp)
        # ety=pts_entropy(ctmp)
        pts_fft_kde(ctmp)
        return 0

ply=o3d.io.read_point_cloud("/home/llg/dataset/4_Mimosa.ply")
ply_xyz=np.asarray(ply.points)
print(ply_xyz.shape)
sigularity(ply_xyz)
