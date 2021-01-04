import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import open3d as o3d
import os
from sklearn.neighbors import KernelDensity
from KDEpy import FFTKDE
from scipy.signal import find_peaks
import sys

def OutlierGrade(input_path,fig_name):
    str_cmd="../build/MGS -i "+ input_path +" -o 1.csv -M OutlierGrade"
    os.system(str_cmd)
    file=pd.read_csv("1.csv")
    og=np.asarray(file)[:,0]
    dty=np.asarray(file)[:,1]
    x, y = FFTKDE(kernel="gaussian", bw=1).fit(og).evaluate()
    # print(np.max(y))
    # print("og=%s" % str(len(y[y<np.max(y)*0.5])/len(x)))    
    st=y>np.max(y)*0.5
    flag=0
    seg_start=0
    seg_length=0
    for i in range(len(st)):
        if st[i]==True and flag==0:
            flag=1
            seg_start=x[i]           
        if st[i]==False and flag==1:
            flag=0
            seg_length=seg_length+(x[i-1]-seg_start)            
    og=seg_length/50    
    # print("og:  %.2f" % (x[y==np.max(y)]/np.max(x)))
    plt.plot(x,y,'r')    
    plt.savefig(fig_name)
    return og

def Homogeneity(input_path):
    str_cmd="../build/MGS -i "+ input_path +" -o 2.csv -M Homogeneity"
    os.system(str_cmd)
    file=pd.read_csv("2.csv")
    uty=np.asarray(file)[:,0]
    x, y = FFTKDE(kernel="gaussian", bw=0.01).fit(uty).evaluate()        
    plt.rcParams.update({'font.size': 22})
    plt.plot(x,y)
    peaks, _ = find_peaks(y, height=0.01)
    plt.plot(np.array([x[peaks]]), y[peaks]+0.3*np.ones([1,peaks.shape[0]]), "rv", markersize=8)  
    fname=input_path.split("/")[-1].split(".")[0]    
    plt.savefig("fig/Homogeneity_"+fname+".png")
    # print("Hm: %d" % (len(peaks)))
    return len(peaks)

def Sigularity(input_path,fig_name="1.png"):
    str_cmd="../build/MGS -i "+ input_path +" -o 3.csv -M Singularity"
    os.system(str_cmd)
    df=pd.read_csv("3.csv")
    df=np.asarray(df)
    print("Sg: %.2f" % (df[0,0]))
    return df[0,0]
    # sg=df.to_numpy()[:,0]
    # x, y = FFTKDE(kernel="gaussian", bw=0.01).fit(sg).evaluate()
    # plt.plot(x,y)
    # plt.savefig(fig_name)
    # print(x[y==np.max(y)])
    # if sg[0,0] > 0.13:
    #     print("sg:  %.2f    --> 4th-style" % sg[0,0])
    # else:
    #     print("sg:  %.2f    --> 3rd-style" % sg[0,0])




# if __name__ == '__main__':
#     input_path=[]
#     output_path=[]
#     mode=[]
#     for i in range(len(sys.argv)):
#         if "-i"==sys.argv[i]:
#             input_path=sys.argv[i+1]
#         elif "-o"==sys.argv[i]:
#             output_path=sys.argv[i+1]
#         elif "-M"==sys.argv[i]:
#             mode=sys.argv[i+1]
    
#     # input_path="/home/llg/Experiment/dtu/m051/stl051_total.ply"
#     # mode="og"

#     if "og"==mode:        
#         og=OutlierGrade(input_path,"og.png")
#         print("%s og: %.2f" % (input_path.split("/")[-1],og))
#     elif "hm"==mode:
#         hm=Homogeneity(input_path)
#         print("%s hm: %.2f" % (input_path.split("/")[-1],hm))
#     elif "sg"==mode:
#         sg=Sigularity(input_path)
#         print("%s sg: %.2f" % (input_path.split("/")[-1],sg))
#     elif "eval3"==mode:
#         rst_sg=Sigularity(input_path)
#         if rst_sg>4:
#             print("4th-style")
#         else:
#             og=OutlierGrade(input_path,"og.png")
#             if og<0.1:
#                 print("1st-style")
#             else:
#                 hm=Homogeneity(input_path)
#                 if hm<=1:
#                     print("2nd-style")
#                 else:
#                     print("3rd-style")
#     elif "eval2"==mode:        
#         og=OutlierGrade(input_path,"og.png")
#         if og<0.1:
#             print("%s: 1st-style, og %.2f" % (input_path.split("/")[-1],og))
#         else:
#             hm=Homogeneity(input_path)
#             if hm<=1:
#                 print("%s: 3nd-style, og %.2f, hm %.2f" % (input_path.split("/")[-1],og,hm))
#             else:
#                 print("%s: 2rd-style, og %.2f, hm %.2f" % (input_path.split("/")[-1],og,hm))

    # Style 1
    # Sigularity("/home/llg/dataset/3Dlib/Ignatius.ply")
    # Sigularity("/home/llg/dataset/3Dlib/Ignatius_COLMAP.ply")
# OutlierGrade("/home/llg/dataset_paper/Ignatius_COLMAP.ply","1.png")
OutlierGrade("/home/llg/dataset_paper/Ignatius.ply","1.png")
    # Homogeneity("/home/llg/dataset/3Dlib/Ignatius.ply")

    # Style 2
    # Sigularity("/home/llg/dataset/3Dlib/Ignatius_COLMAP.ply")
    # OutlierGrade("/home/llg/dataset/3Dlib/Ignatius_COLMAP.ply")
    # Homogeneity("/home/llg/dataset/3Dlib/Ignatius_COLMAP.ply")

    # Style 3
    # Sigularity("/home/llg/dataset/3Dlib/torch_points.ply")
    # OutlierGrade("/home/llg/dataset/3Dlib/torch_points.ply")
    # Homogeneity("/home/llg/dataset/3Dlib/torch_points.ply")

    # Style 4
    # Sigularity("/home/llg/dataset/3Dlib/dog_colmap.ply")
    # OutlierGrade("/home/llg/dataset/3Dlib/dog_colmap.ply")
    # Homogeneity("/home/llg/dataset/3Dlib/dog_colmap.ply")

    # OutlierGrade("/home/llg/Human/H_FRM_0145_clipped.ply","1_1.png")
    # OutlierGrade("/home/llg/Human/fused_20200114_01_FRM_0073_1920.ply","1_2.png")
    # OutlierGrade("/home/llg/dataset/3Dlib/Ignatius_COLMAP.ply","1_2.png")
    # OutlierGrade("/home/llg/dataset/3Dlib/Ignatius_COLMAP_clipped.ply","1_2.png")
    # OutlierGrade("/home/llg/dataset/3Dlib/cat_colmap.ply","1_2.png")
    # OutlierGrade("/home/llg/dataset/3Dlib/Ignatius.ply","1_2.png")
    # OutlierGrade("/home/llg/dataset/3Dlib/Ignatius_COLMAP.ply","1_2.png")
    # Homogeneity("/home/llg/dataset/3Dlib/Ignatius.ply")
    # Homogeneity("/home/llg/dataset/3Dlib/Ignatius_COLMAP.ply")
    # Homogeneity("/home/llg/dataset/3Dlib/Barn_COLMAP.ply")
    # Homogeneity("/home/llg/dataset/3Dlib/torch_points.ply")
    # Sigularity("/home/llg/dataset/3Dlib/Meetingroom_COLMAP.ply")


    # Sigularity("/home/llg/dataset/3Dlib/cat_colmap.ply","3_1.png")
    # Sigularity("/home/llg/dataset/3Dlib/dog_colmap.ply","3_2.png")
    # Sigularity("/home/llg/dataset/3Dlib/Ignatius.ply","3_3.png")
    # Sigularity("/home/llg/dataset/3Dlib/Ignatius_COLMAP.ply","3_4.png")