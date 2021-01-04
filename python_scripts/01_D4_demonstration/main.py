import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import pywt
dist=pd.read_csv("../../Result/dist.txt")
dist=np.asarray(dist)
cA=pd.read_csv("../../Result/cA.txt")
cA=np.asarray(cA)
cD=pd.read_csv("../../Result/cD.txt")
cD=np.asarray(cD)

ridx=201
# my version
plt.subplot(3,1,1)
plt.plot(dist[ridx,:])
plt.subplot(3,1,2)
plt.plot(cA[ridx,:])
plt.subplot(3,1,3)
plt.plot(cD[ridx,:])
plt.savefig("1.png")

# lib version
dat=dist[ridx,:]
cA_lib, cD_lib = pywt.dwt(dat, 'db2')
plt.figure()
plt.subplot(3,1,1)
plt.plot(dat)
plt.subplot(3,1,2)
plt.plot(cA_lib)
plt.subplot(3,1,3)
plt.plot(cD_lib)
plt.savefig("2.png")