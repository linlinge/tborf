import pandas as pd
import matplotlib.pyplot as plt 
dat=pd.read_csv("../Result/PDF.csv")
plt.plot(dat)
# plt.plot(0.1,0.5,'*')
plt.show()
