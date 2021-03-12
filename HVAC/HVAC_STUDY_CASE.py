# Code made for Sergio Andrés Díaz Ariza
# 11 March 2021
# License MIT
# Introduction to Control: Python Program HAVOC Assigment

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

T = np.genfromtxt("Data1", delimiter=';')
P = np.genfromtxt("Data2", delimiter=';')


fig,ax = plt.subplots(2)
ax[0].plot(T[0],T[1],'r',alpha=0.75)
ax[1].plot(P[0],P[1],'b',alpha=0.75)
ax[0].set_ylim(20,30)
#ax[0].set_title("\t$vs$\t Time",fontsize= 16)
#ax[1].set_title("Power \t$vs$\t Time Surface",fontsize= 14)
ax[1].set_xlabel("Time \t$[s]$ ",fontsize= 14)
ax[0].set_ylabel("Temperature \t$[^{\circ}C]$ ",fontsize= 14)
ax[1].set_ylabel("Power\t$[W]$",fontsize= 14)
fig.suptitle('HVAC CONTROL\nTemperature $[^{\circ}C]$ & Power $[W]$', fontsize=16)
plt.show()
