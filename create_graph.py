import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import csv
#import pandas as pd

record_V1 =  np.genfromtxt("res.dat")

data1 = np.genfromtxt("phi_calres_01_00_04_mpi.dat")
data2 = np.genfromtxt("phi_calres_01_01_04_mpi.dat")
data3 = np.genfromtxt("phi_calres_01_02_04_mpi.dat")
data4 = np.genfromtxt("phi_calres_01_03_04_mpi.dat")
data5 = np.genfromtxt("phi_calres_02_00_04_mpi.dat")
data6 = np.genfromtxt("phi_calres_02_01_04_mpi.dat")
data7 = np.genfromtxt("phi_calres_02_02_04_mpi.dat")
data8 = np.genfromtxt("phi_calres_02_03_04_mpi.dat")
data9 = np.genfromtxt("phi_calres_03_00_04_mpi.dat")
data10 = np.genfromtxt("phi_calres_03_01_04_mpi.dat")
data11= np.genfromtxt("phi_calres_03_02_04_mpi.dat")
data12= np.genfromtxt("phi_calres_03_03_04_mpi.dat")
data13 = np.genfromtxt("phi_calres_04_00_04_mpi.dat")
data14 = np.genfromtxt("phi_calres_04_01_04_mpi.dat")
data15 = np.genfromtxt("phi_calres_04_02_04_mpi.dat")
data16= np.genfromtxt("phi_calres_04_03_04_mpi.dat")
graph_t2 = []
dt = 1*pow(10,-8)
n = 1
while(n<=int(1.2e07)):
  graph_t2.append(n*dt)
  n = n + 1
"""
graph_t2.insert(320,3.2*pow(10,-5))
graph_t2.insert(240,2.4*pow(10,-5))
graph_t2.insert(160,1.6*pow(10,-5))
graph_t2.insert(80,0.8*pow(10,-5))
"""
rel_err1 = []
for n in range(1,3000000):
  rel_err1.append( abs( (data1[n]-record_V1[n]) / record_V1[n]) )
for n in range(3000000):
  rel_err1.append( abs( (data2[n]-record_V1[int(n+3000000)]) / record_V1[int(n+3000000)]) )
for n in range(3000000):
  rel_err1.append( abs( (data3[n]-record_V1[int(n+6000000)]) / record_V1[int(n+6000000)]) )
for n in range(3000001):
  rel_err1.append( abs( (data4[n]-record_V1[int(n+9000000)]) / record_V1[int(n+9000000)]) )

rel_err2 = []
for n in range(1,3000000):
  rel_err2.append( abs( (data5[n]-record_V1[n]) / record_V1[n]) )
for n in range(3000000):
  rel_err2.append( abs( (data6[n]-record_V1[int(n+3000000)]) / record_V1[int(n+3000000)]) )
for n in range(3000000):
  rel_err2.append( abs( (data7[n]-record_V1[int(n+6000000)]) / record_V1[int(n+6000000)]) )
for n in range(3000001):
  rel_err2.append( abs( (data8[n]-record_V1[int(n+9000000)]) / record_V1[int(n+9000000)]) )

rel_err3 = []
for n in range(1,3000000):
  rel_err3.append( abs( (data9[n]-record_V1[n]) / record_V1[n]) )
for n in range(3000000):
  rel_err3.append( abs( (data10[n]-record_V1[int(n+3000000)]) / record_V1[int(n+3000000)]) )
for n in range(3000000):
  rel_err3.append( abs( (data11[n]-record_V1[int(n+6000000)]) / record_V1[int(n+6000000)]) )
for n in range(3000001):
  rel_err3.append( abs( (data12[n]-record_V1[int(n+9000000)]) / record_V1[int(n+9000000)]) )

rel_err4 = []
for n in range(1,3000000):
  rel_err4.append( abs( (data13[n]-record_V1[n]) / record_V1[n]) )
for n in range(3000000):
  rel_err4.append( abs( (data14[n]-record_V1[int(n+3000000)]) / record_V1[int(n+3000000)]) )
for n in range(3000000):
  rel_err4.append( abs( (data15[n]-record_V1[int(n+6000000)]) / record_V1[int(n+6000000)]) )
for n in range(3000001):
  rel_err4.append( abs( (data16[n]-record_V1[int(n+9000000)]) / record_V1[int(n+9000000)]) )

#グラフを表示する領域を，figオブジェクトとして作成．
plt.rcParams["font.size"] = 18
#fig, ax = plt.figure(figsize = (12,8))
fig, ax = plt.subplots(2,2,figsize = (12,8))

l1,l2,l3,l4 = "iteration = 1","iteration = 2","iteration = 3","iteration = 4" # 各ラベル

ax[0,0].plot(graph_t2,rel_err1 , color='r',label=l1)
ax[0,1].plot(graph_t2,rel_err2 , color='r',label=l2)
ax[1,0].plot(graph_t2,rel_err3 , color='r',label=l3)
ax[1,1].plot(graph_t2,rel_err4 , color='r',label=l4)
#plt.legend()
fig.supxlabel("Time(sec)")
fig.supylabel("Relative error")

ax[0, 0].set_yscale('log')
ax[0, 1].set_yscale('log')
ax[1, 0].set_yscale('log')
ax[1, 1].set_yscale('log')


ax[0][0].set_ylim(1e-17,1e-7)
ax[0][1].set_ylim(1e-17,1e-7)
ax[1][0].set_ylim(1e-17,1e-7)
ax[1][1].set_ylim(1e-17,1e-7)

#plt.ylim([1e-16, 1e-8]) dt6*10^-7
#plt.ylim([1e-16, 1e-7])
#plt.ylim([1e-16, 3e-3])

plt.savefig("relative_error.png")

plt.show()