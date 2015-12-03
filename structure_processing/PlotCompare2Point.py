import ReadMSFunction as RMS
import MKStemp as mks
import PredictPlasticity as PP
import ParticleDist as PD
import matplotlib.pyplot as plt
import Compare2Point as C2P
import numpy as np
import cPickle
import os

dir_img = "C:\Users\pkern3\Documents\MIC-AL7075-PARTICLES\img\Microstructures_Images"
dir_recon = "C:\Users\pkern3\Documents\MIC-AL7075-PARTICLES\large_predict\\random"
labels = ["Scan 1", "Scan 2", "Scan 3", "Scan 4", "Scan 5", "Scan 6", "Reconstruction"]

corr2d, corr3d = C2P.get2Points(dir_img, dir_recon)
for i in range(corr2d.shape[0]):
    corr2d[i,:,:,0] /= corr2d[i,99,99,0]
    corr2d[i,:,:,1] /= corr2d[i,99,99,1]
corr3d[0,:,:,:,0] /= corr3d[0,99,99,99,0]
corr3d[0,:,:,:,1] /= corr3d[0,99,99,99,1]
x = range(corr2d.shape[1])
handles = []
colors = ['r','g','m','b','c','m']
for i in range(corr2d.shape[0]):
    temp, = plt.plot(x,corr2d[i,99,:,0],colors[i])
    handles.append(temp)
temp, = plt.plot(x,corr3d[0,99,:,99,0])
handles.append(temp)
plt.legend(handles,labels)
plt.ylabel("Normalized Matrix Autocorrelation")
plt.xlabel("X Distance")
plt.show()
handles = []
for i in range(corr2d.shape[0]):
    temp, = plt.plot(x,corr2d[i,99,:,1],colors[i])
    handles.append(temp)
temp, = plt.plot(x,corr3d[0,99,:,99,1])
handles.append(temp)
plt.legend(handles,labels)
plt.ylabel("Normalized Particle Autocorrelation")
plt.xlabel("X Distance")
plt.show()