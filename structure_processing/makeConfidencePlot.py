import matplotlib.pyplot as plt
import numpy as np
import cPickle
import os
import Compare2Point as C2P
from scipy import stats

handles = []

temp_f = open("CompareReconstructions.p", "r")
sample_slices = cPickle.load(temp_f)
temp_f.close()
for i in range(len(sample_slices)):
    sample_slices[i] = np.expand_dims(sample_slices[i], axis=0)
samples = np.concatenate(sample_slices, axis=0)
confidence = np.zeros((2,samples.shape[1]))
for i in range(samples.shape[1]):
    temp = samples[:,i]
    mu, sigma = np.mean(temp), np.std(temp)
    if(sigma>0):
        temp = stats.norm.interval(0.95, loc=mu, scale=sigma)
        confidence[0,i] = temp[0]
        confidence[1,i] = temp[1]
    else:
        confidence[:,i] = mu
for i in range(confidence.shape[0]):
    temp, = plt.plot(confidence[i], 'k--')
    handles.append(temp)
    
dir_img = "C:\Users\pkern3\Documents\MIC-AL7075-PARTICLES\img\Microstructures_Images"
dir_recon = "C:\Users\pkern3\Documents\MIC-AL7075-PARTICLES\large_predict"
labels = ["Reconstruction CI Upper", "Reconstruction CI Lower", "Scan 1", "Scan 2", "Scan 3", "Scan 4", "Scan 5", "Scan 6"]

corr2d, corr3d = C2P.get2Points(dir_img, dir_recon)
for i in range(corr2d.shape[0]):
    corr2d[i,:,:,0] /= corr2d[i,99,99,0]
    corr2d[i,:,:,1] /= corr2d[i,99,99,1]
x = range(corr2d.shape[1])
colors = ['r','g','m','b','c','m']
for i in range(corr2d.shape[0]):
    temp, = plt.plot(x,corr2d[i,99,:,1],colors[i])
    handles.append(temp)
plt.legend(handles,labels)
plt.ylabel("Normalized Particle Autocorrelation")
plt.xlabel("X Distance")
plt.show()
