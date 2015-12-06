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
dir_recon = "C:\Users\pkern3\Documents\MIC-AL7075-PARTICLES\large_predict"
labels = []
handles = []
colors = ["#543005", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#f5f5f5", "#c7eae5", "#80cdc1", "#35978f", "#01665e", "#003c30"]
sample_slices = []
old_dir = os.getcwd()


j = 1
temp_path = os.path.join(dir_recon, str(j))
while(os.path.exists(temp_path)):
    corr2d, corr3d = C2P.get2Points(dir_img, temp_path)
    corr3d[0,:,:,:,1] /= corr3d[0,99,99,99,1]
    x = range(corr2d.shape[1])
    color = colors[j%len(colors)]
    temp, = plt.plot(x,corr3d[0,99,:,99,1],color)
    handles.append(temp)
    labels.append("Reconstruction %d" % j)
    sample_slices.append(corr3d[0,99,:,99,1])
    j += 1
    temp_path = os.path.join(dir_recon, str(j))
os.chdir(old_dir)
temp_f = open("CompareReconstructions.p", "w")
cPickle.dump(sample_slices, temp_f)
temp_f.close()
plt.legend(handles,labels)
plt.ylabel("Normalized Particle Autocorrelation")
plt.xlabel("X Distance")
plt.show()

