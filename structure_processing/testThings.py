import ReadMSFunction as RMS
import MKStemp as mks
import PredictPlasticity as PP
import ParticleDist as PD
import matplotlib.pyplot as plt
import numpy as np
import cPickle
import os
import gc
import scipy.ndimage.filters as filters

colors = ["#543005", "#8c510a", "#bf812d", "#dfc27d", "#80cdc1", "#35978f", "#01665e", "#003c30"]

def plotHistograms(values, filters, num=1):
    filter = 1
    val_filter = filters==filter
    thresh_num = 100
    thresh_filter = 6
    handles = []
    labels = []
    bins = 50
    #range = (np.min(values), np.max(values))
    range = (0, .001)
    while(val_filter.any()):
        temp_vals = values[val_filter]
        samples = len(temp_vals)
        if(samples < thresh_num and filter > 1) or filter > thresh_filter:
            break
        #q75, q25 = np.percentile(temp_vals, [75 ,25])
        #iqr = q75 - q25
        #bins = (np.max(temp_vals) - np.min(temp_vals))/(2*iqr*samples**(-1/3.0))
        hist, bins = np.histogram(temp_vals, bins)
        hist = np.asarray(hist, dtype=np.float_)
        hist /= float(np.sum(hist))
        centers = (bins[:-1] + bins[1:])/2
        
        temp_f = open("histo_%d_%d.p" % (filter, num), "w")
        cPickle.dump((centers, hist), temp_f)
        temp_f.close()
        
        print("finished distance %d" % filter)
        filter += 1
        val_filter = filters==filter
        
        color = colors[(filter-2)%len(colors)]
        handle, = plt.plot(centers, hist, color=color)
        handles.append(handle)
        labels.append("Element Distance (%d,%d]"  % (filter-2,filter-1))
    plt.legend(handles,labels)
    plt.ylabel("Normalized Frequency")
    plt.xlabel("Fatemi-Socie FIP")
    plt.savefig("Plot%d" % num)
    plt.clf()
    #plt.show()
    
loads = [[.002*-.3595,.002,.002*-.3595,0,0,0],[0,0,0,0,0,0.002]]
exists = True
start_dir = os.getcwd()
training_dir = "C:/Users/pkern3/Documents/MIC-AL7075-PARTICLES/training_data/"
#training_dir = os.path.join(start_dir, "mks_training")
i = 1
kernels = [0.25, 0.5, 0.75, 1]
while exists:
    pred_dir = "C:/Users/pkern3/Documents/MIC-AL7075-PARTICLES/large_predict/1"
    #pred_dir = os.path.join(start_dir, "mks_predictions", str(i))
    ms = RMS.readDirectory(pred_dir)
    ms = ms[-1]
    ms = ms[:100,:100,:100]
    distances = PD.getDistances(ms)
    print("Finished Particle Dist")
    for j in range(0,len(loads)):
        os.chdir(start_dir)
        strains = mks.TrainPredict(False, [training_dir, pred_dir], loads[j])
        strains = strains[:,-1,...]
        fips = PP.predFIPs(ms, strains)
        del strains
        gc.collect()
        #plotHistograms(fips, distances,j)
        mask = np.abs(ms.astype(float)-1)
        #print(mask[0,:5,:5])
        #print(fips[0,:5,:5])
        for k in range(len(kernels)):
            blurred_mask = filters.gaussian_filter(mask, kernels[k], mode='wrap')
            blurred_fips = filters.gaussian_filter(fips, kernels[k], mode='wrap')
            blurred_fips /= blurred_mask
            os.chdir(pred_dir)
            plotHistograms(blurred_fips, distances,j*len(kernels)+k)
            del blurred_fips
            del blurred_mask
            gc.collect()
        del fips
        del mask
        gc.collect()
    del ms
    del distances
    i += 1
    pred_dir = os.path.join(start_dir, "mks_predictions", str(i))
    exists = os.path.exists(pred_dir)