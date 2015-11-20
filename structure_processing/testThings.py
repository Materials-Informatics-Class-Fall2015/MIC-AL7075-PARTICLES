import ReadMSFunction as RMS
import MKStemp as mks
import PredictPlasticity as PP
import ParticleDist as PD
import matplotlib.pyplot as plt
import numpy as np
import cPickle

colors = ["#543005", "#8c510a", "#bf812d", "#dfc27d", "#80cdc1", "#35978f", "#01665e", "#003c30"]

def plotHistograms(values, filters):
    filter = 1
    val_filter = filters==filter
    thresh_num = 100
    thresh_filter = 6
    handles = []
    labels = []
    bins = 50
    range = (np.min(values), np.max(values))
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
        
        temp_f = open("histo_%d.p" % filter, "w")
        cPickle.dump((centers, hist), temp_f)
        temp_f.close()
        
        color = colors[(filter-1)%len(colors)]
        handle, = plt.plot(centers, hist, color=color)
        handles.append(handle)
        labels.append("Element Distance <= %d"  % filter)
        print("Finished distance %d" % filter)
        filter += 1
        val_filter = filters==filter
    plt.legend(handles,labels)
    plt.ylabel("Normalized Frequency")
    plt.xlabel("Fatemi-Socie FIP")
    plt.show()

pred_dir = "C:/Users/pkern3/Documents/MIC-AL7075-PARTICLES/large_predict"
ms = RMS.readDirectory(pred_dir)
ms = ms[-1]
ms = ms[:101,:101,:101]
strains = mks.TrainPredict(False, ["C:/Users/pkern3/Documents/MIC-AL7075-PARTICLES/training_data/", pred_dir])
strains = strains[:,-1,...]
fips = PP.predFIPs(ms, strains)
print("Finished FIPs")
distances = PD.getDistances(ms)
print("Finished Particle Dist")
plotHistograms(fips, distances)