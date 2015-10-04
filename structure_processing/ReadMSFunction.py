import os
from pymks import PrimitiveBasis
from pymks.stats import correlate
from pymks.tools import draw_correlations
import numpy as np

def readDirectory(directory):
	prim_basis = PrimitiveBasis(n_states=2, domain=[1, 2])
	filenames=[]
	for f in os.listdir(directory):
    	print(f)
        if(f.find("trial_Phases")!=-1 and f.find("_old")==-1):
            filenames.append(f)
    ms_list = [0]*len(filenames)
    for f in filenames:
        num = int(f[f.rfind("_")+1:f.rfind(".txt")])
        elem_f = os.path.join(directory, "trial_elem_grains_%d.txt" % num)
        elem_f = open(elem_f, "r")
        elem_f.readline()
        line = elem_f.readline()
        shape = map(int, line.split(",")[:3])
        f = os.path.join(directory, f)
        g_f = open(f, "r")
        g_f.readline()
        grains = np.asarray([])
        for line in g_f:
            grains = np.append(grains, int(line))
        ms_f = np.asarray([])
        for line in elem_f:
            ms_f = np.append(ms_f, grains[int(line)-1])
        ms_f = np.reshape(ms_f, [1].extend(shape), order='F')
        elem_f.close()
        g_f.close()
        ms_f = np.expand_dims(ms_f, axis=0)
        X_ = prim_basis.discretize(ms_f)
        print(X_.shape)
        ms_list[num] = X_
	ms_list = np.concatenate(ms_list, axis=0)
	print(ms_list.shape)
    return ms_list
    
def plotCorrelations(ms_list):
    for i, X_ in enumerate(ms_list):
        
        print(X_.shape)
        X_corr = correlate(X_)
        print X_corr[0].shape
        # average over z layers
        avg_corr = np.zeros(X_corr[0].shape[:-2] + (X_corr[0].shape[-1],))
        num_slices = X_corr[0].shape[-2]
        for slice in range(num_slices):
            avg_corr += X_corr[0][:,:,slice]/num_slices
        for corr_num in range(X_corr[0].shape[-1]):
            to_plot = avg_corr[:,:,corr_num]
            to_plot = np.expand_dims(to_plot, axis=2)
            print(to_plot.shape)
            draw_correlations(to_plot)
    

if __name__ == "__main__":
    import sys
    dir = os.getcwd()
    if(os.path.isdir(sys.argv[-1])):
        dir = sys.argv[-1]
    ms_list = readDirectory(dir)
    print(ms_list.shape)
    plotCorrelations(ms_list)
