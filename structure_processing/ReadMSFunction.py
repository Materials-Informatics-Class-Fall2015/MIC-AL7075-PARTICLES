import os
import pymks
from pymks import PrimitiveBasis
from pymks.stats import correlate
from pymks.tools import draw_correlations
import numpy as np

def readDirectory(directory):
    filenames=[]
    for f in os.listdir(directory):
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
        ms_f = np.reshape(ms_f, shape, order='F')
        elem_f.close()
        g_f.close()
        ms_list[num] = ms_f
    return ms_list
    
def plotCorrelations(ms_list):
    prim_basis = PrimitiveBasis(n_states=2, domain=[1, 2])
    for i, ms in enumerate(ms_list):       
        X_ = prim_basis.discretize(ms)
        print(X_.shape)
        X_corr = correlate(X_)
        print X_corr[0].shape
        draw_correlations(X_corr[0])
    

if __name__ == "__main__":
    import sys
    dir = os.getcwd()
    if(os.path.isdir(sys.argv[-1])):
        dir = sys.argv[-1]
    ms_list = readDirectory(dir)
    plotCorrelations(ms_list)
