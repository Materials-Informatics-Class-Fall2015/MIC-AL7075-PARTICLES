import os
import pymks
from pymks import PrimitiveBasis
from pymks.stats import correlate
from pymks.tools import draw_correlations
import numpy as np

def readDirectory(directory):
    filenames=[]
    for f in os.listdir(directory):
        if(f.find("E_el_min1")!=-1):
            filenames.append(f)
    strain_list = [0]*len(filenames)
    voigt = [0,4,8,1,2,5]
    for f in filenames:
        num = int(f[f.rfind("_")+1:f.rfind(".csv")])
        elem_f = os.path.join(directory, "trial_elem_grains_%d.txt" % num)
        elem_f = open(elem_f, "r")
        elem_f.readline()
        line = elem_f.readline()
        shape = map(int, line.split(",")[:3])
        new_shape = [6,1]
        new_shape.extend(shape)
        f = os.path.join(directory, f)
        strains = np.asarray([])
        strain_f = open(f,"r")
        for line in strain_f:
            temp = map(float,line.split(","))
            temp = np.asarray(temp)
            temp = temp[voigt]
            strains = np.append(strains,temp)
        strains = np.reshape(strains, new_shape, order='F')
        elem_f.close()
        strain_list[num] = strains
        
    strain_list = np.concatenate(strain_list, axis=1)
    print("final strains list size %s" % str(strain_list.shape))
    return strain_list


if __name__ == "__main__":
    import sys
    dir = os.getcwd()
    if(os.path.isdir(sys.argv[-1])):
        dir = sys.argv[-1]
    ms_list = readDirectory(dir)
    plotCorrelations(ms_list)