import os
from pymks import PrimitiveBasis
from pymks.stats import correlate
from pymks.tools import draw_correlations
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.io as sio
from scipy import misc
from scipy import interpolate
import ReadMSFunction as RMS
import cPickle

def getErrors(corr1, corr2):
    print(corr1.shape)
    print(corr2.shape)
    ## normalize by the volume fraction
    center_x = (corr1.shape[0]-1)/2
    center_y = (corr1.shape[1]-1)/2
    ## only need one unique correlation
    temp1 = corr1[...,0]/corr1[center_x,center_y,0]
    temp2 = corr2[...,0]/corr2[center_x,center_y,0]
    diff = temp1 - temp2
    raw_se = (diff)**2
    mse = np.average(raw_se)
    raw_ae = np.abs(diff)
    mae = np.average(raw_ae)
    return (mse,mae,raw_se,raw_ae)

if __name__ == "__main__":
    import sys
    dir = os.getcwd()
    dir2 = os.getcwd()
    if(os.path.isdir(sys.argv[-1])):
        dir = sys.argv[-1]
    if(os.path.isdir(sys.argv[-2])):
        dir2 = sys.argv[-2]
    print dir
    print dir2
    ## get our raw data (change the image names and scaling)
    o_sc = 5
    images = RMS.readImages(dir, ["L-T-James-Large-refined.png","L-T-James-refined.png","L-T-James-refined-3.png","L-T-James-refined-4-t.png","L-T-James-refined-5.png","L-T-James-refined-2.png"], [o_sc/1,o_sc/(164/33.75),o_sc/(248/33.75),o_sc/(248/33.75),o_sc/(248/33.75),o_sc/(252/33.75)])
    large_ms = RMS.readDirectory(dir2)
    ## trim large MS to be odd in all spatial directions
    temp = large_ms.shape
    new_size = []
    for i in range(1,len(temp)):
        if(temp[i]%2==0):
            new_size.append(temp[i]-1)
        else:
            new_size.append(temp[i])
    # new_size = [101,101,101]
    large_ms = large_ms[:,:new_size[0],:new_size[1],:new_size[2]]
    prim_basis = PrimitiveBasis(2, [0, 1])
    large_ms = prim_basis.discretize(large_ms)
    ## do 2 point stats on large MS
    corr3d = correlate(large_ms)
    corr3d = corr3d.astype(np.float64)
    ## do 2 point stats on all the image MS
    images = prim_basis.discretize(images)
    corr2d = correlate(images)
    # RMS.plotCorrelations(corr2d)
    center_x = corr2d.shape[1]/2
    center_y = corr2d.shape[2]/2
    p_x = new_size[0]/2+1
    m_x = new_size[0]/2
    p_y = new_size[1]/2+1
    m_y = new_size[1]/2
    corr2d = corr2d[:,center_x-m_x:center_x+p_x,center_y-m_y:center_y+p_y]
    corr2d = corr2d.astype(np.float64)
    # RMS.plotCorrelations(corr3d)
    
    ## quantify the goodness of fit for the large MS in each plane 
    center = (new_size[-1]-1)/2
    ms_slice = corr3d[...,center,:]
    validation_errors = []
    v_err_f = open("Recon_Errors.csv", "w")
    ms_err_f = open("MS_Errors.csv", "w")
    ## write out raw errors to files
    for i in range(len(corr2d)):
        temp = getErrors(corr2d[i], ms_slice[0])
        validation_errors.append(temp[:2])
        f = open("MS%d_recon_2pointdiff.csv" % i, "w")
        temp_1 = temp[2].flat
        temp_2 = temp[3].flat
        for k in range(temp[2].size):
            f.write("%e,%e\n" % (temp_1[k],temp_2[k]))
        f.close()
        v_err_f.write("%e,%e\n" % (temp[0], temp[1]))
    microstructure_errors = []
    
    for i in range(len(corr2d)):
        for j in range(len(corr2d)):
            if(i < j):
                temp = getErrors(corr2d[i], corr2d[j])
                microstructure_errors.append(temp[:2])
                f = open("MS%d_%d_2pointdiff.csv" % (i,j), "w")
                temp_1 = temp[2].flat
                temp_2 = temp[3].flat
                for k in range(temp[2].size):
                    f.write("%e,%e\n" % (temp_1[k],temp_2[k]))
                f.close()
                ms_err_f.write("%d,%d,%e,%e\n" % (i,j,temp[0],temp[1]))
    ## plot the mean errors
    names = ["MSE", "MAE"]
    labels = ["", "Experimental", "Validation", ""]
    print(microstructure_errors)
    print(validation_errors)
    temp = open("errors.p", "w")
    cPickle.dump([microstructure_errors, validation_errors], temp)
    temp.close()
    for i in range(len(names)):
        for e in microstructure_errors:
            plt.plot(1,e[i],'k*')
        for e in validation_errors:
            plt.plot(2,e[i],'k*')
        plt.xlim(0,3)
        plt.xticks(range(len(labels)), labels)
        plt.ylabel(names[i])
        plt.savefig(names[i])
        plt.clf()
    v_err_f.close()
    ms_err_f.close()
    