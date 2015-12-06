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
    print("compare error shapes")
    print(corr1.shape)
    print(corr2.shape)
    ## normalize by the volume fraction
    center_x = (corr1.shape[0]-1)/2
    center_y = (corr1.shape[1]-1)/2
    print("center (%d,%d)" % (center_x,center_y))
    ## only need one unique correlation (use the particle autocorrelation)
    temp1 = corr1[...,1]/corr1[center_x,center_y,1]
    temp2 = corr2[...,1]/corr2[center_x,center_y,1]
    diff = temp1 - temp2
    #diff = diff[(center_x-35):(center_x+35),(center_y-35):(center_y+35)]
    print("max diff %f, max val 1 %f, min val 1 %f, max val 2 %f, min val 2 %f" % (np.max(np.abs(diff)), np.max(temp1), np.min(temp1), np.max(temp2), np.min(temp2)))
    plt.matshow(diff)
    plt.show()
    raw_se = (diff)**2
    mse = np.average(raw_se)
    raw_ae = np.abs(diff)
    mae = np.average(raw_ae)
    return (mse,mae,raw_se,raw_ae)
    
def makeMSodd(ms):
    new_size = []
    temp = ms.shape
    threshold = 799
    for i in range(1,len(temp)):
        if(temp[i]%2==0 and temp[i] < threshold):
            new_size.append(temp[i]-1)
        elif(temp[i] < threshold):
            new_size.append(temp[i])
        else:
            new_size.append(threshold)
    # new_size = [101,101,101]
    if(len(new_size) > 2):
        ms = ms[:,:new_size[0],:new_size[1],:new_size[2]]
    else:
        ms = ms[:,:new_size[0],:new_size[1]]
    return ms
    
def get2Points(dir_img, dir_recon):
    ## get our raw data (change the image names and scaling)
    o_sc = 5/2.0
    images = RMS.readImages(dir_img, ["L-T-James-Large-refined.png","L-T-James-refined.png","L-T-James-refined-3.png","L-T-James-refined-4-t.png","L-T-James-refined-5.png","L-T-James-refined-2.png"], [o_sc/1,o_sc/(164/33.75),o_sc/(248/33.75),o_sc/(248/33.75),o_sc/(248/33.75),o_sc/(252/33.75)])
    #images = RMS.readImages(dir_img, ["L-T-James-refined.png"], [o_sc/(248/33.75)])
    #images = images[::-1]
    image_corr = []
    prim_basis = PrimitiveBasis(2, [0, 1])
    for image in images:
        image = makeMSodd(image)
        image = prim_basis.discretize(image)
        temp = correlate(image)
        print(temp.shape)
        image_corr.append(temp)
    corr2d = RMS.trimMS(image_corr)
    large_ms = RMS.readDirectory(dir_recon)
    ## trim large MS to be odd in all spatial directions
    large_ms = makeMSodd(large_ms)
    large_ms = prim_basis.discretize(large_ms)
    ## do 2 point stats on large MS
    corr3d = correlate(large_ms)
    corr3d = corr3d.astype(np.float64)
    ## do 2 point stats on all the image MS
    #RMS.plotCorrelations(corr3d)
    #RMS.plotCorrelations(corr2d)
    center_x = corr2d.shape[1]/2
    center_y = corr2d.shape[2]/2
    new_size = corr3d.shape
    p_x = new_size[1]/2+1
    m_x = new_size[1]/2
    p_y = new_size[2]/2+1
    m_y = new_size[2]/2
    corr2d = corr2d[:,center_x-m_x:center_x+p_x,center_y-m_y:center_y+p_y]
    corr2d = corr2d.astype(np.float64)
    
    return (corr2d, corr3d)

if __name__ == "__main__":
    import sys
    dir_img = os.getcwd()
    dir_recon = os.getcwd()
    if(os.path.isdir(sys.argv[-1])):
        dir_img = sys.argv[-1]
    if(os.path.isdir(sys.argv[-2])):
        dir_recon = sys.argv[-2]
    print dir_img
    print dir_recon
    
    corr2d, corr3d = get2Points(dir_img, dir_recon)
    
    ## quantify the goodness of fit for the large MS in each plane
    new_size = corr2d.shape
    center = (new_size[-2]-1)/2
    print("center of the microstructure reconstruction z direction %d" % center)
    ms_slice = corr3d[...,center,:]
    validation_errors = []
    v_err_f = open("Recon_Errors.csv", "w")
    ms_err_f = open("MS_Errors.csv", "w")
    ## write out raw errors to files
    for i in range(len(corr2d)):
        print("Compare %d and reconstruction" % i)
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
                print("Compare %d and %d" % (i,j))
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
    