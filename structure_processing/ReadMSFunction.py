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


def readDirectory(directory):
    prim_basis = PrimitiveBasis(n_states=2, domain=[1, 2])
    filenames = []
    mat_files = []
    for f in os.listdir(directory):
        if(f.find("trial_Phases")!=-1 and f.find("_old")==-1):
            filenames.append(f)
        elif(f.find(".mat")==(len(f)-4)):
            mat_files.append(f)
    print(filenames)
    print(mat_files)
    ms_list = [0]*len(filenames)
    for f in filenames:
        num = int(f[f.rfind("_")+1:f.rfind(".txt")])
        elem_f = os.path.join(directory, "trial_elem_grains_%d.txt" % num)
        elem_f = open(elem_f, "r")
        elem_f.readline()
        line = elem_f.readline()
        shape = map(int, line.split(",")[:3])
        new_shape = [1]
        new_shape.extend(shape)
        f = os.path.join(directory, f)
        g_f = open(f, "r")
        g_f.readline()
        phases = np.asarray([])
        for line in g_f:
            phases = np.append(phases, int(line))
        ms_f = np.asarray([])
        for line in elem_f:
            ms_f = np.append(ms_f, phases[int(line)-1]-1)
        
        ms_f = np.reshape(ms_f, new_shape, order='F')
        elem_f.close()
        g_f.close()
        
        ms_list[num] = ms_f
    if(len(filenames)>0):
        ms_list = np.concatenate(ms_list, axis=0)
    temp_list = []
    for f in mat_files:
        f = os.path.join(directory,f)
        temp_list.append(readMatlab(f))
    if(len(temp_list)>0):
        if(len(ms_list)>0):
            ms_list = np.concatenate(temp_list.extend(ms_list), axis=0)
        else:
            ms_list = np.concatenate(temp_list, axis=0)
    print("final MS List shape: %s" % str(ms_list.shape))
    return ms_list
    
def scaleMS(ms, scale):
    x = np.arange(ms.shape[0])
    y = np.arange(ms.shape[1])
    z = ms
    xx = np.linspace(x.min(), x.max(), scale*ms.shape[0])
    yy = np.linspace(y.min(), y.max(), scale*ms.shape[1])
    kernel = interpolate.RectBivariateSpline(x,y,z,kx=2,ky=2)
    ms = kernel(xx,yy)
    threshold = 0.5
    filter = ms < threshold
    ms[filter] = 0
    filter = ms >= threshold
    ms[filter] = 1
    return ms
    
def readImages(directory, filenames, scales):
    """ Return a list of microstructures scaled and trimmed to the minimum value """
    temp_list = []
    threshold = 150
    for f, s in zip(filenames,scales):
        f = os.path.join(directory, f)
        ms = misc.imread(f)
        ms = np.average(ms, axis=-1)
        filter = ms < threshold
        ms[filter] = 1
        filter = ms >= threshold
        ms[filter] = 0
        ms = scaleMS(ms, s)
        #plt.matshow(ms)
        #plt.show()
        ms = np.expand_dims(ms, axis=0)
        print(ms.shape)
        temp_list.append(ms)
    return temp_list
    
def trimMS(temp_list):
    min_x = 10**10
    min_y = 10**10
    
    for ms in temp_list:
        if(ms.shape[1] < min_x):
            min_x = ms.shape[1]
        if(ms.shape[2] < min_y):
            min_y = ms.shape[2]
    ms_list = None
    for ms in temp_list:
        center_x = (ms.shape[1]-1)/2
        center_y = (ms.shape[2]-1)/2
        m_x = center_x - (min_x-1)/2
        m_y = center_y - (min_y-1)/2
        ms = ms[:,m_x:m_x+min_x, m_y:m_y+min_y]
        
        if(ms_list == None):
            ms_list = ms
        else:
            print(ms_list.shape)
            ms_list = np.concatenate((ms_list, ms), axis=0)
    return ms_list
    
def plotCorrelations(correlations):
    print(correlations.shape)
    for i, X_corr in enumerate(correlations):
        # average over z layers
        if(len(X_corr.shape)>=4):
            ## shape is x,y,z,correlation
            center = (X_corr.shape[2]-1)/2
            avg_corr = X_corr[:,:,center,:]
        else:
            avg_corr = X_corr
        # draw_correlations(avg_corr)
        for corr_num in range(X_corr.shape[-1]):
            to_plot = avg_corr[:,:,corr_num]
            # print(to_plot.dtype)
            plt.matshow(to_plot)
            plt.show()
            # to_plot = np.expand_dims(to_plot, axis=2)
            # print(to_plot.shape)
            # draw_correlations(to_plot)
    

def plotMS(ms_list):
    # only plot the first MS for now
    ms = ms_list[0]
    print(np.min(ms))
    print(np.max(ms))
    x = []
    y = []
    z = []
    for i in range(ms.shape[0]):
        for j in range(ms.shape[1]):
            for k in range(ms.shape[2]):
                if(ms[i,j,k]==1):
                    x.append(i)
                    y.append(j)
                    z.append(k)
        # print(i)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x,y,z,s=10)
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.show()
    
def computeStats(ms_list):
    ## particle volume fraction
    for ms in ms_list:
        count = np.sum(ms)
        vf = count/float(ms.size)
        print("Volume Fraction is: %f" % (vf))
        
def readMatlab(file):
    """ Read in a .mat file and return in the normal MS format """
    mat_contents = sio.loadmat(file)
    ## assume will always read the S_star variable with the reconstruction
    ms = mat_contents['S_star']
    ms = np.expand_dims(ms, axis=0)
    ms = 1-ms
    print(ms.shape)
    
    return ms
    
if __name__ == "__main__":
    import sys
    dir = os.getcwd()
    if(os.path.isdir(sys.argv[-1])):
        dir = sys.argv[-1]
    readImages(dir, ["L-S.png", "L-T.png", "L-T-James-refined-5.png"], [1.5, 1.5, 1])
    # ms_list = readDirectory(dir)
    # ms_list = readMatlab(sys.argv[-1])
    # print(ms_list.shape)
    # plotMS(ms_list)
    # computeStats(ms_list)
