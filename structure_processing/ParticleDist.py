import ReadMSFunction as RM
import Queue
import sys
import os
import numpy as np

def getNeighbors(dims):
    """ Return a list of the neighboring elements sizex27 
        remember we use Fortran ordering to wrap/unwrap """
    num_el = dims[0]*dims[1]*dims[2]
    neighbors = np.zeros((num_el,26))
    neighbors -= 1
    for z in range(dims[2]):
        for y in range(dims[1]):
            for x in range(dims[0]):
                index = 0
                for dx in range(-1,2):
                    for dy in range(-1,2):
                        for dz in range(-1,2):
                            n_x = (x+dx+dims[0])%dims[0]
                            n_y = (y+dy+dims[1])%dims[1]
                            n_z = (z+dz+dims[2])%dims[2]
                            if(not (dx==0 and dy==0 and dz==0)):
                                neighbors[x+y*dims[0]+z*dims[0]*dims[1],index] = n_x+n_y*dims[0]+n_z*dims[0]*dims[1]
                                index += 1
    return neighbors
    
def getDistances(ms):
    s = ms.shape
    dist = np.zeros(s)
    dist -= 1
    q = Queue.PriorityQueue()
    for i in range(ms.shape[0]):
        for j in range(ms.shape[1]):
            for k in range(ms.shape[2]):
                if(ms[i,j,k] == 1):
                    dist[i,j,k] = 0
                    q.put((0, i, j, k, i, j, k))
    while (not q.empty()):
        temp = q.get()
        d = temp[0] + 1
        x = temp[1]
        y = temp[2]
        z = temp[3]
        for dx in range(-1,2):
            for dy in range(-1,2):
                for dz in range(-1,2):
                    n_x = x + dx
                    n_y = y + dy
                    n_z = z + dz
                    n_x = (x+dx+s[0])%s[0]
                    n_y = (y+dy+s[1])%s[1]
                    n_z = (z+dz+s[2])%s[2]
                    if(dist[n_x,n_y,n_z] == -1 and not (dx==0 and dy==0 and dz==0)):
                        diff_x = min(abs(n_x-temp[4]), s[0]-abs(n_x-temp[4]))
                        diff_y = min(abs(n_y-temp[5]), s[1]-abs(n_y-temp[5]))
                        diff_z = min(abs(n_z-temp[6]), s[2]-abs(n_z-temp[6]))
                        # print("element %d %d %d" % (n_x, n_y, n_z))
                        # print("distances %d %d %d" % (diff_x, diff_y, diff_z))
                        dist[n_x,n_y,n_z] = (diff_x**2 + diff_y**2 + diff_z**2)**(.5)
                        q.put((d, n_x, n_y, n_z, temp[4], temp[5], temp[6]))
    return dist
    
def write(dir):
    ms_list = RM.readDirectory(dir)
    os.chdir(dir)
    for i in range(ms_list.shape[0]):
        dist = getDistances(ms_list[i])
        # print("Check dist")
        # center = (dist.shape[0]-1)/2
        # print(dist[center,center,center])
        # print(dist[center,center,center+3])
        # print(dist[center,center+3,center+3])
        dist = np.reshape(dist, (dist.size,1), order='F')
        f = open("particle_dist_%d.csv" % i, 'w')
        for j in range(dist.size):
            f.write("%f\n" % dist[j,0])
        f.close()


if(__name__ == "__main__"):
    write(sys.argv[-1])