import matplotlib.pyplot as plt
import cPickle
import os
colors = ["#543005", "#8c510a", "#bf812d", "#dfc27d", "#80cdc1", "#35978f", "#01665e", "#003c30"]

def plotDirs(base_directory):
    for temp in os.listdir(base_directory):
        temp = os.path.join(base_directory, temp)
        if(os.path.isdir(temp)):
            for i in range(2):
                plotNum(temp, i)

def plotNum(directory, num):
    i = 1
    f = os.path.join(directory, "histo_%d_%d.p" % (i, num))
    handles = []
    labels = []
    while(os.path.exists(f)):
        f = open(f)
        temp = cPickle.load(f)
        f.close()
        centers = temp[0]
        hist = temp[1]
        color = colors[(i-1)%len(colors)]
        handle, = plt.plot(centers, hist, color=color)
        handles.append(handle)
        labels.append("Element Distance <= %d"  % i)
        #print("Finished distance %d" % i)
        i += 1
        f = os.path.join(directory, "histo_%d_%d.p" % (i, num))
    plt.legend(handles,labels)
    plt.ylabel("Normalized Frequency")
    plt.xlabel("Fatemi-Socie FIP")
    os.chdir(directory)
    plt.savefig("Plot%d" % num)
    plt.clf()
    #plt.show()
    
if(__name__=="__main__"):
    import sys
    plotDirs(sys.argv[-1])