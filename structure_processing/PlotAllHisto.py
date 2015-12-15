import matplotlib.pyplot as plt
import cPickle
import os
colors = ["#543005", "#8c510a", "#bf812d", "#dfc27d", "#80cdc1", "#35978f", "#01665e", "#003c30"]

def plotDirs(base_directory):
    for temp in os.listdir(base_directory):
        temp = os.path.join(base_directory, temp)
        if(os.path.isdir(temp)):
            for i in range(10):
                plotNum(temp, i)

def plotNum(directory, num):
    
    i = 1
    f = os.path.join(directory, "histo_%d_%d.p" % (i, num))
    if(not os.path.exists(f)):
        return
        
    fig = plt.figure(facecolor="white", figsize=(16,12), dpi=100)
    ax = fig.add_subplot(111)
    handles = []
    labels = []
    while(os.path.exists(f)):
        f = open(f)
        temp = cPickle.load(f)
        f.close()
        centers = temp[0]
        hist = temp[1]
        color = colors[(i-1)%len(colors)]
        handle, = ax.plot(centers*10000, hist, color=color)
        handles.append(handle)
        labels.append("Element Distance (%d,%d]"  % (i-1,i))
        #print("Finished distance %d" % i)
        i += 1
        f = os.path.join(directory, "histo_%d_%d.p" % (i, num))
    plt.legend(handles,labels)
    plt.ylabel("Normalized Frequency")
    plt.xlabel("Fatemi-Socie FIP x 10$^{-4}$")
    plt.ylim([0,0.6])
    plt.xlim([0,1])
    os.chdir(directory)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)
    plt.savefig("Plot%d" % num)
    plt.clf()
    plt.close(fig)
    #plt.show()
    
if(__name__=="__main__"):
    import sys
    plotDirs(sys.argv[-1])