import os
for f in os.listdir(os.getcwd()):
    if(f.find(".inp")!=-1 and f.find("main")==0):
        num = f[f.rfind("_")+1:].split(".")[0]
        break
types = ['.txt', '.csv', '.p']
for f in os.listdir(os.getcwd()):
    for ext in types:
        index = f.rfind(ext)
        if(index==len(f)-len(ext)):
            try:
                os.rename(f, f[:index] + "_" + str(num) + f[index:])
            except:
                print(f + " could not be renamed")