import os 

types = ['.txt', '.csv', '.p']
for f in os.listdir(os.getcwd()):
    for ext in types:
        index = f.rfind(ext)
        index2 = f.rfind("_")
        isNum = False
        try:
            #print(f[index2+1:index])
            int(f[index2+1:index])
            isNum = True
        except:
            pass
        if(index==len(f)-len(ext) and index2>0 and isNum):
            try:
                os.rename(f, f[:index2] +  f[index:])
            except:
                print(f + " could not be renamed, shorter file name likely already exists")