# An example of extracting data from an ABAQUS .odb file
# This code requires the use of abaqus python
# More notes to come.

from sys import *
from string import *
from math import *
from odbAccess import *
from abaqusConstants import *




""" Script file to be run in ABAQUS python environment.
    Expect input file argument first and output csv second
    Currently exports stress/strain data for first step, will later add multiple
    functions to things like FIPs etc given input arguments """
    
#we let it crash out if these aren't supplied
#arguments are all arguments to the abaqus command
#0: abaqus, 1: cae, 2: noGUI 3: script name ... 10: first actual argument

#make these inputs/dynamically linked to the actual changing needs of this in the future
minStep = 7
maxStep = 8

logFile = open("log.txt","w")

fileName = argv[-1]
logFile.write("in: " + fileName + "\n")
outputName = "plasticStrainAmplitudes.csv"
logFile.write("out: " + outputName + "\n")


# Extract the odb, assembly, and part
odb = openOdb(fileName)
assembly = odb.rootAssembly
instance = assembly.instances.keys()[0]


#open all output files before iteration
f = open(outputName, "w")
f.write("Ep11, Ep12, Ep13, Ep21, Ep22, Ep23, Ep31, Ep32, Ep33\n")

numElements = len(assembly.instances[instance].elements)
logFile.write("number of elements " + str(numElements) + "\n")

# Iterate through all steps
stepNum = 0
#perhaps move this to be a parameter eventually
stepsPerCycle = 3
frameCount = 0

#read element volumes from our list instead of the odb file
el_f = open("Element_Volume.txt")
evols = [0]*numElements
for i in range(numElements):
    l = el_f.readline()
    evols[i] = float(l)
    

for step in odb.steps.values():
    firstFrame = odb.steps[step.name].frames[0]

    #number of rows will be number of elements in model

    if (stepNum==minStep or stepNum==maxStep):
        
        avgEp = [0]*9
        for e in range(9):
            totalV = 0
            Ep = firstFrame.fieldOutputs["SDV" + str(85+e)]
            for i in range(numElements):
                v = evols[i]
                totalV += v
                avgEp[e] += Ep.values[i].data*1024*v
            avgEp[e] /= totalV*1024
        if(stepNum==minStep):
            minEp = avgEp
        else:
            maxEp = avgEp
    if(stepNum > maxStep):
        break
    stepNum += 1

for i in range(len(maxEp)):
    f.write(str((maxEp[i]-minEp[i])/2))
    if(i!=len(maxEp)-1):
        f.write(",")
f.write("\n")
    

f.close()
logFile.close()
