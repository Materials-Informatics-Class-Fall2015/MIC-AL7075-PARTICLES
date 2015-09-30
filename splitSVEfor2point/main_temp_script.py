"""
SCRIPT:CG_Block_script.py
AUTHOR:Bill Musinski
Edited by:Gustavo M. Castelluccio
"""

from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import testUtils
import random
import array

##########################
#### Define Variables ####
##########################

CP_xdim=0.028
CP_ydim=0.028
CP_zdim=0.028

##Mesh Size

CP_mesh_size=0.0007

aspect_ratio=1
aspect_ratio_1=1
featureSize=0.1
######################
#### Create Model ####
######################
MyModel = mdb.Model(name="IN_Block_CG")

## Delete Old Model
del mdb.models["Model-1"]

#########################
#### Create Geometry ####
#########################

## Create Overall Geometry ##

s = MyModel.ConstrainedSketch(name="Part1_Sketch", sheetSize=20.0)

g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints

s.rectangle(point1=(0,0),point2=(CP_xdim,CP_ydim))

Part1 = MyModel.Part(name="Part1", dimensionality=THREE_D,
    type=DEFORMABLE_BODY)

Part1.BaseSolidExtrude(sketch=s, depth=CP_zdim)
del mdb.models["IN_Block_CG"].sketches["Part1_Sketch"]

##END GEOMETRY CREATION##

##############################
#### Mesh the CP Geometry ####
##############################

#### Seed Part ####

e = Part1.edges

Part1.seedPart(size=CP_mesh_size, deviationFactor=0.1)

Part1Region=Part1.cells
Part1.setMeshControls(regions=Part1Region, elemShape=HEX, technique=SWEEP)

## Generate meshes

Part1.generateMesh(regions=Part1.cells)

#########################
#### Create Assembly ####
#########################

MyAssembly = MyModel.rootAssembly
AsmInstance = MyAssembly.Instance(name="AsmInst", part=Part1, dependent=ON)

####################
#### Create Job ####
####################

Job="main_temp"

mdb.Job(name=Job, model="IN_Block_CG", type=ANALYSIS, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, description="",
    parallelizationMethodExplicit=DOMAIN, multiprocessingMode=DEFAULT, numDomains=1, userSubroutine="", numCpus=1,
 scratch="", echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF)

mdb.jobs[Job].writeInput(consistencyChecking=OFF)
