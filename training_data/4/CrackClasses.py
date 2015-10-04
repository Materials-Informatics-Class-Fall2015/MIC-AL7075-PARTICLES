## All internal Python classes will use 0 indexing and only change to 1 indexing for export to Fortran
import numpy as np
import Main as m

class Element:
    def __init__(self):
        self.neighbors = []
        self.grain = 0
        self.centroid = np.asarray([0,0,0])
        self.volume = 0

class Point:
    def __init__(self):
        self.id = 0
        self.direction = np.asarray([0,0,0])
        self.point = np.asarray([0,0,0])
        self.element = 0

class PerimeterElement:
    def __init__(self):
        self.id = 0
        self.prop = 0
        self.plane = np.asarray([0,0,0])
        self.points = set()
    
class Grain:
    def __init__(self):
        self.phase = 0
        self.elements = []
        self.orientation = np.asarray([0,0,0])
    
class Crack:
    def __init__(self):
        self.centroid = np.asarray([0,0,0])
        self.length = 0
        self.cycles = 0
        self.normal = np.asarray([0,0,0])
        self.plane = np.asarray([0,0,0])
        self.lengths = {"Projected Area": [], "Maximum Length": [], "Minimum Length": [], "Median Length": []}
        self.cycle_history = []
        self.roughness = []
    def updateCrackProperties(self, points, cycles):
        """ Update crack history tracking based on current point locations and number of cycles at this iso-life contour """
        self.cycle_history.append(cycles)
        self.cycles = cycles
        projected = []
        lengths = []
        for p in points:
            (plane_dist, point) = m.ProjectPoint(self.plane, self.centroid, p.point)
            projected.append(point)
            dist = np.linalg.norm(p.point-self.centroid)
            lengths.append(dist)
        self.lengths["Maximum Length"].append(np.max(lengths))
        self.lengths["Minimum Length"].append(np.min(lengths))
        self.lengths["Median Length"].append(np.median(lengths))
        area = 0
        area += triangleArea(self.centroid, projected[0], projected[-1])
        for i in range(len(projected)-1):
            area += triangleArea(self.centroid, projected[i], projected[i+1])
        self.lengths["Projected Area"].append(area**.5)
        ## TODO update crack roughness as well and actually call this from MSC/NUC code
        
def triangleArea(point1,point2,point3):
    a = np.linalg.norm(point1-point2)
    b = np.linalg.norm(point2-point3)
    c = np.linalg.norm(point3-point1)
    return .25*((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c))**(.5)
        
        
class Edge:
    def __init__(self):
        self.v1 = 0
        self.v2 = 0
        self.cost = 0
    def getCost(self):
        return self.cost
    def __repr__(self):
        return "%d->%d" % (self.v1.id, self.v2.id)
        
class Tree:
    def __init__(self):
        self.v_set = set()
    def __repr__(self):
        return "Tree: " + str(self.v_set)
        
class Vertex:
    def __init__(self):
        self.tree = None
        self.edges = []
        self.id = 0
    def __repr__(self):
        st = "Vertex " + str(self.id) + " in tree: " + str(self.tree) + "\nwith edges: "
        for counter, e in enumerate(self.edges):
            st += str(e) + ", "
            if(counter%15 == 0):
                st += "\n"
            
        return  st