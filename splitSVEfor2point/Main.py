import CrackClasses as CC
import argparse
import logging
import os
import platform
import numpy as np
import math
import glob
import linecache
import cPickle

runout_life = 1E9

def writeElementSet(ith_grain_to_crack, numSet):
    global eliminated_elements, perimeter_elements
    string = '*Elset, elset=Elements_in_crack%d_%d\n' % (ith_grain_to_crack, numSet)
    items_per_line = 15
    for counter , elem in enumerate(eliminated_elements):
        string = string + str(elem+1) + ', '
        if (counter%items_per_line==0 and counter>1):
            string = string + '\n'
    string += '\n'

    string += '*Elset, elset=Perimeter_Elements%d_%d\n' % (ith_grain_to_crack, numSet)
    for counter , elem in enumerate(perimeter_elements):
        string = string + str(elem+1) + ', '
        if (counter%items_per_line==0 and counter>1):
            string = string + '\n'
    string += '\n'
    if(ith_grain_to_crack==0):
        mode = "w"
    else:
        mode = "a"
    f = open(os.path.join(__location__, 'crack_sets.txt'), mode)
    f.write(string)
    f.close()

def Nucleation():
    print("Nothing actually here yet")
    ## will incorporate the incubation (particle cracking) and nucleation (extension into matrix)
    global perimeter_elements, elements, grains, crack, element_fip, eliminated_elements, points
    perimeter_elements = {}
    ## interior elements (will be cracked)
    ids = [1098]
    # plane = np.asarray([2.0**.5,2.0**.5,0])
    plane = np.asarray([0,1,0])
        
    eliminated_elements = {}
    for i in ids:
        eliminated_elements[i] = True
    fid = os.path.join(__location__, "cracked_elem.txt")
    f = open(fid, "w")
    for e in eliminated_elements:
        f.write("%d   %f   %f   %f\n" % (e+1, plane[0], plane[1], plane[2]))
    f.close()
    
    ## get crack centroid by averaging perimeter centroids after nucleation
    crack = CC.Crack()
    crack.centroid = np.asarray([0,0,0])
    volume = 0
    for e in eliminated_elements:
        temp_el = elements[e]
        crack.centroid = crack.centroid + temp_el.volume*temp_el.centroid
        volume += temp_el.volume
    crack.centroid /= volume
    crack.plane = plane
    print("crack centroid: " + str(crack.centroid))
    
    ## set up crack perimeter points
    ## for this set up the radius as 1 times the element size
    points = []
    radius = elements[0].volume**(1/3.0)
    seedNum = 360
    normal_matrix = np.concatenate([np.reshape(plane,[1,3]), np.zeros([2,3])])
    rank, null_space = null(normal_matrix)
    v1 = null_space[:,0]
    v1 /= np.linalg.norm(v1)
    v2 = null_space[:,1]
    v2 /= np.linalg.norm(v2)
    v1 = np.reshape(v1, 3)
    v2 = np.reshape(v2, 3)
    for i in range(seedNum):
        p = CC.Point()
        p.id = i
        p.element = ids[0]
        th = i/(seedNum/2.0)*math.pi
        temp_dir = np.cos(th)*v1*radius + np.sin(th)*v2*radius
        p.direction = temp_dir / np.linalg.norm(temp_dir)
        p.point = crack.centroid + temp_dir
        points.append(p)
    perimeter_elements = {}
    for point in points:
        if(point.element not in perimeter_elements):
            temp = CC.PerimeterElement()
            temp.id = point.element
            temp.points.add(point.id)
            perimeter_elements[point.element] = temp
        perimeter_elements[point.element].points.add(point.id)
    ChangePointElements(ids[0], range(len(points)))
    
def MSC(ith_grain_to_crack):
    global perimeter_elements, elements, grains, element_fip, eliminated_elements, crack_set_counter, crack, points

    avg_prop = 0
    target_distance = elements[0].volume**(1/3.0)*1.1
    # target_distance = .0025 ## mesoscopic propagation distance
    ## create list of perimeter elements
    perimeter_elements = {}
    for point in points:
        if(point.element not in perimeter_elements):
            temp = CC.PerimeterElement()
            temp.id = point.element
            temp.points.add(point.id)
            perimeter_elements[point.element] = temp
        perimeter_elements[point.element].points.add(point.id)
    if(-1 in perimeter_elements):
        temp = perimeter_elements.pop(-1)
        logging.debug("We have lost %d points to propagating past neighbors (out of volume perhaps?)" % len(temp.points))
    ## calculate average propagation rate and use mesoscopic extension length to calculate number of cycles to extend
    ## choosing propagation cycles based on average propagation rate of perimeter is less sensitive to small driving force values
    all_props = []
    for key in perimeter_elements:
        UpdatePropagation(perimeter_elements[key])
        temp_prop = perimeter_elements[key].prop
        all_props.append(temp_prop)

    max_prop = max(all_props)
    avg_prop = np.median(all_props)
    if(avg_prop>0):
        prop_cycles = target_distance/avg_prop
    else:
        prop_cycles = 0
        logging.error("No crack propagation due to no driving force")
        return
    logging.info("max prop: %e" % (max_prop))
    logging.info("prop: %e prop_cycles: %e" %(avg_prop,prop_cycles))
    
    
    ## limit the maximum step so that we don't expand beyond one element neighbor at a time
    ## will only work in cubic meshes, need a way to expand this for non-uniform meshes
    side_len = elements[0].volume**(1/3.0)
    while prop_cycles > 0 :
        step_cycles = min(prop_cycles,side_len/max_prop)
        crack.cycles += step_cycles
        prop_cycles -= step_cycles
        ## propagate outward for each of the points
        oldPerim = set(perimeter_elements.keys())
        for p in points:
            if(p.element!=-1):
                temp_perim = perimeter_elements[p.element]
                plane = temp_perim.plane
                p.direction
                a = -np.dot(plane,p.direction)/np.dot(plane,crack.plane)
                proj_dir = p.direction+crack.plane*a
                proj_dir /= np.linalg.norm(proj_dir)
                ## NOTE: assume that the propgation rate is constant within the current element (is this valid for large element and small grain?)
                p.point += proj_dir*step_cycles*temp_perim.prop
                temp_dist = np.linalg.norm(p.point-crack.centroid)
                if(temp_dist>crack.length):
                    crack.length = temp_dist
        newPerim = set()
        iterkeys = perimeter_elements.keys()
        for base_id in iterkeys:
            newElems = ChangePointElements(base_id, perimeter_elements[base_id].points)
            newPerim = newPerim.union(newElems)
        newPerim = newPerim.difference(oldPerim)
        ## reassess the new perimeter elements for crack propagation
        for e in newPerim:
            UpdatePropagation(perimeter_elements[e])
            temp_prop = perimeter_elements[e].prop
            if(temp_prop > max_prop):
                max_prop = temp_prop
        AppendCrackHistory()
    ## remove old perimeter elements and add to crack and eliminated_elements
    fid = os.path.join(__location__, "cracked_elem.txt")
    f = open(fid, "a")
    iterkeys = perimeter_elements.keys()
    for key in iterkeys:
        if(len(perimeter_elements[key].points)<1):
            e = perimeter_elements.pop(key)
            f.write("%d   %f   %f   %f\n" % (e.id+1, e.plane[0], e.plane[1], e.plane[2]))
            eliminated_elements[key] = True
    f.close()
    
    writeElementSet(ith_grain_to_crack, crack_set_counter)
    crack_set_counter += 1

def null(a, rtol=1e-5):
    u, s, v = np.linalg.svd(a)
    rank = (s > rtol*s[0]).sum()
    return rank, v[rank:].T.copy()    

def AppendCrackHistory():
    """ Copy points to be saved in the crack history variable """
    global crack_history, points
    temp_list = []
    for p in points:
        temp_list.append(np.copy(p.point))
    crack_history.append(temp_list)
    
def ChangePointElements(base_id, temp_points):
    """ for each point in the base_element see which neighbor it has propagated into based on its updated position
        updates the point objects passed in as well as the perimeter_elements who now contain points
        Return new elements expanded to """
    global elements, perimeter_elements, points
    ## don't iterate 26 times over those points that are within this element since this is most likely base case
    extend_element = elements[base_id]
    temp_temp_points = []
    for i in temp_points:
        point = points[i]
        within_element = np.max(np.absolute(point.point-extend_element.centroid)) < (extend_element.volume**(1.0/3))/2
        if(not within_element):
            point.element = -1
            temp_temp_points.append(point.id)
    temp_points = temp_temp_points
    new_elements = set()
    base_element = perimeter_elements[base_id]
    for n in elements[base_id].neighbors:
        extend_element = elements[n]
        temp_temp_points = []
        for i in temp_points:
            point = points[i]
            within_element = np.max(np.absolute(point.point-extend_element.centroid)) < (extend_element.volume**(1.0/3))/2
            if(within_element):
                point.element = n
                if(n not in perimeter_elements):
                    temp_perim = CC.PerimeterElement()
                    temp_perim.id = n
                    perimeter_elements[n] = temp_perim
                    new_elements.add(n)
                ## update perimeter_elements knowledge of point positioning
                base_element.points.remove(point.id)
                perimeter_elements[n].points.add(point.id)
            else:
                temp_temp_points.append(point.id)
        temp_points = temp_temp_points
    ## ignore this for now because elements can be assigned outside of the volume if you crack too close to the edge
    # if(len(temp_points)>0):
        # raise ValueError('Some points were not assigned to neighbor/current element for element %d' % base_id)
    return new_elements
    
def IntermediatePlane(a1,b1,c1,d1,a2,b2,c2,d2,scale1,scale2):

    #Make sure everything is in Hessian Normal form:
    denominator1 = math.sqrt(math.pow(a1, 2) + math.pow(b1, 2) + math.pow(c1, 2))
    n1_x = a1 / denominator1
    n1_y = b1 / denominator1
    n1_z = c1 / denominator1
    p1   = d1 / denominator1

    denominator2 = math.sqrt(math.pow(a2, 2) + math.pow(b2, 2) + math.pow(c2, 2))
    n2_x = a2 / denominator2
    n2_y = b2 / denominator2
    n2_z = c2 / denominator2
    p2   = d2 / denominator2

    # Hessian normal form of plane is n dot x = -p
    # http://mathworld.wolfram.com/Plane.html
    # http://mathworld.wolfram.com/Plane-PlaneIntersection.html
    # solve for point (X1,X2,X3) which lies on line of intersection
    M = np.zeros([2,3])
    b = np.zeros([2,1])

    M[0][0] = n1_x
    M[1][0] = n2_x
    M[0][1] = n1_y
    M[1][1] = n2_y
    M[0][2] = n1_z 
    M[1][2] = n2_z

    b[0][0] = -p1
    b[1][0] = -p2

    X = np.linalg.lstsq(M,b)
    X1 = float(X[0][0])
    X2 = float(X[0][1])
    X3 = float(X[0][2])

    #Check that point lies on planes
    tolerance = 1E-6 #arbitrary tolerance
    if abs(a1*X1 + b1*X2 + c1*X3 + d1) > tolerance:
        print 'Solution 1 is shit'
        print a1*X1 + b1*X2 + c1*X3 + d1
    elif abs(a2*X1 + b2*X2 + c2*X3 + d2) > tolerance:
        print 'Solution 2 is shit'
        print a2*X1 + b2*X2 + c2*X3 + d2

    a3 = a1*scale1 + a2*scale2
    b3 = b1*scale1 + b2*scale2
    c3 = c1*scale1 + c2*scale2
    
    denominator3 = math.sqrt(math.pow(a3, 2) + math.pow(b3, 2) + math.pow(c3, 2))
    try:
        a3 = a3/denominator3
    except ZeroDivisionError:
         a3 = a3
    try:
        b3 = b3/denominator3
    except ZeroDivisionError:
         b3 = b3
    try:
        c3 = c3/denominator3
    except ZeroDivisionError:
         c3 = c3
    d3 = -a3*X1 -b3*X2 - c3*X3
    
    plane = np.asarray([a3,b3,c3])
    plane = np.reshape(plane, [3,1])

    return (plane, d3)

def StorePickles():
    global points, elements, grains, eliminated_elements, crack, crack_history
    fid = os.path.join(__location__, 'grains.p')
    fileObject = open(fid,'wb')
    cPickle.dump(grains, fileObject)
    fileObject.close()
    fid = os.path.join(__location__, 'elements.p')
    fileObject = open(fid,'wb')
    cPickle.dump(elements, fileObject)
    fileObject.close()
    fid = os.path.join(__location__, 'points.p')
    fileObject = open(fid,'wb')
    cPickle.dump(points, fileObject)
    fileObject.close()
    fid = os.path.join(__location__, 'eliminated_elements.p')
    fileObject = open(fid,'wb')
    cPickle.dump(eliminated_elements, fileObject)
    fileObject.close()
    fid = os.path.join(__location__, 'crack.p')
    fileObject = open(fid,'wb')
    cPickle.dump(crack, fileObject)
    fileObject.close()
    fid = os.path.join(__location__, 'crack_history.p')
    fileObject = open(fid,'wb')
    cPickle.dump(crack_history, fileObject)
    fileObject.close()
    
def ProjectPoint(norm, base, point):
    """ Get the point project onto a specific 3D planes, and the distance away from the plane """
    diff = point-base
    dist = np.dot(diff,norm)
    dist = float(dist)
    point = point - norm*(dist)
    dist = abs(dist)
    return (dist, point)
    
def UpdatePropagation(el, SII=True):
    """ Update the plane and life to crack this perimeter element """
    global perimeter_elements, elements, grains, element_fip, crack
    
    ## temp, uncalibrated fatigue parameters
    A_fs = 0.0789
    Phi_irr = 0.35
    DCTD_th = 2.9E-7
    
    temp = element_fip[el.id, :]
    sorted_ids = np.argsort(temp)
    SS_max_1 = sorted_ids[-1]
    plane_1 = SS_max_1/3
    plane_2 = plane_1
    index = -1
    while(plane_2==plane_1):
        index -= 1
        SS_max_2 = sorted_ids[index]
        plane_2 = SS_max_2/3
    max_FIP_1 = temp[sorted_ids[-1]]
    max_FIP_2 = temp[sorted_ids[index]]
    
    
    FIP = max_FIP_1

    # plane = 1 ## TODO temporary remove this XXX errors abound if left in
    g = elements[el.id].grain
    eulers = grains[g].orientation
    p1 = PlaneConstants(plane_1, eulers)
    if(SII):
        direction = elements[el.id].centroid - crack.centroid ## only used in tie breaker situations
        eps = max_FIP_1*(2.0**-8.0)
        p2 = PlaneConstants(plane_2, eulers)
        ## test which propagation normal has highest total propagation
        v1 = max_FIP_2*p2 + max_FIP_1*p1
        v2 = max_FIP_2*p2 - max_FIP_1*p1
        if((np.linalg.norm(v1)-np.linalg.norm(v2))>eps):
            FIP = np.linalg.norm(v1)
            p1 = v1/FIP
        elif((np.linalg.norm(v2)-np.linalg.norm(v1))>eps):
            FIP = np.linalg.norm(v2)
            p1 = v2/FIP
        else:
            logging.debug("Went to a tie breaker for SII plane in element %d" % el.id)
            ## project the direction onto both of the base planes
            (unused, dir_proj1) = ProjectPoint(p1, np.asarray([0,0,0]), direction)
            (unused, dir_proj2) = ProjectPoint(p2, np.asarray([0,0,0]), direction)
            avg_dir = (dir_proj1*max_FIP_1 + dir_proj2*max_FIP_2)/(max_FIP_1+max_FIP_2)
            ## compare the two possible intermediate planes to see which projects the propagation direction along both planes better
            (p1_1, d3) = IntermediatePlane(-float(p1[0]), -float(p1[1]), -float(p1[2]), 0, float(p2[0]), float(p2[1]), float(p2[2]), 0, max_FIP_1, max_FIP_2)
            (p1_2, d3) = IntermediatePlane(float(p1[0]), float(p1[1]), float(p1[2]), 0, float(p2[0]), float(p2[1]), float(p2[2]), 0, max_FIP_1, max_FIP_2)
            (unused, proj_int1) = ProjectPoint(p1_1, np.asarray([0,0,0]), avg_dir)
            (unused, proj_int2) = ProjectPoint(p1_2, np.asarray([0,0,0]), avg_dir)
            len1 = np.linalg.norm(proj_int1)
            len2 = np.linalg.norm(proj_int2)
            if(len1 > len2):
                p1 = p1_1
            else:
                p1 = p1_2
            ## project the FIPs onto the intermediate plane based on the previously projected propagation directions
            dir_proj1 /= np.linalg.norm(dir_proj1)
            dir_proj2 /= np.linalg.norm(dir_proj2)
            (unused, dir_proj1) = ProjectPoint(p1, np.asarray([0,0,0]), dir_proj1)
            (unused, dir_proj2) = ProjectPoint(p1, np.asarray([0,0,0]), dir_proj2)
            FIP = np.linalg.norm(dir_proj1)*max_FIP_1 + np.linalg.norm(dir_proj2)*max_FIP_2
    el.plane = p1
    prop = Phi_irr * (FIP*A_fs - DCTD_th)
    if(prop>0):
        el.prop = prop
    else:
        el.prop = 0    
        
def PlaneConstants(plane, eulers):
    #Define slip plane normals
    if plane == 0:
        n = np.matrix('1;1;1')
    elif plane == 1:
        n = np.matrix('-1;1;1') 
    elif plane == 2:
        n = np.matrix('-1;-1;1')   
    elif plane == 3:
        n = np.matrix('1;-1;1')         
        
    psi_ang = eulers[0]
    theta_ang = eulers[1]
    phi_ang = eulers[2]
    
    #Define direction cosines    
    s1 = math.sin(psi_ang)
    c1 = math.cos(psi_ang)
    s2 = math.sin(theta_ang)
    c2 = math.cos(theta_ang)
    s3 = math.sin(phi_ang)
    c3 = math.cos(phi_ang)

    #Generate rotation matrix
    R = np.zeros([3,3])

    R[0][0] = c1*c3-s1*s3*c2
    R[1][0] = s1*c3+c1*s3*c2
    R[2][0] = s3*s2
    R[0][1] = -c1*s3-s1*c3*c2
    R[1][1] = -s1*s3+c1*c3*c2 
    R[2][1] = c3*s2 
    R[0][2] = s1*s2 
    R[1][2] = -c1*s2 
    R[2][2] = c2     

    R =  np.matrix(R)
    
#    print R

    #Rotate vector normal to slip plane to global cords
    n_rot = np.dot(R,n)

    #Normalize vector to length 1
    n_rot = n_rot/np.linalg.norm(n_rot)
    n_rot = np.asarray(n_rot)
    n_rot = np.squeeze(n_rot)
    return n_rot

def InitVariables():
    global points, elements, grains, eliminated_elements, crack, crack_history
    
    ## if pickled files exist, read them in
    fid = os.path.join(__location__, 'grains.p')

    if os.path.exists(fid) == True:
        fileObject = open(fid,'r')
        grains = cPickle.load(fileObject)
        fileObject.close()
        fid = os.path.join(__location__, 'elements.p')
        fileObject = open(fid,'r')
        elements = cPickle.load(fileObject)
        fileObject.close()
        fid = os.path.join(__location__, 'points.p')
        fileObject = open(fid,'r')
        points = cPickle.load(fileObject)
        fileObject.close()
        fid = os.path.join(__location__, 'eliminated_elements.p')
        fileObject = open(fid,'r')
        eliminated_elements = cPickle.load(fileObject)
        fileObject.close()
        fid = os.path.join(__location__, 'crack.p')
        fileObject = open(fid,'r')
        crack = cPickle.load(fileObject)
        fileObject.close()
        fid = os.path.join(__location__, 'crack_history.p')
        fileObject = open(fid,'r')
        crack_history = cPickle.load(fileObject)
        fileObject.close()
    else:
        ## if not, read the data in from text files
        grains = []
        fid2 = os.path.join(__location__, 'Grains.txt')
        file = open(fid2,'r')
        for line in file:
            data = map(float,line.split(' , '))
            temp = CC.Grain()
            temp.orientation = np.asarray([data[1], data[2], data[3]])
            if(len(data)>4):
                temp.phase = data[4]
            else:
                temp.phase = 1
            grains.append(temp)
        file.close()

        ## init the element list
        elements = []
        fid2 = os.path.join(__location__, 'Neighbors_el.txt')
        fid = os.path.join(__location__, 'El_pos.txt')
        fid3 = os.path.join(__location__, 'Element_Volume.txt')
        f3 = open(fid3, 'r')
        f2 = open(fid2,'r')
        f1 = open(fid,'r')

        for line in f1:
            data = line.split(' ')
            temp = CC.Element()
            temp.centroid = np.asarray([float(data[1]), float(data[2]), float(data[3])])
            line = f2.readline()
            data = line.split(',')
            data = map(int, data)
            for i in range(len(data)):
                data[i] -= 1
            data = set(data)
            if -1 in data:
                data.remove(-1)
            temp.neighbors = data
            temp.volume = float(f3.readline().strip())
            elements.append(temp)
        f1.close()
        f2.close()
        f3.close()

        ## update grain->element and element->grain mappings
        for i in range(len(grains)):
            temp = ListElementsInGrain(i+1)
            for j in range(len(temp)):
                temp[j] -= 1
                elements[temp[j]].grain = i
            grains[i].elements = temp
            
        ## initialize the crack history
        crack_history = []
            
#--------------------------------------------------------------------------------	
# Function to return the elements within a grain
#--------------------------------------------------------------------------------
def ListElementsInGrain(grain):
    input_file = glob.glob(os.path.join(__location__, 'main_*[0-9]*.inp'))[0]         
    lookup = '*Elset, elset=Grain%s_set' % (grain)
    line_num = []
    with open(input_file) as fid:
        for num, line in enumerate(fid, 1):
            if lookup in line:
                line_num.append(num)
#				print 'grain found at line:', num
                                            
    inRange = True
    index = 1
    grain_list = []
    
    while inRange == True:
        #Test if line contains element numbers
        if '*' in linecache.getline(input_file, (line_num[0] + index)):
            inRange = False
            break
        #Add elements in line to list of elements
        list = (linecache.getline(input_file, (line_num[0] + index))).split(' ')
        #clean the list and change it to integers
        for item in list:
            grain_list.append((item.replace("\n", "")).replace(",", ""))
        while "" in grain_list: grain_list.remove("")
        grain_list = map(int, grain_list)	
        #Update index
        index = index + 1

    return [grain_list][0]

#--------------------------------------------------------------------------------	
#Initialize the FIP matrix to be indexed FIP[element][slip system])
#--------------------------------------------------------------------------------
def FIPelUpdate(ith_cracked_grain):
    global elements

    FIP = np.zeros((len(elements),12))

    #Process Data:
    if ith_cracked_grain == 1:
        filename = os.path.join(__location__, 'FIP_Nuc_el.txt')
    else:
        filename = os.path.join(__location__, ('FIP_MSC' + str(ith_cracked_grain) + '_el.txt'))

    #Open the FIP file
    rf = open(filename, 'r')

    #Parse data from FIP file into lists:
    for line in rf:
        str1 = line.split()
        element_num = int(str1[0]) 
        slip_sys_num = int(str1[1])
        grain_num = int(str1[2])
        FIP_val = float(str1[3])
        FIP[element_num-1,slip_sys_num-1] = FIP_val
            
    return FIP


parser = argparse.ArgumentParser(description='Description of your program')
parser.add_argument('-nH','--noHistory', help='run without any history tracking', required=False,action="store_true")
parser.add_argument('-SI','--StageI_only', help='calculate life using origional stage I growth algorithm', required=False,action="store_true")
parser.add_argument('-nABQ','--nonABQUScall', help='Main.py is being called from outside the Abaqus enviroment', required=False,action="store_true")
parser.add_argument('-eN','--exportNormals', help='call should also export crack plane normal values', required=False,action="store_true")
args = parser.parse_args()

crack_set_counter = 1
#Initialize the full path to the folder where the script is being run:
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

#logging.basicConfig(stream=sys.stderr, level=logging.INFO)
fid = os.path.join(__location__, 'PythonLog.txt')
if (os.path.exists(fid) and args.nonABQUScall): #remove old log file if its there
    os.remove(fid)
logging.basicConfig(filename=fid, level=logging.DEBUG)

print platform.python_version()
logging.info(platform.python_version())
logging.info("----------------Run Main Program----------------")

#Determine which call is currently in progress based on presence of FIP files:
grains_to_crack = 1 #Number of grains for which there is FIP data available to crack
num_cracked_grains = 0 #Number of currently cracked grains

while( os.path.exists(os.path.join(__location__, ('FIP_MSC' + str(grains_to_crack+1) + '_el.txt'))) ):
    grains_to_crack += 1

num_cracked_grains = grains_to_crack - 1

if args.nonABQUScall:
    num_cracked_grains = 0
    for f in os.listdir(__location__):
        if(f.find(".p")==len(f)-2):
            os.remove(os.path.join(__location__,f))
    MSC_eval_start = 2
else:
    MSC_eval_start = grains_to_crack
    ## skip the MSC evaluation if this a nucleation step
    if(grains_to_crack == 1):
        MSC_eval_start = 2

InitVariables()

if num_cracked_grains == 0:
    logging.info('----------------------------------EVAL NUC----------------------------------')
    element_fip = FIPelUpdate(1)
    Nucleation()
    fid = os.path.join(__location__, 'CrackGrowth_py.txt')
    f = open(fid,'w')
    f.write(str(crack.cycles) + '    '  + str(crack.length) + '\n')
    f.close()
    writeElementSet(0, 0)

for ith_grain_to_crack in range(MSC_eval_start,(grains_to_crack+1)):
    logging.info('----------------------------------EVAL MSC%d----------------------------------', ith_grain_to_crack)
    element_fip = FIPelUpdate(ith_grain_to_crack)
    MSC(ith_grain_to_crack)
    if(ith_grain_to_crack==grains_to_crack or args.nonABQUScall):
        fid = os.path.join(__location__, 'CrackGrowth_py.txt')
        f = open(fid,'a')
        f.write(str(crack.cycles) + '    '  + str(crack.length) + '\n')
        f.close()
    # raw_input("Enter to continue")
logging.info('----------------------------------LIFE SUMMARY----------------------------------')
logging.info('crack_len: '+str(crack.length))
logging.info('crack_cycles: '+str(crack.cycles))
logging.info('----------------------------------Main.py Complete----------------------------------')

StorePickles()
