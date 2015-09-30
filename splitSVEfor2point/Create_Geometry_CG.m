% This is my first attempt at running Abaqus through Matlab.
% The reason I want to do this is so that I can create and modify voronoi 
% tesellations.

function []=Create_Geometry_CG(job,geom,mesh,options)
%% Write Python Script

%if we have already created a mesh for this folder, do not run again
if(exist( strcat(job.Job_Name_temp, '.inp'), 'file')==2)
   return 
end

fid = fopen([job.Job_Name_temp '_script.py'], 'wt');  % create the script file
% new option to create notched or through-hole meshes
%TODO add job type to the variable that gets passed into this function
if(isfield(mesh,'type'))
    type = mesh.type;
else
    type = 'smooth';
end

%create a rough voxelated hole using our own smarts because it was too much
%effort to integrate that other piece of code
if ~isempty(strfind(type, 'rough'))
    x = geom.x;
    y = geom.y;
    z = geom.z;
    res = mesh.CP_mesh_size;
    fid = fopen([job.Job_Name_temp '.inp'], 'wt');
    featureSize = str2double(mesh.featureSize);
    if(isnan(featureSize))
        featureSize = 0;
    end
    
    nodes = zeros(uint32((x/res+1)*(y/res+1)),3);
    n = 0;
    jumps = zeros(uint32(y/res+1),1);
    index = 0;
    if ~isempty(strfind(type,'hole'))
        %scale hole based on basic pattern by finding even or odd number of
        %elements
        featureSize = floor(featureSize/100*x/res+.5);
        if(mod(featureSize,2)==0)
            offset = 0;
        else
            offset = res/2;
        end
        featureSize = featureSize/2*1.01-.7;
        featureSize = (res*featureSize)^2;
        m_y = floor(y/res/2)*res;
        m_x = floor(x/2/res)*res;
        for j = 0:res:y
            index = index+1;
            jumps(index) = n+1;
            for i = 0:res:x
                if( ((j-m_y-offset)^2 + (i-m_x-offset)^2) > featureSize )
                   n=n+1;
                   nodes(n,:) = [i,j,0];
                end
            end
        end
    elseif ~isempty(strfind(type,'crack'))
        featureSize = featureSize/100*x-res/2;
        elemsRemove = ceil(str2double(options.feature_height_percent)/100*y/res);
        nodesRemove = elemsRemove-1;
        center = (floor(y/res/2)+.25)*res;
        height = nodesRemove*.5*res;
        for j = 0:res:y
            index = index+1;
            jumps(index) = n+1;
            for i = 0:res:x
                if(i > featureSize || abs(j-center)>height)
                   n=n+1;
                   nodes(n,:) = [i,j,0];
                end
            end
        end
    end
    %trim all empty positions, n = number of nodes on an xy plane
    nodes = nodes(1:n,:);
    %sweep all nodes back in the z direction 
    all_nodes = zeros(uint32(n*(z/res+1)),3);
    layer = 0;
    for i = 1:n:length(all_nodes)-1
       all_nodes(i:(i+n-1),:) = [nodes(:,1:2), ones(n,1).*layer];
       layer = layer + res;
    end
    %create the elements from each node
    elements = zeros(uint32(x*y/res^2),8);
    e_count = 0;
    if ~isempty(strfind(type,'hole'))
        for j = 1:(y/res) 
            bottom = jumps(j)+1;
            start = jumps(j+1);
            top = start+1;
            while(bottom < start)
                %if proper element, assign as such
                if(all_nodes(top,1)==all_nodes(bottom,1) && (all_nodes(bottom,1)-all_nodes(bottom-1,1))<res*1.1 && (all_nodes(top,1)-all_nodes(top-1,1))<res*1.1)
                    e_count = e_count+1;
                    elements(e_count,:) = [top, bottom, bottom+n, top+n, top-1, bottom-1, bottom+n-1, top+n-1];
                    bottom = bottom + 1;
                    top = top + 1;
                %if not, increment lower x value node and try again
                elseif(all_nodes(top,1) < all_nodes(bottom,1) || (all_nodes(top,1)- all_nodes(top-1,1)) > res*1.1)
                    top = top + 1;
                else
                    bottom = bottom + 1;
                end
            end
        end
    elseif ~isempty(strfind(type,'crack'))
        featureSize = featureSize + res;
        height = elemsRemove*.5*res;
        for j = 1:(y/res) 
            bottom = jumps(j)+1;
            start = jumps(j+1);
            top = start+1;
            while(bottom < start)
                %if proper element, assign as such
                if(all_nodes(top,1)==all_nodes(bottom,1) && (all_nodes(bottom,1)-all_nodes(bottom-1,1))<res*1.1 && (all_nodes(top,1)-all_nodes(top-1,1))<res*1.1)
                    if(all_nodes(top,1) > featureSize || abs((all_nodes(top,2)+all_nodes(bottom,2))/2-center)>height)
                        e_count = e_count+1;
                        elements(e_count,:) = [top, bottom, bottom+n, top+n, top-1, bottom-1, bottom+n-1, top+n-1];
                    end
                    bottom = bottom + 1;
                    top = top + 1;
                %if not, increment lower x value node and try again
                elseif(all_nodes(top,1) < all_nodes(bottom,1) || (all_nodes(top,1)- all_nodes(top-1,1)) > res*1.1)
                    top = top + 1;
                else
                    bottom = bottom + 1;
                end
            end
        end
    end
    elements = elements(1:e_count,:);
    all_elements = zeros(length(elements)*uint32(z/res),8);
    layer = 0;
    for i = 1:e_count:length(all_elements)-1
       all_elements(i:(i+e_count-1),:) = elements + n*layer;
       layer = layer + 1;
    end
    %space out the header of this inp like an abaqus generated one would so
    %later code works the same
    for i = 1:8
    fprintf(fid, '\n');
    end
    %write out the abaqus node and element definitions
    fprintf(fid, '*Node\n');
    for i = 1:length(all_nodes)
       n = all_nodes(i,:);
       fprintf(fid, '%d,\t\t%d,\t\t%d,\t\t%d\n',i,n(1),n(2),n(3));
    end
    fprintf(fid, '*Element, type=C3D8R\n');
    for i = 1:length(all_elements)
       n = all_elements(i,:);
       fprintf(fid, '%d,%d,%d,%d,%d,%d,%d,%d,%d\n',i,n(1),n(2),n(3),n(4),n(5),n(6),n(7),n(8));
    end
    fclose(fid);

%create geometry using ABAQUS python script
else
    % The following block is used for comments in the header of the python
    % script
    fprintf(fid, ['"""' '\n']);
    fprintf(fid, ['SCRIPT:'     'CG_Block_script.py' '\n']);
    fprintf(fid, ['AUTHOR:'     'Bill Musinski' '\n']);
    fprintf(fid, ['Edited by:'     'Gustavo M. Castelluccio' '\n']);
    fprintf(fid, ['"""' '\n' '\n']);

    fprintf(fid, ['from abaqus import *' '\n']);
    fprintf(fid, ['from abaqusConstants import *' '\n']);
    fprintf(fid, ['import __main__' '\n']);
    fprintf(fid, ['import section' '\n']);
    fprintf(fid, ['import regionToolset' '\n']);
    fprintf(fid, ['import displayGroupMdbToolset as dgm' '\n']);
    fprintf(fid, ['import part' '\n']);
    fprintf(fid, ['import material' '\n']);
    fprintf(fid, ['import assembly' '\n']);
    fprintf(fid, ['import step' '\n']);
    fprintf(fid, ['import interaction' '\n']);
    fprintf(fid, ['import load' '\n']);
    fprintf(fid, ['import mesh' '\n']);
    fprintf(fid, ['import job' '\n']);
    fprintf(fid, ['import sketch' '\n']);
    fprintf(fid, ['import visualization' '\n']);
    fprintf(fid, ['import xyPlot' '\n']);
    fprintf(fid, ['import displayGroupOdbToolset as dgo' '\n']);
    fprintf(fid, ['import connectorBehavior' '\n']);
    fprintf(fid, ['import testUtils' '\n']);
    fprintf(fid, ['import random' '\n']);
    fprintf(fid, ['import array' '\n' '\n']);


    fprintf(fid, ['##########################' '\n']); 
    fprintf(fid, ['#### Define Variables ####' '\n']);
    fprintf(fid, ['##########################' '\n' '\n']);

    fprintf(fid, ['CP_xdim=' num2str(geom.x) '\n']);
    fprintf(fid, ['CP_ydim=' num2str(geom.y) '\n']);
    fprintf(fid, ['CP_zdim=' num2str(geom.z) '\n' '\n']);
    
    % Add crack variables
    if(strcmp(type,'crack') || strcmp(type,'notch'))
        fprintf(fid, ['crack_thickness=' num2str(str2double(options.feature_height_percent)*geom.y/100) '\n']);
        fprintf(fid, ['crack_depth=' num2str(str2double(mesh.featureSize)*geom.x/100) ' - crack_thickness/2\n']);
    end

    fprintf(fid, ['##Mesh Size' '\n' '\n']);

    fprintf(fid, ['CP_mesh_size=' num2str(mesh.CP_mesh_size) '\n' '\n']);
    %TODO add mesh aspect ratio to options
    if(isfield(mesh,'ar'))
        if(strcmp(mesh.ar,''))
            mesh.ar = '1';
        end
        fprintf(fid, ['aspect_ratio=' mesh.ar '\n']);
    else
        fprintf(fid, ['aspect_ratio=1\n']);
    end
    if(isfield(options,'aspect_ratio_1'))
        if(strcmp(options.aspect_ratio_1,''))
            options.aspect_ratio_1 = '1';
        end
        fprintf(fid, ['aspect_ratio_1=' options.aspect_ratio_1 '\n']);
    else
        fprintf(fid, ['aspect_ratio_1=1\n']);
    end
    if(isfield(mesh,'featureSize') && ~strcmp(mesh.featureSize,'') )
        fprintf(fid, ['featureSize=' num2str(str2double(mesh.featureSize)/200) '\n']);
    else
        fprintf(fid, ['featureSize=' num2str(.1) '\n']);
    end


    fprintf(fid, ['######################' '\n']);
    fprintf(fid, ['#### Create Model ####' '\n']);
    fprintf(fid, ['######################' '\n']);

    fprintf(fid, ['MyModel = mdb.Model(name="IN_Block_CG")' '\n' '\n']);

    fprintf(fid, ['## Delete Old Model' '\n']);
    fprintf(fid, ['del mdb.models["Model-1"]' '\n' '\n']);


    fprintf(fid, ['#########################' '\n']);
    fprintf(fid, ['#### Create Geometry ####' '\n']);
    fprintf(fid, ['#########################' '\n' '\n']);

    fprintf(fid, ['## Create Overall Geometry ##' '\n' '\n']);

    fprintf(fid, ['s = MyModel.ConstrainedSketch(name="Part1_Sketch", sheetSize=20.0)' '\n' '\n']);
    fprintf(fid, ['g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints' '\n' '\n']);

    if(strcmp(type,'smooth') || strcmp(type,'hole'))
        fprintf(fid, ['s.rectangle(point1=(0,0),point2=(CP_xdim,CP_ydim))' '\n' '\n']);
    end
    if(strcmp(type,'hole'))
       fprintf(fid, ['s.CircleByCenterPerimeter(center=(CP_xdim/2, CP_xdim/2), point1=(CP_xdim/2, CP_xdim*(.5+featureSize)))' '\n' '\n']);
    elseif(strcmp(type,'notch'))
       fprintf(fid, ['topPoint = (0.0, CP_ydim*0.5+crack_thickness/2)' '\n' '\n']);
       fprintf(fid, ['bottomPoint = (0, CP_xdim*0.5-crack_thickness/2)' '\n' '\n']);
       fprintf(fid, ['s.Arc3Points(point1=topPoint, point2=bottomPoint, point3=(CP_xdim*featureSize, CP_xdim/2))' '\n' '\n']);
       fprintf(fid, ['s.Line(point1=bottomPoint, point2=(0.0, 0.0))' '\n' '\n']);
       fprintf(fid, ['s.Line(point1=(0.0, 0.0), point2=(CP_xdim, 0.0))' '\n' '\n']);
       fprintf(fid, ['s.Line(point1=(CP_xdim, 0.0), point2=(CP_xdim, CP_ydim))' '\n' '\n']);
       fprintf(fid, ['s.Line(point1=(CP_xdim, CP_ydim), point2=(0.0, CP_ydim))' '\n' '\n']);
       fprintf(fid, ['s.Line(point1=(0.0, CP_ydim), point2=topPoint)' '\n' '\n']);
    elseif(strcmp(type,'crack'))
         fprintf(fid, ['s.Line(point1=(0.0,0.0), point2=(0.0,(CP_ydim/2.)-(crack_thickness/2.)))' '\n' '\n']);
         fprintf(fid, ['s.VerticalConstraint(entity=g[2], addUndoState=False)' '\n' '\n']);
         fprintf(fid, ['s.Line(point1=(0.0,(CP_ydim/2.)-(crack_thickness/2.)), point2=(crack_depth,(CP_ydim/2.)-(crack_thickness/2.)))' '\n' '\n']);
         fprintf(fid, ['s.HorizontalConstraint(entity=g[3], addUndoState=False)' '\n' '\n']);
         fprintf(fid, ['s.PerpendicularConstraint(entity1=g[2], entity2=g[3], addUndoState=False)' '\n' '\n']);
         fprintf(fid, ['s.Line(point1=(crack_depth,(CP_ydim/2.)+(crack_thickness/2.)), point2=(0.0,(CP_ydim/2.)+(crack_thickness/2.)))' '\n' '\n']);
         fprintf(fid, ['s.HorizontalConstraint(entity=g[4], addUndoState=False)' '\n' '\n']);
         fprintf(fid, ['s.Line(point1=(0.0,(CP_ydim/2.)+(crack_thickness/2.)), point2=(0.0,CP_ydim))' '\n' '\n']);
         fprintf(fid, ['s.VerticalConstraint(entity=g[5], addUndoState=False)' '\n' '\n']);
         fprintf(fid, ['s.PerpendicularConstraint(entity1=g[4], entity2=g[5], addUndoState=False)' '\n' '\n']);
         fprintf(fid, ['s.Line(point1=(0.0,CP_ydim), point2=(CP_xdim,CP_ydim))' '\n' '\n']);
         fprintf(fid, ['s.HorizontalConstraint(entity=g[6], addUndoState=False)' '\n' '\n']);
         fprintf(fid, ['s.PerpendicularConstraint(entity1=g[5], entity2=g[6], addUndoState=False)' '\n' '\n']);
         fprintf(fid, ['s.Line(point1=(CP_xdim, CP_ydim), point2=(CP_xdim, 0.0))' '\n' '\n']);
         fprintf(fid, ['s.VerticalConstraint(entity=g[7], addUndoState=False)' '\n' '\n']);
         fprintf(fid, ['s.PerpendicularConstraint(entity1=g[6], entity2=g[7], addUndoState=False)' '\n' '\n']);
         fprintf(fid, ['s.Line(point1=(CP_xdim, 0.0), point2=(0.0, 0.0))' '\n' '\n']);
         fprintf(fid, ['s.HorizontalConstraint(entity=g[8], addUndoState=False)' '\n' '\n']);
         fprintf(fid, ['s.PerpendicularConstraint(entity1=g[7], entity2=g[8], addUndoState=False)' '\n' '\n']);
         fprintf(fid, ['s.Arc3Points(point1=(crack_depth,(CP_ydim/2.)-(crack_thickness/2.)), point2=(crack_depth,(CP_ydim/2.)+(crack_thickness/2.)), point3=(crack_depth + (crack_thickness/2.), CP_ydim/2.))' '\n' '\n']);
    end

    fprintf(fid, ['Part1 = MyModel.Part(name="Part1", dimensionality=THREE_D,' '\n']);
    fprintf(fid, ['    type=DEFORMABLE_BODY)' '\n' '\n']);

    fprintf(fid, ['Part1.BaseSolidExtrude(sketch=s, depth=CP_zdim)' '\n']);

    fprintf(fid, ['del mdb.models["IN_Block_CG"].sketches["Part1_Sketch"]' '\n' '\n']);

    fprintf(fid, ['##END GEOMETRY CREATION##' '\n' '\n']);

    fprintf(fid, ['##############################' '\n']);
    fprintf(fid, ['#### Mesh the CP Geometry ####' '\n']);
    fprintf(fid, ['##############################' '\n' '\n']);

    fprintf(fid, ['#### Seed Part ####' '\n' '\n']);

    fprintf(fid, ['e = Part1.edges' '\n' '\n']);

    % set edge seeds based on type of generated structure
    if(strcmp(type,'smooth'))
        fprintf(fid, ['Part1.seedPart(size=CP_mesh_size, deviationFactor=0.1)' '\n' '\n']);
    elseif(strcmp(type,'hole'))
        fprintf(fid, ['pickedEdges = e.getSequenceFromMask(mask=(''[#10 ]'', ), )' '\n' '\n']);
        fprintf(fid, ['Part1.seedEdgeBySize(edges=pickedEdges, size=CP_mesh_size, deviationFactor=0.1, constraint=FINER)' '\n' '\n']);
        fprintf(fid, ['#seed in z direction' '\n' '\n']);
        fprintf(fid, ['pickedEdges = e.getSequenceFromMask(mask=(''[#8 ]'', ), )' '\n' '\n']);
        fprintf(fid, ['Part1.seedEdgeBySize(edges=pickedEdges, size=CP_mesh_size/aspect_ratio_1, deviationFactor=0.1, constraint=FINER)' '\n' '\n']);
        fprintf(fid, ['#seed hole edge' '\n' '\n']);
        fprintf(fid, ['pickedEdges = e.getSequenceFromMask(mask=(''[#1000 ]'', ), )' '\n' '\n']);
        fprintf(fid, ['Part1.seedEdgeBySize(edges=pickedEdges, size=CP_mesh_size/aspect_ratio, deviationFactor=0.1, constraint=FINER)' '\n' '\n']);
    elseif(strcmp(type,'notch'))
        fprintf(fid, ['pickedEdges = e.getSequenceFromMask(mask=(''[#1 ]'', ), )' '\n' '\n']);
        fprintf(fid, ['Part1.seedEdgeBySize(edges=pickedEdges, size=CP_mesh_size, deviationFactor=0.1)' '\n' '\n']);
        fprintf(fid, ['#seed in z direction' '\n' '\n']);
        fprintf(fid, ['pickedEdges = e.getSequenceFromMask(mask=(''[#8 ]'', ), )' '\n' '\n']);
        fprintf(fid, ['Part1.seedEdgeBySize(edges=pickedEdges, size=CP_mesh_size/aspect_ratio_1, deviationFactor=0.1, constraint=FINER)' '\n' '\n']);
        fprintf(fid, ['#seed in notch edge' '\n' '\n']);
        fprintf(fid, 'if(abs(crack_thickness-CP_ydim)>CP_ydim/1000):\n');
        fprintf(fid, ['\tpickedEdges = e.getSequenceFromMask(mask=(''[#2000 ]'', ), )' '\n' '\n']);
        fprintf(fid, ['\tPart1.seedEdgeBySize(edges=pickedEdges, size=CP_mesh_size/aspect_ratio, deviationFactor=0.1, constraint=FINER)' '\n' '\n']);
    elseif(strcmp(type,'crack'))
         fprintf(fid, ['# Seed bottom side' '\n' '\n']);
         fprintf(fid, ['pickedEdges = e.getSequenceFromMask(mask=(''[#400000 ]'', ), )' '\n' '\n']);
         fprintf(fid, ['Part1.seedEdgeBySize(edges=pickedEdges, size=CP_mesh_size, deviationFactor=0.1, constraint=FINER)' '\n' '\n']);
         fprintf(fid, ['# Seed right side' '\n' '\n']);
         fprintf(fid, ['pickedEdges = e.getSequenceFromMask(mask=(''[#80000 ]'', ), )' '\n' '\n']);
         fprintf(fid, ['Part1.seedEdgeBySize(edges=pickedEdges, size=CP_mesh_size, deviationFactor=0.1, constraint=FINER)' '\n' '\n']);
         fprintf(fid, ['# Seed top side' '\n' '\n']);
         fprintf(fid, ['pickedEdges = e.getSequenceFromMask(mask=(''[#10000 ]'', ), )' '\n' '\n']);
         fprintf(fid, ['Part1.seedEdgeBySize(edges=pickedEdges, size=CP_mesh_size, deviationFactor=0.1, constraint=FINER)' '\n' '\n']);
         fprintf(fid, ['# Seed in z direction' '\n' '\n']);
         fprintf(fid, ['pickedEdges = e.getSequenceFromMask(mask=(''[#20000 ]'', ), )' '\n' '\n']);
         fprintf(fid, ['Part1.seedEdgeBySize(edges=pickedEdges, size=CP_mesh_size/aspect_ratio_1, deviationFactor=0.1, constraint=FINER)' '\n' '\n']);
         fprintf(fid, ['# Seed top of crack' '\n' '\n']);
         fprintf(fid, ['pickedEdges = e.getSequenceFromMask(mask=(''[#400 ]'', ), )' '\n' '\n']);
         fprintf(fid, ['Part1.seedEdgeByBias(biasMethod=SINGLE, constraint=FINER, end1Edges=pickedEdges, maxSize=CP_mesh_size, minSize=CP_mesh_size/aspect_ratio)' '\n' '\n']);
         fprintf(fid, ['# Seed bottom of crack' '\n' '\n']);
         fprintf(fid, ['pickedEdges = e.getSequenceFromMask(mask=(''[#10 ]'', ), )' '\n' '\n']);
         fprintf(fid, ['Part1.seedEdgeByBias(biasMethod=SINGLE, constraint=FINER, end2Edges=pickedEdges, maxSize=CP_mesh_size, minSize=CP_mesh_size/aspect_ratio)' '\n' '\n']);
         fprintf(fid, ['# Seed crack tip' '\n' '\n']);
         fprintf(fid, ['pickedEdges = e.getSequenceFromMask(mask=(''[#80 ]'', ), )' '\n' '\n']);
         fprintf(fid, ['Part1.seedEdgeBySize(edges=pickedEdges, size=CP_mesh_size/aspect_ratio, deviationFactor=0.1, constraint=FINER)' '\n' '\n']);
         
    end


    fprintf(fid, ['Part1Region=Part1.cells' '\n']);
    fprintf(fid, ['Part1.setMeshControls(regions=Part1Region, elemShape=HEX, technique=SWEEP)' '\n' '\n']);

    fprintf(fid, ['## Generate meshes' '\n' '\n']);

    fprintf(fid, ['Part1.generateMesh(regions=Part1.cells)' '\n' '\n']);

    fprintf(fid, ['#########################' '\n']);
    fprintf(fid, ['#### Create Assembly ####' '\n']);
    fprintf(fid, ['#########################' '\n' '\n']);

    fprintf(fid, ['MyAssembly = MyModel.rootAssembly' '\n']);
    fprintf(fid, ['AsmInstance = MyAssembly.Instance(name="AsmInst", part=Part1, dependent=ON)' '\n' '\n']);

    fprintf(fid, ['####################' '\n']);
    fprintf(fid, ['#### Create Job ####' '\n']);
    fprintf(fid, ['####################' '\n' '\n']);

    fprintf(fid, ['Job=' '"' job.Job_Name_temp '"' '\n' '\n']);

    fprintf(fid, ['mdb.Job(name=Job, model="IN_Block_CG", type=ANALYSIS, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, description="",' '\n']);
    fprintf(fid, ['    parallelizationMethodExplicit=DOMAIN, multiprocessingMode=DEFAULT, numDomains=1, userSubroutine="", numCpus=1,' '\n']);
    fprintf(fid, [' scratch="", echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF)' '\n' '\n']);

    fprintf(fid, ['mdb.jobs[Job].writeInput(consistencyChecking=OFF)' '\n']);

    fclose(fid);

    % Run python script to create temporary input file with mesh info only
    system(['abaqus cae noGUI=' job.Job_Name_temp '_script.py']);
end

return
end
