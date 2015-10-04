%File:          Nodal_Values_preprocess.m
%Description:   This finds nodal xyz and connectivity values.
%
%Author:        Bill Musinski
%Date:          6/15/09

function [nxyz,nc,n_Nodes,n_El,el_centroid,V_el,Element_Neighbors, Element_Neighbors_common_nodes, Node_El]=Nodal_Values_preprocess(job,n_header_lines, options)
%% Get Node and Element Data from Input File

% Open input file
fid = fopen(job.temp_fid);

% Bypass the header lines in the input file:
% n_header_lines=number of lines in input file up til first node xyz
% description

for yy=1:n_header_lines
    fgetl(fid);
end

% Get node xyz coordinate values
% Node number, x, y, z
node_xyz=textscan(fid,'%f%f%f%f','Delimiter',',');

% Skip line in input file (=*Element, type=C3D8R (HEX))
fgetl(fid);

% Get nodal connectivity values (Only applies to tet elements. In order to
% apply to other element types, have to change number of '%f') for C3D4(TET)
% %f%f%f%f%f and for C3D8R(HEX) %f%f%f%f%f%f%f%f%f
node_con=textscan(fid,'%f%f%f%f%f%f%f%f%f','Delimiter',',');

fclose(fid);

% Input cells into arrays
% x, y, z indexed by node number
nxyz=[node_xyz{1,2}(:,1), node_xyz{1,3}(:,1), node_xyz{1,4}(:,1)];
% node numbers for each element (brick) indexed by element number
nc=[node_con{1,2}(:,1), node_con{1,3}(:,1), node_con{1,4}(:,1), node_con{1,5}(:,1), node_con{1,6}(:,1), node_con{1,7}(:,1), node_con{1,8}(:,1), node_con{1,9}(:,1)];

% Scale the model

nxyz=nxyz;

% Get number of total nodes/elements:
n_Nodes=size(nxyz,1);
n_El=size(nc,1);

%% Matlab Computations - Element Centroids and Volumes

% Element centroids for all elements
el_centroid=zeros(n_El,4);

for jj=1:n_El
    el_centroid(jj,1)=1/8*(nxyz(nc(jj,1),1)+nxyz(nc(jj,2),1)+nxyz(nc(jj,3),1)+nxyz(nc(jj,4),1)+nxyz(nc(jj,5),1)+nxyz(nc(jj,6),1)+nxyz(nc(jj,7),1)+nxyz(nc(jj,8),1));
    el_centroid(jj,2)=1/8*(nxyz(nc(jj,1),2)+nxyz(nc(jj,2),2)+nxyz(nc(jj,3),2)+nxyz(nc(jj,4),2)+nxyz(nc(jj,5),2)+nxyz(nc(jj,6),2)+nxyz(nc(jj,7),2)+nxyz(nc(jj,8),2));
    el_centroid(jj,3)=1/8*(nxyz(nc(jj,1),3)+nxyz(nc(jj,2),3)+nxyz(nc(jj,3),3)+nxyz(nc(jj,4),3)+nxyz(nc(jj,5),3)+nxyz(nc(jj,6),3)+nxyz(nc(jj,7),3)+nxyz(nc(jj,8),3));
    el_centroid(jj,4)=jj;
end

% Volume of tetrahedron:
% V=1/6*abs[(a-d).((b-d)x(c-d))]
% where '.' = dot product
%   and 'x' = cross product

% Volume of all elements
V_el=zeros(n_El,1);
%Vtot = 0;
%Volume calculation only valid for linearly swept meshes along z-direction
for jj=1:n_El
    a=nxyz(nc(jj,1),:);
    b=nxyz(nc(jj,2),:);
    c=nxyz(nc(jj,3),:);
    d=nxyz(nc(jj,4),:);
    e=nxyz(nc(jj,5),:);
    %f=nxyz(nc(jj,6),:);
    %g=nxyz(nc(jj,7),:);
    %h=nxyz(nc(jj,8),:);
    
    %V_el(jj,1)=abs(dot(a-b,cross(a-e,a-d))); % Calculates the volume of each element
    V_el(jj,1)=.5*norm(cross(d-b,c-a))*norm(a-e);
    %Vtot = Vtot + V_el(jj,1);
end

fid = fopen(['Element_Volume_' num2str(options.MS_number) '.txt'],'W');
for jj=1:n_El
    fprintf(fid,[sprintf('%d',V_el(jj)) '\n']);
end
fclose(fid);
%Vtot

%% Find Neighbors via generic connections of node numbers-----------------------------
%only valid for hexahedral elements

Node_El = zeros(n_Nodes,8); % elements attached to each node

Element_Neighbors = zeros(n_El,27); %index 27 is number of neighbors
Element_Neighbors_common_nodes = zeros(n_El,26);

%create reverse mapping of nodes to elements that it makes up
for i=1:n_El
    %for the 8 nodes this element is used in,
    for j=1:8
        %for all possible spaces in array
        for k=1:8
        %find empty space and place element value
        if Node_El(nc(i,j),k) == 0
            Node_El(nc(i,j),k) = i;
            break
        end
        end
    end
end

%get element connections and count number of nodes connecting them
%4 = face, 2 = edge, 1 = single node... 

for i = 1:n_Nodes;
    % for each element attached to node
    for j=1:8
        if Node_El(i,j) == 0 
            break
        end
        % for each element attached to node
        for k = 1:8
            if Node_El(i,k) == 0
                break
            end
            if j~=k
            curEl = Node_El(i,j);
            % for each item in neighbor list
            for ii = 1:26
            % if this element already in list, increment count    
            if Element_Neighbors(curEl,ii) ==  Node_El(i,k)
                Element_Neighbors_common_nodes(curEl,ii) = Element_Neighbors_common_nodes(curEl,ii) + 1;
                break
            % otherwise add to last open location in list and set count to 1
            elseif Element_Neighbors(curEl,ii) == 0
                Element_Neighbors(curEl,ii) = Node_El(i,k);
                Element_Neighbors_common_nodes(curEl,ii) = 1;
                Element_Neighbors(curEl,27) = Element_Neighbors(Node_El(i,j),27) + 1;               
                break
            end
			
            end
            end
        end
    end
end

%FINISH TEST STUFFFF-----------------------------------------
return
end
