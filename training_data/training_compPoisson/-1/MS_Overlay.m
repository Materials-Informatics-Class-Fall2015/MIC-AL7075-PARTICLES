% This function will return descriptors for use in later functions based on
% the overlay of a Dream3D generated structured mesh of grains over an
% unstructured mesh generated using the earlier part of the matlab mesh
% generator

function [grain,ElemGrainNo,RelLength, V_grn, new_to_old_grains]=MS_Overlay(el_centroid,n_El,V_el, grain_file, geom, options)
%% read in Dream3D exported structured mesh elements to grain number mapping

%read mesh dimensions in terms of number of elements
try
    dims = dlmread(grain_file, ',', [1 0 1 5]);
    x_elems = dims(1);
    y_elems = dims(2);
    z_elems = dims(3);
    old_x = double(dims(4));
    old_y = double(dims(5));
    old_z = double(dims(6));
    res_x = old_x / x_elems;
    res_y = old_y / y_elems;
    res_z = old_z / z_elems;
    offset_y = (old_y - geom.y)/2;
    offset_z = (old_z - geom.z)/2;
    offset_x = (old_x - geom.x)/2;
    if(isfield(options,'free_surface'))
        if(strcmp(options.free_surface,'x'))
            offset_x = 0;
        elseif(strcmp(options.free_surface,'y'))
            offset_y = 0;
        elseif(strcmp(options.free_surface,'z'))
            offset_z = 0;
        end
    end
catch ME
    dims = dlmread(grain_file, ',', [1 0 1 2]);
    x_elems = dims(1);
    y_elems = dims(2);
    z_elems = dims(3);
    res_x = geom.x / x_elems;
    res_y = geom.y / y_elems;
    res_z = geom.z / z_elems;
    offset_y = 0;
    offset_z = 0;
    offset_x = 0;
end

if(offset_y < 0 || offset_z < 0)
    error('Attempt to overlay smaller MS onto larger geometry')
end

layer = x_elems*y_elems;

fid = fopen(grain_file);
%strip header
fgets(fid);
%strip geometry line
fgets(fid);
%textscan remainder of file for the mapping of element# to grain number
structured_el_grain = cell2mat(textscan(fid, '%d'));
fclose(fid);

%% map unstructured elements to structured mesh grains based on centroid location
n_grains = max(structured_el_grain);
V_grn = zeros(n_grains,1);
grain_centroids = zeros(n_grains,3);

ElemGrainNo = zeros(n_El,1);

for elem = 1:n_El
    %find associated grain from voxelated mesh
    x_ind = floor((el_centroid(elem,1) + offset_x) / res_x);
    y_ind = floor((el_centroid(elem,2) + offset_y) / res_y);
    z_ind = floor((el_centroid(elem,3) + offset_z) / res_z);
    g_No = structured_el_grain( z_ind * layer + y_ind * x_elems + x_ind + 1);
    %assign grain to element
    
    % if this doesn't make nice looking grains, use neighboring elements to
    % vote on grain
    ElemGrainNo(elem) = g_No;
    %assign element information to grain
    v = V_el(elem);
    V_grn(g_No) = V_grn(g_No) +  v;
    grain_centroids(g_No,1) = grain_centroids(g_No,1) + el_centroid(elem,1)*v;
    grain_centroids(g_No,2) = grain_centroids(g_No,2) + el_centroid(elem,2)*v;
    grain_centroids(g_No,3) = grain_centroids(g_No,3) + el_centroid(elem,3)*v;
    
end

%% account for grains that have 0 volume by creating a mapping between old and new grains
% relate trimmed grain numbers from 'trial_elem_grains' file to final grain
% sets in ABAQUS
new_to_old_grains = [];
old_to_new_grains = zeros(n_grains);
new_num = 1;
for i = 1:n_grains
    if(V_grn(i) > 0)
        new_to_old_grains(new_num) = i;
        old_to_new_grains(i) = new_num;
        new_num = new_num+1;
    end
    
end
n_grains = new_num - 1;


% fix grain centroid and volume numbering
grain_centroids = grain_centroids(new_to_old_grains,:);
V_grn = V_grn(new_to_old_grains);



%% find centroids for grains
for g = 1:n_grains
    for j = 1:3
       grain_centroids(g,j) =  grain_centroids(g,j)/V_grn(g);
    end
end

%% calculate relative position of element to grain centroid

RelLength=zeros(n_El,3);

for e=1:n_El
        % also fix the element numbering to new grain numbers
        g = old_to_new_grains(ElemGrainNo(e,1));
        ElemGrainNo(e,1) = g;
        
        %Vector from grain centroid to each element centroid
        RelLength(e,1)=(el_centroid(e,1)-grain_centroids(g,1));
        RelLength(e,2)=(el_centroid(e,2)-grain_centroids(g,2));
        RelLength(e,3)=(el_centroid(e,3)-grain_centroids(g,3)); 
end

%% write out grain centers and display finish
fid = fopen(['Grain_Centers_' num2str(options.MS_number) '.txt'], 'w'); 
for jj=1:n_grains
    fprintf(fid,'%d %d %d\n', grain_centroids(jj,1), grain_centroids(jj,2), grain_centroids(jj,3)); 
end
fclose(fid);


disp(' ')
disp('Finished mapping grains')

grain = zeros(n_grains,1);

end
