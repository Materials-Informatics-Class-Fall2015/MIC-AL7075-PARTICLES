function [Max_num_layers,num_lines_Min_dist,NumLayer_max,NumLayer_min]=Finish_input_file_CG(job,geom,RelLength,n_Orient,nxyz,n_Nodes,n_El,Grain_Array,n_grains,ElemGrainNo,NodeGrain,...
El_pos,mesh,SBperiod,Nom_strain,Nom_strain_min,Strain_rate,band_jump,grain,scaling,n_grain_fail,t2,...
Mantle_slip_sys,Core_slip_sys,Element_Neighbors, Element_Neighbors_common_nodes, options, V_el, V_grn, new_to_old_grains, Node_El)

%El_pos(:,2:4);
num_lines_Min_dist = 0;
NumLayer_max = 0;
NumLayer_min = 0;
DUMMY = zeros(27*n_El,2);

% normalized distance values
% face neighbor elements, edge neighbor elements, single node neighbors
dist_vals = [1 1.414213 1.732508];

%Preallocated dist array:
%dist = zeros(n_El, n_El);
dist = sparse(n_El, n_El,1);
disp(' ')
disp('Finish input file started')

Neighbor_g=int8(zeros(n_grains,n_grains));

%Neighbors_el=int32(zeros(n_El,26));

N_Elem_Bound=int16(zeros(n_grains,n_grains));

El_GB=int16(zeros(n_El,6));
El_GB_un=int16(zeros(n_El,6));

GB_dirX=int8(zeros(n_El,2));
GB_dirY=int8(zeros(n_El,2));
GB_dirZ=int8(zeros(n_El,2));


%What is going on here?
%Appears to init all possible x,y,z difference vectors for neighboring
%elements in a structured, voxelated mesh
% pp=1;
% for i=[0 1 -1 ]
% for j=[0 1 -1 ]
% for k=[0 1 -1]
% Neigh_vec(pp,1:3)=[i j k];
% pp=pp+1;
% if (i==0 && j==0 && k==0)
% pp=1;
% end
% end
% end
% end



if(strcmp(options.banded, 'TRUE'))

    [mmm,nnn] = size(DUMMY);
    DUMMY;
    counter=0;

    for mm=1:n_El
    temp1= unique(El_GB(mm,:));
    for nn=1:length(temp1)
    El_GB_un(mm,nn)=temp1(nn);
    end
    end
end



%-----------------------------------------------------------

% only write out the necessary files
if(strcmp(options.banded, 'TRUE'))
    
    fid = fopen(['Neighbors_el_' num2str(options.MS_number) '.txt'], 'W');  % Read temp input file
    for jj=1:n_El
    for pp=1:25
    %fprintf(fid,[sprintf('%d',Neighbors_el(jj,pp)) ', ']);
    fprintf(fid,[sprintf('%d',Element_Neighbors(jj,pp)) ', ']);

    end
    %fprintf(fid,[sprintf('%d',Neighbors_el(jj,26)) '\n']);
    fprintf(fid,[sprintf('%d',Element_Neighbors(jj,26)) '\n']);
    end
    fclose(fid);

    disp(' ')
    disp('    write Boundary_el.txt')
    fid = fopen(['Boundary_el_' num2str(options.MS_number) '.txt'], 'W');  % Read temp input file

    for jj=1:n_El

    if El_GB_un(jj,1)==0
       fprintf(fid,[ sprintf('%d',El_GB_un(jj,2)) ' , ' sprintf('%d',El_GB_un(jj,3))  ' , ' sprintf('%d',El_GB_un(jj,4)) ' , ' sprintf('%d',El_GB_un(jj,5)) ' , ' sprintf('%d',El_GB_un(jj,6)) ' , ' sprintf('%d',El_GB_un(jj,1))  '\n']);
    else
       fprintf(fid,[ sprintf('%d',El_GB_un(jj,1)) ' , ' sprintf('%d',El_GB_un(jj,2))  ' , ' sprintf('%d',El_GB_un(jj,3)) ' , ' sprintf('%d',El_GB_un(jj,4)) ' , ' sprintf('%d',El_GB_un(jj,5)) ' , ' sprintf('%d',El_GB_un(jj,6))  '\n']);
    end
    end
    fclose(fid);


    disp(' ')
    disp('    write Boundary_dir.txt')

    fid = fopen(['Boundary_dir_' num2str(options.MS_number) '.txt'], 'W');  % Read temp input file

    for jj=1:n_El

    %fprintf(fid,[ sprintf('%d',GB_dirX(jj,1)) ' , ' sprintf('%d',GB_dirX(jj,2))  ' , ' sprintf('%d',GB_dirY(jj,1)) ' , ' sprintf('%d',GB_dirY(jj,2)) ' , ' sprintf('%d',GB_dirZ(jj,1)) ' , ' sprintf('%d',GB_dirZ(jj,2))  '\n']);

    fprintf(fid,[ sprintf('%d',GB_dirX(jj,1)) ' , ' sprintf('%d',GB_dirX(jj,2)) ' , ' sprintf('%d',GB_dirY(jj,1))  ' , ' sprintf('%d',GB_dirY(jj,2)) ' , ' sprintf('%d',GB_dirZ(jj,1)) ' , ' sprintf('%d',GB_dirZ(jj,2))  '\n']);


    end
    fclose(fid);
end



fid = fopen([job.Job_Name_temp '.inp'], 'r');  % Read temp input file

% Get important lines from temp input file
A=cell(n_El+n_Nodes+2,1);

%skip the header and start of *Part
for ii = 1:8
    fgetl(fid);
end

for ii=1:(n_El+n_Nodes+2);
    A{ii,1}=fgetl(fid);
end

fclose(fid);

fid = fopen(job.inp_fid,'W');  %Open input file to write to
% copy lines to new file
for ii=1:length(A);
    fprintf(fid,[A{ii,1} '\n']);
end

% write out reference node to use in applying load
% fprintf(fid,['*Node\n999999,0,0,0\n']);

% Element set for CP region
fprintf(fid,['*Elset, elset=CP_set, Generate' '\n']);
fprintf(fid,['1, ' sprintf('%d',n_El) ', 1' '\n']);


% Create CP materials

%% Chose and Assing orientation to elements
% needed here because Matlab code handles all the slicing of grains into
% layers and stuff

%check presence of orientation file
orientation_file = ['trial_EulerAngles_' num2str(options.MS_number) '.txt'];
files = dir(orientation_file);
of_exists = ~isempty(files);

elem_orient = zeros(n_El,3);

%store randomly assigned Euler angles for each element
if(~of_exists)
    g_orient = zeros(n_grains, 4);
    
    for jj=1:n_grains
        phi1= rand*4*pi;
        PHI= rand*2*pi;
        phi2= rand*8*pi;
        
        g_orient(jj,:) = [double(jj), phi1, PHI, phi2];
        
        
        for ii=1:numel(Grain_Array{jj,1})
                elem_orient(Grain_Array{jj,1}(ii,1),1)=phi1;
                elem_orient(Grain_Array{jj,1}(ii,1),2)=PHI;
                elem_orient(Grain_Array{jj,1}(ii,1),3)=phi2;
    %            elem_grain(Grain_Array{jj,1}(ii,1))=jj;
        end
    end
%store Dream3D assigned Bunge-Euler angles for each element
else
    or_f = fopen(orientation_file);
    
    'found orientation file and using to assign grain orientation'
    
    %parse number of grains line
    fgets(or_f);
    %parse time line?
    fgets(or_f);
    %parse header line
    fgets(or_f);
    %textscan remainder of file for the mapping of element# to grain number
    g_orient = cell2mat(textscan(or_f, '%f  %f  %f  %f'));
    g_orient = g_orient(new_to_old_grains, :);
    g_orient = [g_orient(:,1) , g_orient(:,2)*pi/180 , g_orient(:,3)*pi/180 , g_orient(:,4)*pi/180];
    fclose(or_f);
    for jj=1:n_grains
        for ii=1:numel(Grain_Array{jj,1})
                %no conversion necessary, both Bunge-Euler
                e = Grain_Array{jj,1}(ii,1);
                elem_orient(e,1)=g_orient(jj,2);
                elem_orient(e,2)=g_orient(jj,3);
                elem_orient(e,3)=g_orient(jj,4);
        end
    end
end

%Create CP materials
  y(1:12,1:3)=1;
  y(4:9,1)=-1;
  y(7:12,2)=-1;
  

RLength = zeros(n_El, 4);
NumLayer = zeros(n_El, 4);
for ii=1:n_El
    %CALCULATE SLIP DIRECTIONS
    s1 = sin(elem_orient(ii,1));
    c1 = cos(elem_orient(ii,1));
    s2 = sin(elem_orient(ii,2));
    c2 = cos(elem_orient(ii,2));
    s3 = sin(elem_orient(ii,3));
    c3 = cos(elem_orient(ii,3));
    %BUNGE EULER NOTATION  
    dir_cos(1,1) = c1*c3-s1*s3*c2;
    dir_cos(1,2) = s1*c3+c1*s3*c2;
    dir_cos(1,3) =  s3*s2;
    dir_cos(2,1) = -c1*s3-s1*c3*c2;
    dir_cos(2,2) = -s1*s3+c1*c3*c2;
    dir_cos(2,3) = c3*s2;
    dir_cos(3,1) = s1*s2;
    dir_cos(3,2) = -c1*s2;
    dir_cos(3,3) = c2;

    %Project slip normal direction onto directionality of grain

    xm0=zeros(3,1);

    for n = 1:4
     xm0(:) = dir_cos(1,:)*y(n*3,1)+dir_cos(2,:)*y(n*3,2)+dir_cos(3,:)*y(n*3,3);
     xm0(:) = xm0(:)/norm(xm0);
    %Distance from the Element Center to the Grain center 
        RLength(ii,n)=  xm0(1) * RelLength(ii,1) + xm0(2) * RelLength(ii,2) + xm0(3) * RelLength(ii,3) ;
        NumLayer(ii,n)=int8(ceil(RLength(ii,n)/SBperiod));
    end
end


NumLayer_max=max(max(NumLayer));

NumLayer_min=min(min(NumLayer));

%Max_num_layers=abs(NumLayer_max)+abs(NumLayer_min)+1
Max_num_layers=NumLayer_max-NumLayer_min+1;

NumLayer(:,:)=NumLayer(:,:)-NumLayer_min+1;

%redefine the grain number so that they are in order

% num_g=unique(ElemGrainNo(:,1));
% for jj=1:size(num_g)
% for ll=1:n_El
% 	if ElemGrainNo(ll,1)==num_g(jj)
% 	ElemGrainNo(ll,1)=jj;
% 	end
% end
% end

if(strcmp(options.banded, 'TRUE'))
    Grain_layers(1:n_grains,1:Max_num_layers,1:4)=int16(0);


    for ii=1:n_El
    for jj=1:4
    Grain_layers(ElemGrainNo(ii,1),NumLayer(ii,jj),jj)=Grain_layers(ElemGrainNo(ii,1),NumLayer(ii,jj),jj) +1;
    end
    end
    

    max_elem_in_grain = int32(5*(n_El/n_grains)); %some safety factor for the largest grains
    ElemGrainNo;
    Elem_in_Grain = zeros(n_grains,max_elem_in_grain);

    for i=1:n_grains
    count = 2;
    for j=1:n_El
        if ElemGrainNo(j,1) == i;
            Elem_in_Grain(i,1) = Elem_in_Grain(i,1) + 1;
            Elem_in_Grain(i,count) = j;
            count = count + 1;
        end	
    end
    end

    Elem_in_Grain;

    % Create Layer_Grain element sets
    for ii=1:n_grains
    for jj=1:Max_num_layers
    for kk=1:4

    if Grain_layers(ii,jj,kk)>0

     fprintf(fid,['*Elset, elset=Bands_Grain' sprintf('%d',ii) '_layer' sprintf('%d',jj)  '_plane' sprintf('%d',kk)  '\n']);

       a=1;

        for KK = 1:Elem_in_Grain(ii,1)
            ll = Elem_in_Grain(ii,KK+1);
    %    for ll=1:n_El
     %        if ii==ElemGrainNo(ll,1)
             if (NumLayer(ll,kk))==jj

                if a==Grain_layers(ElemGrainNo(ll,1),NumLayer(ll,kk),kk)
    %			   fprintf(fid,sprintf('%d',ll));
                   fprintf(fid,num2str(ll));			   
                else
    %				fprintf(fid,[sprintf('%d',ll) ', ']);
                    fprintf(fid,[num2str(ll) ', ']);
                end

                a=a+1; 

                if a==17; %at the end of line
                   fprintf(fid,'\n');
                   a=1;
                end


             end



        end
        fprintf(fid,'\n');	
    end

    end
    end
    end
end



% Create Grain element sets
for jj=1:n_grains
    fprintf(fid,['*Elset, elset=Grain' sprintf('%d',jj) '_set' '\n']);
    a=1;
    for ii=1:numel(Grain_Array{jj,1})
        if ii==numel(Grain_Array{jj,1})
            fprintf(fid,sprintf('%d',Grain_Array{jj,1}(ii,1)));
        else
            fprintf(fid,[sprintf('%d',Grain_Array{jj,1}(ii,1)) ', ']);
            a=a+1;
            if a==17;
                fprintf(fid,'\n');
                a=1;
            end
        end
    end
    fprintf(fid,'\n');
end



% Create element sets that are neighboring 
if(strcmp(options.banded, 'TRUE'))
    for jj=1:n_grains
    for ii=1:n_grains

    if N_Elem_Bound(jj,ii)~=0

        fprintf(fid,['*Elset, elset=Boundary_Grain' sprintf('%d',(jj)) '_&_' sprintf('%d',(ii)) '\n']);
        a=1;
    for kk=1:N_Elem_Bound(jj,ii)

              if kk==N_Elem_Bound(jj,ii)
                fprintf(fid,[sprintf('%d',(GB(jj,ii,kk)))]);
              else
                fprintf(fid,[sprintf('%d',(GB(jj,ii,kk))) ', ']);
              end
                a=a+1;
                if a==17;
                    fprintf(fid,'\n');
                    a=1;
                end
    end
        fprintf(fid,'\n');
    end
    end
    end
end

% Create element sets that are neighboring 
if(strcmp(options.banded, 'TRUE'))
    % Create Element element sets
    for jj=1:n_El
        fprintf(fid,['*Elset, elset=Elem' sprintf('%d',(jj)) '_set' '\n']);
        fprintf(fid,[ sprintf('%d',(jj)) '\n']);
    end
end

%fprintf(fid,['*End Part' '\n' '**' '\n' '**' '\n' '** ASSEMBLY' '\n' '**' '\n']);
%fprintf(fid,['*Assembly, name=Assembly' '\n' '**' '\n']);
%fprintf(fid,['*Instance, name=AsmInst, part=Part1' '\n']);
%fprintf(fid,['*End Instance' '\n' '**' '\n']);


%Bottom face nodes
fprintf(fid,['*Nset, nset=Bottom_Face_Geometry' '\n']);

if(isfield(mesh,'type'))
    type = mesh.type;
else
    type = 'smooth';
end

%if the mesh was generated by ABAQUS, sort nodes in increasing x,y,z order
%so that the node sets match up for use in periodic boundary conditions
if isempty(strfind(type, 'rough'))
    [nxyz,I] = sortrows(nxyz,[1,2,3]);
else
   I = 1:length(nxyz);
end

Bot_nodes=I(nxyz(:,2)==0);
a=1;
for ii=1:numel(Bot_nodes)
    if ii==numel(Bot_nodes)
        fprintf(fid,num2str(Bot_nodes(ii)));
    else
        fprintf(fid,[num2str(Bot_nodes(ii)) ', ']);
        a=a+1;
        if a==17;
            fprintf(fid,'\n');
            a=1;
        end
    end
end

fprintf(fid,'\n');

%Right face nodes
fprintf(fid,['*Nset, nset=Left_Face_Geometry' '\n']);

Left_nodes=I(nxyz(:,1)==0);
a=1;
for ii=1:numel(Left_nodes)
    if ii==numel(Left_nodes)
        fprintf(fid,num2str(Left_nodes(ii)));
    else
        fprintf(fid,[num2str(Left_nodes(ii)) ', ']);
        a=a+1;
        if a==17;
            fprintf(fid,'\n');
            a=1;
        end
    end
end

fprintf(fid,'\n');

%Top face nodes
fprintf(fid,['*Nset, nset=Top_Face_Geometry' '\n']);

Top_nodes=I(nxyz(:,2)>geom.y-0.001);
a=1;
for ii=1:numel(Top_nodes)
    if ii==numel(Top_nodes)
        fprintf(fid,num2str(Top_nodes(ii)));
    else
        fprintf(fid,[num2str(Top_nodes(ii)) ', ']);
        a=a+1;
        if a==17;
            fprintf(fid,'\n');
            a=1;
        end
    end
end

fprintf(fid,'\n');


%Front face nodes
fprintf(fid,['*Nset, nset=Front_Face_Geometry' '\n']);

Front_nodes=I(nxyz(:,1)>geom.x-0.001);
a=1;
for ii=1:numel(Front_nodes)
    if ii==numel(Front_nodes)
        fprintf(fid,num2str(Front_nodes(ii)));
    else
        fprintf(fid,[num2str(Front_nodes(ii)) ', ']);
        a=a+1;
        if a==17;
            fprintf(fid,'\n');
            a=1;
        end
    end
end

fprintf(fid,'\n');


%Back face nodes
fprintf(fid,['*Nset, nset=Back_Face_Geometry' '\n']);

Back_nodes=I(nxyz(:,3)==0);
a=1;
for ii=1:numel(Back_nodes)
    if ii==numel(Back_nodes)
        fprintf(fid,num2str(Back_nodes(ii)));
    else
        fprintf(fid,[num2str(Back_nodes(ii)) ', ']);
        a=a+1;
        if a==17;
            fprintf(fid,'\n');
            a=1;
        end
    end
end

fprintf(fid,'\n');


Geom_Dif=min([geom.x,geom.y,geom.z])/1000;

% Vertices (Corners) of Geometry

V000=I(nxyz(:,1)>-Geom_Dif & nxyz(:,1)<Geom_Dif...
    & nxyz(:,2)>-Geom_Dif & nxyz(:,2)<Geom_Dif...
    & nxyz(:,3)>-Geom_Dif & nxyz(:,3)<Geom_Dif);

V001=I(nxyz(:,1)>-Geom_Dif & nxyz(:,1)<Geom_Dif...
    & nxyz(:,2)>-Geom_Dif & nxyz(:,2)<Geom_Dif...
    & nxyz(:,3)>geom.z-Geom_Dif & nxyz(:,3)<geom.z+Geom_Dif);

V010=I(nxyz(:,1)>-Geom_Dif & nxyz(:,1)<Geom_Dif... 
    & nxyz(:,2)>geom.y-Geom_Dif & nxyz(:,2)<geom.y+Geom_Dif...
    & nxyz(:,3)>-Geom_Dif & nxyz(:,3)<Geom_Dif);

V011=I(nxyz(:,1)>-Geom_Dif & nxyz(:,1)<Geom_Dif...
    & nxyz(:,2)>geom.y-Geom_Dif & nxyz(:,2)<geom.y+Geom_Dif...
    & nxyz(:,3)>geom.z-Geom_Dif & nxyz(:,3)<geom.z+Geom_Dif);

V100=I(nxyz(:,1)>geom.x-Geom_Dif & nxyz(:,1)<geom.x+Geom_Dif...
    & nxyz(:,2)>-Geom_Dif & nxyz(:,2)<Geom_Dif...
    & nxyz(:,3)>-Geom_Dif & nxyz(:,3)<Geom_Dif);

V101=I(nxyz(:,1)>geom.x-Geom_Dif & nxyz(:,1)<geom.x+Geom_Dif...
    & nxyz(:,2)>-Geom_Dif & nxyz(:,2)<Geom_Dif...
    & nxyz(:,3)>geom.z-Geom_Dif & nxyz(:,3)<geom.z+Geom_Dif);

V110=I(nxyz(:,1)>geom.x-Geom_Dif & nxyz(:,1)<geom.x+Geom_Dif...
    & nxyz(:,2)>geom.y-Geom_Dif & nxyz(:,2)<geom.y+Geom_Dif...
    & nxyz(:,3)>-Geom_Dif & nxyz(:,3)<Geom_Dif);

V111=I(nxyz(:,1)>geom.x-Geom_Dif & nxyz(:,1)<geom.x+Geom_Dif...
    & nxyz(:,2)>geom.y-Geom_Dif & nxyz(:,2)<geom.y+Geom_Dif...
    & nxyz(:,3)>geom.z-Geom_Dif & nxyz(:,3)<geom.z+Geom_Dif);
    


fprintf(fid,['*Nset, nset=V000' '\n']);
fprintf(fid,[num2str(V000)  '\n']);
fprintf(fid,['*Nset, nset=V001' '\n']);
fprintf(fid,[num2str(V001)  '\n']);
fprintf(fid,['*Nset, nset=V010' '\n']);
fprintf(fid,[num2str(V010)  '\n']);
fprintf(fid,['*Nset, nset=V011' '\n']);
fprintf(fid,[num2str(V011)  '\n']);
fprintf(fid,['*Nset, nset=V100' '\n']);
fprintf(fid,[num2str(V100)  '\n']);
fprintf(fid,['*Nset, nset=V101' '\n']);
fprintf(fid,[num2str(V101)  '\n']);
fprintf(fid,['*Nset, nset=V110' '\n']);
fprintf(fid,[num2str(V110)  '\n']);
fprintf(fid,['*Nset, nset=V111' '\n']);
fprintf(fid,[num2str(V111)  '\n']);


% Edges of Geometry

EX000=I(nxyz(:,1)>Geom_Dif & nxyz(:,1)<geom.x-Geom_Dif...
    & nxyz(:,2)>-Geom_Dif & nxyz(:,2)<Geom_Dif...
    & nxyz(:,3)>-Geom_Dif & nxyz(:,3)<Geom_Dif);

EX001=I(nxyz(:,1)>Geom_Dif & nxyz(:,1)<geom.x-Geom_Dif...
    & nxyz(:,2)>-Geom_Dif & nxyz(:,2)<Geom_Dif...
    & nxyz(:,3)>geom.z-Geom_Dif & nxyz(:,3)<geom.z+Geom_Dif);

EX010=I(nxyz(:,1)>Geom_Dif & nxyz(:,1)<geom.x-Geom_Dif...
    & nxyz(:,2)>geom.y-Geom_Dif & nxyz(:,2)<geom.y+Geom_Dif...
    & nxyz(:,3)>-Geom_Dif & nxyz(:,3)<Geom_Dif);

EX011=I(nxyz(:,1)>Geom_Dif & nxyz(:,1)<geom.x-Geom_Dif...
    & nxyz(:,2)>geom.y-Geom_Dif & nxyz(:,2)<geom.y+Geom_Dif...
    & nxyz(:,3)>geom.z-Geom_Dif & nxyz(:,3)<geom.z+Geom_Dif);

EY000=I( nxyz(:,1)>-Geom_Dif & nxyz(:,1)<Geom_Dif...
    & nxyz(:,2)>Geom_Dif & nxyz(:,2)<geom.y-Geom_Dif...
    & nxyz(:,3)>-Geom_Dif & nxyz(:,3)<Geom_Dif);

EY001=I( nxyz(:,1)>-Geom_Dif & nxyz(:,1)<Geom_Dif...
    & nxyz(:,2)>Geom_Dif & nxyz(:,2)<geom.y-Geom_Dif...
    & nxyz(:,3)>geom.z-Geom_Dif & nxyz(:,3)<geom.z+Geom_Dif);

EY100=I(nxyz(:,1)>geom.x-Geom_Dif & nxyz(:,1)<geom.x+Geom_Dif...
    & nxyz(:,2)>Geom_Dif & nxyz(:,2)<geom.y-Geom_Dif...
    & nxyz(:,3)>-Geom_Dif & nxyz(:,3)<Geom_Dif);

EY101=I(nxyz(:,1)>geom.x-Geom_Dif & nxyz(:,1)<geom.x+Geom_Dif...
    & nxyz(:,2)>Geom_Dif & nxyz(:,2)<geom.y-Geom_Dif...
    & nxyz(:,3)>geom.z-Geom_Dif & nxyz(:,3)<geom.z+Geom_Dif);
    
EZ000=I(nxyz(:,1)>-Geom_Dif & nxyz(:,1)<Geom_Dif...
    & nxyz(:,2)>-Geom_Dif & nxyz(:,2)<Geom_Dif...
    & nxyz(:,3)>Geom_Dif & nxyz(:,3)<geom.z-Geom_Dif);

EZ010=I(nxyz(:,1)>-Geom_Dif & nxyz(:,1)<Geom_Dif...
    & nxyz(:,2)>geom.y-Geom_Dif & nxyz(:,2)<geom.y+Geom_Dif...
    & nxyz(:,3)>Geom_Dif & nxyz(:,3)<geom.z-Geom_Dif );

EZ100=I(nxyz(:,1)>geom.x-Geom_Dif & nxyz(:,1)<geom.x+Geom_Dif...
    & nxyz(:,2)>-Geom_Dif & nxyz(:,2)<Geom_Dif...
    & nxyz(:,3)>Geom_Dif & nxyz(:,3)<geom.z-Geom_Dif );

EZ110=I(nxyz(:,1)>geom.x-Geom_Dif & nxyz(:,1)<geom.x+Geom_Dif...
    & nxyz(:,2)>geom.y-Geom_Dif & nxyz(:,2)<geom.y+Geom_Dif...
    & nxyz(:,3)>Geom_Dif & nxyz(:,3)<geom.z-Geom_Dif );






 fprintf(fid,['*Nset, nset=EX000' '\n']);
    a=1;
    for kk=1:(length(EX000))
      if kk==length(EX000)
         fprintf(fid,[num2str(EX000(kk))]);
      else
        fprintf(fid,[num2str(EX000(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;
          end
      end
    end
    fprintf(fid,'\n');

fprintf(fid,['*Nset, nset=EX001' '\n']);
    a=1;
    for kk=1:(length(EX001))
      if kk==length(EX001)
         fprintf(fid,[num2str(EX001(kk))]);
      else
        fprintf(fid,[num2str(EX001(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;
          end
      end
    end
    fprintf(fid,'\n');

fprintf(fid,['*Nset, nset=EX010' '\n']);
    a=1;
    for kk=1:(length(EX010))
      if kk==length(EX010)
         fprintf(fid,[num2str(EX010(kk))]);
      else
        fprintf(fid,[num2str(EX010(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;
          end
      end
    end
    fprintf(fid,'\n');






fprintf(fid,['*Nset, nset=EX011' '\n']);
    a=1;
    for kk=1:(length(EX011))
      if kk==length(EX011)
         fprintf(fid,[num2str(EX011(kk))]);
      else
        fprintf(fid,[num2str(EX011(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;
          end
      end
    end
    fprintf(fid,'\n');

 fprintf(fid,['*Nset, nset=EY000' '\n']);
    a=1;
    for kk=1:(length(EY000))
      if kk==length(EY000)
         fprintf(fid,[num2str(EY000(kk))]);
      else
        fprintf(fid,[num2str(EY000(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;   
          end
      end
    end
    fprintf(fid,'\n');

fprintf(fid,['*Nset, nset=EY001' '\n']);
    a=1;
    for kk=1:(length(EY001))
      if kk==length(EY001)
         fprintf(fid,[num2str(EY001(kk))]);
      else
        fprintf(fid,[num2str(EY001(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;   
          end
      end
    end
    fprintf(fid,'\n');


     
 fprintf(fid,['*Nset, nset=EY100' '\n']);
    a=1;
    for kk=1:(length(EY100))
      if kk==length(EY100)
         fprintf(fid,[num2str(EY100(kk))]);
      else
        fprintf(fid,[num2str(EY100(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;   
          end
      end
    end
    fprintf(fid,'\n');

fprintf(fid,['*Nset, nset=EY101' '\n']);
    a=1;
    for kk=1:(length(EY101))
      if kk==length(EY101)
         fprintf(fid,[num2str(EY101(kk))]);
      else
        fprintf(fid,[num2str(EY101(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;   
          end
      end
    end
    fprintf(fid,'\n');

fprintf(fid,['*Nset, nset=EZ000' '\n']);
    a=1;
    for kk=1:(length(EZ000))
      if kk==length(EZ000)
         fprintf(fid,[num2str(EZ000(kk))]);
      else
        fprintf(fid,[num2str(EZ000(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;
          end
      end
    end
    fprintf(fid,'\n');


fprintf(fid,['*Nset, nset=EZ010' '\n']);
    a=1;
    for kk=1:(length(EZ010))
      if kk==length(EZ010)
         fprintf(fid,[num2str(EZ010(kk))]);
      else
        fprintf(fid,[num2str(EZ010(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;
          end
      end
    end
    fprintf(fid,'\n');

fprintf(fid,['*Nset, nset=EZ100' '\n']);
    a=1;
    for kk=1:(length(EZ100))
      if kk==length(EZ100)
         fprintf(fid,[num2str(EZ100(kk))]);
      else
        fprintf(fid,[num2str(EZ100(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;
          end
      end
    end
    fprintf(fid,'\n');

fprintf(fid,['*Nset, nset=EZ110' '\n']);
    a=1;
    for kk=1:(length(EZ110))
      if kk==length(EZ110)
         fprintf(fid,[num2str(EZ110(kk))]);
      else
        fprintf(fid,[num2str(EZ110(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;
          end
      end
    end
    fprintf(fid,'\n');


     
    % faces of Geometry

FXN=I(nxyz(:,1)>-Geom_Dif & nxyz(:,1)<Geom_Dif...
    & nxyz(:,2)>Geom_Dif & nxyz(:,2)<geom.y-Geom_Dif...
    & nxyz(:,3)>Geom_Dif & nxyz(:,3)<geom.z-Geom_Dif);

FXP=I(nxyz(:,1)>geom.x-Geom_Dif & nxyz(:,1)<geom.x+Geom_Dif...
    & nxyz(:,2)>Geom_Dif & nxyz(:,2)<geom.y-Geom_Dif...
    & nxyz(:,3)>Geom_Dif & nxyz(:,3)<geom.z-Geom_Dif);

FYN=I(nxyz(:,1)>Geom_Dif & nxyz(:,1)<geom.x-Geom_Dif...
    & nxyz(:,2)>-Geom_Dif & nxyz(:,2)<Geom_Dif...
    & nxyz(:,3)>Geom_Dif & nxyz(:,3)<geom.z-Geom_Dif);

FYP=I(nxyz(:,1)>Geom_Dif & nxyz(:,1)<geom.x-Geom_Dif...
    & nxyz(:,2)>geom.y-Geom_Dif & nxyz(:,2)<geom.y+Geom_Dif...
    & nxyz(:,3)>Geom_Dif & nxyz(:,3)<geom.z-Geom_Dif);

FZN=I(nxyz(:,1)>Geom_Dif & nxyz(:,1)<geom.x-Geom_Dif...
    & nxyz(:,2)>Geom_Dif & nxyz(:,2)<geom.y-Geom_Dif...
    & nxyz(:,3)>-Geom_Dif & nxyz(:,3)<Geom_Dif);

FZP=I(nxyz(:,1)>Geom_Dif & nxyz(:,1)<geom.x-Geom_Dif...
    & nxyz(:,2)>Geom_Dif & nxyz(:,2)<geom.y-Geom_Dif...
    & nxyz(:,3)>geom.z-Geom_Dif & nxyz(:,3)<geom.z+Geom_Dif);


Write_Surface_Sets(Node_El, FXN, 'FXN', 'S2', fid)
Write_Surface_Sets(Node_El, FXP, 'FXP', 'S1', fid)
Write_Surface_Sets(Node_El, FYN, 'FYN', 'S4', fid)
Write_Surface_Sets(Node_El, FYP, 'FYP', 'S6', fid)
Write_Surface_Sets(Node_El, FZN, 'FZN', 'S3', fid)
Write_Surface_Sets(Node_El, FZP, 'FZP', 'S5', fid)

fprintf(fid,['*Nset, nset=FXN' '\n']);
    a=1;
    for kk=1:(length(FXN))
      if kk==length(FXN)
         fprintf(fid,[num2str(FXN(kk))]);
      else
        fprintf(fid,[num2str(FXN(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;
          end
      end
    end
    fprintf(fid,'\n');
     


fprintf(fid,['*Nset, nset=FXP' '\n']);
    a=1;
    for kk=1:(length(FXP))
      if kk==length(FXP)
         fprintf(fid,[num2str(FXP(kk))]);
      else
        fprintf(fid,[num2str(FXP(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;
          end
      end
    end
    fprintf(fid,'\n');

fprintf(fid,['*Nset, nset=FYN' '\n']);
    a=1;
    for kk=1:(length(FYN))
      if kk==length(FYN)
         fprintf(fid,[num2str(FYN(kk))]);
      else
        fprintf(fid,[num2str(FYN(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;
          end
      end
    end
    fprintf(fid,'\n');

fprintf(fid,['*Nset, nset=FYP' '\n']);
    a=1;
    for kk=1:(length(FYP))
      if kk==length(FYP)
         fprintf(fid,[num2str(FYP(kk)) '\n']);
      else
        fprintf(fid,[num2str(FYP(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;
          end
      end
    end
    fprintf(fid,'\n');



fprintf(fid,['*Nset, nset=FZN' '\n']);
    a=1;
    for kk=1:(length(FZN))
      if kk==length(FZN)
         fprintf(fid,[num2str(FZN(kk))]);
      else
       fprintf(fid,[num2str(FZN(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;
          end
      end
    end
    fprintf(fid,'\n');

fprintf(fid,['*Nset, nset=FZP' '\n']);
    a=1;
   for kk=1:(length(FZP))
      if kk==length(FZP)
         fprintf(fid,[num2str(FZP(kk))]);
      else
        fprintf(fid,[num2str(FZP(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;
          end
      end
    end
    fprintf(fid,'\n');
    
% assign node sets to use when periodicity is broken in the x direction
% vectors in z direction at varying heights
x = geom.x;
y = geom.y;
res = mesh.CP_mesh_size;
for i=1:1:y/res-1
    temp_nset = I(nxyz(:,1)>x-Geom_Dif & nxyz(:,1)<x+Geom_Dif...
    & nxyz(:,2)>i*res-Geom_Dif & nxyz(:,2)<i*res+Geom_Dif);
    fprintf(fid,['*Nset, nset=VZ100_' num2str(i) '\n']);
    a=1;
    for kk=1:(length(temp_nset))
      if kk==length(temp_nset)
         fprintf(fid,num2str(temp_nset(kk)));
      else
        fprintf(fid,[num2str(temp_nset(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;
          end
      end
    end
    fprintf(fid,'\n');
end

% assign node sets for crack tip top and bototm
if ~isempty(strfind(type, 'crack'))
    x = geom.x;
    y = geom.y;
    res = mesh.CP_mesh_size;
    featureSize = str2double(mesh.featureSize)*x/100;
    height = str2double(options.crack_height)/100*y;
    
    if ~isempty(strfind(type, 'rough'))
        elemsRemove = ceil(height/res);
        center = y/2+(-mod(y/res,2)/2+mod(elemsRemove,2)/2)*res;
        crack_top = center + elemsRemove/2*res;
        crack_bottom = center - elemsRemove/2*res;
        crack_in = (floor(featureSize/res-.5)+1)*res;
        
        crack_top_nodes=I(nxyz(:,1)>crack_in-res-Geom_Dif & nxyz(:,1)<crack_in-res+Geom_Dif...
        & nxyz(:,2)>crack_top-Geom_Dif & nxyz(:,2)<crack_top+Geom_Dif);

        fprintf(fid,['*Nset, nset=crack_tip_top_1' '\n']);
        a=1;
        for kk=1:(length(crack_top_nodes))
          if kk==length(crack_top_nodes)
             fprintf(fid,[num2str(crack_top_nodes(kk))]);
          else
            fprintf(fid,[num2str(crack_top_nodes(kk)) ', ' ]);
              a=a+1;
              if a==17;
                  fprintf(fid,'\n');
                  a=1;
              end
          end
        end
        fprintf(fid,'\n');

        crack_bottom_nodes=I(nxyz(:,1)>crack_in-res-Geom_Dif & nxyz(:,1)<crack_in-res+Geom_Dif...
        & nxyz(:,2)>crack_bottom-Geom_Dif & nxyz(:,2)<crack_bottom+Geom_Dif);


        fprintf(fid,['*Nset, nset=crack_tip_bottom_1' '\n']);
        a=1;
        for kk=1:(length(crack_bottom_nodes))
          if kk==length(crack_bottom_nodes)
             fprintf(fid,[num2str(crack_bottom_nodes(kk))]);
          else
            fprintf(fid,[num2str(crack_bottom_nodes(kk)) ', ' ]);
              a=a+1;
              if a==17;
                  fprintf(fid,'\n');
                  a=1;
              end
          end
        end
        fprintf(fid,'\n');
    
    else
        center = y/2;
        crack_top = center + height/2;
        crack_bottom = center - height/2;
        crack_in = featureSize - height/2;
    end
    
    crack_top_nodes=I(nxyz(:,1)>crack_in-Geom_Dif & nxyz(:,1)<crack_in+Geom_Dif...
    & nxyz(:,2)>crack_top-Geom_Dif & nxyz(:,2)<crack_top+Geom_Dif);

    fprintf(fid,['*Nset, nset=crack_tip_top' '\n']);
    a=1;
    for kk=1:(length(crack_top_nodes))
      if kk==length(crack_top_nodes)
         fprintf(fid,[num2str(crack_top_nodes(kk))]);
      else
        fprintf(fid,[num2str(crack_top_nodes(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;
          end
      end
    end
    fprintf(fid,'\n');
    
    crack_bottom_nodes=I(nxyz(:,1)>crack_in-Geom_Dif & nxyz(:,1)<crack_in+Geom_Dif...
    & nxyz(:,2)>crack_bottom-Geom_Dif & nxyz(:,2)<crack_bottom+Geom_Dif);


    fprintf(fid,['*Nset, nset=crack_tip_bottom' '\n']);
    a=1;
    for kk=1:(length(crack_bottom_nodes))
      if kk==length(crack_bottom_nodes)
         fprintf(fid,[num2str(crack_bottom_nodes(kk))]);
      else
        fprintf(fid,[num2str(crack_bottom_nodes(kk)) ', ' ]);
          a=a+1;
          if a==17;
              fprintf(fid,'\n');
              a=1;
          end
      end
    end
    fprintf(fid,'\n');
end
       

%fprintf(fid,['*End Assembly' '\n' '**' '\n' '** MATERIALS' '\n' '**' '\n']);

% Create CP materials





if(strcmp(options.banded, 'TRUE'))
    Cnum=zeros(n_grains,Max_num_layers,4);

    for i=1:n_El
    for j=1:4
    Cnum(ElemGrainNo(i,1),NumLayer(i,j),j) = Cnum(ElemGrainNo(i,1),NumLayer(i,j),j)+1.0;
    end
    end

    Cnum_max=max(max(max(Cnum)));
end


%fprintf(fid,['**CHECK' '\n']);
%for ii=1:n_Nodes

%fprintf(fid,['** ' num2str(NodeGrain(ii,1)) ' , ' num2str(NodeGrain(ii,2)) ' , ' num2str(NodeGrain(ii,3)) ' , ' num2str(NodeGrain(ii,4)) ' , ' num2str(NodeGrain(ii,5)) ' , ' num2str(NodeGrain(ii,6)) ' , ' num2str(NodeGrain(ii,7)) ' , ' num2str(NodeGrain(ii,8)) '\n']);

%end



fclose(fid);

%% Write out auxiliary files
'Writing auxiliary files'

if(strcmp(options.banded, 'TRUE'))
    Min_dist_bands(1:n_grains,1:n_grains,1:Max_num_layers,1:Max_num_layers,1:4,1:4)=false;

    
    Grain_angle=zeros(n_grains,3);


    disp('    Min_dist progress:')
    percent = 1; %increment to display min_dist progress

    for Z=1:(counter-1);

        ii = DUMMY(Z,1);
        jj = DUMMY(Z,2);



        %Print current progress on Min_dist:
        if (Z/(counter-1))*100 >= percent    
        PROG = ['     ', sprintf('%d',percent)];
        disp(PROG)
        percent = percent + 1;
        end

    %	for jj=(ii+1):n_El


            if dist(ii,jj)*mesh.CP_mesh_size<band_jump+0.01 && dist(ii,jj) > 0
            if ElemGrainNo(ii,1)~=ElemGrainNo(jj,1)		
    %		if Neighbor_g(ElemGrainNo(ii,1),ElemGrainNo(jj,1))==1		
    %        full(dist(ii,jj))
                for mm=1:4
                for nn=1:4

%                Min_dist_bands(ElemGrainNo(ii,1),ElemGrainNo(jj,1),NumLayer(ii,mm),NumLayer(jj,nn),mm,nn) =true;

                end
                end

    %		end
            end 
            end

    %	end



    end
end


% only write out the necessary files
if(strcmp(options.banded, 'TRUE'))
    fid = fopen(['Num_layer_' num2str(options.MS_number) '.txt'], 'W');  % Read temp input file
    for ii=1:n_El
        fprintf(fid,[  sprintf('%i', ii)  ' , ' sprintf('%i', ElemGrainNo(ii,1))  ' , '   sprintf('%i', NumLayer(ii,1))  ' , ' sprintf('%i', NumLayer(ii,2))  ' , ' sprintf('%i', NumLayer(ii,3))  ' , ' sprintf('%i',NumLayer(ii,4))  ' , '  sprintf('%.4f',elem_orient(ii,1)) ',  ' sprintf('%.4f',elem_orient(ii,2)) ' , '  sprintf('%.4f',elem_orient(ii,3))  ' \n'  ]);

        
    end
    fclose(fid);
end

% only write out if the diameters have changed from padding
if(strcmp(options.banded, 'TRUE') || (isfield(options,'free_surface_padding') && (str2double(options.free_surface_padding) > 0)) || (isfield(options,'repeat') && strcmp(options.repeat,'TRUE')) )
    
    Grain_diam=2*((3/(4*pi)*V_grn).^(1/3));

    fid = fopen(['Grains_' num2str(options.MS_number) '.txt'], 'W');  % Read temp input file
    for ii=1:n_grains
        fprintf(fid,[ sprintf('%0.7f', Grain_diam(ii) ) ' , ' sprintf('%0.4f',g_orient(ii,2) ) ' , ' sprintf('%0.4f', g_orient(ii,3) )  ' , ' sprintf('%0.4f', g_orient(ii,4))    ' \n'  ]);
    end
    fclose(fid);
end


%fid = fopen(['Dist.txt'], 'W');  % Read temp input file
%for ii=1:n_El
%for jj=1:n_El
%if dist(ii,jj)<1.5
%fprintf(fid,[ num2str(ii)  ' , '  num2str(jj)  ' , ' num2str( dist(ii,jj))  ' \n'  ]);
%end
%end
%end
%fclose(fid);



% only write out the necessary files
if(strcmp(options.banded, 'TRUE'))

    num_lines_Min_dist=0;
	
	Min_dist_bands_written(1:n_grains,1:n_grains,1:Max_num_layers,1:Max_num_layers,1:4,1:4)=false;
	
    fid = fopen(['Min_dist_' num2str(options.MS_number) '.txt'], 'W');  % Read temp input file

	disp(' ')
	disp('    Dist loop progress:')
	percent  = 1;

	counter = 1;

	%Make the distance for non-neighbor elements twice the mesh size
	%Mesh_Len_2 = 2*mesh.CP_mesh_size;
	%Mesh_Len_2*dist;

	for jj=1:n_El

		if (jj/n_El)*100 >= percent    
		PROG = ['       ', sprintf('%d', percent)];
		disp(PROG)
		percent = percent + 1;
		end

		dummy1=int32(1);
		dummyX1=int32(1);
		dummyY1=int32(1);
		dummyZ1=int32(1);


		% we already have neighboring element lists, so find neighboring grains
		% list
		for kk = 1:Element_Neighbors(jj,27)
			ii = Element_Neighbors(jj,kk);
		% 	dist1=(El_pos(jj,2)-El_pos(ii,2))/mesh.CP_mesh_size;
		% 	dist2=(El_pos(jj,3)-El_pos(ii,3))/mesh.CP_mesh_size;
		% 	dist3=(El_pos(jj,4)-El_pos(ii,4))/mesh.CP_mesh_size;
		% 	
		% 	dist4 = abs(dist1)+abs(dist2)+abs(dist3);
		% 	
		% %     if (dist4 < 1 || jj < 0 || ii < 0)
		% %        disp(' Uh oh ') 
		% %     end
		%     
		% 	dist(jj,ii) = dist_vals(dist4);
		% 	dist(jj,ii);
		% 
		% 		
		% %	dist(jj,ii)=((dist1)^2+(dist2)^2+(dist3)^2)^0.5;
		% 	
		% 	
		% 	if (dist(jj,ii)<1.9 && jj~=ii) %they are neighbors
		% 	%pp=1;
		% 	%while (Neighbors_el(jj,pp)~=0 )
		% 	%pp=pp+1;
		% 	%end
		% 	%Neighbors_el(jj,pp)=ii;
		% 
		% 	[chk1, chk2]=min(sum(abs(Neigh_vec(:,1:3)+ [ones(26,1)*dist1 ones(26,1)*dist2 ones(26,1)*dist3] ),2));
		% 	Neighbors_el(jj,chk2)=ii;
		% 	end


			%if (ElemGrainNo(jj,1)~=ElemGrainNo(ii,1) && dist(jj,ii)<1.9) %they are neighbors
			if (ElemGrainNo(jj,1)~=ElemGrainNo(ii,1)) && Element_Neighbors_common_nodes(jj,kk)==4 %face-connevtivity
				Neighbor_g(ElemGrainNo(jj,1),ElemGrainNo(ii,1))=1; 
                if (ElemGrainNo(jj,1)<ElemGrainNo(ii,1)) %so Min_dist only written once 

                    for mm=1:4
                        for nn=1:4
                                
                            if Min_dist_bands_written(ElemGrainNo(jj,1),ElemGrainNo(ii,1),NumLayer(jj,mm),NumLayer(ii,nn),mm,nn)==false                                

                                fprintf(fid,[ sprintf('%d',ElemGrainNo(jj,1)) ' , '  sprintf('%d',ElemGrainNo(ii,1)) ' , ' sprintf('%d',NumLayer(jj,mm)) ' , ' sprintf('%d',NumLayer(ii,nn)) ' , ' sprintf('%d',mm)    ' , '   sprintf('%d',nn)    ' \n'  ]);

                                num_lines_Min_dist=num_lines_Min_dist+1;

                                Min_dist_bands_written(ElemGrainNo(jj,1),ElemGrainNo(ii,1),NumLayer(jj,mm),NumLayer(ii,nn),mm,nn)=true;

                            end
                        end
                    end
                end
			end
			



		
		end


		%for ii=1:n_El

		%dist1=(El_pos(jj,2)-El_pos(ii,2))/mesh.CP_mesh_size;
		%dist2=(El_pos(jj,3)-El_pos(ii,3))/mesh.CP_mesh_size;
		%dist3=(El_pos(jj,4)-El_pos(ii,4))/mesh.CP_mesh_size;

		%dist(jj,ii)=((dist1)^2+(dist2)^2+(dist3)^2)^0.5;



		% if (dist(jj,ii)<1.9 && jj~=ii) %they are neighbors
		% %pp=1;
		% %while (Neighbors_el(jj,pp)~=0 )
		% %pp=pp+1;
		% %end
		% %Neighbors_el(jj,pp)=ii;

		% [chk1, chk2]=min(sum(abs(Neigh_vec(:,1:3)+ [ones(26,1)*dist1 ones(26,1)*dist2 ones(26,1)*dist3] ),2));
		% Neighbors_el(jj,chk2)=ii;
		% end


		% if (ElemGrainNo(jj,1)~=ElemGrainNo(ii,1) && dist(jj,ii)<1.9) %they are neighbors
			% Neighbor_g(ElemGrainNo(jj,1),ElemGrainNo(ii,1))=1;
			% DUMMY(counter,1) = ii;
			% DUMMY(counter,2) = jj;
			% counter = counter + 1;
		% end


		% if (ElemGrainNo(jj,1)~=ElemGrainNo(ii,1) && abs(dist1)<1.01 && abs(dist2)<0.01 && abs(dist3)<0.01) %they are neighbors
		% N_Elem_Bound(ElemGrainNo(jj,1),ElemGrainNo(ii,1))=N_Elem_Bound(ElemGrainNo(jj,1),ElemGrainNo(ii,1))+1; % number of elements on the bounday set
		% GB(ElemGrainNo(jj,1),ElemGrainNo(ii,1),N_Elem_Bound(ElemGrainNo(jj,1),ElemGrainNo(ii,1)))=jj;
		% El_GB(jj,dummy1)=ElemGrainNo(ii,1) ;
		% dummy1=dummy1+1;
		% GB_dirX(jj,dummyX1)=sign(dist1)*1;  %ii;
		% dummyX1=dummyX1+1;

		% elseif (ElemGrainNo(jj,1)~=ElemGrainNo(ii,1) && abs(dist2)<1.01 && abs(dist1)<0.01 && abs(dist3)<0.01 ) %they are neighbors
		% N_Elem_Bound(ElemGrainNo(jj,1),ElemGrainNo(ii,1))=N_Elem_Bound(ElemGrainNo(jj,1),ElemGrainNo(ii,1))+1; % number of elements on the bounday set
		% GB(ElemGrainNo(jj,1),ElemGrainNo(ii,1),N_Elem_Bound(ElemGrainNo(jj,1),ElemGrainNo(ii,1)))=jj;
		% El_GB(jj,dummy1)=ElemGrainNo(ii,1) ;
		% dummy1=dummy1+1;
		% GB_dirY(jj,dummyY1)=sign(dist2)*1; %ii;
		% dummyY1=dummyY1+1;

		% elseif (ElemGrainNo(jj,1)~=ElemGrainNo(ii,1) && abs(dist3)<1.01 && abs(dist2)<0.01 && abs(dist1)<0.01 ) %they are neighbors
		% N_Elem_Bound(ElemGrainNo(jj,1),ElemGrainNo(ii,1))=N_Elem_Bound(ElemGrainNo(jj,1),ElemGrainNo(ii,1))+1; % number of elements on the bounday set
		% GB(ElemGrainNo(jj,1),ElemGrainNo(ii,1),N_Elem_Bound(ElemGrainNo(jj,1),ElemGrainNo(ii,1)))=jj;
		% El_GB(jj,dummy1)=ElemGrainNo(ii,1) ;
		% dummy1=dummy1+1;
		% GB_dirZ(jj,dummyZ1)=sign(dist3)*1; %ii;
		% dummyZ1=dummyZ1+1;
		% end

	end % loop over all elements

	fclose(fid);   

end

Neighbor_num=squeeze(sum(int16(Neighbor_g),1));

fid = fopen(['Neighbor_grains_' num2str(options.MS_number) '.txt'], 'W');  % Read temp input file

for jj=1:n_grains
fprintf(fid,[sprintf('%d',Neighbor_num(jj)) '\n']);
%fprintf(fid,[num2str(jj)  ', ']);
    for ii=1:n_grains
        if (jj~=ii &&  Neighbor_g(jj,ii)==1)
        fprintf(fid,[sprintf('%d',ii)  ', ']);
        end
    end
fprintf(fid,['\n']);
end

fclose(fid);

return
end

function Write_Surface_Sets(Node_El, Face_set, name, face, fid)
    fprintf(fid,['*Elset, elset=' name '_surf_elset' '\n']);
    a=1;
    for kk=1:(length(Face_set))
      for jj=1:8
          if (kk==length(Face_set) && (jj==8 || Node_El(Face_set(kk),jj+1)==0))
             fprintf(fid,num2str(Node_El(Face_set(kk),jj)));
             break
          elseif(Node_El(Face_set(kk),jj)~=0)
            fprintf(fid,[num2str(Node_El(Face_set(kk),jj)) ', ' ]);
              a=a+1;
              if a==17;
                  fprintf(fid,'\n');
                  a=1;
              end
          end
      end
    end
    fprintf(fid,'\n');
    fprintf(fid, ['*Surface, type=ELEMENT, name=Surf_' name '\n' name '_surf_elset, ' face '\n']);
end
