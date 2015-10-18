function[NodeGrain,El_pos]=Grain_Boundary(geom,mesh,nc,n_El,n_grains,n_Nodes,ElemGrainNo,el_centroid, options)
NodeGrain =[];
% if(strcmp(options.banded, 'TRUE'))
% 
%     NodeGrain=zeros(n_Nodes,8);
% 
%     %find the eight possible grains that this node can border/be contained by
%     for ii=1:n_El
%     % for each node making up the element
%     for jj=1:8
%       % find location in matrix to store value
%       for kk=1:8       
%          if NodeGrain(nc(ii,jj),kk)==0,break,end         
%       end
%          NodeGrain(nc(ii,jj),kk)=ElemGrainNo(ii,1);            
%     end   
%     end
% else
%     NodeGrain = [0];
% end

El_pos(:,1)=el_centroid(:,4);  % This is the element number
El_pos(:,2)=el_centroid(:,1);  % X coordinate of the centroid 
El_pos(:,3)=el_centroid(:,2);  % Y coordinate of the centroid
El_pos(:,4)=el_centroid(:,3);  % Z coordinate of the centroid

Sort_el_pos=sortrows(El_pos,[2 3 4]);

fid = fopen(['El_pos_' num2str(options.MS_number) '.txt'], 'w');  % Write temp input file
for jj=1:n_El
  fprintf(fid,['Elem' num2str(El_pos(jj,1)) ' ' num2str(El_pos(jj,2)) ' ' num2str(El_pos(jj,3)) ' ' num2str(El_pos(jj,4)) ' ' '\n']);
end
fclose(fid);

fid = fopen(['Sort_el_pos_' num2str(options.MS_number) '.txt'], 'w');  % Write temp input file
for jj=1:n_El
  fprintf(fid,['Elem' num2str(Sort_el_pos(jj,1)) ' ' num2str(Sort_el_pos(jj,2)) ' ' num2str(Sort_el_pos(jj,3)) ' ' num2str(Sort_el_pos(jj,4)) ' ' '\n']);
end
fclose(fid);


% only write out the necessary files
% if(strcmp(options.banded, 'TRUE'))
%     x0=Sort_el_pos(1,2);
%     y0=Sort_el_pos(1,3);
%     z0=Sort_el_pos(1,4);
% 
%     nx=0;
%     ny=0;
%     nz=0;
% 
% 
%     for ii=1:(n_El-1)
%         if (ElemGrainNo(ii,1) ~= ElemGrainNo(ii+1,1)) 
%             if  ((Sort_el_pos(ii,2)-Sort_el_pos(1,2))< geom.x-1.5*mesh.CP_mesh_size)
%                 nx=nx+1;
%                 Boundary_x(nx)=ii;
%             end 
% 
%             if  ((Sort_el_pos(ii,3)-Sort_el_pos(1,3))< geom.y-1.5*mesh.CP_mesh_size)
%                 ny=ny+1;
%                 Boundary_y(ny)=ii;
%             end
% 
%             if  ((Sort_el_pos(ii,4)-Sort_el_pos(1,4))< geom.z-1.5*mesh.CP_mesh_size)
%                 nz=nz+1;
%                 Boundary_z(nz)=ii;
%             end
% 
%         end
%     end
% 
%     fid = fopen(['Boundary_x' num2str(options.MS_number) '.txt'], 'w');  % Write temp input file
%     for jj=1:nx
%       fprintf(fid,[num2str(Boundary_x(jj))  '\n']);
%     end
% 
%     fid = fopen(['Boundary_y' num2str(options.MS_number) '.txt'], 'w');  % Write temp input file
%     for jj=1:ny
%       fprintf(fid,[num2str(Boundary_y(jj))  '\n']);
%     end
% 
%     fid = fopen(['Boundary_z' num2str(options.MS_number) '.txt'], 'w');  % Write temp input file
%     for jj=1:nz
%       fprintf(fid,[num2str(Boundary_z(jj))  '\n']);
%     end
%     fclose(fid);
% end


%This part creates a file that has the neighbouring grains for each grain


%for jj=1:n_El
%for ii=1:n_El
%dist=((El_pos(jj,2)-El_pos(ii,2))^2+(El_pos(jj,3)-El_pos(ii,3))^2+(El_pos(jj,4)-El_pos(ii,4))^2)^0.5/mesh.CP_mesh_size;
%if (ElemGrainNo(jj,1)~=ElemGrainNo(ii,1) & dist<1.9) %they are neighbors
%Neighbor_g(ElemGrainNo(jj,1),ElemGrainNo(ii,1))=1;
%end

%if (ElemGrainNo(jj,1)~=ElemGrainNo(ii,1) & dist<1.01) %they are neighbors
;
%end

%end
%end

%Neighbor_num=squeeze(sum(Neighbor_g,1));

%fid = fopen(['Neighbor_grains.txt'], 'w');  % Read temp input file

%for jj=1:n_grains
%fprintf(fid,[num2str(Neighbor_num(jj)) '\n']);

%for ii=1:n_grains
%if (jj~=ii &  Neighbor_g(jj,ii)==1)
%fprintf(fid,[num2str(ii)  ', ']);
%end
%end
%fprintf(fid,['\n']);
%end

%fclose(fid);


 return
 end

