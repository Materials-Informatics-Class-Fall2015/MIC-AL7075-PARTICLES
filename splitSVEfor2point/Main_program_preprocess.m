% Author: Bill Musinski in "Novel methods for microstructure-sensitive
% probabilistic fatigue notch factor", Ms Thesis Georgia Institute of technology, 2010.
% Adapted and modified by: Gustavo M. Castelluccio 2010-2011
% Date: 9/18/09
% This program does all preprocessing work to make the coarse
% grain geometry that assigns grains according to combined
% circular packing and voronoi algorithm.  It creates the 
% input file to run.

function[]=Main_program_preprocess(Input_file_name,n_header_lines,geom,mesh,band_jump,Big_grain1,Big_grain2,Big_grain3,...
d_gr,sigma_fit,mu_fit,SBperiod,Nom_strain,Nom_strain_min,Strain_rate,Sigma_y,alpha_g,Cont_FS_trans,...
Phi_irr,DCTD_th,A_fs,refine,scaling,n_grain_fail,d_gr_ref,Eval_Nuc,Eval_MSC,...
Pg,Phi_irr_int,C_int_0,P_beta,Q_ox,B_int,Cepsilon_int,Clambda_int,t2,Nmulti_begins,Mantle_slip_sys,Core_slip_sys,Stage_II,n_elem_fail, options)

%% Main Program
A=cell(1,1);
job.Job_Name=A{1,1};

%% Input values for script/model

tic
SimTime.clock1=clock;
job.Job_Name=Input_file_name;
job.Job_Name_temp='main_temp';

job.inp_fid=strcat(job.Job_Name, '.inp');
job.temp_fid=strcat(job.Job_Name_temp, '.inp');

mesh.type = options.geometry;
if(isfield(options,'aspect_ratio'))
    mesh.ar = options.aspect_ratio;
end
if(isfield(options,'feature_percent'))
    mesh.featureSize = options.feature_percent;
end

% Create geometry mesh only 
disp(' ')
disp('Call CreateGeom')
Create_Geometry_CG(job,geom,mesh,options);
disp(' ')
disp('Return CreateGeom')

% Find all interesting data of mesh
disp(' ')
disp('Call Nodal_Values_preprocess')
[nxyz,nc,n_Nodes,n_El,el_centroid,V_el,Element_Neighbors, Element_Neighbors_common_nodes]=Nodal_Values_preprocess(job,n_header_lines, options);
disp(' ')
disp('Return Nodal_Values_preprocess')

% Create neighbor list
%[neigh]=Neighbor_List(nc,n_El);

%check if D3D file exists with same microstructure number, if so use this
%to determine grains, else do grain packing
grain_file = ['trial_elem_grains_' num2str(options.MS_number) '.txt'];
files = dir(grain_file);
D3D_exists = ~isempty(files);

% Circular Packing
if(~D3D_exists)

    n_grains=int32(sum(V_el)/(4/3*pi*d_gr^3));
    if(n_grains < 1)
        n_grains = 1;
    end
    %n_grains=1; % Uncomment this line to use a fixed number of grains
    n_Orient = n_grains;
    max_Iter=n_grains*5; 

    [grain,ElemGrainNo,V_grn,RelLength]=Circular_Packing_CG(d_gr,sigma_fit,mu_fit,geom,n_grains, el_centroid, max_Iter,n_El,V_el,Big_grain1,Big_grain2,Big_grain3, options);
    new_to_old_grains = zeros(n_grains,1);
    for i=1:n_grains
        new_to_old_grains(i) = (i);
    end
else
   % new function to calculate element grains based on Dream3D microstructure
   [grain,ElemGrainNo,RelLength,V_grn,new_to_old_grains]=MS_Overlay(el_centroid,n_El,V_el,grain_file, geom, options);
   n_grains = length(grain);
   n_Orient = n_grains;
end

NodeGrain = [];
El_pos = [];
if(strcmp(options.banded, 'TRUE'))
    disp(' ')
    disp('Call Grain_Boundary')
    [NodeGrain,El_pos]=Grain_Boundary(geom,mesh,nc,n_El,n_grains,n_Nodes,ElemGrainNo,el_centroid, options);
    disp(' ')
    disp('Return Grain_Boundary')
end

% Assign grains
Grain_Array=cell(n_grains,1);

for jj=1:n_grains
    el=find(ElemGrainNo==jj);
    %Transform local element numbers (CP) to global element numbers
    Grain_Array{jj,1}=el;
end

 % Append rest of information on input file
 [Max_num_layers,num_lines_Min_dist,NumLayer_max,NumLayer_min]=Finish_input_file_CG(job,geom,RelLength,n_Orient,nxyz,n_Nodes,n_El,Grain_Array,n_grains,ElemGrainNo,...
NodeGrain,El_pos,mesh,SBperiod,Nom_strain,Nom_strain_min,Strain_rate,...
band_jump,grain,scaling,n_grain_fail,t2,Mantle_slip_sys,Core_slip_sys,Element_Neighbors, Element_Neighbors_common_nodes, options, V_el, V_grn, new_to_old_grains);
 
SimTime.clock2=clock;

% Save information
save(job.Job_Name);

% Delete files that are not needed
delete ('abaqus.rpy');
%delete ([job.Job_Name '_temp.inp']);
delete ([job.Job_Name '_temp_script.py']);

disp(' ')
disp('write Geom_Def')


    fid = fopen(['Geom_Def_' num2str(options.MS_number) '.txt'], 'w');  % Read temp input file
    fprintf(fid,['        parameter(           '  '\n']);
    fprintf(fid,['     &            num_slip_sys    = 12, ! Number of slip planes' '\n']);
    fprintf(fid,['     &            Mantle_slip_sys    = ' num2str(Mantle_slip_sys) ', ! Number of slip planes' '\n']);
    fprintf(fid,['     &            Core_slip_sys    = ' num2str(Core_slip_sys) ', ! Number of slip planes' '\n']);
    fprintf(fid,['     &            max_loops       = 10, ! Number of loops to iterate in the line search ' '\n']);
    fprintf(fid,['     &            Tolerance       = 0.000001, ! ' '\n']);
    fprintf(fid,['     &            num_elem       = ' num2str(n_El) ', ! Number of elements' '\n']);
    fprintf(fid,['     &            num_grains       = ' num2str(n_grains) ',  ! Number of grains' '\n']);
    fprintf(fid,['     &            max_num_layers  = ' num2str(Max_num_layers+1) ', ! ! Maximum number of layers' '\n']);
	fprintf(fid,['     &            NumLayer_max  = ' num2str(NumLayer_max) ', ! ! Maximum number of layers' '\n']);	
	fprintf(fid,['     &            NumLayer_min  = ' num2str(NumLayer_min) ', ! ! Minimum number of layers' '\n']);
    fprintf(fid,['     &            num_grain_fail = ' num2str(n_grain_fail)  ', ! ! Max number of grains to fail' '\n']);
    fprintf(fid,['     &            num_elem_fail = ' num2str(n_elem_fail)  ', ! ! Max number of elements to fail' '\n']);
    fprintf(fid,['     &            mesh_size  = ' num2str(mesh.CP_mesh_size)  ', ! ! mesh_size ' '\n']);
    fprintf(fid,['     &            band_width  = ' num2str(SBperiod)  ', ! ! width of the bands ' '\n']);
    fprintf(fid,['     &            num_lines_Min_dist  = ' num2str(num_lines_Min_dist)  ' ) ! ! Number of lines to read in Min_dis.txt ' '\n']);
    fclose(fid);    

%     fid = fopen(['Definitions.txt'], 'w'); 
%     fprintf(fid,['c*********************************************************' '\n']);
%     fprintf(fid,['c Definition of loading conditions' '\n']);
%     fprintf(fid,['c*********************************************************' '\n']);
%     fprintf(fid,['          H_cycle=' num2str(Nom_strain/Strain_rate)  '\n' '\n']);
%     fprintf(fid,['          ramp_time=' num2str(Nom_strain/Strain_rate)  '\n' '\n']);
%     fprintf(fid,['          hcycle_time=' num2str(Nom_strain/Strain_rate)  '\n' '\n']);
%     fprintf(fid,['          Nmulti_begins=' num2str(Nmulti_begins)  '\n']);
%     fprintf(fid,['          NStage_II=' num2str(Stage_II)  '\n']);
%     fprintf(fid,['          Eval_Nuc= ' num2str(Eval_Nuc) ' \n' '\n']);
%     fprintf(fid,['          Eval_MSC= ' num2str(Eval_MSC) '  \n' '\n']);
%     fprintf(fid,['c*********************************************************' '\n']);
%     fprintf(fid,['c Definition of global mechanical properties' '\n']);
%     fprintf(fid,['c*********************************************************' '\n']);
%     fprintf(fid,['          Sigma_y=' num2str(Sigma_y)  '\n']);
%     fprintf(fid,['          Cont_FS_trans=' num2str(Cont_FS_trans)  '\n' '\n']);
%     fprintf(fid,['c*********************************************************' '\n']);
%     fprintf(fid,['c  Definition of nucleation fatigue life parameters' '\n']);
%     fprintf(fid,['c*********************************************************' '\n']);
%     fprintf(fid,['          alpha_g=' num2str(alpha_g) '  ![cycles*mm] \n' '\n']);
%     fprintf(fid,['c*********************************************************' '\n']);
%     fprintf(fid,['c   Definition of transgranular  fatigue life parameters' '\n']);
%     fprintf(fid,['c*********************************************************' '\n']);
%     fprintf(fid,['          Phi_irr='  num2str(Phi_irr) '            !Irreversibility factor  [cycles^-1]' '\n']);
%     fprintf(fid,['          DCTD_th='  num2str(DCTD_th) '          !Delta CTD threshold  [mm]' '\n']);
%     fprintf(fid,['          A_fs='  num2str(A_fs) '           ! Proportionality factor  [mm]' '\n']);
%     fprintf(fid,['          d_gr_ref=' num2str(d_gr_ref) '             !diameter of coarse grain  [mm]'  '\n' '\n']);
%     fprintf(fid,['          Pg= '  num2str(Pg) '              !controls evolution of FIP as band cracks [unitless]'  '\n' '\n']);
%     fprintf(fid,['c*********************************************************' '\n']);
%     fprintf(fid,['c   Definition of intergranular fatigue life parameters' '\n']);
%     fprintf(fid,['c*********************************************************' '\n']);
%     fprintf(fid,['          Phi_irr_int='  num2str(Phi_irr_int) '            !Irreversibility factor' '\n']);
%     fprintf(fid,['          C_int_0='  num2str(C_int_0) '            !Irreversibility factor' '\n']);
%     fprintf(fid,['          P_beta='  num2str(P_beta) '          !Delta CTD threshold' '\n']);
%     fprintf(fid,['          Q_ox='  num2str(Q_ox) '           ! Proporcionality factor ' '\n']);
%     fprintf(fid,['          B_int='  num2str(B_int) '           ! Proporcionality factor ' '\n']);
%     fprintf(fid,['          Cepsilon_int='  num2str(Cepsilon_int) '           ! Proporcionality factor ' '\n']);
%     fprintf(fid,['          Clambda_int=' num2str(Clambda_int) '             !diameter of coarse grain in microns'  '\n' '\n']);


%     fclose(fid); 

%     fid = fopen(['Common_block_SS304v04.txt'], 'w');
%     fprintf(fid,['         real*8 d_gamma_dot,gamma_dot,gamma_try,g0,' '\n']);
%     fprintf(fid,['     &   g,a0,a,tau,sigma,dam_elem'  '\n']);
%     fprintf(fid,['         real*8 gamma_cum_element_max,gamma_cum_element_min,' '\n']);
%     fprintf(fid,['     &   gamma_cum_init_even_step' '\n']);
%     fprintf(fid,['         real*4 Cnum,Cnum_fail' '\n']);
%     fprintf(fid,['         integer*4 Elem_pos' '\n']);
%     fprintf(fid,['         real*8 sigma_gl,disAngle,d_gr_nd' '\n']);
%     fprintf(fid,['         integer*8 N_msc,N_history,N_msc_min,NGrain_failed' '\n']);
%     fprintf(fid,['         integer*8 Nstep_failed,N_nuc_min,Nfailed,Eval_Nuc' '\n']);
%     fprintf(fid,['         integer*8 N_El_GB,N_int_local,N_int,Num_int' '\n']);
%     fprintf(fid,['         real*8 Grains_info,GB_dir' '\n']);
%     fprintf(fid,['         integer*8 NTemp2,NTemp3' '\n']);
%     fprintf(fid,['         integer*8 Neighbor1,Neighbor_grains,Neighbor_num' '\n']);
%     fprintf(fid,['         integer*8 N_nuc' '\n']);
%     fprintf(fid,['         integer*8 Nelem_fail, Num_failed_el' '\n']);
%     fprintf(fid,['         integer*8 N_history_el,N_history_el_old' '\n']);
%     fprintf(fid,['         integer*8 NStage_II' '\n']);
%     fprintf(fid,['         logical Min_dist,Flg' '\n']);
%     fprintf(fid,['         real*8 CGR_msc,Crack_plane_temp,sigma_crack_plane' '\n']);
%     fprintf(fid,['         real*8 Crack_plane_norm,dir_cos_gl'  '\n' '\n' '\n']);
% 
%     fprintf(fid,['       common/KUEXT0/ Elem_pos(num_elem,5) ' '\n']);
% 
%     fprintf(fid,['       common/KUEXT1/ Cnum(num_grains,max_num_layers,4) ' '\n']);
%     fprintf(fid,['       common/KUEXT1b/ Cnum_fail(num_grain_fail+1) ' ' \n' '\n']);
% 
%     fprintf(fid,['       common/KUEXT2/ gamma_dot_element(num_elem,num_slip_sys)  ' '\n']);
%     fprintf(fid,['       common/KUEXT2b/ gamma_cum_element(num_elem,num_slip_sys)' '\n']);
%     fprintf(fid,['       common/KUEXT2c/ gamma_cum_element_max(num_elem,num_slip_sys)' '\n']);
%     fprintf(fid,['       common/KUEXT2d/ gamma_cum_element_min(num_elem,num_slip_sys)' '\n']);
%     fprintf(fid,['       common/KUEXT2e/ gamma_cum_init_even_step(num_elem,num_slip_sys)' '\n']);
%     fprintf(fid,['       common/KUEXT2f/ sigma_gl(num_elem,num_slip_sys)' '\n' '\n']);
%     fprintf(fid,['       common/KUEXT2g/ dam_elem(num_elem)' '\n' '\n']);
% 
%     fprintf(fid,['       common/KUEXT3/ N_nuc(num_grains,max_num_layers,num_slip_sys)' '\n']);
%     fprintf(fid,['       common/KUEXT3b/ N_msc(num_grains,max_num_layers,num_slip_sys)' '\n']);
%     fprintf(fid,['       common/KUEXT3c/ N_history(num_grains,num_grain_fail)' '\n']);
%     fprintf(fid,['       common/KUEXT3d/ N_msc_min(4,num_grain_fail+1)' '\n']);
%     fprintf(fid,['       common/KUEXT3e/ NGrain_failed(4,num_grain_fail+1)' '\n']);
%     fprintf(fid,['       common/KUEXT3f/ Nstep_failed(num_grain_fail+1)' '\n']);
%     fprintf(fid,['       common/KUEXT3g/ N_nuc_min(4)' '\n']);
%     fprintf(fid,['       common/KUEXT3h/ Nfailed' '\n' '\n' '\n']);
% 
%     fprintf(fid,['       common/KUEXT4/ Neighbor1(num_grains*2*num_grains)' '\n']);
%     fprintf(fid,['       common/KUEXT4b/ Neighbor_grains(num_grains,num_grains)' '\n']);
%     fprintf(fid,['       common/KUEXT4c/ Neighbor_num(num_grains)' '\n' '\n']);
% 
%     fprintf(fid,['       common/KUEXT5/ N_El_GB(num_elem,6)' '\n' ]);
%     fprintf(fid,['       common/KUEXT5b/ N_int(num_grains,num_grains)' '\n']);
%     fprintf(fid,['       common/KUEXT5c/ N_int_local(num_elem)' '\n']);
%     fprintf(fid,['       common/KUEXT5d/ Num_int(num_grains,num_grains)' '\n' '\n']);
% 
%     fprintf(fid,['       common/KUEXT7/ Spk2_gl(num_elem,3,3) ! 2nd Piola-Kirchhoff global var' '\n']);
%     fprintf(fid,['       common/KUEXT7b/ GB_dir(num_elem,6)'  '\n']);
%     fprintf(fid,['       common/KUEXT7c/ S_prin_max(num_elem)'  '\n' '\n']);
% 
%     fprintf(fid,['       common/KUEXT8/' '\n']);
%     fprintf(fid,['     &   Min_dist(num_grains,num_grains,max_num_layers,max_num_layers,4,4)' '\n']);
%     fprintf(fid,['       common/KUEXT8b/ Flg' '\n' '\n']);
% 
%     fprintf(fid,['       common/KUEXT9/ disAngle(num_grains,num_grains)' '\n']);
%     fprintf(fid,['       common/KUEXT9b/ disAngleMax(num_grains,num_grains)' '\n']);
%     fprintf(fid,['       common/KUEXT9c/ Grs(num_grains,4)' '\n'  '\n']);
% 
%     fprintf(fid,['       common/KUEXT10/ d_gr_nd(num_grains,max_num_layers,num_slip_sys) ! This is the grain size contribution from the first neighbor in grain (k)' '\n' '\n']);
% 
%     fprintf(fid,['       common/KUEXT11/ Neighbors_el(num_elem,26)' '\n']);
%     fprintf(fid,['       common/KUEXT11b/ Nelem_fail(100) ! ID number for those elements damaged' '\n']);
%     fprintf(fid,['       common/KUEXT11c/ Num_failed_el ! number of elements cracked' '\n'  '\n']);
% 
%     fprintf(fid,['       common/KUEXT12/ N_history_el(num_elem),N_history_el_old(num_elem) '  '\n']);
%     fprintf(fid,['       common/KUEXT12b/ Crack_plane_temp(3,num_elem) '  '\n']);
%     fprintf(fid,['       common/KUEXT12c/ sigma_crack_plane(num_elem) ' '\n'  '\n']);
%     fprintf(fid,['       common/KUEXT12d/  Crack_plane_norm(3,num_elem)' '\n'  '\n']);
%     fprintf(fid,['       common/KUEXT12e/ dir_cos_gl(3,3,num_elem) ' '\n'  '\n']);
% 
% %    fprintf(fid,['        common/CONSTANTS/' '\n']);
% %    fprintf(fid,['        SAVE /KUEXT0/ , /KUEXT1/, /KUEXT2/, /KUEXT3/, ' '\n']);
% %    fprintf(fid,['     &  /KUEXT4/, /KUEXT5/, /KUEXT5b/, /KUEXT5c/, ' '\n']);
% %    fprintf(fid,['     &  /KUEXT7/ , /KUEXT8/ , /KUEXT9/ , /KUEXT10/' '\n']);
% %    fprintf(fid,['     &  /KUEXT7/ , /KUEXT8/ , /KUEXT9/ , /KUEXT10/' '\n']);
%     fclose(fid);


%clc;   
    n_El % Number of elements
    n_grains  % Number of grains
    Max_num_layers  % Maximum number of layers
%clear; 

toc

return
end
