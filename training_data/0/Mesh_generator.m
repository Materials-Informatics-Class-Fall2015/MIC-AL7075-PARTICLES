% Author: Bill Musinski in "Novel methods for microstructure-sensitive
% probabilistic fatigue notch factor", Ms Thesis Georgia Institute of technology, 2010.
% Adapted and modified by: Gustavo M. Castelluccio 2010-2011, Conor Hennessey 2013-2015
% Date: 9/18/09
% Driver to generate a mesh with grain with periodic boundary conditions


function[]=Mesh_generator(options_filename)
%clear all

%Determine number of this microstructural instantiation
s =  options_filename;
MS_try = regexp(s, '(?<=_)[0-9]+','match');
if ~isempty(MS_try)
    [num , status] = str2num(MS_try{1});
    if status ==1
        MS_number = num; %Instantiation number
    else
        MS_number = 1; 
    end
else
    MS_number = 1;
end

%*********************************************************
% Read additional parameters into struct from file
%*********************************************************
options = struct;
if(nargin > 0 && ~strcmp(options_filename,''))
    o = tdfread(options_filename);
    for i=1:size(o.option, 1)
       options.(strtrim(o.option(i,:))) = strtrim(o.value(i,:));
    end
end
options.MS_number = MS_number;


%Seed the 'random' number generator:
temp = clock;
temp2 = num2str(temp(1,6)*1000000); %takes second value and x by 1000000 for Unix based OS
seedVal = str2num(strcat(temp2(5:end),num2str(MS_number))); %append MS_instantation to last digits of seconds
%seedVal = 1
rand('twister',seedVal);

%**********************************************************************************

% Definition of the parameters to create the mesh

% Name of the input file
Input_file_name=strcat('main_',num2str(MS_number));
n_header_lines=9; % Number of lines in the header of the input file. This can change with abaqus version. For abaqus 6.9 it is 9, for abaqus 6.7 it is 8

%*********************************************************
% Definition of mesh attibute
%*********************************************************
scaling=1;
n_grain_fail=20;
geom.x=str2double(options.x)*scaling; % Dimensions of the cube to model in microns along the X axis
geom.y=str2double(options.y)*scaling; % Dimensions of the cube to model in microns along the Y axis
geom.z=str2double(options.z)*scaling; % Dimensions of the cube to model in microns along the Z axis
refine=95*scaling;
mesh.CP_mesh_size=str2double(options.res)*scaling;% Size to seed the part
n_elem_fail=uint8(0.3*ceil((geom.x/mesh.CP_mesh_size)*(geom.y/mesh.CP_mesh_size)*(geom.z/mesh.CP_mesh_size)));
SBperiod=str2double(options.band_thickness)*mesh.CP_mesh_size;% Bands thickness
band_jump=str2double(options.res)*scaling;% Minimum number of elements to consider bands as connected

%*********************************************************
% lognormal Distribution of grains to be created
%*********************************************************
d_gr=7.0/1000;  %RADIUS of grains in mm!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mu_fit=-0.1;% Mean grain size distribution pdf parameters
sigma_fit=0.4;% Grain size distribution pdf parameters
Big_grain1=0; % DIAMETER of grain at the center.
Big_grain2=0; % Size of the grain to the bottom
Big_grain3=0; % Size of the grain to the right.

%*********************************************************
% Loading
%*********************************************************
Nom_strain = 1.4; % Nominal peak strain applied in percentage
Nom_strain_min = 0.0; % Nominal minimum strain applied in percentage
Strain_rate=0.05; % Strain rate applied
Mantle_slip_sys=5;
Core_slip_sys=2;
t1=Nom_strain/Strain_rate;
t2=0;
t3=Nom_strain/Strain_rate;
t4=t2;
%Nmulti_begins=(n_grain_fail+1);
Nmulti_begins=0;
Stage_II=1;
if t2==0
Eval_Nuc=2;
Eval_MSC=6;
else
Eval_Nuc=9;
Eval_MSC=8;
end

%*********************************************************
% Definition of global mechanical properties
%*********************************************************
Sigma_y=517.;  % Cyclic yield strength [MPa]

%*********************************************************
% Definition of nucleation fatigue life parameters
%*********************************************************
alpha_g= 2.9e-7; %[cycles*mm]
Cont_FS_trans=0.5;  % FIP_trans constant multiplying the normal stress for transgranular failure

%*********************************************************
% Definition of transgranular fatigue life parameters
%*********************************************************

Phi_irr=0.35;    %Irreversibility factor [cycles^-1]
DCTD_th=2.86e-7;        %Delta CTD threshold [mm]
A_fs=.0331;           %Proportionality factor [mm]
d_gr_ref=.014; % !diameter of coarse grain [mm]
Pg=0.5;
%*********************************************************
% Definition of intergranular fatigue life parameters
%*********************************************************
irr_exp_int=0.05;
C_int_0=0.002077;
Phi_irr_int=(0.9*t1+t2+0.1*t3+0*t4)^(0.5-irr_exp_int);    %Irreversibility factor
P_beta=1.;
Q_ox=241;
B_int=20;
Cepsilon_int=0.1;
Clambda_int=2;






%**********************************************************************************

% Run the generator
Main_program_preprocess(Input_file_name,n_header_lines,geom,mesh,band_jump,Big_grain1,Big_grain2,Big_grain3,...
d_gr,sigma_fit,mu_fit,SBperiod,Nom_strain,Nom_strain_min,Strain_rate,Sigma_y,alpha_g,Cont_FS_trans,...
Phi_irr,DCTD_th,A_fs,refine,scaling,n_grain_fail,d_gr_ref,Eval_Nuc,Eval_MSC,...
Pg,Phi_irr_int,C_int_0,P_beta,Q_ox,B_int,Cepsilon_int,Clambda_int,t2,Nmulti_begins,...
Mantle_slip_sys,Core_slip_sys,Stage_II,n_elem_fail,options)
%**********************************************************************************


quit
return
end

