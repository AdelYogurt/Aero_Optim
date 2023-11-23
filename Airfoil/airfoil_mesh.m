clc;
clear;
close all hidden;

% cd ..
% addpath cgns4m;
% addpath coord;
% addpath mesh;
% addpath geometry;
% addpath model;
% addpath optimizer;
% startup_cgns4m;
% cd Airfoil;

%% read CGNS file

mesh_filestr='SU2/NACA0012.cgns';
marker_name_list={'AIRFOIL_UP','AIRFOIL_LOW'};
[point_list,marker_index_list]=readBCCGNS(mesh_filestr,marker_name_list);
save('mesh_data_airfoil.mat','point_list','marker_index_list')
mesh_point=calMesh(point_list,marker_index_list);
save('mesh_data_airfoil.mat','mesh_point','-append')

%% pre calculate curve coord

% gen object
control_point_up=importdata('geom/NACA0012_CSTshape_up.txt');
control_point_low=importdata('geom/NACA0012_CSTshape_low.txt');
airfoil=AirfoilGeom(control_point_low,control_point_up);

% load data
load('mesh_data_airfoil.mat','point_list','marker_index_list');
mesh_coord=airfoil.calCoord(point_list,marker_index_list);
save('mesh_data_airfoil.mat','mesh_coord','-append')

% % write coord
% writeCoord(mesh_coord,1,'airfoil');

%% calculate deform point coordinate

% % base parameter
% control_point_up=importdata('geom/Clark_Y_CSTshape_up.txt');
% control_point_low=importdata('geom/Clark_Y_CSTshape_low.txt');
% 
% control_point_up=importdata('geom/NACA0012_CSTshape_up.txt');
% control_point_low=importdata('geom/NACA0012_CSTshape_low.txt');
% 
% control_point_up=importdata('geom/NACA4412_CSTshape_up.txt');
% control_point_low=importdata('geom/NACA4412_CSTshape_low.txt');
% 
% control_point_up=importdata('geom/RAE2822_CSTshape_up.txt');
% control_point_low=importdata('geom/RAE2822_CSTshape_low.txt');
% 
% % gen object
% airfoil=AirfoilGeom(control_point_low,control_point_up);
% 
% % load data and calculate deform mesh point
% load('mesh_data_airfoil.mat','mesh_coord')
% mesh_point=airfoil.calMeshPoint(mesh_coord);
% save('mesh_data_airfoil.mat','mesh_point','-append')
% 
% % % write mesh point
% % writePoint(mesh_point,'airfoil_deform.dat',2)
% 
% % draw mesh point
% drawPoint(mesh_point,2)

