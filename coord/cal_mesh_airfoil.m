clc;
clear;
% close all hidden;

%% read CGNS file

mesh_filestr='NACA0012.cgns';
marker_name_list={'AIRFOIL_UP','AIRFOIL_LOW'};

[point_list,marker_index_list]=readBCCGNS(mesh_filestr,marker_name_list);

save('mesh_data_airfoil.mat','point_list','marker_index_list')

mesh_point=calMesh(point_list,marker_index_list);
save('mesh_data_airfoil.mat','mesh_point','-append')

%% read SU2 file

% mesh_filestr='mesh_RAE2822_turb.su2';
% 
% [point_list,element_list,marker_list]=readMeshSU2(cgns_filestr);
% marker_element=marker_list.('AIRFOIL');
% marker_index=unique(marker_element(2:end,:));
% Bool=point_list(marker_index,2) >= 0;
% marker_index_list.AIRFOIL_UP=marker_index(Bool);
% marker_index_list.AIRFOIL_LOW=marker_index(~Bool);
% 
% save('mesh_data_airfoil.mat','point_list','marker_index_list')

%% pre calculate curve coord

% gen object
control_point_up=importdata('airfoil_NACA0012_shape_up.txt');
control_point_low=importdata('airfoil_NACA0012_shape_low.txt');
airfoil=AirfoilCST(control_point_low,control_point_up);

% load data
load('mesh_data_airfoil.mat','point_list','marker_index_list');
mesh_coord=airfoil.calCoord(point_list,marker_index_list);
save('mesh_data_airfoil.mat','mesh_coord','-append')

% % write coord
% writeCoord(mesh_coord,1,'airfoil');

%% calculate deform point coordinate

% % base parameter
% control_point_up=importdata('airfoil_NACA0012_shape_up.txt');
% control_point_low=importdata('airfoil_NACA0012_shape_low.txt');
% 
% % control_point_up=importdata('airfoil_RAE2822_shape_up.txt');
% % control_point_low=importdata('airfoil_RAE2822_shape_low.txt');
% 
% % gen object
% airfoil=AirfoilCST(control_point_low,control_point_up);
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

