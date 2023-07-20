clc;
clear;
close all hidden;

%% read CGNS file

cgns_filestr='..\Airfoil\SU2\NACA0012.cgns';
marker_name_list={'AIRFOIL_UP','AIRFOIL_LOW'};

[point_list,marker_index_list]=readCGNSMesh(cgns_filestr,marker_name_list);

save('mesh_data_airfoil.mat','point_list','marker_index_list')

%% pre calculate curve coord

% gen object
control_point_up=importdata('airfoil_NACA0012_shape_up.txt');
control_point_low=importdata('airfoil_NACA0012_shape_low.txt');
airfoil=Airfoil2DCST(control_point_up,control_point_low);

% load data
load('mesh_data_airfoil.mat','point_list','marker_index_list');
airfoil_coord=airfoil.calCoord(point_list,marker_index_list);
save('mesh_data_airfoil.mat','airfoil_coord','-append')

% % write coord
% writeCoord(airfoil_coord,1);

%% calculate deform point coordinate

% base parameter
% control_point_up=importdata('airfoil_NACA0012_shape_up.txt');
% control_point_low=importdata('airfoil_NACA0012_shape_low.txt');

control_point_up=importdata('airfoil_RAE2822_shape_up.txt');
control_point_low=importdata('airfoil_RAE2822_shape_low.txt');

% optimal parmeter
mid=[control_point_up(2:end-1,2);control_point_low(2:end-1,2)]';

variable_num=10;
low_bou=mid-0.05;
up_bou=mid+0.05;
% x=rand(1,variable_num).*(up_bou-low_bou)+low_bou;
x=(up_bou+low_bou)/2;
% x=up_bou;
% x=low_bou;
% x=[ones(1,5)*0,ones(1,5)*0.1];

% gen object
[control_point_up,control_point_low]=Airfoil2DCST.decode(control_point_up,control_point_low,x);
airfoil=Airfoil2DCST(control_point_up,control_point_low);

% load data and calculate deform mesh point
load('mesh_data_airfoil.mat','airfoil_coord')
mesh_point=airfoil.calMeshPoint(airfoil_coord);
save('mesh_data_airfoil.mat','mesh_point','-append')

% % write mesh point
% writePoint(mesh_point,'airfoil_deform.dat',2)

% draw mesh point
drawPoint(mesh_point,2)


