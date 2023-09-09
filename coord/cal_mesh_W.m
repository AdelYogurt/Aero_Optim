clc;
clear;
close all hidden;

% optimal parmeter
total_length=4; % total length
low_bou=[0.7,0.7, 2.4,0.4,0.5,0.5, 0.1,0.1,0.01, -0.5,0.5];
up_bou=[1.0,1.0, 3.2,0.8,2.0,2.0, 0.5,0.5,0.04, 0.5,2.0];

%% read CGNS file

mesh_filestr='W.cgns';
marker_name_list={'head_up','head_low','head_side','back'};
[point_list,marker_index_list]=readBCCGNS(mesh_filestr,marker_name_list);
save('mesh_data_W.mat','point_list','marker_index_list')

%% pre calculate waverider coord

% gen object
x=[0.8500,0.8500, 2.8000,0.6000,1.2500,1.2500, 0.3000,0.3000,0.0250, 0,1.2500];
waverider=Waverider(total_length,x,'Dia');

% load data
load('mesh_data_W.mat','point_list','marker_index_list');

% divide marker into surface
back_index=marker_index_list.back;
back_point=point_list(back_index,:);

% analysis point
[marker_index_list.('back_up'),marker_index_list.('back_low')]=divideUL(back_index,back_point);
marker_index_list=rmfield(marker_index_list,'back');

mesh_coord=waverider.calCoord(point_list,marker_index_list);
save('mesh_data_W.mat','mesh_coord','marker_index_list','-append')

% % write coord
% writeCoord(mesh_coord,2,'WWD');

%% calculate deform point coordinate

x=[0.8500,0.8500, 2.8000,0.6000,1.2500,1.2500, 0.3000,0.3000,0.0250, 0,1.2500];

% gen object
waverider=Waverider(total_length,x);

% load data and calculate deform mesh point
% load('mesh_data_W.mat','mesh_coord')
mesh_data=waverider.calMeshPoint(mesh_coord);
save('mesh_data_W.mat','mesh_data','-append')

% % write mesh point
% writePoint(mesh_point,'WWD_deform.dat',3)

% draw mesh point
drawPoint(mesh_data,3)
