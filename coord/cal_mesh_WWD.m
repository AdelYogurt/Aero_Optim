clc;
clear;
close all hidden;

% optimal parmeter
total_length=4; % total length
low_bou=[0.7,0.7, ...
    2.4,0.4,0.5,0.5, ...
    0.1,0.1,0.01, ...
    0.4,0.4,0.4,0.6,0.5,0.01,0.01];
up_bou=[1.0,1.0, ...
    3.2,0.8,2.0,2.0, ...
    0.5,0.5,0.04, ...
    0.6,0.6,0.6,0.8,0.7,0.05,0.05];

%% read CGNS file

mesh_filestr='WWD.cgns';

marker_name_list={'WAVERIDER','wing_front','tri_wing_front','head_side',...
    'WING_SIDE','WING','TRI_WING','symmetry'};

[point_list,marker_index_list]=readBCCGNS(mesh_filestr,marker_name_list);

save('mesh_data_WWD.mat','point_list','marker_index_list')

%% pre calculate waverider coord

% gen object
x=[0.85,0.85,2.8,0.6,1.25,1.25,0.3,0.3,0.025,0.5,0.5,0.5,0.7,0.6,0.03,0.03];
% x=[0.85,0.85,2.8,0.6,1.25,1.25,0.3,0.3,0.0125,0.5,0.5,0.5,0.7,0.6,0.03,0.03];
% x=[0.85,0.85,2.8,0.7,1.25,1.25,0.3,0.3,0.0125,0.3,0.7,1,0.7,0.6,0.03,0.03];

waverider_wing=WaveriderWingDia(total_length,x);
[par_width,par_hight_up,par_hight_low,...
    par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
    par_rho1,par_rho12,par_rho23,par_WS1,par_WS2,par_WT_up,par_WT_low]=WaveriderWingDia.decode(x);
y_cut=waverider_wing.y_cut;
head_length=waverider_wing.head_length;

% load data
load('mesh_data_WWD.mat','point_list','marker_index_list');

% divide marker into surface
waverider_index=marker_index_list.WAVERIDER;
tri_wing_index=marker_index_list.TRI_WING;
wing_index=marker_index_list.WING;
wing_side_index=marker_index_list.WING_SIDE;

waverider_point=point_list(waverider_index,:);
tri_wing_point=point_list(tri_wing_index,:);
wing_point=point_list(wing_index,:);
wing_side_point=point_list(wing_side_index,:);

% analysis point
[wing_side_index,marker_index_list.('wing_side_front'),wing_side_point]=divideFB(wing_side_index,wing_side_point,4-4*par_rho1*par_rho12*par_rho23-1e-9);
[marker_index_list.('wing_side_up'),marker_index_list.('wing_side_low')]=divideUL(wing_side_index,wing_side_point);

[tri_wing_back_index,tri_wing_index,tri_wing_back_point,tri_wing_point]=divideFB(tri_wing_index,tri_wing_point,4-1e-9);
[marker_index_list.('tri_wing_up'),marker_index_list.('tri_wing_low')]=divideUL(tri_wing_index,tri_wing_point);
[marker_index_list.('tri_wing_back_up'),marker_index_list.('tri_wing_back_low')]=divideUL(tri_wing_back_index,tri_wing_back_point);

[wing_back_index,wing_index,wing_back_point,wing_point]=divideFB(wing_index,wing_point,4-1e-9);
[marker_index_list.('wing_up'),marker_index_list.('wing_low')]=divideUL(wing_index,wing_point);
[marker_index_list.('wing_back_up'),marker_index_list.('wing_back_low')]=divideUL(wing_back_index,wing_back_point);

[waverider_back_index,head_index,waverider_back_point,head_point]=divideFB(waverider_index,waverider_point,head_length+1e-9);
[body_back_index,body_index,body_back_point,body_point]=divideFB(waverider_back_index,waverider_back_point,4-1e-9);
[marker_index_list.('head_up'),marker_index_list.('head_low')]=divideUL(head_index,head_point);
[marker_index_list.('body_up'),marker_index_list.('body_low')]=divideUL(body_index,body_point);
[marker_index_list.('body_back_up'),marker_index_list.('body_back_low')]=divideUL(body_back_index,body_back_point);

mesh_coord=waverider_wing.calCoord(point_list,marker_index_list);
save('mesh_data_WWD.mat','mesh_coord','marker_index_list','-append')

% % write coord
% writeCoord(mesh_coord,2,'WWD');

%% calculate deform point coordinate

x=[0.85,0.85,2.8,0.6,1.25,1.25,0.3,0.3,0.0125,0.5,0.5,0.5,0.7,0.6,0.03,0.03];
% x=[0.85,0.85,2.8,0.6,1.25,1.25,0.3,0.3,0.0125,0.5,0.5,0.5,0.7,0.6,0.03,0.03];
% x=[0.85,0.85,2.8,0.7,1.25,1.25,0.3,0.3,0.0125,0.3,0.7,1,0.7,0.6,0.03,0.03];

% gen object
waverider_wing=WaveriderWingDia(total_length,x);

% load data and calculate deform mesh point
load('mesh_data_WWD.mat','mesh_coord')
mesh_data=waverider_wing.calMeshPoint(mesh_coord);
save('mesh_data_WWD.mat','mesh_data','-append')

% % write mesh point
% writePoint(mesh_point,'WWD_deform.dat',3)

% draw mesh point
drawPoint(mesh_data,3)

% stag_deform_Y_up=@(Y,Z) (min((Y/(3*par_R)).^2+((Z-par_R)/(par_R)).^2,1)).^(1-par_T/0.7);
% stag_deform_Y_low=@(Y,Z) (min((Y/(3*par_R)).^2+((Z+par_R)/(par_R)).^2,1)).^(1-par_T/0.7);
% stag_deform_Y_mid=@(Y,Z) (min((Y/(3*par_R)).^2+((Z)/(par_R)).^2,1)).^(1-par_T*1/0.7);

