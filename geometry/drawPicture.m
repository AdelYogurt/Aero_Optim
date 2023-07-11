clc;
clear;
close all hidden;

cd ..

addpath([pwd,'\lib'])
addpath([pwd,'\geometry'])

total_length=4; % total length

vari_num=15;
low_bou=[0.7,0.7, ...
    2.4,0.6,0.5,0.5, ...
    0.1,0.1,0.005, ...
    0.25,0.6,0.6,0.5,0.01,0.01];
up_bou=[1.0,1.0, ...
    3.2,0.8,2.0,2.0, ...
    0.5,0.5,0.02, ...
    0.35,0.8,0.8,0.7,0.05,0.05];
% x=rand(1,vari_num).*(up_bou-low_bou)+low_bou;
% x=(up_bou+low_bou)/2;
% x=up_bou;
% x=low_bou;

xi_grid_num_head=50; % head x direction gird num
eta_grid_num_head=20; % head and body y direction gird num
xi_grid_num_body=20; % body x direction gird num
eta_grid_num_wing=6; % wing y direction gird num
edge_gird_num=12; % edge gird num

%% write to INP file
% [par_width,par_hight_up,par_hight_low,...
%     par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
%     par_rho1,par_rho12,par_WS1,par_WS2,par_TWU,par_TWL]=decode(x);
% waverider_wing=WaveriderWingDia...
%     (total_length,par_width,par_hight_up,par_hight_low,...
%     par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
%     par_rho1,par_rho12,par_WS1,par_WS2,par_TWU,par_TWL);

% [X,Y,Z]=waverider_wing.calSurfaceMatrix(xi_grid_num_head,eta_grid_num_head,xi_grid_num_body,eta_grid_num_wing,edge_gird_num);

% part_list=cell(1,3);
% 
% network_list={'head_up','head_low','body_up','body_low','body_back_up','body_back_low'};
% mesh_list=createMesh(X,Y,Z,network_list);
% part.name='main_body';
% part.mesh_list=mesh_list;
% part_list{1}=part;
% 
% network_list={'tri_wing_up','tri_wing_front','tri_wing_back_up','tri_wing_low','tri_wing_back_low'};
% mesh_list=createMesh(X,Y,Z,network_list);
% part.name='tri_wing';
% part.mesh_list=mesh_list;
% part_list{2}=part;
% 
% network_list={'wing_up','wing_front','wing_back_up','wing_low','wing_back_low','wing_side_up','wing_side_low'};
% mesh_list=createMesh(X,Y,Z,network_list);
% part.name='wing';
% part.mesh_list=mesh_list;
% part_list{3}=part;
% 
% [part_1,point_list]=convertWGSToMesh(part_list{1},false(1));
% writeMeshINP([part_1.name,'.inp'],part_1,point_list)
% 
% [part_2,point_list]=convertWGSToMesh(part_list{2},false(1));
% writeMeshINP([part_2.name,'.inp'],part_2,point_list)
% 
% [part_3,point_list]=convertWGSToMesh(part_list{3},false(1));
% writeMeshINP([part_3.name,'.inp'],part_3,point_list)

%% demo cutting

% par_R=0;
% 
% waverider_wing=WaveriderWingDia...
%     (total_length+0.001,par_width,par_hight_up,par_hight_low,...
%     par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R+0.001,...
%     par_rho1,par_rho12,par_WS1,par_WS2,par_TWU,par_TWL);
% 
% waverider_base=WaveriderBase...
%     (total_length,par_width,par_hight_up,par_hight_low,...
%     par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R);
% 
% part=waverider_base.getWGSMesh('waverider_base',xi_grid_num_head,eta_grid_num_head,edge_gird_num);
% [part,point_list]=convertWGSToMesh(part,false(1));
% writeMeshINP('waverider_base',part,point_list);
% 
% [X,Y,Z]=waverider_wing.calSurfaceMatrix(xi_grid_num_head,eta_grid_num_head,xi_grid_num_body,eta_grid_num_wing,edge_gird_num);
% 
% network_list={'body_up','body_low','body_back_up','body_back_low'};
% mesh_list=createMesh(X,Y,Z,network_list);
% part.name='waverider_wing_cut';
% part.mesh_list=mesh_list;
% [part,point_list]=convertWGSToMesh(part,false(1));
% writeMeshINP('waverider_wing_cut',part,point_list);

%% demo multi shape

% X=[0.800,0.850,2.500,0.600,2.000,1.300,0.300,0.200,0.013,0.250,0.800,0.780,0.550,0.040,0.0250
%     0.600,0.800,2.500,0.700,0.600,0.650,0.400,0.450,0.008,0.300,0.650,0.700,0.650,0.030,0.0200
%     0.750,0.900,3.000,0.800,1.000,1.800,0.150,0.300,0.0170,0.350,0.770,0.600,0.550,0.013,0.040];
% 
% x=X(1,:);
% [par_width,par_hight_up,par_hight_low,...
%     par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
%     par_rho1,par_rho12,par_WS1,par_WS2,par_TWU,par_TWL]=decode(x);
% waverider_wing=WaveriderWingDia...
%     (total_length,par_width,par_hight_up,par_hight_low,...
%     par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
%     par_rho1,par_rho12,par_WS1,par_WS2,par_TWU,par_TWL);
% part=waverider_wing.getWGSMesh('waverider_wing_diaI',xi_grid_num_head,eta_grid_num_head,xi_grid_num_body,eta_grid_num_wing,edge_gird_num);
% writeMeshSTL('waverider_wing_diaI',convertWGSToSTL(part));
% 
% x=X(2,:);
% [par_width,par_hight_up,par_hight_low,...
%     par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
%     par_rho1,par_rho12,par_WS1,par_WS2,par_TWU,par_TWL]=decode(x);
% waverider_wing=WaveriderWingDia...
%     (total_length,par_width,par_hight_up,par_hight_low,...
%     par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
%     par_rho1,par_rho12,par_WS1,par_WS2,par_TWU,par_TWL);
% part=waverider_wing.getWGSMesh('waverider_wing_diaII',xi_grid_num_head,eta_grid_num_head,xi_grid_num_body,eta_grid_num_wing,edge_gird_num);
% writeMeshSTL('waverider_wing_diaII',convertWGSToSTL(part));
% 
% x=X(3,:);
% [par_width,par_hight_up,par_hight_low,...
%     par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
%     par_rho1,par_rho12,par_WS1,par_WS2,par_TWU,par_TWL]=decode(x);
% waverider_wing=WaveriderWingDia...
%     (total_length,par_width,par_hight_up,par_hight_low,...
%     par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
%     par_rho1,par_rho12,par_WS1,par_WS2,par_TWU,par_TWL);
% part=waverider_wing.getWGSMesh('waverider_wing_diaIII',xi_grid_num_head,eta_grid_num_head,xi_grid_num_body,eta_grid_num_wing,edge_gird_num);
% writeMeshSTL('waverider_wing_diaIII',convertWGSToSTL(part));
