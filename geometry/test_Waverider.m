clc;
clear;
close all hidden;

%% base parameter

total_length=4; % total length

vari_num=16;
low_bou=[0.7,0.7, ...
    2.4,0.4,0.5,0.5, ...
    0.1,0.1,0.01, ...
    0.4,0.4,0.4,0.6,0.5,0.01,0.01];
up_bou=[1.0,1.0, ...
    3.2,0.8,2.0,2.0, ...
    0.5,0.5,0.04, ...
    0.6,0.6,0.6,0.8,0.7,0.05,0.05];

% x=rand(1,vari_num).*(up_bou-low_bou)+low_bou;
% x=(up_bou+low_bou)/2;
% x=[0.85,0.85,2.8,0.7,1.25,1.25,0.3,0.3,0.0125,0.3,0.7,1,0.7,0.6,0.03,0.03];
x=[0.778,0.988,3.112,0.6,1.538,1.720,0.409,0.225,0.0125,0.339,0.6,1,0.661,0.642,0.01,0.01];
% x=up_bou;
% x=low_bou;

% load('optres.mat')
% x=x_MF_best;

%% create waverider

% waverider_wing=WaveriderWingFei...
%     (total_length,par_W,par_H_up,par_H_low,...
%     par_T,par_MS_up,par_MS_low,par_NB_up,par_NB_low,par_R,...
%     par_rho1,par_rho12,par_WS1,par_WS2);

% waverider_wing=WaveriderWingTri...
%     (total_length,par_W,par_H_up,par_H_low,...
%     par_T,par_MS_up,par_MS_low,par_NB_up,par_NB_low,par_R,...
%     par_rho1,par_rho12,par_WS1,par_WS2,par_TWU,par_TWL);

waverider_wing=WaveriderWingDia(total_length,x);

% waverider_base=WaveriderBase...
%     (total_length,par_width,par_hight_up,par_hight_low,...
%     par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R);

%% draw

% fine grid
xi_grid_num_head=40; % head x direction gird num
eta_grid_num_head=20; % head and body y direction gird num
xi_grid_num_body=20; % body x direction gird num
eta_grid_num_wing=8; % wing y direction gird num
edge_gird_num=6; % edge gird num

% % rough grid
% xi_grid_num_head=20; % head x direction gird num
% eta_grid_num_head=10; % head and body y direction gird num
% xi_grid_num_body=10; % body x direction gird num
% eta_grid_num_wing=4; % wing y direction gird num
% edge_gird_num=3; % edge gird num

% waverider_wing.drawBody(xi_grid_num_head,eta_grid_num_head,xi_grid_num_body,eta_grid_num_wing,edge_gird_num);
waverider_wing.drawBody()
% waverider_base.drawBody(xi_grid_num_head,eta_grid_num_head,edge_gird_num);

%% write to wgs

% part=waverider_wing.getWGSMesh('waverider_wing_dia',xi_grid_num_head,eta_grid_num_head,xi_grid_num_body,eta_grid_num_wing,edge_gird_num);
% writeMeshWGS('waverider_wing_dia',part);

%% write to STL

% part=waverider_wing.getWGSMesh('waverider_wing_dia',xi_grid_num_head,eta_grid_num_head,xi_grid_num_body,eta_grid_num_wing,edge_gird_num);
% writeMeshSTL('waverider_wing_dia',convertWGSToSTL(part));

% part=waverider_base.getWGSMesh('waverider_base',xi_grid_num_head,eta_grid_num_head,edge_gird_num);
% writeMeshSTL('waverider_base',convertWGSToSTL(part));

%% write to STEP
% waverider_wing.writeStepOpenShell('waverider_wing_dia',1e-4);

%% function

function mesh_list=createMesh(X_total,Y_total,Z_total,network_list)
mesh_list=cell(1,length(network_list));
for network_index=1:length(network_list)
    network_name=network_list{network_index};
    mesh.X=X_total.(network_name);
    mesh.Y=Y_total.(network_name);
    mesh.Z=Z_total.(network_name);
    mesh.element_type='wgs';
    mesh_list{network_index}=mesh;
end
end

% % calculate ratio
% y_R=(par_R/total_length)^par_T*par_W*total_width;
% yd_R=(par_R/total_length)^(par_T-1)*par_W*total_width/total_length;
% 
% a=fsolve(@(a) (a-par_R)^2/a/a-y_R*(a-par_R)/a/a/yd_R,2,...
%     optimoptions('fsolve','Display','none'));
% b_sq=y_R*y_R/((a-par_R)^2/a/a);
% Ro=b_sq^2/a;
% 
% if Ro > par_R
%     D_x=par_R;
%     D_z=Ro;
% else
%     D_x=Ro;
%     D_z=par_R;
% end
