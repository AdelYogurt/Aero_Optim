clc;
clear;
close all hidden;

%% base parameter

total_length=4; % total length

vari_num=16;
low_bou=[0.7,0.7, ...
    2.4,0.4,0.5,0.5, ...
    0.1,0.1,0.005, ...
    0.4,0.4,0.4,0.6,0.5,0.01,0.01];
up_bou=[1.0,1.0, ...
    3.2,0.8,2.0,2.0, ...
    0.5,0.5,0.02, ...
    0.6,0.6,0.6,0.8,0.7,0.05,0.05];

% x=rand(1,vari_num).*(up_bou-low_bou)+low_bou;
x=(up_bou+low_bou)/2;
% x=up_bou;
% x=low_bou;

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

% load('optres.mat')
% x=x_MF_best;

[par_width,par_hight_up,par_hight_low,...
    par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
    par_rho1,par_rho12,par_rho23,par_WS1,par_WS2,par_TWU,par_TWL]=decode(x);

%% create waverider

% waverider_wing=WaveriderWingFei...
%     (total_length,par_W,par_H_up,par_H_low,...
%     par_T,par_MS_up,par_MS_low,par_NB_up,par_NB_low,par_R,...
%     par_rho1,par_rho12,par_WS1,par_WS2);

% waverider_wing=WaveriderWingTri...
%     (total_length,par_W,par_H_up,par_H_low,...
%     par_T,par_MS_up,par_MS_low,par_NB_up,par_NB_low,par_R,...
%     par_rho1,par_rho12,par_WS1,par_WS2,par_TWU,par_TWL);

waverider_wing=WaveriderWingDia...
    (total_length,par_width,par_hight_up,par_hight_low,...
    par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
    par_rho1,par_rho12,par_rho23,par_WS1,par_WS2,par_TWU,par_TWL);

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
waverider_wing.writeStepShell('waverider_wing_dia');

%% function

function [par_width,par_hight_up,par_hight_low,...
    par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
    par_rho1,par_rho12,par_rho23,par_WS1,par_WS2,par_TWU,par_TWL]=decode(x)
% length parameter
par_M_up=x(1);
par_M_low=x(2);
% width parameter
par_width=x(3);
par_T=x(4);
par_N_up=x(5);
par_N_low=x(6);
% height parameter
par_hight_up=x(7);
par_hight_low=x(8);
par_R=x(9);
% wing parameter
par_rho1=x(10);
par_rho12=x(11);
par_rho23=x(12);
par_WS1=x(13);
par_WS2=x(14);
par_TWU=x(15);
par_TWL=x(16);

end

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
