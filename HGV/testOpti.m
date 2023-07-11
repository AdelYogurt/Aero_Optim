clc;
clear;
close all hidden;

addpath([pwd,'\lib'])
addpath([pwd,'\geometry'])

total_length=4; % total length

% length parameter
par_MS_up=0.752; % 0.3-3.0 1
par_MS_low=0.973; % 0.3-3.0 1
% width parameter
par_W=2*1.634; % 0.5-1.5 1 *total_width
par_T=1/(0.805*2); % 0.3-0.9 0.6
par_NB_up=0.796; % 0.3-3.0 1
par_NB_low=1.948; % 0.3-3.0 1
% height parameter
par_H_up=0.446; % 0.05-0.2 0.1 *total_height
par_H_low=0.165; % 0.05-0.2 0.1 *total_height
par_R=0.005; % 0.005-0.02 0.01
% wing parameter
par_rho1=1-0.793; % 0.2-0.6 0.4
par_rho12=0.773; % 0.4-1 0.7
par_WS1=0.778; % 0.5-1 0.7
par_WS2=0.409; % 0.5-1 0.7
par_TWU=0.01; % 0.05-0.2 0.1
par_TWL=0.02; % 0.05-0.2 0.1

xi_grid_num_head=30; % head x direction gird num
eta_grid_num_head=20; % head and body y direction gird num
xi_grid_num_body=10; % body x direction gird num
eta_grid_num_wing=10; % wing y direction gird num
edge_gird_num=ceil((par_R/0.0015)/2)*2; % edge gird num

waverider_fei=WaveriderWingFei...
    (total_length,par_W,par_H_up,par_H_low,...
    par_T,par_MS_up,par_MS_low,par_NB_up,par_NB_low,par_R,...
    par_rho1,par_rho12,par_WS1,par_WS2);

[X,Y,Z]=waverider_fei.calSurfaceMatrix(xi_grid_num_head,eta_grid_num_head,xi_grid_num_body,eta_grid_num_wing,edge_gird_num);

% waverider_fei.drawBody(xi_grid_num_head,eta_grid_num_head,xi_grid_num_body,eta_grid_num_wing,edge_gird_num);

part=waverider_fei.getWGSPart('waverider_wing_fei',xi_grid_num_head,eta_grid_num_head,xi_grid_num_body,eta_grid_num_wing,edge_gird_num);
writeMeshSTL('waverider_wing_fei',convertWGSToSTL(part));
