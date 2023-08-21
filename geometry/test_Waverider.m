clc;
clear;
close all hidden;

%% base parameter

total_length=4; % total length

% vari_num=16;
% low_bou=[0.7,0.7, ...
%     2.4,0.4,0.5,0.5, ...
%     0.1,0.1,0.01, ...
%     0.4,0.4,0.4,0.6,0.5,0.01,0.01];
% up_bou=[1.0,1.0, ...
%     3.2,0.8,2.0,2.0, ...
%     0.5,0.5,0.04, ...
%     0.6,0.6,0.6,0.8,0.7,0.05,0.05];

% x=rand(1,vari_num).*(up_bou-low_bou)+low_bou;
% x=(up_bou+low_bou)/2;
% x=[0.85,0.85,2.8,0.7,1.25,1.25,0.3,0.3,0.0125,0.3,0.7,1,0.7,0.6,0.03,0.03];
% x=[0.778,0.988,3.112,0.6,1.538,1.720,0.409,0.225,0.0125,0.339,0.6,1,0.661,0.642,0.01,0.01];

% x=[0.752,0.973,2*1.634,1/(0.805*2),0.796,1.948,0.446,0.165,0.005,1-0.793,0.773,1,0.778,0.409,0.04,0.01]; % optimal result of WaveriderWingFei in fei thesis
x=[0.750,0.750,2*1.400,1/(1.000*2),1.250,1.250,0.300,0.300,0.005,1-0.700,0.700,1,0.750,0.600,0.04,0.01]; % origin shape of WaveriderWingFei in fei thesis

% x=up_bou;
% x=low_bou;

% load('optres.mat')
% x=x_MF_best;

%% create waverider

% waverider_wing=WaveriderWing(total_length,x,'Fei');
% waverider_wing=WaveriderWing(total_length,x,'Dia');
% waverider_wing=WaveriderWing(total_length,x,'Tri');

% waverider_base=WaveriderBase...
%     (total_length,par_width,par_hight_up,par_hight_low,...
%     par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R);

%% draw

% fine grid
% xi_grid_num_head=40; % head x direction gird num
% eta_grid_num_head=20; % head and body y direction gird num
% xi_grid_num_body=20; % body x direction gird num
% eta_grid_num_wing=8; % wing y direction gird num
% edge_gird_num=6; % edge gird num

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

% part=waverider_wing.getWGSMesh('waverider_wing_dia');
% writeMeshSTL('waverider_wing_dia',convertWGSToSTL(part));

% part=waverider_base.getWGSMesh('waverider_base',xi_grid_num_head,eta_grid_num_head,edge_gird_num);
% writeMeshSTL('waverider_base',convertWGSToSTL(part));

%% write to INP

% part=waverider_wing.getWGSMesh('waverider_wing_dia');
% [part_list,point_list]=convertWGSToMesh(part,false(1));
% writeMeshINP('waverider_wing_dia',part_list,point_list);

%% write to STEP
% waverider_wing.writeStepOpenShell('waverider_wing_dia',1e-4);

%% cone waverider

% Ma_in=13.8;
% T_in=89.3;
% P_in=951.5;
% beta=10/180*pi; % oblique shock wave angle
% 
% R_0=0.05;
% L_total=0.4;
% W_total=0.16;
% 
% x_0=R_0/tan(beta);
% R_total_sq=(tan(beta)*(x_0+L_total))^2;
% H_total=(sqrt(R_total_sq-(W_total/2)^2)-R_0);
% lead_edge_fcn=@(z) -R_0+H_total*(cos(pi*z/W_total)-1);
% 
% % drawFunction(lead_edge_fcn,0,W_total/2);
% 
% waverider=WaveriderCone(Ma_in,T_in,P_in,beta,...
%     lead_edge_fcn,R_0,L_total,W_total);
% 
% U_num=11;
% V_num=11;
% W_num=11;
% 
% waverider.drawBody(U_num,[],W_num)

% surf_total=waverider.calSurfaceMatrix(U_num,[],W_num);
% part.mesh_list=surf_total;
% writeMeshSTL('waverider_cone',convertWGSToSTL(part));

% waverider.writeStepOpenShell('waverider_cone',U_num,[],W_num);

%%

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
