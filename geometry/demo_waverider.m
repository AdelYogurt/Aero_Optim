clc;
clear;
close all hidden;

%% cone waverider

Ma_in=13.8;
T_in=89.3;
P_in=951.5;
beta=10/180*pi; % oblique shock wave angle

R_0=0.05;
L_total=0.4;
W_total=0.16;

x_0=R_0/tan(beta);
R_total_sq=(tan(beta)*(x_0+L_total))^2;
H_total=(sqrt(R_total_sq-(W_total/2)^2)-R_0);
lead_edge_fcn=@(z) -R_0+H_total*(cos(pi*z/W_total)-1);
line((0:0.1:1)*W_total/2,lead_edge_fcn((0:0.1:1)*W_total/2));
waverider=WaveriderCone('',Ma_in,T_in,P_in,beta,...
    lead_edge_fcn,R_0,L_total,W_total);
U_num=21;
V_num=21;
W_num=21;
waverider.drawBody([],U_num,[],W_num)
% part.mesh_list=waverider.calSurfaceMatrix(U_num,[],W_num);
% writeMeshSTL('waverider_cone',convertWGSToSTL(part));
waverider.writeStepOpenShell('waverider_cone.step',U_num,[],W_num);
% 
% surf_total=waverider.calSurfaceMatrix(U_num,[],W_num);
% surf_total=waverider.reverseUV(surf_total);
% X_up=surf_total{2}.X;Y_up=surf_total{2}.Y;Z_up=surf_total{2}.Z;
% X_low=surf_total{1}.X;Y_low=surf_total{1}.Y;Z_low=surf_total{1}.Z;
% save('wing/cone_waverider.mat','X_up','Y_up','Z_up','X_low','Y_low','Z_low')
