clc;
clear;
close all hidden;

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
% line((0:0.1:1)*W_total/2,lead_edge_fcn((0:0.1:1)*W_total/2));
% waverider=WaveriderCone('',Ma_in,T_in,P_in,beta,...
%     lead_edge_fcn,R_0,L_total,W_total);
% U_num=21;
% V_num=21;
% W_num=21;
% waverider.drawBody([],U_num,[],W_num);axis equal;
% waverider.writeStepOpenShell('waverider_cone.step',U_num,[],W_num);
% 
% surf_total=waverider.calSurfaceMatrix(U_num,[],W_num);
% surf_total=waverider.reverseUV(surf_total);
% X_up=surf_total{2}.X;Y_up=surf_total{2}.Y;Z_up=surf_total{2}.Z;
% X_low=surf_total{1}.X;Y_low=surf_total{1}.Y;Z_low=surf_total{1}.Z;
% save('wing/cone_waverider.mat','X_up','Y_up','Z_up','X_low','Y_low','Z_low')

%% Fit cone waverider

% load('wing/cone_waverider.mat');
% X_up=X_up-X_up(1,1);Z_up=Z_up-Z_up(1,1);
% X_low=X_low-X_low(end,1);Z_low=Z_low-Z_low(end,1);
% X_low=flipud(X_low);Y_low=flipud(Y_low);Z_low=flipud(Z_low);
% X_LE=X_up(end,:);Y_LE=Y_up(end,:);Z_LE=Z_up(end,:);
% LX_LE=X_LE(end)-X_LE(1);LY_LE=Y_LE(end)-Y_LE(1);LZ_LE=Z_LE(end)-Z_LE(1);
% U_LE=(X_LE-X_LE(1))/LX_LE;V_LE=(Y_LE-Y_LE(1))/LY_LE;W_LE=(Z_LE-Z_LE(1))/LZ_LE;
% par_F=log(U_LE(2:end-1)')\log(W_LE(2:end-1)');par_G=LZ_LE; % fit lead edge
% deform_Z=@(U,V) U.^par_F*par_G;
% par_T=log(U_LE(2:end-1)')\log(V_LE(2:end-1)'); % fit lead edge
% 
% u_degree=2;v_degree=3;u_ctrl_num=3;v_ctrl_num=4;
% C_par_X=[0,0];C_par_Y=[par_T,0];C_par_ZV=[0.5,0.5];C_par_ZU=[1,0];sym_x=false;sym_y=true;
% surf_up=SurfaceCST('',C_par_X,C_par_Y,C_par_ZV,C_par_ZU,sym_x,sym_y);
% surf_low=SurfaceCST('',C_par_X,C_par_Y,C_par_ZV,C_par_ZU,sym_x,sym_y);
% 
% surf_up.addDeform([],[],deform_Z)
% surf_low.addDeform([],[],deform_Z)
% 
% U=(X_up(1,:)-min(X_up(1,:)))/(max(X_up(1,:))-min(X_up(1,:)));
% V=(Y_up(:,end)-min(Y_up(:,end)))/(max(Y_up(:,end))-min(Y_up(:,end)));
% 
% surf_up.addShapeBSpline(X_up,Y_up,Z_up,true,u_degree,v_degree,[],[],[],[],u_ctrl_num,v_ctrl_num,U,V);
% surf_low.addShapeBSpline(X_low,Y_low,Z_low,true,u_degree,v_degree,[],[],[],[],u_ctrl_num,v_ctrl_num,U,V);
% 
% axe_hdl=axes(figure());
% surf_up.drawSurface(axe_hdl,100,10);
% surf_low.drawSurface(axe_hdl,100,10);
% axis equal
