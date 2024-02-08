clc;
clear;
close all hidden;

%% Edge

% Degree=3;group_num=3;
% Points=[0,0,0;(0.5:0.5:group_num)',repmat([1;0],group_num,1),repmat([1;0],group_num,1)];
% crv=Curve(Points,Degree);
% edg=Edge('',crv);
% 
% axe_hdl=axes(figure());
% edg.displayModel(axe_hdl);
% axis equal;

%% Wire

% Degree=3;group_num=3;
% Points_1=[0,0,0;(0.5:0.5:group_num)',repmat([1;0],group_num,1),repmat([1;0],group_num,1)];
% Points_2=[0,0,0;repmat([-1;0],group_num,1),(0.5:0.5:group_num)',repmat([1;0],group_num,1)];
% edg_1=Edge('',Curve(Points_1,Degree));
% edg_2=Edge('',Curve(Points_2,Degree));
% wir=Wire('',[edg_1;edg_2]);
% 
% axe_hdl=axes(figure());
% wir.displayModel(axe_hdl);
% axis equal;

%% Face

% UDegree=3;VDegree=3;point_num=5;
% [point_X,point_Y]=meshgrid(linspace(0,1,point_num),linspace(0,1,point_num));point_Z=rands(point_num,point_num)*0.1+0.5;
% srf=Surface(cat(3,point_X,point_Y,point_Z),UDegree,VDegree);
% fce=Face('',srf);
% 
% axe_hdl=axes(figure());
% fce.displayModel(axe_hdl);
% axis equal;

%% Shell

% UDegree=3;VDegree=3;point_num=5;
% [point_X,point_Y]=meshgrid(linspace(0,2,point_num*2-1),linspace(0,1,point_num));point_Z=rands(point_num,point_num*2-1)*0.1+0.5;
% fce_1=Face('',Surface(cat(3,point_X(:,1:point_num),point_Y(:,1:point_num),point_Z(:,1:point_num)),UDegree,VDegree));
% fce_2=Face('',Surface(cat(3,point_X(:,point_num:end),point_Y(:,point_num:end),point_Z(:,point_num:end)),UDegree,VDegree));
% shl=Shell('',[fce_1;fce_2]);
% 
% axe_hdl=axes(figure());
% shl.displayModel(axe_hdl);
% axis equal;

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
% waverider.displayGeom([],U_num,[],W_num);axis equal;
% waverider.generateShell();
% writeGeomStep(waverider,'waverider_cone.step');
% 
% srf_total=waverider.calGeom(U_num,[],W_num);
% srf_total=waverider.reverseUV(srf_total);
% X_up=srf_total{2}.X;Y_up=srf_total{2}.Y;Z_up=srf_total{2}.Z;
% X_low=srf_total{1}.X;Y_low=srf_total{1}.Y;Z_low=srf_total{1}.Z;
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
% UDegree=2;VDegree=3;u_pole_num=3;v_pole_num=4;
% C_par_X=[0,0];C_par_Y=[par_T,0];C_par_ZV=[0.5,0.5];C_par_ZU=[1,0];sym_x=false;sym_y=true;
% srf_up=SurfaceCST(C_par_X,C_par_Y,C_par_ZV,C_par_ZU,sym_x,sym_y);
% srf_low=SurfaceCST(C_par_X,C_par_Y,C_par_ZV,C_par_ZU,sym_x,sym_y);
% 
% srf_up.addDeform([],[],deform_Z);
% srf_low.addDeform([],[],deform_Z);
% 
% U_node=(X_up(1,:)-min(X_up(1,:)))/(max(X_up(1,:))-min(X_up(1,:)));
% V_node=(Y_up(:,end)-min(Y_up(:,end)))/(max(Y_up(:,end))-min(Y_up(:,end)));
% 
% srf_up.fitSpline(cat(3,X_up,Y_up,Z_up),UDegree,VDegree,u_pole_num,v_pole_num,U_node,V_node);
% srf_low.fitSpline(cat(3,X_low,Y_low,Z_low),UDegree,VDegree,u_pole_num,v_pole_num,U_node,V_node);
% 
% axe_hdl=axes(figure());
% srf_up.displayGeom(axe_hdl);
% srf_up.displayPoles(axe_hdl);
% srf_low.displayGeom(axe_hdl);
% srf_low.displayPoles(axe_hdl);
% axis equal
