clc;
clear;
close all hidden;

%% Spline Curve

% Degree=3;group_num=3;
% Points=[0,0,0;(0.5:0.5:group_num)',repmat([1;0],group_num,1),repmat([1;0],group_num,1)];
% axe_hdl=axes(figure());
% scatter3(axe_hdl,Points(:,1),Points(:,2),Points(:,3),'Marker','*','MarkerEdgeColor',[0.9290 0.6940 0.1250]);
% crv=Curve(Points,Degree);
% crv=GeomApp.interpPointToCurve(Points,Degree);
% crv=Curve([1,0;1,1;0,1],[],[],[],[1,sqrt(2)/2,1]); % arc
% crv.displayGeom(axe_hdl);
% crv.displayPoles(axe_hdl);
% axis equal;

%% CST2D Curve

% LX=1;LY=0.06;C_par_Y=[0.5,1];sym=false;
% Degree=3;pole_num=5;
% ctrl_X=linspace(0,1,pole_num)';ctrl_Y=rand(pole_num,1);
% axe_hdl=axes(figure());
% crv=CurveCST(C_par_Y,sym,LX,LY);
% crv.addSpline(Curve([ctrl_X,ctrl_Y],Degree));
% 
% crv.displayGeom(axe_hdl);
% crv.displayPoles(axe_hdl);
% 
% crv=crv.convertSpline();
% crv.displayGeom(axe_hdl,struct('Color','g','LineStyle','--'));
% crv.displayPoles(axe_hdl,struct('Color','k','LineStyle',':','Marker','s'));
% axis equal;

%% fit airfoil with CST2D

% airfoil_data=importdata('airfoil\NACA0012.txt');
% airfoil_up=airfoil_data(1:67,:);airfoil_low=airfoil_data(68:end,:);
% 
% airfoil_data=importdata('airfoil\NACA4412.txt');
% airfoil_up=airfoil_data(1:18,:);airfoil_low=airfoil_data(19:end,:);
% 
% airfoil_data=importdata('airfoil\RAE2822.txt');
% airfoil_up=airfoil_data(1:65,:);airfoil_low=airfoil_data(66:end,:);
% 
% airfoil_data=importdata('airfoil\Clark_Y.txt');
% airfoil_up=airfoil_data(1:61,:);airfoil_low=airfoil_data(62:end,:);
% 
% LX=1;LY=1;C_par=[0.5,1];sym=false;
% Degree=4;pole_num=5;
% 
% crv_up=CurveCST(C_par,sym,LX,LY);
% crv_low=CurveCST(C_par,sym,LX,-LY);
% crv_up.fitSpline(airfoil_up,Degree,pole_num,airfoil_up(:,1));
% crv_low.fitSpline(airfoil_low,Degree,pole_num,airfoil_low(:,1));
% 
% crv_up.optimClass();
% crv_low.optimClass();
% 
% axe_hdl=axes(figure());
% line(axe_hdl,airfoil_up(:,1),airfoil_up(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.8500 0.3250 0.0980])
% line(axe_hdl,airfoil_low(:,1),airfoil_low(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.8500 0.3250 0.0980])
% crv_up.displayGeom(axe_hdl);
% crv_low.displayGeom(axe_hdl);
% axis equal;
% 
% writematrix(crv_up.Poles,'CSTshape_up.txt')
% writematrix(crv_low.Poles,'CSTshape_low.txt')

%% Spline Surface

% UDegree=3;VDegree=3;point_num=5;
% [point_X,point_Y]=meshgrid(linspace(0,1,point_num),linspace(0,1,point_num));point_Z=rands(point_num,point_num)*0.1+0.5;
% axe_hdl=axes(figure());
% surface(axe_hdl,point_X,point_Y,point_Z,'Marker','*','MarkerEdgeColor','r','LineStyle','none','FaceAlpha',0);
% srf=Surface(cat(3,point_X,point_Y,point_Z),UDegree,VDegree);
% srf=GeomApp.interpPointToSurface(cat(3,point_X,point_Y,point_Z),UDegree,VDegree);
% srf=Surface(cat(3,[1,1,0;1,1,0;0,0,0],[0,1,1;0,1,1;0,0,0],[0,0,0;1,1,1;1,1,1]),[],[],[],[],[],[],[1,sqrt(2)/2,1;sqrt(2)/2,0.5,sqrt(2)/2;1,sqrt(2)/2,1]); % sphere
% srf.displayGeom(axe_hdl);
% srf.displayPoles(axe_hdl);
% axis equal;

%% CST Surface

% LX=2;LY=1/2;LZ=0.3;
% C_par_X=[0,0];C_par_Y=[0.5,0];C_par_ZV=[15,15];C_par_ZU=[0.5,0];sym_x=false;sym_y=true;
% srf_up=SurfaceCST(C_par_X,C_par_Y,C_par_ZV,C_par_ZU,sym_x,sym_y,LX,LY,LZ);
% srf_low=SurfaceCST(C_par_X,C_par_Y,C_par_ZV,C_par_ZU,sym_x,sym_y,LX,-LY,LZ);
% axe_hdl=axes(figure());
% srf_up.displayGeom(axe_hdl);
% srf_low.displayGeom(axe_hdl);
% srf_up=srf_up.convertSpline();
% srf_low=srf_low.convertSpline();
% srf_up.displayGeom(axe_hdl,struct('FaceColor','g','FaceAlpha',0.5,'LineStyle','none'));
% srf_up.displayPoles(axe_hdl);
% srf_low.displayGeom(axe_hdl,struct('FaceColor','g','FaceAlpha',0.5,'LineStyle','none'));
% srf_low.displayPoles(axe_hdl);
% axis equal

%% mapping generate Surface

% axe_hdl=axes(figure());
% axis equal;view(3);
% 
% line_u0=[0,0,0.5;0,1,2;0,0,0]';
% line_u1=[0,0.5,1;0,1,2.5;2.5,2,2]';
% line_v0=[0,0,0;0,1,0;0,1,2.5]';
% line_v1=[0.5,0.5,1;2,1.5,2.5;0,1,2]';
% edge_u0=Curve(line_u0);
% edge_u1=Curve(line_u1);
% edge_0v=Curve(line_v0);
% edge_1v=Curve(line_v1);
% edge_u0.displayGeom(axe_hdl);
% edge_1v.displayGeom(axe_hdl);
% edge_u1.displayGeom(axe_hdl);
% edge_0v.displayGeom(axe_hdl);
% 
% Poles=GeomApp.MapGrid(line_u0,line_u1,line_v0,line_v1);
% srf=Surface(Poles);
% srf.displayGeom(axe_hdl,struct('LineStyle','none','FaceAlpha',0.5));
% srf.displayPoles(axe_hdl);

%% Fit 3D wing

% load('wing/wing_RAE2822_NACA4412.mat');
% UDegree=3;VDegree=3;u_pole_num=4;v_pole_num=8;
% C_par_X=[0,0];C_par_Y=[0,0];C_par_ZV=[0,0];C_par_ZU=[0.5,1];sym_x=false;sym_y=false;
% srf_up=SurfaceCST(C_par_X,C_par_Y,C_par_ZV,C_par_ZU,sym_x,sym_y);
% srf_low=SurfaceCST(C_par_X,C_par_Y,C_par_ZV,C_par_ZU,sym_x,sym_y);
% 
% U_node_up=X_up(1,:)/max(X_up(1,:));
% V_node_up=Y_up(:,1)/max(Y_up(:,1));
% U_node_low=X_low(1,:)/max(X_low(1,:));
% V_node_low=Y_low(:,1)/max(Y_low(:,1));
% 
% srf_up.fitSpline(cat(3,X_up,Y_up,Z_up),UDegree,VDegree,u_pole_num,v_pole_num,U_node_up,V_node_up);
% srf_low.fitSpline(cat(3,X_low,Y_low,Z_low),UDegree,VDegree,u_pole_num,v_pole_num,U_node_low,V_node_low);
% 
% axe_hdl=axes(figure());
% srf_up.displayGeom(axe_hdl);
% srf_up.displayPoles(axe_hdl);
% srf_low.displayGeom(axe_hdl);
% srf_low.displayPoles(axe_hdl);
% axis equal

%% generate 3D wing

% % ----generate 3D wing by interplate airfoil between two section----
% ang_W_LE=27.1; % deg
% span_W=1171.3;
% chord_W_root=240.1;
% chord_W_tip=60.6;
% span_mid=235.1*2;
% 
% [edge_up_1,edge_low_1]=airfoilFit(importdata('airfoil\airfoil_RAE2822.txt'));
% [edge_up_2,edge_low_2]=airfoilFit(importdata('airfoil\airfoil_NACA4412.txt'));
% 
% [U,V]=meshgrid(sin(linspace(0,1,51)*pi/2).^2,linspace(0,1,51));
% 
% [X_up,Y_up,Z_up]=calWing(U,V,chord_W_root,chord_W_tip,ang_W_LE,span_W,span_mid,...
%     edge_up_1,edge_up_1,edge_up_2);
% [X_low,Y_low,Z_low]=calWing(U,V,chord_W_root,chord_W_tip,ang_W_LE,span_W,span_mid,...
%     edge_low_1,edge_low_1,edge_low_2);
% 
% surface(X_up,Y_up,Z_up);
% surface(X_low,Y_low,Z_low);
% axis equal;view(3);
% 
% save('wing/wing_RAE2822_NACA4412.mat','X_up','Y_up','Z_up','X_low','Y_low','Z_low');
% 
% function Points=calWing(U,V,chord_W_root,chord_W_tip,ang_W_LE,span_W,span_mid,...
%     airfoil_root,airfoil_mid,airfoil_tip)
% v_mid=span_mid/span_W;
% 
% X_deform_fcn=@(U,V) V*tan(ang_W_LE/180*pi)*span_W/2;
% 
% tan_ang_W_TE=(tan(ang_W_LE/180*pi)*span_W/2+chord_W_tip-chord_W_root)/(span_W-span_mid)*2;
% point_TE_x=chord_W_root-tan_ang_W_TE*span_mid/2;
% 
% chord_W_in_fcn=@(U,V) chord_W_root-V*span_W/2*tan(ang_W_LE/180*pi);
% chord_W_out_fcn=@(U,V) point_TE_x-V*(point_TE_x-chord_W_tip);
% 
% X=zeros(size(U));
% Y=zeros(size(U));
% Z=zeros(size(U));
% 
% airfoil_in_fcn=@(U,V) airfoilChord(U,V/v_mid,airfoil_root,airfoil_mid);
% airfoil_out_fcn=@(U,V) airfoilChord(U,(V-v_mid)/(1-v_mid),airfoil_mid,airfoil_tip);
% 
% idx=find(V(:,1)> v_mid, 1);
% X(1:idx-1,:)=U(1:idx-1,:).*chord_W_in_fcn(U(1:idx-1,:),V(1:idx-1,:))+X_deform_fcn(U(1:idx-1,:),V(1:idx-1,:));
% Y(1:idx-1,:)=V(1:idx-1,:).*span_W/2;
% Z(1:idx-1,:)=airfoil_in_fcn(U(1:idx-1,:),V(1:idx-1,:)).*chord_W_in_fcn(U(1:idx-1,:),V(1:idx-1,:));
% 
% X(idx:end,:)=U(idx:end,:).*chord_W_out_fcn(U(idx:end,:),V(idx:end,:))+X_deform_fcn(U(idx:end,:),V(idx:end,:));
% Y(idx:end,:)=V(idx:end,:).*span_W/2;
% Z(idx:end,:)=airfoil_out_fcn(U(idx:end,:),V(idx:end,:)).*chord_W_out_fcn(U(idx:end,:),V(idx:end,:));
% 
%     function Z=airfoilChord(U,V,airfoil_1,airfoil_2)
%         [~,Y_1]=airfoil_1.calPoints(U(1,:));
%         [~,Y_2]=airfoil_2.calPoints(U(1,:));
%         Z=Y_1'.*(1-V)+Y_2'.*V;
%     end
% end
% 
% function [crv_up,crv_low]=airfoilFit(airfoil_data)
% LX=1;LY=1;C_par=[0.5,1];sym=false;
% Degree=5;pole_num=6;
% 
% u_num=size(airfoil_data,1)/2;
% airfoil_up=airfoil_data(1:u_num,:);airfoil_low=airfoil_data(u_num+1:end,:);
% 
% crv_up=CurveCST(LX,LY,C_par,sym);
% crv_low=CurveCST(LX,LY,C_par,sym);
% crv_up.addShapeBSpline(airfoil_up(:,1),airfoil_up(:,2),Degree,[],[],pole_num,airfoil_up(:,1));
% crv_low.addShapeBSpline(airfoil_low(:,1),airfoil_low(:,2),Degree,[],[],pole_num,airfoil_low(:,1));
% end
% % ----end----