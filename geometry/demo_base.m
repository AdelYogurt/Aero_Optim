clc;
clear;
close all hidden;

%% BSpline curve

% degree=3;point_num=10;
% point_X=rand(point_num,1);point_Y=rand(point_num,1);point_Z=rand(point_num,1);
% axe_hdl=axes(figure());
% scatter3(axe_hdl,point_X,point_Y,point_Z,'Marker','*','MarkerEdgeColor',[0.9290 0.6940 0.1250]);
% curve=CurveBSpline('',point_X,point_Y,point_Z,[],degree);
% curve=CurveBSpline('',point_X,point_Y,point_Z,true,degree);
% curve.drawCurve(axe_hdl);
% axis equal

%% CST2D curve

% LX=1;LY=1;C_par_Y=[0.5,1];sym=false;
% degree=3;ctrl_num=5;
% ctrl_X=linspace(0,1,ctrl_num);ctrl_Y=rand(ctrl_num,1);
% axe_hdl=axes(figure());
% curve=CurveCST2D('',C_par_Y,sym,LX,LY);
% curve.addShapeBSpline(ctrl_X,ctrl_Y,[],degree);
% curve.drawCurve(axe_hdl);
% axis equal;

%% BCST2D curve

% LX=1;LY=0.5;C_par_Y=[1,1];sym=true;
% degree=3;ctrl_num=5;
% ctrl_X=linspace(0,1,ctrl_num);ctrl_Y=rand(ctrl_num,1)*0.5;
% axe_hdl=axes(figure());
% curve=CurveBCST2D('',C_par_Y,sym,LX,LY);
% curve.addShapeBSpline(ctrl_X,ctrl_Y,[],degree);
% curve.drawCurve(axe_hdl);
% curve.addBlunt(1,0.01,-0.1);
% curve.drawCurve(axe_hdl);
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
% degree=4;ctrl_num=5;
% 
% curve_up=CurveCST2D('',C_par,sym,LX,LY);
% curve_low=CurveCST2D('',C_par,sym,LX,-LY);
% curve_up.addShapeBSpline(airfoil_up(:,1),airfoil_up(:,2),true,degree,[],[],ctrl_num,airfoil_up(:,1));
% curve_low.addShapeBSpline(airfoil_low(:,1),airfoil_low(:,2),true,degree,[],[],ctrl_num,airfoil_low(:,1));
% 
% curve_up.optimClass();
% curve_low.optimClass();
% 
% axe_hdl=axes(figure());
% line(axe_hdl,airfoil_up(:,1),airfoil_up(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.8500 0.3250 0.0980])
% line(axe_hdl,airfoil_low(:,1),airfoil_low(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor',[0.8500 0.3250 0.0980])
% curve_up.drawCurve(axe_hdl);
% curve_low.drawCurve(axe_hdl);
% axis equal;
% 
% writematrix([curve_up.ctrl_X,curve_up.ctrl_Y,curve_up.ctrl_Z],'CSTshape_up.txt')
% writematrix([curve_low.ctrl_X,curve_low.ctrl_Y,curve_low.ctrl_Z],'CSTshape_low.txt')

%% BSpline surface

% u_degree=3;v_degree=3;point_num=10;
% [point_X,point_Y]=meshgrid(linspace(0,1,point_num),linspace(0,1,point_num));point_Z=rands(point_num,point_num)*0.1+0.5;
% axe_hdl=axes(figure());
% surface(axe_hdl,point_X,point_Y,point_Z,'Marker','*','MarkerEdgeColor','r','LineStyle','none','FaceAlpha',0);
% surf=SurfaceBSpline('',point_X,point_Y,point_Z,[],u_degree,v_degree);
% surf=SurfaceBSpline('',point_X,point_Y,point_Z,true,u_degree,v_degree);
% surf.drawSurface(axe_hdl);
% axis equal;

%% CST surface

% LX=2;LY=1/2;LZ=0.3;
% C_par_X=[0,0];C_par_Y=[0.5,0];C_par_ZV=[15,15];C_par_ZU=[0.5,0];sym_x=false;sym_y=true;
% surf_up=SurfaceCST('',C_par_X,C_par_Y,C_par_ZV,C_par_ZU,sym_x,sym_y,LX,LY,LZ);
% surf_low=SurfaceCST('',C_par_X,C_par_Y,C_par_ZV,C_par_ZU,sym_x,sym_y,LX,-LY,LZ);
% axe_hdl=axes(figure());
% surf_up.drawSurface(axe_hdl);
% surf_low.drawSurface(axe_hdl);
% axis equal

%% Coons surface

% curve_u0=CurveBSpline('',[0,0,0.5],[0,1,2],[0,0,0]);
% curve_u1=CurveBSpline('',[0,0.5,1],[0,1,2.5],[2.5,2,2]);
% curve_0v=CurveBSpline('',[0,0,0],[0,1,0],[0,1,2.5]);
% curve_1v=CurveBSpline('',[0.5,0.5,1],[2,1.5,2.5],[0,1,2]);
% 
% axe_hdl=axes(figure());
% axis equal;view(3);
% 
% curve_u0.drawCurve(axe_hdl);
% curve_1v.drawCurve(axe_hdl);
% curve_u1.drawCurve(axe_hdl);
% curve_0v.drawCurve(axe_hdl);
% 
% surf=SurfaceCoons('',@(u) curve_u0.calPoint(u),@(u) curve_u1.calPoint(u),...
%     @(v) curve_0v.calPoint(v),@(v) curve_1v.calPoint(v));
% [U,V]=meshgrid(0:0.1:1,0:0.1:1);
% [X,Y,Z]=surf.calPoint(U,V);
% surface(X,Y,Z);

%% mapping generate surface

% axe_hdl=axes(figure());
% axis equal;view(3);
% 
% line_u0=[0,0,0.5;0,1,2;0,0,0]';
% line_u1=[0,0.5,1;0,1,2.5;2.5,2,2]';
% line_v0=[0,0,0;0,1,0;0,1,2.5]';
% line_v1=[0.5,0.5,1;2,1.5,2.5;0,1,2]';
% curve_u0=CurveBSpline('',[0,0,0.5],[0,1,2],[0,0,0]);
% curve_u1=CurveBSpline('',[0,0.5,1],[0,1,2.5],[2.5,2,2]);
% curve_0v=CurveBSpline('',[0,0,0],[0,1,0],[0,1,2.5]);
% curve_1v=CurveBSpline('',[0.5,0.5,1],[2,1.5,2.5],[0,1,2]);
% curve_u0.drawCurve(axe_hdl);
% curve_1v.drawCurve(axe_hdl);
% curve_u1.drawCurve(axe_hdl);
% curve_0v.drawCurve(axe_hdl);
% 
% [point_X,point_Y,point_Z]=geomMapGrid(line_u0,line_u1,line_v0,line_v1);
% surf=SurfaceBSpline('',point_X,point_Y,point_Z);
% surf.drawSurface(axe_hdl,[],[],struct('LineStyle','none','FaceAlpha',0.5));

%% Fit 3D wing

% load('wing/wing_RAE2822_NACA4412.mat');
% u_degree=3;v_degree=3;u_ctrl_num=4;v_ctrl_num=8;
% C_par_X=[0,0];C_par_Y=[0,0];C_par_ZV=[0,0];C_par_ZU=[0.5,1];sym_x=false;sym_y=false;
% surf_up=SurfaceCST('',C_par_X,C_par_Y,C_par_ZV,C_par_ZU,sym_x,sym_y);
% surf_low=SurfaceCST('',C_par_X,C_par_Y,C_par_ZV,C_par_ZU,sym_x,sym_y);
% surf_up.addShapeBSpline(X_up,Y_up,Z_up,true,u_degree,v_degree,[],[],[],[],u_ctrl_num,v_ctrl_num,X_up(1,:)/max(X_up(1,:)),Y_up(:,1)/max(Y_up(:,1)));
% surf_low.addShapeBSpline(X_low,Y_low,Z_low,true,u_degree,v_degree,[],[],[],[],u_ctrl_num,v_ctrl_num,X_low(1,:)/max(X_low(1,:)),Y_low(:,1)/max(Y_low(:,1)));
% 
% axe_hdl=axes(figure());
% surf_up.drawSurface(axe_hdl,20,20);
% surf_low.drawSurface(axe_hdl,20,20);
% axis equal

%% generate 3D wing

% % ----generate 3D wing by interplate airfoil between two section----
% ang_W_LE=27.1; % deg
% span_W=1171.3;
% chord_W_root=240.1;
% chord_W_tip=60.6;
% span_mid=235.1*2;
% 
% [curve_up_1,curve_low_1]=airfoilFit(importdata('airfoil\airfoil_RAE2822.txt'));
% [curve_up_2,curve_low_2]=airfoilFit(importdata('airfoil\airfoil_NACA4412.txt'));
% 
% [U,V]=meshgrid(sin(linspace(0,1,51)*pi/2).^2,linspace(0,1,51));
% 
% [X_up,Y_up,Z_up]=calWing(U,V,chord_W_root,chord_W_tip,ang_W_LE,span_W,span_mid,...
%     curve_up_1,curve_up_1,curve_up_2);
% [X_low,Y_low,Z_low]=calWing(U,V,chord_W_root,chord_W_tip,ang_W_LE,span_W,span_mid,...
%     curve_low_1,curve_low_1,curve_low_2);
% 
% surface(X_up,Y_up,Z_up);
% surface(X_low,Y_low,Z_low);
% axis equal;view(3);
% 
% save('wing/wing_RAE2822_NACA4412.mat','X_up','Y_up','Z_up','X_low','Y_low','Z_low');
% 
% function [X,Y,Z]=calWing(U,V,chord_W_root,chord_W_tip,ang_W_LE,span_W,span_mid,...
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
%         [~,Y_1]=airfoil_1.calPoint(U(1,:));
%         [~,Y_2]=airfoil_2.calPoint(U(1,:));
%         Z=Y_1'.*(1-V)+Y_2'.*V;
%     end
% end
% 
% function [curve_up,curve_low]=airfoilFit(airfoil_data)
% LX=1;LY=1;C_par=[0.5,1];sym=false;
% degree=5;ctrl_num=6;
% 
% u_num=size(airfoil_data,1)/2;
% airfoil_up=airfoil_data(1:u_num,:);airfoil_low=airfoil_data(u_num+1:end,:);
% 
% curve_up=CurveCST2D('',LX,LY,C_par,sym);
% curve_low=CurveCST2D('',LX,LY,C_par,sym);
% curve_up.addShapeBSpline(airfoil_up(:,1),airfoil_up(:,2),degree,[],[],ctrl_num,airfoil_up(:,1));
% curve_low.addShapeBSpline(airfoil_low(:,1),airfoil_low(:,2),degree,[],[],ctrl_num,airfoil_low(:,1));
% end
% % ----end----