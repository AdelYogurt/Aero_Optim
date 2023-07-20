clc;
clear;
close all hidden;

% % point_list=[
% %     0,0;
% %     0,1;
% %     1,1;];
% point_list=rand(5,2);
% 
% % u_list=[0,0,0.2,0.33333,0.66666,0.7,1,1];
% % u_list=[0,0,0,0,1,1,1,1];
% % u_list=[];
% % curve=CurveBSpline('',point_list,[],3);
% % curve=CurveBSpline('',[],point_list,3,[]);
% curve=CurveBSpline('',[],point_list,3);
% curve.drawCurve();
% axis equal;

% surf_CST_B=SurfaceCST3D('',1,1,-0.5,[],[],[0,0],[0,0]);
% 
% % [X,Y,Z,XI,PSI]=surf_CST_B.calSurface(1,2);
% % surface(X,Y,Z);
% % view(3);
% 
% surf_B=surf_CST_B.getSurfaceBSpline(1,2);
% 
% body=Body({surf_B});
% body.writeStepOpenShell('generate.step');

% total_length=4; % total length

% vari_num=16;
% low_bou=[0.7,0.7, ...
%     2.4,0.4,0.5,0.5, ...
%     0.1,0.1,0.005, ...
%     0.4,0.4,0.4,0.6,0.5,0.01,0.01];
% up_bou=[1.0,1.0, ...
%     3.2,0.8,2.0,2.0, ...
%     0.5,0.5,0.02, ...
%     0.6,0.6,0.6,0.8,0.7,0.05,0.05];
% x=(up_bou+low_bou)/2;
% 
% waverider_wing=WaveriderWingDia(total_length,x);
% 
% body=Body({waverider_wing.getSurface('wing_up'),waverider_wing.getSurface('wing_low')});
% body.writeStepOpenShell('generate.step');
