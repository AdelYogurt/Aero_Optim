clc;
clear;
close all hidden;

% point_list=[
%     0,0;
%     1,0;
%     1,1;
%     0,1;
%     0,0;]+1;
% 
% u_list=[0,0,0.2,0.33333,0.66666,0.7,1,1];
% % u_list=[0,0,0.0,1/3,2/3,1.0,1,1];
% % u_list=[];
% % curve=CurveBSpline(point_list,[],3);
% curve=CurveBSpline([],point_list,3,[]);
% curve.drawCurve();

surf_CST_A=SurfaceCST3D(1,1,0.5,[],[0.5,0],[0.5,0.5],[0.5,0.5]);
surf_CST_B=SurfaceCST3D(1,1,-0.5,[],[],[0.5,0.5],[0.5,0.5]);

[X,Y,Z]=surf_CST_A.calSurface();
surface(X,Y,Z);

% body=Body([surf_CST_A,surf_CST_B]);
% body.writeStepShell('generate.step');

