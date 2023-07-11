clc;
clear;
close all hidden;

addpath geometry;
addpath model;
addpath optimal;
addpath SU2;

% my_logger=Logger('optNACA0012.log');
% problem=NACA0012AirfoilProblem(128,2.5,0.0,0.8,273.15,101325.0,273.15,...
%     'optNACA0012','SU2','model/SU2_temp/opt',[],[]);
problem=RAE2822AirfoilProblem(128,2.5,0.0,0.8,273.15,101325.0,273.15,...
    'opt_RAE2822','SU2','model/SU2_temp/opt',[],[]);

% max_NFE=200;
% optimizer=OptSACORS(max_NFE,300);
% disp('optimization begin')
% [x_best,obj_best,NFE,output]=optimizer.optimize(problem);
% save('optres_SACO_RS.mat','x_best','obj_best','NFE','output');
% disp('optimization all done')

x=[0.0207,0.0604,0.0551,0.0571,0.0728,...
    0.0946,0.0275,0.1057,0.0469,0.0930];
[obj,con,coneq]=problem.objcon_fcn(x)


