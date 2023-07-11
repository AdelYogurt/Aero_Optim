clc;
clear;
close all hidden;

cd ..

addpath geometry;
addpath model;
addpath optimal;
addpath SU2;

my_logger=Logger('multiprocess/opt_SU2_SACO_RS.log');
problem=RAE2822AirfoilProblem(128,2.5,0.0,0.8,273.15,101325.0,273.15,...
    'opt_SACO_RS','SU2','model/SU2_temp/opt_SACO_RS',[],my_logger);

max_NFE=200;
optimizer=OptSACORS(max_NFE,300);
disp('optimization begin')
[x_best,obj_best,NFE,output]=optimizer.optimize(problem);
save('optres_SACO_RS.mat','x_best','obj_best','NFE','output');
disp('optimization all done')




