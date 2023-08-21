clc;
clear;
close all hidden;

cd ..

addpath cgns4m;
addpath coord;
addpath mesh;
addpath geometry;
addpath model;
addpath optimal;

startup_cgns4m;

%% define problem

% basical parameter
if ispc()
    partitions=4;
else
    partitions=128;
end
run_description=[];
out_logger=[];

% geometry module
geo_param.solver='HICKS_HENNE';
geo_param.mesh_point=load('Airfoil/Fluent/mesh_data_airfoil.mat','mesh_point').mesh_point;
hicks_hemme_define.MARKER=[repmat({'AIRFOIL_LOW'},5,1);repmat({'AIRFOIL_UP'},5,1)];
hicks_hemme_define.PARAM=[[zeros(5,1),(1/6:1/6:5/6)'];[ones(5,1),(1/6:1/6:5/6)']];
geo_param.hicks_hemme_define=hicks_hemme_define;

% geo_param.solver='CST';
% geo_param.mesh_coord=load('Airfoil/Fluent/mesh_data_airfoil.mat','mesh_coord').mesh_coord;

% mesh module
mesh_param.solver='SU2_DEF';
mesh_param.initial_mesh_filestr='Airfoil/Fluent/NACA0012.cgns';
mesh_param.SU2_DEF_param='Airfoil/SU2/NACA0012_deform.cfg';
mesh_param.dat_filestr='model/SU2_temp/optiaml/NACA0012_deform.dat';
mesh_param.mesh_filestr='model/SU2_temp/optiaml/airfoil.su2';

% CFD module
CFD_param.solver='Fluent';
CFD_param.fluent_jou_filestr='Airfoil/Fluent/airfoil.jou';
if ispc()
    CFD_param.fluent_dir='D:/Program_Files/ANSYS Inc/v222/fluent/ntbin/win64';
else
    CFD_param.fluent_dir='/public3/wshome/ws173/scg7758/software-scg7758/install231/install231/ansys_inc/v231/fluent/bin';
end
CFD_param.solver_dimension=2;
CFD_param.out_filename='fluent_history.out';

% CFD_param.solver='SU2_CFD';
% CFD_param.SU2_CFD_param='Airfoil/SU2/airfoil.cfg';

problem=ProblemAirfoil(partitions,geo_param,mesh_param,CFD_param,[],[],run_description,out_logger);

%% define optimize

max_NFE=200;

optimizer=OptSACORS(max_NFE,300,1e-6,0,[]);
% optimizer=OptSADE(max_NFE,300,1e-6,0,[]);

disp('optimization begin')
[x_best,obj_best,NFE,output]=optimizer.optimize(problem);
save('Airfoil/optres_SACORS.mat','x_best','obj_best','NFE','output');
% save('Airfoil/optres_SADE.mat','x_best','obj_best','NFE','output');
disp('optimization all done')
