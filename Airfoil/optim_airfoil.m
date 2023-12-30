clc;
clear;
close all hidden;

% cd ..
% addpath cgns4m;
% addpath coord;
% addpath mesh;
% addpath geometry;
% addpath model;
% addpath optimizer;
% startup_cgns4m;
% cd Airfoil;

%% define problem

% basical parameter
if ispc(),partitions=1;
else,partitions=128;end
run_desc=[];
out_logger=[];

% geometry module
% geom_param.solver='HICKS_HENNE';
% geom_param.mesh_point=load('mesh_data_airfoil.mat','mesh_point').mesh_point;
% hicks_hemme_define.MARKER=[repmat({'AIRFOIL_LOW'},5,1);repmat({'AIRFOIL_UP'},5,1)];
% hicks_hemme_define.PARAM=[[zeros(5,1),(1/6:1/6:5/6)'];[ones(5,1),(1/6:1/6:5/6)']];
% geom_param.hicks_hemme_define=hicks_hemme_define;

geom_param.solver='CST';
geom_param.mesh_coord=load('mesh_data_airfoil.mat','mesh_coord').mesh_coord;

% mesh module
mesh_param.solver='SU2_DEF';
mesh_param.initial_mesh_filestr='Fluent/NACA0012.cgns';
mesh_param.SU2_DEF_param='SU2/NACA0012_deform.cfg';
mesh_param.dat_filestr='geom/NACA0012_deform.dat';
mesh_param.mesh_filestr='geom/airfoil.su2';

% CFD module
CFD_param.solver='Fluent';
CFD_param.fluent_jou_filestr='Fluent/airfoil.jou';
CFD_param.fluent_dir='';
CFD_param.solver_dimension=2;
CFD_param.out_filename='fluent_history.out';

% CFD_param.solver='SU2_CFD';
% CFD_param.SU2_CFD_param='SU2/airfoil.cfg';

% model=AirfoilModel(partitions,geom_param,mesh_param,CFD_param,[],[],run_desc,out_logger);
% geo_in.up=[linspace(0,1,5)',rand(5,1)*0.1];
% geo_in.low=[linspace(0,1,5)',-rand(5,1)*0.1];
% [geo_out,mesh_out,CFD_out]=model.solveAirfoil(geo_in);

problem=AirfoilProblem(partitions,geom_param,mesh_param,CFD_param,[],[],run_desc,out_logger);
x=(problem.low_bou+problem.up_bou)/2;
problem.drawX(x)
[obj,con,coneq]=problem.objcon_fcn(x);

%% define optimize

% max_NFE=200;
% 
% optimizer=OptimRBFCDE(max_NFE,300,1e-6,0);
% 
% disp('optimization begin')
% [x_best,obj_best,NFE,output,con_best,coneq_best,vio_best]=optimizer.optimize(problem);
% save('optimres.mat','x_best','obj_best','NFE','output','con_best','coneq_best','vio_best');
% disp('optimization all done')

%% parallel run

% repeat_num=10;
% data_list=cell(1,repeat_num);
% parfor par_idx=1:repeat_num
%     run_desc=num2str(par_idx);
% 
%     % mesh module
%     mesh_param=struct();
%     mesh_param.solver='SU2_DEF';
%     mesh_param.initial_mesh_filestr='Fluent/NACA0012.cgns';
%     mesh_param.SU2_DEF_param=SU2Config('SU2/NACA0012_deform.cfg');
%     mesh_param.dat_filestr=['geom/',num2str(par_idx),'/NACA0012_deform.dat'];
%     mesh_param.mesh_filestr=['geom/',num2str(par_idx),'/airfoil.su2'];
% 
%     problem=AirfoilProblem(partitions,geom_param,mesh_param,CFD_param,[],[],run_desc,out_logger);
% 
%     x=rand(1,problem.vari_num).*(problem.up_bou-problem.low_bou)+problem.low_bou;
%     [obj,con,coneq]=problem.objcon_fcn(x);
%     data_list{par_idx}={obj,con,coneq};
%     disp(obj);disp(con);disp(coneq);
% end
