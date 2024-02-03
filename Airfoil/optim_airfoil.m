clc;
clear;
close all hidden;

%% PATH

cd ..
addpath cgns4m;
addpath coord;
addpath mesh;
addpath geometry;
addpath model;
startup_cgns4m;
cd Airfoil;

%% define problem

% basical parameter
if ispc(),partitions=1;
else,partitions=8;end
run_desc=[];
out_logger=[];

% geometry module
% geom_param.solver='HICKS_HENNE';
% geom_param.mesh_point=load('mesh/mesh_data_airfoil.mat','mesh_point').mesh_point;
% hicks_hemme_define.MARKER=[repmat({'AIRFOIL_LOW'},5,1);repmat({'AIRFOIL_UP'},5,1)];
% hicks_hemme_define.PARAM=[[zeros(5,1),(1/6:1/6:5/6)'];[ones(5,1),(1/6:1/6:5/6)']];
% geom_param.hicks_hemme_define=hicks_hemme_define;

geom_param.solver='CST';
geom_param.mesh_coord=load('mesh/mesh_data_airfoil_26.mat','mesh_coord').mesh_coord;

% mesh module
mesh_param.solver='SU2_DEF';
mesh_param.initial_mesh_filestr='mesh/NACA0012_26.cgns';
mesh_param.SU2_DEF_param='SU2/NACA0012_deform.cfg';
mesh_param.dat_filestr='geom/NACA0012_deform.dat';
mesh_param.mesh_filestr='geom/airfoil.su2';

% CFD module
CFD_param.solver='Fluent';
CFD_param.fluent_jou_filestr='Fluent/airfoil_80_4.jou';
CFD_param.fluent_dir='';
CFD_param.solver_dimension=2;
CFD_param.out_filename='fluent_history.out';

% CFD_param.solver='SU2_CFD';
% CFD_param.SU2_CFD_param='SU2/airfoil.cfg';

model=AirfoilModel(partitions,geom_param,mesh_param,CFD_param,[],[],run_desc,out_logger);
geo_in.C_par_low=[0.5,1.0];
geo_in.Poles_low=importdata('geom/RAE2822_CSTshape_low.txt');
geo_in.C_par_up=[0.5,1.0];
geo_in.Poles_up=importdata('geom/RAE2822_CSTshape_up.txt');
[geo_out,mesh_out,CFD_out]=model.solveModel(geo_in);

% problem=AirfoilProblem(partitions,geom_param,mesh_param,CFD_param,[],[],run_desc,out_logger);
% x=(problem.low_bou+problem.up_bou)/2;
% problem.drawX(x)
% [obj,con,coneq]=problem.objcon_fcn(x);

% problem=PhysicsProblem(partitions,geom_param,mesh_param,CFD_param,[],[],run_desc,out_logger);
% x=(problem.low_bou+problem.up_bou)/2;
% problem.drawX(x)
% [obj,con,coneq]=problem.objcon_fcn(x);

%% optimize PATH

cd ../..
addpath(genpath('Srgt_Optim'));
cd Aero_Optim/Airfoil;

%% define optimize

% problem=AirfoilProblem(partitions,geom_param,mesh_param,CFD_param,[],[],run_desc,out_logger);
% problem=PhysicsProblem(partitions,geom_param,mesh_param,CFD_param,[],[],run_desc,out_logger);
% NFE_max=100;iter_max=300;obj_torl=1e-6;con_torl=0;

% optimizer=OptimRBFCDE(NFE_max,iter_max,obj_torl,con_torl);
% optimizer=OptimSKO(NFE_max,iter_max,obj_torl,con_torl);

% disp('optimization begin')
% [x_best,obj_best,NFE,output,con_best,coneq_best,vio_best]=optimizer.optimize(problem);
% disp('optimization all done')
% save('optres.mat','x_best','obj_best','NFE','output','con_best','coneq_best','vio_best');

%% parallel run

% NFE_max=100;iter_max=300;obj_torl=1e-6;con_torl=0;
% repeat_num=12;
% 
% data_list=cell(1,repeat_num);
% parfor par_idx=1:repeat_num
%     run_desc=['optim',num2str(par_idx)];
% 
%     % mesh module
%     mesh_param=struct();
%     mesh_param.solver='SU2_DEF';
%     mesh_param.initial_mesh_filestr='mesh/NACA0012.cgns';
%     mesh_param.SU2_DEF_param=SU2Config('SU2/NACA0012_deform.cfg');
%     mesh_param.dat_filestr=['optim/',num2str(par_idx),'/NACA0012_deform.dat'];
%     mesh_param.mesh_filestr=['optim/',num2str(par_idx),'/airfoil.su2'];
% 
%     % problem=AirfoilProblem(partitions,geom_param,mesh_param,CFD_param,[],[],run_desc,out_logger);
%     problem=PhysicsProblem(partitions,geom_param,mesh_param,CFD_param,[],[],run_desc,out_logger);
% 
%     optimizer=OptimSKO(NFE_max,iter_max,obj_torl,con_torl);
% 
%     % x=rand(1,problem.vari_num).*(problem.up_bou-problem.low_bou)+problem.low_bou;
%     % [obj,con,coneq]=problem.objcon_fcn(x);
%     % data_list{par_idx}={obj,con,coneq};
%     % disp(obj);disp(con);disp(coneq);
% 
%     disp('optimization begin')
%     [x_best,obj_best,NFE,output,con_best,coneq_best,vio_best]=optimizer.optimize(problem);
%     disp('optimization all done')
%     data_list{par_idx}=output;
% end
% 
% save('optim.mat','data_list');
