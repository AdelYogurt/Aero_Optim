clc;
clear;
close all hidden;

load('model_LDratio.mat');

problem=AirfoilProblem(1);

x=X(randi(100),:);
control_point_up=[linspace(0,1,7)',x(1:7)'];
control_point_low=[linspace(0,1,7)',x(8:14)'];
total_point_list=problem.prePoint(control_point_up,control_point_low,'coord_data');
hold on;
scatter(total_point_list{1}(:,2),total_point_list{1}(:,3));
scatter(total_point_list{2}(:,2),total_point_list{2}(:,3));
