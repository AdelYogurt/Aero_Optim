clc;
clear;
close all hidden;

addpath mesh\

fcn=@(x) x^0.7;
low_bou=0;
up_bou=1;
torl=1e-3;
max_level=50;

[X,Fval,node_list]=girdAdapt1D(fcn,low_bou,up_bou,torl,max_level);
line(X,Fval,'marker','o');

% fcn=@(x) sum(x.^0.5);
% low_bou=[0,0];
% up_bou=[1,1];
% torl=1e-3;
% max_level=50;
% 
% [X,Fval,node_list]=girdAdapt2DM(fcn,low_bou,up_bou,torl,max_level);
% line(X,Fval,'marker','o');
