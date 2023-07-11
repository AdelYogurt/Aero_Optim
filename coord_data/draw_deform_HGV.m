clc;
clear;
close all hidden;


data=importdata('WDB_deform.dat');

scatter3(data(:,2),data(:,3),data(:,4))
view(3);
axis equal;

