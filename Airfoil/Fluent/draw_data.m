clc
clear;
close all hidden;

data=importdata('exp_data/NACA0012_surface_flow.csv');
X=data(:,1);
CP=data(:,2);
num=length(X)/2;
XL=X(1:num);CPL=CP(1:num);
XU=X(num+1:end);CPU=CP(num+1:end);

% [XL_opt,CPL_opt,XU_opt,CPU_opt]=getCp('surface_flow.csv');

load('exp_data/NACA0012_experiment_data.mat');

line(XL,CPL,'Color','r','LineWidth',1);
line(XU,CPU,'Color','r','LineWidth',1);
line(NACA0012Cp(:,1),NACA0012Cp(:,2),'Color','k','MarkerFaceColor','k','LineStyle','none','Marker','^');

% line(XL_opt,CPL_opt,'Color','g','LineWidth',1);
% line(XU_opt,CPU_opt,'Color','g','LineWidth',1);

gca().set('YDir','reverse');
legend({'','仿真数据','实验数据'});
box on;xlabel('\it\fontname{times new roman}x');ylabel('\fontname{times new roman}{\itC}_P');

print(gcf(),'NACA0012_CP.png', '-dpng','-r1200');

function [XL,CPL,XU,CPU]=getCpSU2(data_filestr)
% load pressure coefficient
%
SU2_surface=readSU2CSV(data_filestr);
CP=SU2_surface.Pressure_Coefficient;
X=SU2_surface.x;
Y=SU2_surface.y;
Bool=Y>0;
XU=X(Bool);
CPU=CP(Bool);
XL=X(~Bool);
CPL=CP(~Bool);
[XU,idx]=sort(XU);
CPU=CPU(idx);
[XL,idx]=sort(XL);
CPL=CPL(idx);

end


