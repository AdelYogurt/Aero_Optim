clc
clear;
close all hidden;

%% load data

exp_Cp=importdata('RAE2822_Cp_exp.csv');
data=importdata('RAE2822_Cp_Fluent.csv');
[XL,CPL,XU,CPU]=splitCp(data);

%% draw data

fig_hdl=figure();
axe_hdl=axes(fig_hdl);

line(XL,CPL,'Color','r','LineWidth',1);
line(XU,CPU,'Color','r','LineWidth',1);
line(exp_Cp(:,1),exp_Cp(:,2),'Color','k','MarkerFaceColor','k','LineStyle','none','Marker','^');

axe_hdl.set('XTick',[0,0.2,0.4,0.6,0.8,1],'FontName','times new romax')
axe_hdl.set('YDir','reverse','FontSize',10.5);
legend({'','仿真数据','实验数据'});
box on;grid on;xlabel('\it\fontname{times new roman}x');ylabel('\fontname{times new roman}{\itC}_P');
fig_hdl.set('Position',[200,350,380,320]);

%% print

% print(fig_hdl,'NACA0012_CP.emf', '-dmeta','-r600');
% print(fig_hdl,'NACA0012_CP.png', '-dpng','-r1200');

%% function

function [XL,CPL,XU,CPU]=splitCp(data)
X=data(:,1);
CP=data(:,2);
num=length(X)/2;
XL=X(1:num);CPL=CP(1:num);
XU=X(num+1:end);CPU=CP(num+1:end);
end
