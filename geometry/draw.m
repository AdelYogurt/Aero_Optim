clc;
clear;
close all hidden;

LX=2;
torl=1e-3;
N_list=[0.5,0.8,1,1.25,2.0];
color=[
    [0 0.4470 0.7410];
    [0.8500 0.3250 0.0980];
    [0.9290 0.6940 0.1250];
    [0.4940 0.1840 0.5560];
    [0.4660 0.6740 0.1880];
    [0.3010 0.7450 0.9330];
    [0.6350 0.0780 0.1840];
    ];


fig_hdl=figure();
fig_hdl.Position=[200,400,500,300];

axe_hdl_1=axes(fig_hdl);
hold on;
for idx=[1,5,2,4,3]
    N=N_list(idx);
    fcn=@(X) baseFcnClass(X/LX,N,N);
    fplot(fcn,[0,LX],'Color',color(idx,:),'LineWidth',1);
    if N < 1,u=0;v=-0.5;
    elseif N == 1,u=-1/sqrt(5)*0.5;v=-2/sqrt(5)*0.5;
    else,u=-0.5;v=0;end
    quiver(0,0,u,v,'Color',color(idx,:),'LineWidth',1,'MaxHeadSize',0.5);
end
hold off;
box on;grid on;
xlim([-0.6,0.6]);ylim([-0.6,0.6]);
axe_hdl_1.Position=[0.11,0.24,0.38,0.60];
axe_hdl_1.FontSize=10.5;
axe_hdl_1.YTick=[-0.5,0,0.5];
xlabel('\fontname{times new roman}{\itX}');
ylabel('\fontname{times new roman}{\itY}');

axe_hdl_2=axes(fig_hdl);
hold on;
for idx=[1,5,2,4,3]
    N=N_list(idx);
    fcn=@(X) baseFcnClass(X/LX,N,N);
    f_hdl(idx)=fplot(fcn,[0,LX],'Color',color(idx,:),'LineWidth',1);
    norm_par=(N./(N+N)).^N.*(N./(N+N)).^N;
    % [k,x]=calTangTorl(N,torl/N,LX,1/norm_par);
    crv=EdgeCST2D('',[N,N],false,LX);
    [k,~,x]=crv.calTangTorl(torl/N);
    u=-1/sqrt(1+k^2)*0.5;v=-k/sqrt(1+k^2)*0.5;
    quiver(0,0,u,v,'Color',color(idx,:),'LineWidth',1,'MaxHeadSize',0.5);
end
hold off;
box on;grid on;

xlim([-0.6,0.6]);ylim([-0.6,0.6]);
axe_hdl_2.Position=[0.60,0.24,0.38,0.60];
axe_hdl_2.FontSize=10.5;
axe_hdl_2.YTick=[-0.5,0,0.5];
xlabel('\fontname{times new roman}{\itX}');
ylabel('\fontname{times new roman}{\itY}');

txt_hdl_1=text(axe_hdl_1,-0.3,-1.0,0,'\fontname{Times new roman}(a)  \fontname{宋体}解析导矢','FontSize',10.5);
txt_hdl_2=text(axe_hdl_2,-0.3,-1.0,0,'\fontname{Times new roman}(b)  \fontname{宋体}差分导矢','FontSize',10.5);

lgd_hdl=legend(f_hdl,...
    {'{\itN}_1=0.5','{\itN}_1=0.8','{\itN}_1=1.0','{\itN}_1=1.25','{\itN}_1=2.0'});
lgd_hdl.set('Position',[0.01,0.87,0.98,0.1],'Orientation','horizontal','FontSize',10.5)
lgd_hdl.FontName='Times new roman';
lgd_hdl.Box='off';

% print(fig_hdl,'diff_CST_grad.emf','-dmeta');

function [k,x]=calTangTorl(N,torl,KX,KY)
if N == 1,k=N*KY/KX;x=0;
else,k=N*KY/KX*(torl/KY/abs(1-N))^((N-1)/N);x=KX*(k*KX/N/KY)^(1/(N-1));end
end
