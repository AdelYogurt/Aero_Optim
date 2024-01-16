clc;
% clear;
close all hidden;

load('optim.mat')

obj_list=[];
for idx=1:12
    datalib=data_list{idx}.datalib;
    obj_list=[obj_list,datalib.Obj(datalib.Best_idx)];
end

obj_list_phy=[];
for idx=12:24
    datalib=data_list{idx}.datalib;
    obj_list_phy=[obj_list_phy,datalib.Obj(datalib.Best_idx)];
end

aver=mean(obj_list,2);
aver_phy=mean(obj_list_phy,2);

NFE=(1:50)';

%% draw data

fig_hdl=figure();
axe=axes(fig_hdl);

SNI=10;

hold on;box on;grid on;
X_fill=[NFE(SNI:end); flipud(NFE(SNI:end))];
Y_fill=[max(obj_list(SNI:end,:),[],2); flipud(min(obj_list(SNI:end,:),[],2))];

fill_best=fill(X_fill,Y_fill,[0.9290 0.6940 0.1250],'edgealpha', 0, 'facealpha', 0.2);
line_best=line(NFE(SNI:end),aver(SNI:end),'Color',[0.9290 0.6940 0.1250],'LineWidth',1);

X_fill=[NFE(SNI:end); flipud(NFE(SNI:end))];
Y_fill=[max(obj_list_phy(SNI:end,:),[],2); flipud(min(obj_list_phy(SNI:end,:),[],2))];

fill_phy=fill(X_fill,Y_fill,[0.6350 0.0780 0.1840],'edgealpha', 0, 'facealpha', 0.2);
line_phy=line(NFE(SNI:end),aver_phy(SNI:end),'Color',[0.6350 0.0780 0.1840],'LineWidth',1);

line_base=line([0,50],[0.0580,0.0580],'Color',[0.3789,0.4531,0.4961],'lineStyle','--','LineWidth',1);

ylim([0.012,0.024]);
ylabel('\fontname{宋体}阻力系数\fontname{times new roman}C_{D}')

patch_lhd=patch('XData',[0,SNI,SNI,0],'YData',[-10,-10,20,20],'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.5,'LineStyle','none');

axe.set('Parent',fig_hdl,'position',[0.18,0.17,0.78,0.76],'FontSize',10.5,'GridColor',[0 0.4470 0.7410],'GridAlpha',0.8);
axe.GridAlpha=0.2;

lgd_hdl=legend([line_best,line_phy,patch_lhd],...
    {'\fontname{times new roman}EGO',...
    '\fontname{times new roman}PA-EGO','\fontname{宋体}初始采样'},'Location','northeast');
lgd_hdl.set('Parent',fig_hdl,'FontSize',10.5,'Box','on');
xlabel('\fontname{宋体}模型调用次数');

xlim([0,50]);
box on;
xlabel('模型调用次数')

fig_hdl.set('Position',[200,400,420,320]);

% print(fig_hdl,'optim_conv.emf','-dmeta','-r600');
% print(fig_hdl,'optim_conv.png','-dpng','-r1200');
