clc
% clear;
close all hidden;

cut_idx=200;

load('data_lib_5.mat')

% load('Airfoil\optres_SADE.mat','output')
% load('Airfoil\optres_SACORS.mat','output')
% load('Airfoil\optres_FSRBF.mat','output')
% data_lib=output.data_lib;

Obj=data_lib.Obj;
Vio=data_lib.Vio;

cut_idx=min(cut_idx,length(Obj));

result_best_idx=[];
for x_idx=1:length(Obj)
    if isempty(result_best_idx)
        result_best_idx=x_idx;
    else
        vio=Vio(x_idx);
        obj=Obj(x_idx);
        if vio == 0
            if obj < Obj(result_best_idx(end))
                result_best_idx=[result_best_idx;x_idx];
            else
                result_best_idx=[result_best_idx;result_best_idx(end)];
            end
        else
            if vio < Vio(result_best_idx(end))
                result_best_idx=[result_best_idx;x_idx];
            else
                result_best_idx=[result_best_idx;result_best_idx(end)];
            end
        end
    end
end

Obj_best=data_lib.Obj(result_best_idx(1:cut_idx));
Vio_best=data_lib.Vio(result_best_idx(1:cut_idx));
Con_best=data_lib.Con(result_best_idx(1:cut_idx),:);
X_best=data_lib.X(result_best_idx(1:cut_idx),:);

hold on 
yyaxis left
line(1:length(Obj_best),-Obj_best,'linestyle','--','marker','*','MarkerIndices',5:10:195)
line(1:length(Obj),-Obj,'linestyle','--','marker','.','Color','g')
% ylim([1.6,3.5]);
ylabel('升阻比L/D')

yyaxis right;
line(1:length(Vio_best),Vio_best,'linestyle','--','marker','*','Color','r','MarkerIndices',5:10:195)
line(1:length(Vio),Vio,'linestyle','--','marker','.','Color','k')
% ylim([-0.01,0.1]);
ylabel('约束违背度')

hold off

grid on;
xlabel('模型调用次数')

legend('SACO-RS目标函数值','SACO-RS约束违背度','Location','northwest')

% exportgraphics(gcf,'opt_res.emf','Resolution',600)

