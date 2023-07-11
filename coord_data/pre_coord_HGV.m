clc;
clear;
close all hidden

%% base parameter

total_length=4; % total length

vari_num=16;
low_bou=[0.7,0.7, ...
    2.4,0.4,0.5,0.5, ...
    0.1,0.1,0.005, ...
    0.4,0.4,0.4,0.6,0.5,0.01,0.01];
up_bou=[1.0,1.0, ...
    3.2,0.8,2.0,2.0, ...
    0.5,0.5,0.02, ...
    0.6,0.6,0.6,0.8,0.7,0.05,0.05];

x=(up_bou+low_bou)/2;

% length parameter
par_M_up=x(1);
par_M_low=x(2);
% width parameter
par_width=x(3);
par_T=x(4);
par_N_up=x(5);
par_N_low=x(6);
% height parameter
par_hight_up=x(7);
par_hight_low=x(8);
par_R=x(9);
% wing parameter
par_rho1=x(10);
par_rho12=x(11);
par_rho23=x(12);
par_WS1=x(13);
par_WS2=x(14);
par_TWU=x(15);
par_TWL=x(16);

%% gen object

waverider_wing=WaveriderWingDia...
    (total_length,par_width,par_hight_up,par_hight_low,...
    par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
    par_rho1,par_rho12,par_rho23,par_WS1,par_WS2,par_TWU,par_TWL);

y_cut=waverider_wing.y_cut;
head_length=waverider_wing.head_length;

%% load

load('mesh_data.mat');
head_side_index=marker_index_list.HEAD_SIDE;
waverider_index=marker_index_list.WAVERIDER;
tri_wing_front_index=marker_index_list.TRI_WING_FRONT;
tri_wing_index=marker_index_list.TRI_WING;
wing_index=marker_index_list.WING;
wing_front_index=marker_index_list.WING_FRONT;
wing_side_index=marker_index_list.WING_SIDE;

head_side_point=point_list(head_side_index,:);
waverider_point=point_list(waverider_index,:);
tri_wing_front_point=point_list(tri_wing_front_index,:);
tri_wing_point=point_list(tri_wing_index,:);
wing_point=point_list(wing_index,:);
wing_front_point=point_list(wing_front_index,:);
wing_side_point=point_list(wing_side_index,:);

% fine grid
xi_grid_num_head=40; % head x direction gird num
eta_grid_num_head=20; % head and body y direction gird num
xi_grid_num_body=20; % body x direction gird num
eta_grid_num_wing=8; % wing y direction gird num
edge_gird_num=6; % edge gird num

%% analysis point
[wing_side_index,wing_side_point,wing_side_front_index,wing_side_front_point]=getFrobac(wing_side_index,wing_side_point,4-4*par_rho1*par_rho12*par_rho23-1e-9);
[wing_side_up_index,wing_side_up_point,wing_side_low_index,wing_side_low_point]=getUpLow(wing_side_index,wing_side_point);

[tri_wing_back_index,tri_wing_back_point,tri_wing_index,tri_wing_point]=getFrobac(tri_wing_index,tri_wing_point,4-1e-9);
[tri_wing_up_index,tri_wing_up_point,tri_wing_low_index,tri_wing_low_point]=getUpLow(tri_wing_index,tri_wing_point);
[tri_wing_back_up_index,tri_wing_back_up_point,tri_wing_back_low_index,tri_wing_back_low_point]=getUpLow(tri_wing_back_index,tri_wing_back_point);

[wing_back_index,wing_back_point,wing_index,wing_point]=getFrobac(wing_index,wing_point,4-1e-9);
[wing_up_index,wing_up_point,wing_low_index,wing_low_point]=getUpLow(wing_index,wing_point);
[wing_back_up_index,wing_back_up_point,wing_back_low_index,wing_back_low_point]=getUpLow(wing_back_index,wing_back_point);

[waverider_back_index,waverider_back_point,head_index,head_point]=getFrobac(waverider_index,waverider_point,head_length+1e-9);
[body_back_index,body_back_point,body_index,body_point]=getFrobac(waverider_back_index,waverider_back_point,4-1e-9);
[head_up_index,head_up_point,head_low_index,head_low_point]=getUpLow(head_index,head_point);
[body_up_index,body_up_point,body_low_index,body_low_point]=getUpLow(body_index,body_point);
[body_back_up_index,body_back_up_point,body_back_low_index,body_back_low_point]=getUpLow(body_back_index,body_back_point);

% hold on;
% scatter3(head_up_point(:,1),head_up_point(:,2),head_up_point(:,3))
% scatter3(head_low_point(:,1),head_low_point(:,2),head_low_point(:,3))
% scatter3(head_side_point(:,1),head_side_point(:,2),head_side_point(:,3))
% scatter3(body_up_point(:,1),body_up_point(:,2),body_up_point(:,3))
% scatter3(body_low_point(:,1),body_low_point(:,2),body_low_point(:,3))
% hold off;
% axis equal;
% view(3);

%% rebuild coord

coord_data=struct();
head_xi_list=linspace(0,1,xi_grid_num_head+1);
head_xi_list=head_xi_list.^(1/par_T);
coord_data=calCoord(coord_data,waverider_wing,'head_up',head_up_point,head_up_index,head_xi_list,eta_grid_num_head);
coord_data=calCoord(coord_data,waverider_wing,'head_low',head_low_point,head_low_index,head_xi_list,eta_grid_num_head);
coord_data=calCoord(coord_data,waverider_wing,'head_side',head_side_point,head_side_index,fliplr(1-head_xi_list),2*edge_gird_num);
coord_data=calCoord(coord_data,waverider_wing,'body_up',body_up_point,body_up_index,xi_grid_num_body,eta_grid_num_head);
coord_data=calCoord(coord_data,waverider_wing,'body_low',body_low_point,body_low_index,xi_grid_num_body,eta_grid_num_head);
coord_data=calCoord(coord_data,waverider_wing,'body_back_up',body_back_up_point,body_back_up_index,eta_grid_num_head,edge_gird_num);
coord_data=calCoord(coord_data,waverider_wing,'body_back_low',body_back_low_point,body_back_low_index,eta_grid_num_head,edge_gird_num);

coord_data=calCoord(coord_data,waverider_wing,'tri_wing_up',tri_wing_up_point,tri_wing_up_index,eta_grid_num_wing,xi_grid_num_body);
coord_data=calCoord(coord_data,waverider_wing,'tri_wing_low',tri_wing_low_point,tri_wing_low_index,eta_grid_num_wing,xi_grid_num_body);
coord_data=calCoord(coord_data,waverider_wing,'tri_wing_back_up',tri_wing_back_up_point,tri_wing_back_up_index,eta_grid_num_wing,edge_gird_num);
coord_data=calCoord(coord_data,waverider_wing,'tri_wing_back_low',tri_wing_back_low_point,tri_wing_back_low_index,eta_grid_num_wing,edge_gird_num);
coord_data=calCoord(coord_data,waverider_wing,'tri_wing_front',tri_wing_front_point,tri_wing_front_index,eta_grid_num_wing,2*edge_gird_num);

coord_data=calCoord(coord_data,waverider_wing,'wing_up',wing_up_point,wing_up_index,eta_grid_num_wing,xi_grid_num_body);
coord_data=calCoord(coord_data,waverider_wing,'wing_low',wing_low_point,wing_low_index,eta_grid_num_wing,xi_grid_num_body);
coord_data=calCoord(coord_data,waverider_wing,'wing_back_up',wing_back_up_point,wing_back_up_index,eta_grid_num_wing,edge_gird_num);
coord_data=calCoord(coord_data,waverider_wing,'wing_back_low',wing_back_low_point,wing_back_low_index,eta_grid_num_wing,edge_gird_num);
coord_data=calCoord(coord_data,waverider_wing,'wing_front',wing_front_point,wing_front_index,eta_grid_num_wing,2*edge_gird_num);
coord_data=calCoord(coord_data,waverider_wing,'wing_side_up',wing_side_up_point,wing_side_up_index,eta_grid_num_wing,edge_gird_num);
coord_data=calCoord(coord_data,waverider_wing,'wing_side_low',wing_side_low_point,wing_side_low_index,eta_grid_num_wing,edge_gird_num);
coord_data=calCoord(coord_data,waverider_wing,'wing_side_front',wing_side_front_point,wing_side_front_index,edge_gird_num,edge_gird_num);

% % different process of stag and head_side
% [XI_stag,PSI_stag]=waverider_wing.stag.calCoordinate(stag_point(:,1),stag_point(:,2),stag_point(:,3));
% [XI_head_side,PSI_head_side]=waverider_wing.surface_list.head_side.calCoordinate(head_side_point(:,1),head_side_point(:,2),head_side_point(:,3));
% index=[stag_index;head_side_index];
% XI=[XI_head_side;(0.01-(XI_stag*waverider_wing.stag.LX*2/par_wight).^(1/par_T)/(1-par_rho1))+0.99];
% PSI=[PSI_head_side;PSI_stag];
% 
% fid=fopen(['coord_data/','head_side','_local_coord.txt'],'w');
% for point_index=1:size(index,1)
%     fprintf(fid,'%d %12f %12f\n',index(point_index),XI(point_index),PSI(point_index));
% end
% fclose(fid);

save('coord_data.mat','coord_data');
% writeCoord(coord_data);

%% function 

function [index_up,point_up,index_low,point_low]=getUpLow(index,point)
Bool=point(:,3)>1e-9;

index_up=index(Bool,:);
point_up=point(Bool,:);
index_low=index(~Bool,:);
point_low=point(~Bool,:);

end

function [index_front,point_front,index_back,point_back]=getFrobac(index,point,x)
Bool=point(:,1)>x;

index_front=index(Bool,:);
point_front=point(Bool,:);
index_back=index(~Bool,:);
point_back=point(~Bool,:);

end

function writeCoord(coord_data)
%
%
surface_name_list=fieldnames(coord_data);

for surface_idx=1:length(surface_name_list)
    surface_name=surface_name_list{surface_idx};
    index=coord_data.(surface_name).index;
    XI=coord_data.(surface_name).XI;
    PSI=coord_data.(surface_name).PSI;
    fid=fopen([surface_name,'_local_coord.txt'],'w');
    for point_index=1:length(index)
        fprintf(fid,'%d %.12f %.12f\n',index(point_index),XI(point_index),PSI(point_index));
    end
    fclose(fid);
end
end

function coord_data=calCoord(coord_data,body,surface_name,point,index,XI_bou,PSI_bou)
%
%
if strcmp(surface_name,'head_side')
    [XI,PSI]=body.surface_list.(surface_name).calCoordinate(point(:,1),point(:,2),point(:,3));
    %     large_index=find(XI>1);
    %     xi_large=XI(large_index);
    %     xi_large=xi_large-1;
    %     psi_large=PSI(large_index);
    %     scatter(xi_large,psi_large)
    %     text(xi_large,psi_large,num2str(large_index),'FontSize',6)

    %     sum_z_sq=sum(xi_large.^2);
    %     sum_z_sq2=sum(xi_large.^4);
    %     fit_equal=[1,sum_z_sq;sum_z_sq,sum_z_sq2];
    %     c=fit_equal\[sum(xi_large.^2);sum(psi_large.^2.*xi_large.^2)];
    %     xi_edge_func=@(psi) c(1)+c(2)*psi.^2;

    %     symmetry_index=importdata("symmetry.txt");
    %     [sam_index,large_index]=intersect(index,symmetry_index);
    %     xi_large=XI(large_index);
    %     xi_large=xi_large-1;
    %     psi_large=PSI(large_index);
    %     xi_edge_func=@(Z) interpLinear(Z,psi_large,xi_large);
    %
    %     XI(large_index)=1;

    %     XI=XI-XI.*xi_edge_func(PSI);

    scatter(XI,PSI);
else
    [XI,PSI]=body.surface_list.(surface_name).calCoordinate(point(:,1),point(:,2),point(:,3));
end

coord_data.(surface_name).index=index;
coord_data.(surface_name).XI=XI;
coord_data.(surface_name).PSI=PSI;
end

function Y_pred=interpLinear(X_pred,X,Y)
[X,index]=sort(X);
Y=Y(index);
Y_pred=zeros(length(X_pred),1);
for x_index=1:length(X_pred)
    x_pred=X_pred(x_index);
    num=length(X);
    index=num; % search start from last one, find out X samll than x
    while ((index > 1) && (X(index) > x_pred))
        index=index-1;
    end

    if (index == num)
        % out interp
        Y_pred(x_index)=(Y(end)-Y(end-1))/(X(end)-X(end-1))*(x_pred-X(end))+Y(end);
    elseif (index == 0)
        Y_pred(x_index)=(Y(2)-Y(1))/(X(2)-X(1))*(x_pred-X(1))+Y(1);
    else
        % linear interpolation
        Y_pred(x_index)=Y(index)+...
            (Y(index+1)-Y(index))*...
            (x_pred-X(index))/...
            (X(index+1)-X(index));
    end
end

end

