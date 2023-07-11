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

% x=rand(1,variable_num).*(up_bou-low_bou)+low_bou;
x=(up_bou+low_bou)/2;
% x=up_bou;
% x=low_bou;

% load('LHD.mat')
% x=X_HF(1,:);

coord_data=proPoint(x);

%% draw coord

surface_name_list=fieldnames(coord_data);
hold on;
for surface_index=1:length(surface_name_list)
    surface_name=surface_name_list{surface_index};
    scatter3(coord_data.(surface_name).X,coord_data.(surface_name).Y,coord_data.(surface_name).Z)
end
hold off
view(3);
axis equal

%% write coord

% surface_name_list=fieldnames(coord_data);
% fid=fopen('WWDB_deform.dat','w');
% for surface_index=1:length(surface_name_list)
%     surface_name=surface_name_list{surface_index};
% 
%     index=coord_data.(surface_name).index;
%     X=coord_data.(surface_name).X;
%     Y=coord_data.(surface_name).Y;
%     Z=coord_data.(surface_name).Z;
% 
%     for point_index=1:size(index,1)
%         fprintf(fid,'%d %f %f %f\n',index(point_index)-1,X(point_index),Y(point_index),Z(point_index));
%     end
% end
% fclose(fid);


%% function 

function coord_data=proPoint(x)
total_length=4; % total length

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

waverider_wing=WaveriderWingDia...
    (total_length,par_width,par_hight_up,par_hight_low,...
    par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
    par_rho1,par_rho12,par_rho23,par_WS1,par_WS2,par_TWU,par_TWL);

surface_name_list=fieldnames(waverider_wing.surface_list);

stag_deform_Y_up=@(Y,Z) (min((Y/(3*par_R)).^2+((Z-par_R)/(par_R)).^2,1)).^(1-par_T/0.7);
stag_deform_Y_low=@(Y,Z) (min((Y/(3*par_R)).^2+((Z+par_R)/(par_R)).^2,1)).^(1-par_T/0.7);
stag_deform_Y_mid=@(Y,Z) (min((Y/(3*par_R)).^2+((Z)/(par_R)).^2,1)).^(1-par_T*1/0.7);

coord_data=struct();

load("coord_data.mat",'coord_data');
for surface_index=1:length(surface_name_list)
    surface_name=surface_name_list{surface_index};
    %     data_filename=[surface_name,'_local_coord.txt'];
    %     data=importdata(data_filename);
    %     index=data(:,1);XI=data(:,2);PSI=data(:,3);
    %     coord_data.(surface_name).index=index;
    %     coord_data.(surface_name).XI=XI;
    %     coord_data.(surface_name).PSI=PSI;

    XI=coord_data.(surface_name).XI;
    PSI=coord_data.(surface_name).PSI;
    [X,Y,Z]=waverider_wing.surface_list.(surface_name).calPoint(XI,PSI);
    
    if strcmp(surface_name,'head_side')

    end
%     if strcmp(surface_name,'head_up') || strcmp(surface_name,'head_low')
%         % deform of stag
%         Y=Y.*stag_deform_Y_up(Y,Z);
%         Y=Y.*stag_deform_Y_low(Y,Z);
%         Y=Y.*stag_deform_Y_mid(Y,Z);
%     end

    coord_data.(surface_name).X=X;
    coord_data.(surface_name).Y=Y;
    coord_data.(surface_name).Z=Z;
end

% % calculate error
% [index_equal,posi_base,posi_check]=intersect(point_stag(:,1),point_head_up(:,1));
% 
% x=point_head_up(posi_check,2);
% dx=point_head_up(posi_check,2)-point_stag(posi_base,2);
% y=point_head_up(posi_check,3);
% dy=point_head_up(posi_check,3)-point_stag(posi_base,3);
% 
% [x,index]=sort(x);dx=dx(index);
% [y,index]=sort(y);dy=dy(index);
% 
% for point_index=1:size(point_stag,1)
%     y_stag=coord_stag(point_index,1)*waverider_wing.stag_length;
%     x_stag=total_length*(coord_stag(point_index,1)*waverider_wing.stag_length/(par_W/2)).^(1/par_T);
%     point_stag(point_index,2)=point_stag(point_index,2)+interpLinear(x_stag,x,dx);
%     point_stag(point_index,3)=point_stag(point_index,3)+interpLinear(y_stag,y,dy);
% end
% 
% total_point_list{index_stag}=point_stag;


% % process repeat
% for base_index=1:length(surface_name_list)
%     point_base=total_point_list{base_index};
%     index_base=point_base(:,1);
%     for check_index=base_index+1:length(surface_name_list)
%         point_check=total_point_list{check_index};
%         index_check=point_check(:,1);
% 
%         [index_equal,posi_base,posi_check]=intersect(index_base,index_check);
%         if ~isempty(index_equal)
%             point_equal=(point_base(posi_base,2:4)+point_check(posi_check,2:4))/2;
%             point_base(posi_base,2:4)=point_equal;
%             point_check(posi_check,2:4)=point_equal;
%         end
%         total_point_list{check_index}=point_check;
%     end
%     total_point_list{base_index}=point_base;
% end

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
