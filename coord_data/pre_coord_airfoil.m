clc;
clear;
close all hidden

class_par_Y=[0.5,1];

%% gen object
control_point_up=importdata('airfoil_NACA0012_shape_up.txt');
control_point_low=importdata('airfoil_NACA0012_shape_low.txt');

bezier_up=CurveBezier(control_point_up);
bezier_low=CurveBezier(control_point_low);
shape_fun_up=@(x) shapeFunctionCurve(bezier_up,x);
shape_fun_low=@(x) shapeFunctionCurve(bezier_low,x);

airfoil_up=CST2DCurve(1,1,shape_fun_up,class_par_Y);
airfoil_low=CST2DCurve(1,-1,shape_fun_low,class_par_Y);

%% load data
load('mesh_data.mat');
point_idx_up=marker_index_list.AIRFOIL_UP;
point_idx_low=marker_index_list.AIRFOIL_LOW;

airfoil_mesh_up=point_list(point_idx_up,1:2);
airfoil_mesh_low=point_list(point_idx_low,1:2);

%% rebuild coord
writeCalCoord(airfoil_up,'airfoil_up',airfoil_mesh_up,point_idx_up);
writeCalCoord(airfoil_low,'airfoil_low',airfoil_mesh_low,point_idx_low);

%% function 
function XI=writeCalCoord(curve,curve_name,point,index)
XI=curve.calCoordinate(point(:,1));

fid=fopen([curve_name,'_local_coord.txt'],'w');
for point_index=1:size(point,1)
    fprintf(fid,'%d,%.12f\n',index(point_index),XI(point_index));
end
fclose(fid);
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

function s=shapeFunctionCurve(curve,u)
s=curve.interpPoint(u);
s=s(:,2)';
end

function result=divedeSingular(y,c,x)
num__=length(y);
if num__ ~= length(c)
    error('divedeSingular: A, B size error');
end

boolean=c == 0;
c(boolean)=1;

result=y./c;

% interpolation by close point if B is zero
for i=1:num__
    if boolean(i)
        if i==1
            result(i)=result(i+1)-(x(i+1)-x(i))*(result(i+2)-result(i+1))/(x(i+2)-x(i+1));
        elseif i == num__
            result(i)=result(i-1)+(x(i)-x(i-1))*(result(i-1)-result(i-2))/(x(i-1)-x(i-2));
        else
            result(i)=result(i-1)+(x(i)-x(i-1))*(result(i+1)-result(i-1))/(x(i+1)-x(i-1));
        end
    end
end

end

function c=typicalFunction(x,N1,N2,NC)
c=x.^N1.*(1-x).^N2/NC;
end

function nomlz_par=calNormPar(N1,N2)
% calculate normailize class function parameter by N1, N2
%
if N1 == 0 && N2 == 0
    nomlz_par=1;
else
    nomlz_par=(N1/(N1+N2))^N1*(N2/(N1+N2))^N2;
end
end
