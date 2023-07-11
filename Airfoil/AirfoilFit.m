clc;
clear;
close all hidden;

addpath('data\')
addpath('geometry\')

class_par_Y=[0.5,1];
NC=calNormPar(class_par_Y(1),class_par_Y(2));
class_fun=@(x) typicalFunction(x,class_par_Y(1),class_par_Y(2),NC);

max_order=6;

%% initialize, calculate node_point to fit shape_fun
figure(1);
airfoil_data=importdata('NACA0012.txt');
airfoil_data_up=airfoil_data(1:67,1:2);
airfoil_data_low=airfoil_data(68:end,1:2);

node_point_up=airfoil_data_up;
node_point_low=airfoil_data_low;
node_point_up(:,2)=divedeSingular(airfoil_data_up(:,2),class_fun(airfoil_data_up(:,1)),airfoil_data_up(:,1));
node_point_low(:,2)=divedeSingular(airfoil_data_low(:,2),-class_fun(airfoil_data_low(:,1)),airfoil_data_low(:,1));
node_point_up(end,2)=0.0552;
node_point_low(end,2)=0.0552;

line(node_point_up(:,1),node_point_up(:,2),'Color','r');
line(node_point_low(:,1),node_point_low(:,2),'Color','b');

% airfoil_data=importdata('RAE2822.txt');
% airfoil_data_up=airfoil_data(1:65,1:2);
% airfoil_data_low=airfoil_data(66:end,1:2);
% 
% node_point_up=airfoil_data_up;
% node_point_low=airfoil_data_low;
% node_point_up(:,2)=divedeSingular(airfoil_data_up(:,2),class_fun(airfoil_data_up(:,1)),airfoil_data_up(:,1));
% node_point_low(:,2)=divedeSingular(airfoil_data_low(:,2),-class_fun(airfoil_data_low(:,1)),airfoil_data_low(:,1));
% 
% line(node_point_up(:,1),node_point_up(:,2),'Color','r','LineStyle','--');
% line(node_point_low(:,1),node_point_low(:,2),'Color','b','LineStyle','--');

%% fit with Bezier

bezier_up=CurveBezier(node_point_up,max_order,node_point_up(:,1));
bezier_low=CurveBezier(node_point_low,max_order,node_point_low(:,1));

% up_bou=ones(1,8)*0.15;
% low_bou=ones(1,8)*-0.05;
% control_point_up=rand(1,8).*(up_bou-low_bou)+low_bou;
% control_point_up=[linspace(0,1,8)',control_point_up'];
% control_point_low=rand(1,8).*(up_bou-low_bou)+low_bou;
% control_point_low=[linspace(0,1,8)',control_point_low'];
% bezier_up=CurveBezier(control_point_up);
% bezier_low=CurveBezier(control_point_low);
% 
% shape_fun_up=@(x) shapeFunctionCurve(bezier_up,x);
% shape_fun_low=@(x) shapeFunctionCurve(bezier_low,x);

%% write control point coordinate
writeCoord(bezier_up,'airfoil_shape_up');
writeCoord(bezier_low,'airfoil_shape_low');

%% fit with BSpline

bezier_up=CurveBSpline(node_point_up);
bezier_low=CurveBSpline(node_point_low);
shape_fun_up=@(x) shapeFunctionCurve(bezier_up,x);
shape_fun_low=@(x) shapeFunctionCurve(bezier_low,x);

%% draw
CST2D_curve_up=CST2DCurve(1,1,shape_fun_up,class_par_Y);
CST2D_curve_low=CST2DCurve(1,-1,shape_fun_low,class_par_Y);

[draw_X,draw_Y_up]=CST2D_curve_up.calCurve(50);
[draw_X,draw_Y_low]=CST2D_curve_low.calCurve(50);

figure(2);
line(airfoil_data_up(:,1),airfoil_data_up(:,2),'Color','r','LineStyle','--');
line(airfoil_data_low(:,1),airfoil_data_low(:,2),'Color','b','LineStyle','--');

line(draw_X,draw_Y_up,'Color','r');
line(draw_X,draw_Y_low,'Color','b');

%% function
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

function writeCoord(curve,curve_name)
control_point=curve.control_list;

fid=fopen(['coord_data/',curve_name,'.txt'],'w');
for point_index=1:size(control_point,1)
    fprintf(fid,'%.12f,%.12f\n',control_point(point_index,1),control_point(point_index,2));
end
fclose(fid);
end
