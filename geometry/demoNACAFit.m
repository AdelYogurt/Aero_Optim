clc;
clear;
close all hidden;

addpath('data\')
addpath('geometry\')

%% initialize

class_par_Y=[0.5,1];
NC=calNormPar(class_par_Y(1),class_par_Y(2));
class_fun=@(x) typicalFunction(x,class_par_Y(1),class_par_Y(2),NC);

airfoil_data=importdata('RAE2822.txt');
airfoil_data_up=airfoil_data(1:65,:);
airfoil_data_low=airfoil_data(66:end,:);

%% calculate node_point to fit shape_fun
node_point_up=airfoil_data_up;
node_point_low=airfoil_data_low;
node_point_up(:,2)=divedeSingular(airfoil_data_up(:,2),class_fun(airfoil_data_up(:,1)),airfoil_data_up(:,1));
node_point_low(:,2)=divedeSingular(airfoil_data_low(:,2),-class_fun(airfoil_data_low(:,1)),airfoil_data_low(:,1));

% line(node_point_up(:,1),node_point_up(:,2),'Color','r');
% line(node_point_low(:,1),node_point_low(:,2),'Color','b');

%% fit with Bezier

max_order=5;
bezier_up=CurveBezier(node_point_up,max_order,node_point_up(:,1));
bezier_low=CurveBezier(node_point_low,max_order,node_point_low(:,1));
shape_fun_up=@(x) shapeFunctionCurve(bezier_up,x);
shape_fun_low=@(x) shapeFunctionCurve(bezier_low,x);

%% fit with BSpline

% bezier_up=CurveBSpline(node_point_up);
% bezier_low=CurveBSpline(node_point_low);
% shape_fun_up=@(x) shapeFunctionCurve(bezier_up,x);
% shape_fun_low=@(x) shapeFunctionCurve(bezier_low,x);

%% draw
CST2D_curve_up=CST2DCurve(1,1,shape_fun_up,class_par_Y);
CST2D_curve_low=CST2DCurve(1,-1,shape_fun_low,class_par_Y);

[draw_X,draw_Y_up]=CST2D_curve_up.calCurve(50);
[draw_X,draw_Y_low]=CST2D_curve_low.calCurve(50);

line(airfoil_data_up(:,1),airfoil_data_up(:,2),'Color','r');
line(airfoil_data_low(:,1),airfoil_data_low(:,2),'Color','r');
line(draw_X,draw_Y_up,'Color','b');
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