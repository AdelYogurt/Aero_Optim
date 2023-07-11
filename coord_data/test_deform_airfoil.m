clc;
clear;
close all hidden

%% base parameter
% control_point_up=importdata('airfoil_NACA0012_shape_low.txt');
% control_point_low=importdata('airfoil_NACA0012_shape_up.txt');

control_point_up=importdata('airfoil_RAE2822_shape_low.txt');
control_point_low=importdata('airfoil_RAE2822_shape_up.txt');
mid=[control_point_up(2:end-1,2);control_point_low(2:end-1,2)]';

variable_num=10;
low_bou=mid-0.05;
up_bou=mid+0.05;
% x=rand(1,variable_num).*(up_bou-low_bou)+low_bou;
% x=(up_bou+low_bou)/2;
% x=up_bou;
% x=low_bou;
x=[ones(1,5)*0,ones(1,5)*0.1];

total_point_list=prePoint(x,'');

%% draw coord

for surface_index=1:length(total_point_list)
    point_list=total_point_list{surface_index};
    line(point_list(:,2),point_list(:,3),'marker','.','linestyle','none')
end
axis equal

%% write coord

file_out=fopen('airfoil_deform.dat','w');
for surface_index=1:length(total_point_list)
    point_list=total_point_list{surface_index};

    index=point_list(:,1);
    X=point_list(:,2);
    Y=point_list(:,3);

    for point_index=1:size(index,1)
        fprintf(file_out,'%d %f %f\n',index(point_index)-1,X(point_index),Y(point_index));
    end
end
fclose(file_out);


%% function 
function total_point_list=prePoint(x,coord_dir)
control_number=length(x)/2;
class_par_Y=[0.5,1];

%% gen object
control_point_up=importdata(fullfile(coord_dir,'airfoil_NACA0012_shape_low.txt'));
control_point_low=importdata(fullfile(coord_dir,'airfoil_NACA0012_shape_up.txt'));

control_point_up(2:end-1,2)=x(1:control_number)';
control_point_low(2:end-1,2)=x(control_number+1:end)';

bezier_up=CurveBezier(control_point_up);
bezier_low=CurveBezier(control_point_low);
shape_fun_up=@(x) shapeFunctionCurve(bezier_up,x);
shape_fun_low=@(x) shapeFunctionCurve(bezier_low,x);

airfoil_up=CST2DCurve(1,1,shape_fun_up,class_par_Y);
airfoil_low=CST2DCurve(1,-1,shape_fun_low,class_par_Y);

airfoil_coord_up=importdata(fullfile(coord_dir,'airfoil_up_local_coord.txt'));
airfoil_coord_low=importdata(fullfile(coord_dir,'airfoil_low_local_coord.txt'));

total_point_list=cell(2,1);
[X,Y]=airfoil_up.calPoint(airfoil_coord_up(:,2)');
total_point_list{1}=[airfoil_coord_up(:,1),X',Y'];
[X,Y]=airfoil_low.calPoint(airfoil_coord_low(:,2)');
total_point_list{2}=[airfoil_coord_low(:,1),X',Y'];
end

function s=shapeFunctionCurve(curve,u)
s=curve.interpPoint(u);
s=s(:,2)';
end

