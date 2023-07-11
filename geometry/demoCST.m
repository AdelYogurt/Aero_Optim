clc;
clear;
close all hidden;

%% test 3D surface

LX=2;
LY=1;
LZ=0.3;

% shape_par_Y=[0,0];
shape_par_Y=@(XI) 1-XI*(1-0.4);
shape_par_Z=[0,0];
% class_par_Z={@(XI) sin(XI/2*pi)*15+0.001,@(XI) sin(XI/2*pi)*15+0.001};
class_par_Z={[1,1],[0.001,0.001]};

symmetry_y=false(1);
xi_gird_num=10;
psi_gird_num=10;

[XI,PSI]=meshgrid(linspace(0,1,xi_gird_num+1),linspace(0,1,psi_gird_num+1));

CST3D_surface=CST3DSurface(LX,LY,LZ,shape_par_Y,[],shape_par_Z,class_par_Z,symmetry_y);
CST3D_surface.addDeform(shape_par_Y,[],[])
CST3D_surface.addRotation(90,0,0);
CST3D_surface.addTranslation(5,6,9);
[X,Y,Z]=CST3D_surface.calSurface(XI,PSI);
[inv_XI,inv_PSI]=CST3D_surface.calCoordinate(X,Y,Z);
if any(any((inv_XI-XI) > 1e-9)) || any(any((inv_PSI-PSI) > 1e-9))
    d_XI=inv_XI-XI;
    d_PSI=inv_PSI-PSI;
end

hold on;
surf(X,Y,Z);
hold off;
xlabel('x');
ylabel('y');
zlabel('z');
view(3);
axis equal

%% test 2D surface

% LX=1;
% LY=1;
% 
% xi_gird_num=50;
% psi_gird_num=50;
% 
% [XI,PSI]=meshgrid(linspace(0,1,xi_gird_num+1),linspace(0,1,psi_gird_num+1));
% 
% class_par_X=[1,1];
% class_par_Y=[];
% 
% CST2D_surface=CST2DSurface(LX,LY,[],class_par_X,[],class_par_Y);
% [X,Y]=CST2D_surface.calSurface(XI,PSI);
% Z=zeros(size(X));
% 
% [XI_inv,PSI_inv]=CST2D_surface.calCoordinate(X,Y);
% 
% hold on;
% surf(X,Y,Z);
% hold off;
% xlabel('x');
% ylabel('y');
% zlabel('z');
% view(3);
% axis equal

%% test 2D curve

% LX=1;
% LY=1;
% 
% xi_gird_num=50;
% 
% class_par_Y=[0.3,0.1];
% 
% CST2D_curve=CST2DCurve(LX,LY,[],class_par_Y);
% [X,Y]=CST2D_curve.calCurve(xi_gird_num);
% Z=zeros(size(X));
% 
% [XI_inv]=CST2D_curve.calCoordinate(X);
% 
% line(X,Y);
% xlabel('x');
% ylabel('y');
% axis equal