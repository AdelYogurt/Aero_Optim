clc;
clear;
close all hidden;

% low_bou=0.01*ones(1,10);
% low_bou=[0.01,0.01,0.01,0.01,-0.01,0.01,0.01,0.01,0.01,-0.01];
% up_bou=0.1*ones(1,10);
% x=rand(1,10).*(up_bou-low_bou)+low_bou;
% x=low_bou;
% 
% control_point_low=[linspace(0,1,5);x(1:5)]';
% control_point_up=[linspace(0,1,5);x(6:10)]';
% 
% airfoil=AirfoilCST(control_point_low,control_point_up);
% airfoil.drawCurve();
% axis equal;


%% test 3D surface

LX=2;
LY=1;
LZ=0.3;

shape_par_X=[0,0];
shape_par_X=@(V) 1+V*(1-0.4);
% shape_par_Y=[0,0];
shape_par_Y=@(U) 1-U*(1-0.4);
shape_par_Z=[0,0];
% class_par_Z={@(U) sin(U/2*pi)*15+0.001,@(U) sin(U/2*pi)*15+0.001};
% class_par_Z={[1,1],[0.001,0.001]};
class_par_Z=[0,0];

symmetry_y=false(1);
xi_gird_num=10;
psi_gird_num=10;

[U,V]=meshgrid(linspace(0,1,xi_gird_num+1),linspace(0,1,psi_gird_num+1));

CST3D_surface=SurfaceCST3D('',LX,LY,LZ,shape_par_X,shape_par_Y,shape_par_Z,class_par_Z,symmetry_y);
% CST3D_surface.addDeform(shape_par_Y,[],[])
% CST3D_surface.addRotation(90,0,0);
CST3D_surface.addTranslation(5,6,9);
[X,Y,Z]=CST3D_surface.calSurface(U,V);
[inv_U,inv_V]=CST3D_surface.calCoordinate(X,Y,Z);
if any(any((inv_U-U) > 1e-9)) || any(any((inv_V-V) > 1e-9))
    d_U=inv_U-U;
    d_V=inv_V-V;
end
[inv_X,inv_Y,inv_Z]=CST3D_surface.calSurface(inv_U,inv_V);

hold on;
surf(X,Y,Z);
surf(inv_X,inv_Y,inv_Z-1);
hold off;
xlabel('x');
ylabel('y');
zlabel('z');
view(3);
axis equal


%% test 2D curve

% LX=1;
% LY=1;
% 
% xi_gird_num=50;
% 
% class_par_Y=[0.3,0.1];
% 
% curve=CurveCST2D('',LX,LY,[],class_par_Y);
% [X,Y]=curve.calCurve();
% line(X,Y);
% Z=zeros(size(X));
% 
% [U_inv]=curve.calCoordinate(X);
% 
% line(X,Y);
% xlabel('x');
% ylabel('y');
% axis equal