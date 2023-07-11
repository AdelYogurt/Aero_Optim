clc;
clear;
close all hidden;

R=8314/28.959;

beta=10/180*pi;

R_0=0.5;
L_total=4;

Ma_1=10;
T_1=89.3;
gama=1.4;
p_1=951.5;

a_1=sqrt(R*gama*T_1);
V_1=a_1*Ma_1;
T_0=T_1*(1+(gama-1)/2*Ma_1*Ma_1);
p_01=p_1*(1+(gama-1)/2*Ma_1*Ma_1)^(gama/(gama-1));

a_cr_sq=(2*R*gama*T_0/(gama+1));
a_cr=sqrt(a_cr_sq);

[beta,theta,Ma_2,...
    p_2__p_1,rou_2__rou_1,T_2__T_1,p_02__p_01]=aerodynamicObliqueShock...
    (gama,Ma_1,beta);
T_2=T_2__T_1*T_1;
a_2=sqrt(R*gama*T_2);
V_2=a_2*Ma_2;

differ_conical_flow_velocity=@(theta,V) differConicalFlowVelocity...
    (theta,V,a_cr_sq,gama-1,gama+1);

V_0=[V_2*cos(beta-theta);-V_2*sin(beta-theta)];

[theta_list,V_list]=ode45(differ_conical_flow_velocity,[beta:-1e-3:1e-3],V_0);
V_r_list=V_list(:,1);
V_theta_list=V_list(:,2);

% find V_theta equal to zero to get delta
index_zero=1;
while V_theta_list(index_zero) < 0
    index_zero=index_zero+1;
end
delta=(theta_list(index_zero-1)+theta_list(index_zero))/2;
theta_list((index_zero+1):end)=[];
V_r_list((index_zero+1):end)=[];
V_theta_list((index_zero+1):end)=[];

% plot(theta_list,V_r_list);

axis equal;

FCC_function=@(y) -R_0-y^2+y^4/2;
% drawFunction(FCC_function,0,1);

x_0=R_0/tan(beta);
r_0_sq=(tan(beta)*(x_0+L_total))^2;
W_total=2*fsolve(@(y) (FCC_function(y))^2+y^2-r_0_sq,L_total/2);

sin_beta_sq=sin(beta)^2;
cos_beta=cos(beta);

Y=[0:0.05:W_total/2];
Z=zeros(1,length(Y));
X=zeros(1,length(Y));
R=zeros(1,length(Y));

for index=1:length(Y)
    Z(index)=FCC_function(Y(index));
    R(index)=sqrt((Z(index)^2+Y(index)^2)/sin_beta_sq);
    X(index)=R(index)*cos_beta;
end

% track stream line to get low surface
LS_Y=zeros(1,length(Y));
LS_Z=zeros(1,length(Y));
for index=1:length(LS_Z)-1
    x=X(index);
    y=Y(index);
    z=Z(index);
    r=R(index);
    
    cos_phi=z/r/sin(beta);
    sin_phi=y/r/sin(beta);
    
    [r_list,theta_t_list]=getStreamline...
        (2*((L_total-(x-X(1)))/V_2),r,beta,theta_list,V_r_list,V_theta_list);
    
    % find intersecting line
    index_equal=2;
    while r_list(index_equal)*cos(theta_t_list(index_equal)) < (L_total+X(1))
        index_equal=index_equal+1;
    end
    r_equal=(r_list(index_equal-1)+r_list(index_equal))/2;
    theta_equal=(theta_t_list(index_equal-1)+theta_t_list(index_equal))/2;
    r_sin_theta=r_equal*sin(theta_equal);
    
    LS_Y(index)=r_sin_theta*sin_phi;
    LS_Z(index)=r_sin_theta*cos_phi;
end
LS_Y(end)=Y(end);
LS_Z(end)=Z(end);

line(X,Y,Z);
line((X(1)+L_total)*ones(1,length(Y)),Y,Z);
line((X(1)+L_total)*ones(1,length(Y)),LS_Y,LS_Z);
axis equal;

function [r_list,theta_list]=getStreamline...
    (t_end,r_0,theta_0,theta_list,V_r_list,V_theta_list)
% get streamline of conical flow
%
differ_conical_flow=@(t,X) differConicalFlow...
    (t,X,theta_list,V_r_list,V_theta_list);

X_0=[r_0;theta_0];

[t_list,X_list]=ode45(differ_conical_flow,[0 t_end],X_0);
r_list=X_list(:,1);
theta_list=X_list(:,2);

% plot(r_list.*cos(theta_list),-r_list.*sin(theta_list));
end

function dX=differConicalFlow...
    (t,X,theta_list,V_r_list,V_theta_list)
% equation of conical flow field
% V=[r;theta;V_r;V_theta]
%
r=X(1);
theta=X(2);

index_equal=2;
while index_equal < length(theta_list) && theta_list(index_equal) > theta
    index_equal=index_equal+1;
end

if(index_equal > length(V_r_list))
   disp('?'); 
end

dr=(V_r_list(index_equal-1)+V_r_list(index_equal))/2;
dtheta=(V_theta_list(index_equal-1)+V_theta_list(index_equal))/2/r;
dX=[dr;dtheta];
end

function dV=differConicalFlowVelocity...
    (theta,V,a_cr_sq,gama_sub,gama_plus)
% equation of conical flow field
% V=[V_r;V_theta]
%
V_r=V(1);
V_theta=V(2);
dV_r=V_theta;
a_sq=gama_sub/2*(gama_plus/gama_sub*a_cr_sq-V_r*V_r-V_theta*V_theta);
dV_theta=a_sq*(V_r+V_theta*cot(theta))/(V_theta*V_theta-a_sq)-V_r;
dV=[dV_r;dV_theta];
end

function [beta,theta,Ma_2,...
    p_2__p_1,rou_2__rou_1,T_2__T_1,p_02__p_01]=aerodynamicObliqueShock...
    (gama,Ma_1,beta,theta)
% function to calculate aerodynamic parameter after oblique shock wave
% beta is oblique shock wave angle, theta is pointed wedge angle
% gama is fluid specific heat ratio, Ma_1 is fluid mach number
% p01, p02 is total pressure
%
if isempty(beta)
    beta=asin(sqrt(functionBetaTheta(gama,Ma_1,theta)));
    beta=beta(2);
else
    theta=atan(functionThetaBeta(gama,Ma_1,beta));
end

gama_sub=gama-1;
gama_plus=gama+1;
sin_beta_sq=(sin(beta))^2;
Ma_1_sq=Ma_1*Ma_1;

p_2__p_1=2*gama/gama_plus*Ma_1_sq*sin_beta_sq-gama_sub/gama_plus;
rou_2__rou_1=gama_plus*Ma_1_sq*sin_beta_sq/(gama_sub*Ma_1_sq*sin_beta_sq+2);
T_2__T_1=(gama_sub/gama_plus)^2*(2*gama*Ma_1_sq*sin_beta_sq/gama_sub-1)*...
    (2/gama_sub/Ma_1_sq/sin_beta_sq+1);
p_02__p_01=(2*gama*Ma_1_sq*sin_beta_sq/gama_plus-gama_sub/gama_plus)^(-1/gama_sub)*...
    (gama_plus*Ma_1_sq*sin_beta_sq/(gama_sub*Ma_1_sq*sin_beta_sq+2))^(gama/gama_sub);
Ma_2_sq=(Ma_1_sq+2/gama_sub)/(2*gama*Ma_1_sq*sin_beta_sq/gama_sub-1)+...
    2/gama_sub*Ma_1_sq*cos(beta)^2/(Ma_1_sq*sin_beta_sq+2/gama_sub);
Ma_2=sqrt(Ma_2_sq);

    function sin_beta_sq=functionBetaTheta(gama,Ma_1,theta)
        % function to get sin(beta)^2 by theta
        %
        tan_theta_sq__=tan(theta)^2;
        Ma_1_sq__=Ma_1*Ma_1;
        Ma_1_qu__=Ma_1_sq__*Ma_1_sq__;
        gama_plus__=gama+1;
        c0=1;
        c1=-(tan_theta_sq__*(Ma_1_qu__*gama_plus__*gama_plus__/4+Ma_1_sq__*gama_plus__+1)+...
            (2*Ma_1_sq__+1));
        c2=(tan_theta_sq__*(Ma_1_qu__*gama_plus__+2*Ma_1_sq__)+(Ma_1_qu__+2*Ma_1_sq__));
        c3=-(tan_theta_sq__*Ma_1_qu__+Ma_1_qu__);
        sin_beta_sq=roots([c3 c2 c1 c0]);
        sin_beta_sq=sin_beta_sq(1:2);
    end
    function tan_theta=functionThetaBeta(gama,Ma_1,beta)
        % function to get tan(theta) by beta
        %
        sin_beta__=sin(beta);
        sin_beta_sq__=sin_beta__*sin_beta__;
        tan_beta__=tan(beta);
        Ma_1_sq__=Ma_1*Ma_1;
        tan_theta=(Ma_1_sq__*sin_beta_sq__-1)/...
            (Ma_1_sq__*((gama+1)/2-sin_beta_sq__)+1)/tan_beta__;
    end
end

function drawFunction(draw_function,low_bou,up_bou,...
    grid_number,Y_min,Y_max,figure_handle)
% function to draw one dimension/ two dimension function
%
if nargin < 7
    figure_handle=figure(10);
    if nargin < 6
        Y_max=inf;
        if nargin < 5
            Y_min=-inf;
            if nargin < 4
                grid_number=100;
            end
        end
    end
end
axes_handle=figure_handle.CurrentAxes;
if isempty(axes_handle)
    axes_handle=axes(figure_handle);
end
axes_context=axes_handle.Children;
dimension=size(low_bou,1);

switch dimension
    case 1
        d_bou=(up_bou-low_bou)/grid_number;
        X__=low_bou:d_bou:up_bou;
        fval__=zeros(grid_number+1,1);
        for x_index__=1:(grid_number+1)
            fval__(x_index__)=draw_function(X__(x_index__));
        end
        line(axes_handle,X__,fval__);
        xlabel('X');
        ylabel('value');
        
    case 2
        d_bou=(up_bou-low_bou)/grid_number;
        [X__,Y]=meshgrid(low_bou(1):d_bou(1):up_bou(1),low_bou(2):d_bou(2):up_bou(2));
        fval__=zeros(grid_number+1);
        for x_index__=1:grid_number+1
            for y_index__=1:grid_number+1
                predict_x=([x_index__;y_index__]-1).*d_bou+low_bou;
                fval__(x_index__,y_index__)=draw_function(predict_x);
            end
        end
        fval__(find(fval__ > Y_max))=Y_max;
        fval__(find(fval__ < Y_min))=Y_min;
        axes_context=[axes_context;surface(X__',Y',fval__,'FaceAlpha',0.5,'EdgeColor','none')];
        axes_handle.set('Children',axes_context);
        xlabel('X');
        ylabel('Y');
        zlabel('value');
end
end