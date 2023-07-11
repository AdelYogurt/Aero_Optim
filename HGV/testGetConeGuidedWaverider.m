clc;
clear;
close all hidden;

R=8314/28.959;

beta=10/180*pi;
R_0=0.05;
L_total=0.4;
W_total=0.16;
node_interval=0.01;

Ma_1=15;
gama=1.4;

T_1=89.3;
p_1=951.5;

rou_1=p_1/T_1/R;
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

% judge base on beta can get delta or not
[beta_0,~,~,...
    ~,~,~,~]=aerodynamicObliqueShock...
    (gama,Ma_1,[],0);
if beta < beta_0
    error('getConeGuidedWaverider: beta is too small');
end

% get delta
differ_conical_flow_velocity=@(theta,V) differConicalFlowVelocity...
    (theta,V,a_cr_sq,gama-1,gama+1);
V_0=[V_2*cos(beta-theta);-V_2*sin(beta-theta)];
[theta_list,V_list]=ode45(differ_conical_flow_velocity,[beta:-1e-3:1e-3],V_0);
V_r_list=V_list(:,1);
V_theta_list=V_list(:,2);

% plot(theta_list,V_theta_list);

% find V_theta equal to zero to get delta
[delta,index_zero,out_index]=interpolateLinear(V_theta_list,theta_list,0,1);
theta_list((index_zero+1):end)=[];
V_r_list((index_zero+1):end)=[];
V_theta_list((index_zero+1):end)=[];

x_0=R_0/tan(beta);
R_total_sq=(tan(beta)*(x_0+L_total))^2;
H_total=(sqrt(R_total_sq-(W_total/2)^2)-R_0);
FCC_function=@(z) -R_0+H_total*(cos(pi*z/W_total)-1);
% drawFunction(FCC_function,0,W_total/2);

sin_beta=sin(beta);
sin_beta_sq=sin_beta^2;
cos_beta=cos(beta);

Y=[0:node_interval:W_total/2]';
node_number=length(Y);
if abs((node_number-1)*node_interval-W_total/2) > 1e-3
    Y=[Y;W_total/2];
    node_number=node_number+1;
end
X=zeros(node_number,1);
Z=zeros(node_number,1);
R=zeros(node_number,1);

for node_index=1:node_number
    Z(node_index)=FCC_function(Y(node_index));
    R(node_index)=sqrt(Z(node_index)^2+(Y(node_index)^2))/sin_beta;
    X(node_index)=R(node_index)*cos_beta;
end

% track stream line to get low surface
LS_X=ones(node_number,1)*(x_0+L_total);
LS_Y=zeros(node_number,1);
LS_Z=zeros(node_number,1);

for node_index=1:node_number-1
    x=X(node_index);
    z=Z(node_index);
    y=Y(node_index);
    r=R(node_index);
    
    cos_phi=z/r/sin(beta);
    sin_phi=y/r/sin(beta);
    
    [r_list,theta_t_list]=getStreamline...
        (2*((L_total-(x-x_0))/V_2),r,beta,theta_list,V_r_list,V_theta_list);
    
    % find intersecting line
    r_cos_theta_list=r_list.*cos(theta_t_list);
    [Equal,index_equal,out_index]=interpolateLinear...
        (r_cos_theta_list,[r_list,theta_t_list],(L_total+x_0),1);
    r_equal=Equal(1);
    theta_equal=Equal(2);
    
    r_list(index_equal+1:end)=[];
    theta_t_list(index_equal+1:end)=[];
    line(r_list.*cos(theta_t_list),...
        r_list.*sin(theta_t_list)*sin_phi,...
        r_list.*sin(theta_t_list)*cos_phi);
    
    r_sin_theta=r_equal*sin(theta_equal);
    
    LS_Y(node_index)=r_sin_theta*cos_phi;
    LS_Z(node_index)=r_sin_theta*sin_phi;
end
LS_Y(end)=Z(end);
LS_Z(end)=Y(end);

view(3);
line(X,Y,Z);
line(LS_X,Y,Z);
line(LS_X,LS_Z,LS_Y);
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;


function [r_list,theta_list]=getStreamline...
    (t_end,r_0,theta_0,theta_list,V_r_list,V_theta_list)
% get streamline
%
differ_conical_flow=@(t,X) differConicalFlow...
    (t,X,theta_list,V_r_list,V_theta_list);

X_0=[r_0;theta_0];

[t_list,X_list]=ode45(differ_conical_flow,[0:(t_end)*1e-3:t_end],X_0);
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

[y_equal__,~,~]=interpolateLinear...
    (theta_list,[V_r_list,V_theta_list],theta,-1);

V_r=y_equal__(1);
V_theta=y_equal__(2);

dr=V_r;
dtheta=V_theta/r;
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
% calculate oblique shock parameters
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

function [y,x_index,out_index]=interpolateLinear...
    (x_list,y_list,x,increment)
% interpolate to get data
% increment is 1, decrease is -1
% x_index is small or large than x
% y_list can be mxn, n is variable
%
if length(x_list) ~= length(y_list)
    error('interpolateFirst: length do not match');
end
out_index=0;

if increment*x < increment*x_list(1)
    out_index=-1*increment;
    x_index=0;
    y=y_list(1,:);
    return;
end
if increment*x > increment*x_list(end)
    out_index=1*increment;
    x_index=length(x_list)+1;
    y=y_list(end,:);
    return;
end

% find x in where
x_index=1;
while (increment*x_list(x_index) < increment*x) &&...
        (x_index < length(x_list)-1)
    if x == x_list(x_index)
        y=y_list(x_index,:);
        return;
    end
    x_index=x_index+1;
end
if x == x_list(end)
    x_index=length(x_list);
    y=y_list(end,:);
    return;
end

% first interpolation
d_x__=increment*(x_list(x_index+1)-x_list(x_index));
d_y__=y_list(x_index+1,:)-y_list(x_index,:);

y=increment*(x-x_list(x_index))/d_x__*d_y__+y_list(x_index,:);
end

function drawFunction(object_function,low_bou,up_bou,...
    grid_number,Y_min,Y_max,figure_handle)
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
dimension=size(low_bou,1);

switch dimension
    case 1
        d_bou=(up_bou-low_bou)/grid_number;
        X__=low_bou:d_bou:up_bou;
        fval__=zeros(grid_number+1,1);
        for x_index=1:(grid_number+1)
            fval__(x_index)=object_function(X__(x_index));
        end
        line(axes_handle,X__,fval__);
        xlabel('X');
        ylabel('value');
        
    case 2
        d_bou=(up_bou-low_bou)/grid_number;
        [X__,Y]=meshgrid(low_bou(1):d_bou(1):up_bou(1),low_bou(2):d_bou(2):up_bou(2));
        fval__=zeros(grid_number+1);
        for x_index=1:grid_number+1
            for y_index=1:grid_number+1
                predict_x=([x_index;y_index]-1).*d_bou+low_bou;
                fval__(x_index,y_index)=object_function(predict_x);
            end
        end
        fval__(find(fval__ > Y_max))=Y_max;
        fval__(find(fval__ < Y_min))=Y_min;
        surf(axes_handle,X__',Y',fval__,'FaceAlpha',0.5,'EdgeColor','none');
        xlabel('X');
        ylabel('Y');
        zlabel('value');
end

end