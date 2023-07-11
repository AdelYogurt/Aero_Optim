clc;
clear;
close all hidden;

gama=1.4;
Ma_1=10;

beta=20/180*pi;
theta=[];
% beta=[];
% theta=15/180*pi;

[beta,theta,Ma_2,...
    p_2__p_1,rou_2__rou_1,T_2__T_1,p_02__p_01]=aerodynamicShockWave...
    (gama,Ma_1,beta,theta)

function [beta,theta,Ma_2,...
    p_2__p_1,rou_2__rou_1,T_2__T_1,p_02__p_01]=aerodynamicShockWave...
    (gama,Ma_1,beta,theta)
% function to calculate aerodynamic parameter after oblique shock wave
% beta is oblique shock wave angle, theta is pointed wedge angle
% gama is fluid specific heat ratio, Ma_1 is fluid mach number
% p01, p02 is total pressure
%
if nargin < 4
    theta=[];
    if nargin < 3
        error('aerodynamicShockWave: lack input angle');
    end
end

if Ma_1 < 1
    error('aerodynamicShockWave: Ma_1 less than 1');
end

if isempty(beta)
    % input theta obtain beta
    if abs(theta-pi/2) < 1e-12
        % normal shock
        beta=pi/2;
    else
        % oblique shock wave
        beta=asin(sqrt(functionBetaTheta(gama,Ma_1,theta)));
        beta=beta(2);
    end
end

if isempty(theta)
    % input beta obtain theta
    if abs(beta-pi/2) < 1e-12
        theta=pi/2;
    else
        tan_theta=functionThetaBeta(gama,Ma_1,beta);
        theta=atan(tan_theta);
    end
end

gama_sub=gama-1;
gama_plus=gama+1;
sin_beta_sq=(sin(beta))^2;
Ma_1_sq=Ma_1*Ma_1;

% calculate parameter
if abs(beta-pi/2) < 1e-12
    % normal shock
    p_2__p_1=2*gama/gama_plus*Ma_1_sq-gama_sub/gama_plus;
    rou_2__rou_1=gama_plus*Ma_1_sq/(gama_sub*Ma_1_sq+2);
    T_2__T_1=(gama_sub/gama_plus)^2*(2*gama*Ma_1_sq/gama_sub-1)*...
        (2/gama_sub/Ma_1_sq+1);
    p_02__p_01=(2*gama*Ma_1_sq/gama_plus-gama_sub/gama_plus)^(-1/gama_sub)*...
        (gama_plus*Ma_1_sq/(gama_sub*Ma_1_sq+2))^(gama/gama_sub);
    Ma_2_sq=(Ma_1_sq+2/gama_sub)/(2*gama*Ma_1_sq/gama_sub-1);
    Ma_2=sqrt(Ma_2_sq);
else
    % oblique shock wave    
    p_2__p_1=2*gama/gama_plus*Ma_1_sq*sin_beta_sq-gama_sub/gama_plus;
    rou_2__rou_1=gama_plus*Ma_1_sq*sin_beta_sq/(gama_sub*Ma_1_sq*sin_beta_sq+2);
    T_2__T_1=(gama_sub/gama_plus)^2*(2*gama*Ma_1_sq*sin_beta_sq/gama_sub-1)*...
        (2/gama_sub/Ma_1_sq/sin_beta_sq+1);
    p_02__p_01=(2*gama*Ma_1_sq*sin_beta_sq/gama_plus-gama_sub/gama_plus)^(-1/gama_sub)*...
        (gama_plus*Ma_1_sq*sin_beta_sq/(gama_sub*Ma_1_sq*sin_beta_sq+2))^(gama/gama_sub);
    Ma_2_sq=(Ma_1_sq+2/gama_sub)/(2*gama*Ma_1_sq*sin_beta_sq/gama_sub-1)+...
        2/gama_sub*Ma_1_sq*cos(beta)^2/(Ma_1_sq*sin_beta_sq+2/gama_sub);
    Ma_2=sqrt(Ma_2_sq);
end

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