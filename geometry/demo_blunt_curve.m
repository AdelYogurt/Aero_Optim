clc;
clear;
close all hidden;

%% quadratic blunt

% vctr_up=[1,-0.1];
% par_R=0.2;
% P_up=[0,par_R];
% 
% k=vctr_up(2)/vctr_up(1);
% ploy=[2*k^2,-(4*k+1)*par_R,2*par_R^2];
% d=roots(ploy);
% d=d(end);
% 
% P0=[-d,0];
% P1=[-d,par_R-k*d];
% P2=[0,par_R];
% 
% % P0=[-1.1,0];
% % P1=[-1.1,0.4];
% % P2=[0,0.2];
% 
% crv=Curve([P0;P1;P2]);
% crv.displayGeom(gca());view(2);axis equal;
% 
% dP=-2*P0+2*P1;ddP=2*P0-4*P1+2*P2;
% R=1/((dP(1)*ddP(2)-dP(2)*ddP(1))/norm(dP)^3)
% 
% hold on;
% quiver(P_up(1),P_up(2),vctr_up(1),vctr_up(2));
% crv.displayGeom(gca());view(2);axis equal;
% hold off;

%% cubic blunt

% vctr_low=[1,-0.1];
% vctr_up=[1,1];
% par_R=0.2;
% P_low=[0,-par_R];
% P_up=[0,par_R];
% 
% crv=getBluntEdge(vctr_low,vctr_up,P_low,P_up);
% 
% hold on;
% quiver(P_up(1),P_up(2),vctr_up(1),vctr_up(2));
% quiver(P_low(1),P_low(2),vctr_low(1),vctr_low(2));
% crv.displayGeom(gca());view(2);axis equal;
% hold off;
% 
% function crv=getBluntEdge(vctr_low,vctr_up,P_low,P_up)
% % connect up and low curve by BSpline
% %
% vctr_low=vctr_low/norm(vctr_low);
% vctr_up=vctr_up/norm(vctr_up);
% k_low=vctr_low(2)/vctr_low(1);
% k_up=vctr_up(2)/vctr_up(1);
% theta_low=-atan(k_low);
% theta_up=atan(k_up);
% theta=pi/2+(theta_low+theta_up)/2;
% 
% if theta_low < theta_up
%     P3_low=P_low;
%     k=tan(theta_up-theta+pi);
%     P3_up=[P_up(2)-k_up*P_up(1),P_low(2)-k*P_low(1)]/[-k_up,1;-k,1]';
% elseif theta_low > theta_up
%     P3_up=P_up;
%     k=tan(theta-theta_low);
%     P3_low=[P_low(2)-k_low*P_low(1),P_up(2)-k*P_up(1)]/[-k_low,1;-k,1]';
% else
%     P3_low=P_low;
%     P3_up=P_up;
% end
% 
% P_center=(P3_low+P3_up)/2;
% vctr_center=(vctr_low+vctr_up)/2;
% vctr_center=vctr_center/norm(vctr_center);
% 
% k_1=tan(theta-pi/2);
% R=norm(P3_low-P3_up)/2;
% k_0=-1/R;
% 
% P_up_list=calBluntPoint(k_1,R);
% P_low_list=flipud(P_up_list);
% P_low_list(:,2)=-P_low_list(:,2);
% rotate_matrix=[
%     vctr_center(1),-vctr_center(2),0;
%     vctr_center(2),vctr_center(1),0;
%     0,0,1];
% P_up_list=(P_up_list*rotate_matrix'+[P_center,0]);
% P_low_list=(P_low_list*rotate_matrix'+[P_center,0]);
% crv_up=Curve(P_up_list);
% crv_low=Curve(P_low_list);
% 
% if theta_low < theta_up
%     P_ext=[P3_up,0;P_up,0];
%     crv=GeomApp.JointCurve([crv_low,crv_up,Curve(P_ext)]);
% elseif theta_low > theta_up
%     P_ext=[P_low,0;P3_low,0];
%     crv=GeomApp.JointCurve([Curve(P_ext),crv_low,crv_up]);
% else
%     crv=GeomApp.JointCurve([crv_low,crv_up,[]]);
% end
% 
% end
% 
% function P=calBluntPoint(k_1,R)
% K_0=-1/R;
% if abs(k_1) < sqrt(eps)
%     d=3*R/2;
% else
%     delta=(2*k_1*R+2*R/3)^2-4*k_1^2*R^2;
%     d=(R/k_1+R/3/k_1^2-sqrt(delta)/2/k_1^2);
% end
% 
% P_0=[-d,0,0];
% P_3=[0,R,0];
% v_0=[0,1,0];
% v_1=[1,k_1,0];
% 
% D=P_3-P_0;
% 
% delta_0=cross(v_1,D)/cross(v_1,v_0);
% delta_1=(1.5*K_0*delta_0^2-cross(v_0,D))/cross(v_0,v_1);
% 
% P_1=P_0+delta_0*v_0;
% P_2=P_3+delta_1*v_1;
% 
% P=[P_0;P_1;P_2;P_3];
% end

%% arc blunt

R=0.05;torl=1e-3;N=0.8;LX=1;LY=0.2;
crv=CurveCST([N,N],[],LX,LY);
crv.addSpline(Curve(cat(2,(0:0.2:1)',rand(6,1)-0.5)))
fig_hdl=figure();
axe_hdl=axes(fig_hdl);axis equal;
crv.displayGeom(axe_hdl);
[k1,k2]=crv.calTangTorl(torl/N);
[h1,d1]=calBluntHD(R,k1);
[h2,d2]=calBluntHD(R,k2);
usb=d1/LX;ueb=1-d2/LX;

crv_blunt=CurveCST();
crv_blunt.shape_fcn=@(U) calPoint(U,crv,R,LX,usb,ueb,d1,d2,h1,h2);
crv_blunt.displayGeom(axe_hdl);

function Point=calPoint(U,crv,R,LX,usb,ueb,d1,d2,h1,h2)
Point=zeros(length(U),2);
Point(:,1)=U*LX;
idx=find(U < usb);
if ~isempty(idx),Point(idx,2)=sqrt(R^2-(Point(idx,1)-R).^2+eps);end
idx=find(ueb < U);
if ~isempty(idx),Point(idx,2)=sqrt(R^2-(Point(idx,1)-(LX-R)).^2+eps);end
idx=find(usb <= U & U <= ueb);
if ~isempty(idx)
    u_loc=(U(idx)-usb)/(ueb-usb);L_loc=LX-d1-d2;
    Point(idx,:)=crv.calPoint(u_loc);
    Point(idx,1)=Point(idx,1)*L_loc+d1;
    Point(idx,2)=Point(idx,2)+h1*(1-u_loc)+h2*u_loc;
end
end

function [h,d]=calBluntHD(R,k)
h=R./sqrt(1+k.^2);
d=R-k.*R./sqrt(1+k.^2);
end
