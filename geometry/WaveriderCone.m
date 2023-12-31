classdef WaveriderCone < Shell
    % generate a cone guide waverider by streamline trace
    %

    % aerodynamic properties
    properties
        beta;

        lead_edge_fcn;
        R_0;
        L_total;
        W_total;
    end

    % Conical axial symmetry flow data
    properties
        eta_list;
        V_r_list;
        V_eta_list;
    end

    % main function
    methods
        function self=WaveriderCone(name,Ma_in,T_in,P_in,beta,...
                lead_edge_fcn,R_0,L_total,W_total)
            % lead_edge_fcn is YOZ plane function
            %
            self=self@Shell(name);

            self.beta=beta;
            self.lead_edge_fcn=lead_edge_fcn;
            self.R_0=R_0;
            self.L_total=L_total;
            self.W_total=W_total;

            % air parameter
            R=287.107;
            gama=1.4;

            rho_in=P_in/T_in/R;
            a_in=sqrt(R*gama*T_in);
            V_in=a_in*Ma_in;
            T_0=T_in*(1+(gama-1)/2*Ma_in*Ma_in);
            P_01=P_in*(1+(gama-1)/2*Ma_in*Ma_in)^(gama/(gama-1));

            a_cr_sq=(2*R*gama*T_0/(gama+1));

            [beta,eta,Ma_1,T_2__T_1]=aerodynamicShockWave(gama,Ma_in,beta);
            T_1=T_2__T_1*T_in;
            a_1=sqrt(R*gama*T_1);
            V_1=a_1*Ma_1;

            % check base on beta can get delta or not
            beta_0=aerodynamicShockWave(gama,Ma_in,[],0);
            if beta < beta_0,error('WaveriderCone: beta is too small');end

            % notice: the characteristic of cone shock flow is the physical
            % quantity is the same on the same radial line, so we only need
            % to calculate the different V_r and V_eta on different eta
            differ_fcn=@(eta,V) differConicalFlowVelocity...
                (eta,V,a_cr_sq,gama-1,gama+1);
            [eta_list,V_list]=ode45(differ_fcn,beta:-1e-3:eta,[V_1*cos(beta-eta);-V_1*sin(beta-eta)]);
            V_r_list=V_list(:,1);
            V_eta_list=V_list(:,2);

            % move out V_eta large than zero value
            [~,idx_zero]=interpLinear(0,V_eta_list,eta_list);
            eta_list((idx_zero+1):end)=[];
            V_r_list((idx_zero+1):end)=[];
            V_eta_list((idx_zero+1):end)=[];

            self.eta_list=eta_list;
            self.V_r_list=V_r_list;
            self.V_eta_list=V_eta_list;

        end

        function srf_list=calShell(self,U_num,V_par,W_num)
            % calculate all face point
            %
            sin_beta=sin(self.beta);
            cos_beta=cos(self.beta);

            if isempty(V_par) || V_par ~= fix(V_par)
                % input is V_torlance
                if isempty(V_par)
                    V_par=1e-4*self.W_total;
                end
                min_level=4;
                max_level=10;
                lead_edge=@(y) lineLeadEdge(self.lead_edge_fcn,y,sin_beta,cos_beta);

                [Y,F,~]=GeomApp.meshAdapt1D(lead_edge,0,self.W_total/2,V_par,min_level,max_level,3);
                V_num=length(Y);

                X=F(:,1);Z=F(:,2);R=F(:,3);
            else
                % input is V_num
                V_num=V_par;

                Y=linspace(0,self.W_total/2,V_num);
                X=zeros(V_num,1);
                Z=zeros(V_num,1);
                R=zeros(V_num,1);

                for node_idx=1:V_num
                    Z(node_idx)=self.lead_edge_fcn(Y(node_idx));
                    R(node_idx)=sqrt(Z(node_idx)^2+(Y(node_idx)^2))/sin_beta;
                    X(node_idx)=R(node_idx)*cos_beta;
                end
            end

            % track streamline to get low surface
            vertex_x=self.R_0/tan(self.beta);
            V_r_max=max(self.V_r_list);

            surf_low.X=zeros(V_num,U_num);
            surf_low.Y=zeros(V_num,U_num);
            surf_low.Z=zeros(V_num,U_num);

            for node_idx=1:V_num-1
                x=X(node_idx);
                z=Z(node_idx);
                y=Y(node_idx);
                r=R(node_idx);

                % angle between plane streamline and ZOX place
                cos_psi=z/r/sin(self.beta);
                sin_psi=y/r/sin(self.beta);

                t_end=2*((self.L_total-(x-vertex_x))/V_r_max);
                [SL_t_list,SL_r_list,SL_eta_list]=self.calStreamline...
                    (t_end,r,self.beta,self.eta_list,self.V_r_list,self.V_eta_list);

                % find intersecting line
                r_cos_eta_list=SL_r_list.*cos(SL_eta_list);
                [Equal,index_equal]=interpLinear...
                    ((self.L_total+vertex_x),r_cos_eta_list,[SL_t_list,SL_r_list,SL_eta_list]);
                t_equal=Equal(1);
                r_equal=Equal(2);
                eta_equal=Equal(3);

                SL_t_list(index_equal+2:end)=[];SL_t_list(end)=t_equal;
                SL_r_list(index_equal+2:end)=[];SL_r_list(end)=r_equal;
                SL_eta_list(index_equal+2:end)=[];SL_eta_list(end)=eta_equal;

                t_list=linspace(0,t_equal,U_num);
                [SL,~]=interpLinear...
                    (t_list,SL_t_list,[SL_r_list,SL_eta_list]);
                SL_r_list=SL(:,1);
                SL_eta_list=SL(:,2);

                surf_low.X(V_num-node_idx+1,:)=SL_r_list.*cos(SL_eta_list);
                surf_low.Y(V_num-node_idx+1,:)=SL_r_list.*sin(SL_eta_list)*sin_psi;
                surf_low.Z(V_num-node_idx+1,:)=SL_r_list.*sin(SL_eta_list)*cos_psi;
            end
            surf_low.X(1,:)=self.L_total+vertex_x;
            surf_low.Y(1,:)=self.W_total/2;
            surf_low.Z(1,:)=Z(end);
            surf_low.name='low';
            surf_low.element_type='wgs';

            % extend to get surface up
            surf_up.X=flipud(surf_low.X);
            surf_up.Y=flipud(surf_low.Y);
            surf_up.Z=zeros(V_num,U_num);
            for v_idx=1:V_num
                for u_idx=1:U_num
                    surf_up.Z(v_idx,u_idx)=self.lead_edge_fcn(surf_low.Y(V_num-v_idx+1,u_idx));
                end
            end
            surf_up.name='up';
            surf_up.element_type='wgs';

            % calculate back surface
            surf_back.X=repmat(surf_low.X(:,end),1,W_num);
            surf_back.Y=repmat(surf_low.Y(:,end),1,W_num);
            W=repmat(linspace(0,1,W_num),V_num,1);
            surf_back.Z=surf_low.Z(:,end)+W.*(flipud(surf_up.Z(:,end))-surf_low.Z(:,end));
            surf_back.name='back';
            surf_back.element_type='wgs';

            srf_list={surf_low,surf_up,surf_back};

            function f=lineLeadEdge(lead_edge_fcn,fy,sin_beta,cos_beta)
                fz=lead_edge_fcn(fy);
                fr=sqrt(fz^2+fy^2)/sin_beta;
                fx=fr*cos_beta;
                f=[fx,fz,fr];
            end
        end

        function surf_total=reverseUV(~,surf_total)
            % reverse UV of low and up
            %
            for surf_idx=1:length(surf_total)
                surf=surf_total{surf_idx};
                if strcmp(surf.name,'up')
                    X_old=surf.X;Y_old=surf.Y;Z_old=surf.Z;
                    X=repmat(X_old(1,:),size(X_old,1),1);
                    V=Y_old(:,end)/Y_old(end,end);
                    Y=interp1(X_old(:,1),Y_old(:,1),X(1,:)')';
                    Y=V.*Y;

                    % base on old X, old Y, old Z to interp Z
                    Z=zeros(size(X));
                    Z(:,1)=Z_old(1,1);
                    for u_idx=2:size(X,2)-1
                        X_line=X(:,u_idx);
                        Y_line=[];Z_line=[];
                        v_idx=1;
                        while X_old(v_idx,1) < X_line(1) && v_idx < size(X_old,1)+1
                            Y_line=[Y_line,interp1(X_old(v_idx,:),Y_old(v_idx,:),X_line(1))];
                            Z_line=[Z_line,interp1(X_old(v_idx,:),Z_old(v_idx,:),X_line(1))];
                            v_idx=v_idx+1;
                        end
                        if length(Y_line) < 3
                            error('too less lead edge gird');
                        end
                        Z_line=interp1(Y_line,Z_line,Y(:,u_idx),"spline");
                        Z(:,u_idx)=Z_line;
                    end
                    Z(:,end)=Z_old(:,end);

                    surf.X=X;surf.Y=Y;surf.Z=Z;
                elseif strcmp(surf.name,'low')
                    X_old=surf.X;Y_old=surf.Y;Z_old=surf.Z;

                    X=repmat(X_old(end,:),size(X_old,1),1);
                    V=Y_old(:,end)/Y_old(1,end);
                    Y=interp1(X_old(:,1),Y_old(:,1),X(1,:)')';
                    Y=V.*Y;

                    % base on old X, old Y, old Z to interp Z
                    Z=zeros(size(X));
                    Z(:,1)=Z_old(end,1);
                    for u_idx=2:size(X,2)-1
                        X_line=X(:,u_idx);
                        Y_line=[];Z_line=[];
                        v_idx=size(X_old,1);
                        while X_old(v_idx,1) < X_line(1) && v_idx > 0
                            Y_line=[Y_line,interp1(X_old(v_idx,:),Y_old(v_idx,:),X_line(1))];
                            Z_line=[Z_line,interp1(X_old(v_idx,:),Z_old(v_idx,:),X_line(1))];
                            v_idx=v_idx-1;
                        end
                        if length(Y_line) < 3
                            error('too less lead edge gird');
                        end
                        Z_line=interp1(Y_line,Z_line,Y(:,u_idx),"spline");
                        Z(:,u_idx)=Z_line;
                    end
                    Z(:,end)=Z_old(:,end);

                    surf.X=X;surf.Y=Y;surf.Z=Z;
                end
                surf_total{surf_idx}=surf;
            end
        end
    end

    % visualizate function
    methods
        function writeStepOpenShell(self,step_filestr,U_num,V_num,W_num)
            % write surface into step file
            %
            [~,step_filename,~]=fileparts(step_filestr);

            % write head
            step_file=fopen(step_filestr,'w');
            obj_idx=1;
            obj_idx=writeStepHead(self,step_file,obj_idx,step_filename);

            % write surface
            ADVANCED_FACE_index_list=zeros(1,length(self.face_list));

            surf_total=calShell(self,U_num,V_num,W_num);
            surf_num=length(surf_total);

            for surf_idx=1:surf_num
                fce=surf_total{surf_idx};
                fce_name=fce.name;UDegree=min(size(fce.X,2)-1,3);VDegree=min(size(fce.X,1)-1,3);
                Nodes=cat(3,fce.X,fce.Y,fce.Z);
                fce=GeomApp.VertexToFace(fce_name,Nodes,UDegree,VDegree);
                [step_str,obj_idx,ADVANCED_FACE_index_list(surf_idx)]=fce.getStep(obj_idx);
                fprintf(step_file,step_str);
                fprintf(step_file,'\n');
            end

            % generate OPEN_SHELL
            OPEN_SHELL_index=obj_idx;
            step_str=[num2str(obj_idx,'#%d'),' = OPEN_SHELL ',...
                '( ''NONE'', ',...
                '( ',num2str(ADVANCED_FACE_index_list(1:end-1),'#%d, '),' ',num2str(ADVANCED_FACE_index_list(end),'#%d'),' )',...
                ' );\n'];obj_idx=obj_idx+1;
            fprintf(step_file,step_str);
            fprintf(step_file,'\n');

            % write model
            SHELL_BASED_SURFACE_MODEL_index=obj_idx;
            step_str=[num2str(obj_idx,'#%d'),' = SHELL_BASED_SURFACE_MODEL ',...
                '( ''NONE'', ( ',...
                num2str(OPEN_SHELL_index,'#%d'),' )',...
                ' );\n'];obj_idx=obj_idx+1;
            fprintf(step_file,step_str);
            fprintf(step_file,'\n');

            % write end of step file
            writeStepEnd(self,step_file,obj_idx,step_filename,SHELL_BASED_SURFACE_MODEL_index);

            fclose(step_file);
            clear('step_file');

        end

        function drawShell(self,axes_handle,U_num,V_par,W_num)
            % show all surface of shell
            %
            if isempty(axes_handle),axes_handle=axes(figure());end

            surf_total=self.calShell(U_num,V_par,W_num);

            for surf_index=1:length(surf_total)
                surf=surf_total{surf_index};
                surface(axes_handle,surf.X,surf.Y,surf.Z,'FaceAlpha',0.5);
            end

            xlabel('x');
            ylabel('y');
            zlabel('z');
            view(3);

            % x_range=xlim();
            % y_range=ylim();
            % z_range=zlim();
            % center=[mean(x_range),mean(y_range),mean(z_range)];
            % range=max([x_range(2)-x_range(1),y_range(2)-y_range(1),z_range(2)-z_range(1)])/2;
            % xlim([center(1)-range,center(1)+range]);
            % ylim([center(2)-range,center(2)+range]);
            % zlim([center(3)-range,center(3)+range]);
        end
    end

    % aerodynamic function
    methods(Static)
        function [SL_t_list,SL_r_list,SL_eta_list]=calStreamline...
                (t_end,r_0,eta_0,eta_list,V_r_list,V_eta_list)
            % get streamline
            %
            differ_conical_flow=@(t,X) differConicalFlow...
                (t,X,eta_list,V_r_list,V_eta_list);

            X_0=[r_0,eta_0];

            [SL_t_list,X_list]=ode45(differ_conical_flow,[0:t_end*1e-3:t_end],X_0);
            SL_r_list=X_list(:,1);
            SL_eta_list=X_list(:,2);

        end

    end
end

%% aerodynamci function

function [beta,eta,Ma_2,...
    T_2__T_1,P_2__P_1,rho_2__rho_1,P_02__P_01]=aerodynamicShockWave...
    (gamma,Ma_1,beta,eta)
% function to calculate aerodynamic parameter after oblique shock wave
% beta is oblique shock wave angle, eta is pointed wedge angle
% gamma is fluid specific heat ratio, Ma_1 is fluid mach number
% p01, p02 is total pressure
%
if nargin < 4
    eta=[];
    if nargin < 3
        error('aerodynamicShockWave: lack input angle');
    end
end

if Ma_1 < 1
    error('aerodynamicShockWave: Ma_1 less than 1');
end

if isempty(beta)
    % input eta obtain beta
    if abs(eta-pi/2) < 1e-12
        % normal shock
        beta=pi/2;
    else
        % oblique shock wave
        beta=asin(sqrt(functionBetaeta(gamma,Ma_1,eta)));
        beta=beta(2);
    end
end

if isempty(eta)
    % input beta obtain eta
    if abs(beta-pi/2) < 1e-12
        eta=pi/2;
    else
        tan_eta=functionetaBeta(gamma,Ma_1,beta);
        eta=atan(tan_eta);
    end
end

gamma_sub=gamma-1;
gamma_plus=gamma+1;
sin_beta_sq=(sin(beta))^2;
Ma_1_sq=Ma_1*Ma_1;

% calculate parameter
if abs(beta-pi/2) < 1e-12
    % normal shock
    T_2__T_1=(gamma_sub/gamma_plus)^2*(2*gamma*Ma_1_sq/gamma_sub-1)*...
        (2/gamma_sub/Ma_1_sq+1);
    P_2__P_1=2*gamma/gamma_plus*Ma_1_sq-gamma_sub/gamma_plus;
    rho_2__rho_1=gamma_plus*Ma_1_sq/(gamma_sub*Ma_1_sq+2);
    P_02__P_01=(2*gamma*Ma_1_sq/gamma_plus-gamma_sub/gamma_plus)^(-1/gamma_sub)*...
        (gamma_plus*Ma_1_sq/(gamma_sub*Ma_1_sq+2))^(gamma/gamma_sub);
    Ma_2_sq=(Ma_1_sq+2/gamma_sub)/(2*gamma*Ma_1_sq/gamma_sub-1);
    Ma_2=sqrt(Ma_2_sq);
else
    % oblique shock wave
    T_2__T_1=(gamma_sub/gamma_plus)^2*(2*gamma*Ma_1_sq*sin_beta_sq/gamma_sub-1)*...
        (2/gamma_sub/Ma_1_sq/sin_beta_sq+1);
    P_2__P_1=2*gamma/gamma_plus*Ma_1_sq*sin_beta_sq-gamma_sub/gamma_plus;
    rho_2__rho_1=gamma_plus*Ma_1_sq*sin_beta_sq/(gamma_sub*Ma_1_sq*sin_beta_sq+2);
    P_02__P_01=(2*gamma*Ma_1_sq*sin_beta_sq/gamma_plus-gamma_sub/gamma_plus)^(-1/gamma_sub)*...
        (gamma_plus*Ma_1_sq*sin_beta_sq/(gamma_sub*Ma_1_sq*sin_beta_sq+2))^(gamma/gamma_sub);
    Ma_2_sq=(Ma_1_sq+2/gamma_sub)/(2*gamma*Ma_1_sq*sin_beta_sq/gamma_sub-1)+...
        2/gamma_sub*Ma_1_sq*cos(beta)^2/(Ma_1_sq*sin_beta_sq+2/gamma_sub);
    Ma_2=sqrt(Ma_2_sq);
end

    function sin_beta_sq=functionBetaeta(gamma,Ma_1,eta)
        % function to get sin(beta)^2 by eta
        %
        tan_eta_sq__=tan(eta)^2;
        Ma_1_sq__=Ma_1*Ma_1;
        Ma_1_qu__=Ma_1_sq__*Ma_1_sq__;
        gama_plus__=gamma+1;
        c0=1;
        c1=-(tan_eta_sq__*(Ma_1_qu__*gama_plus__*gama_plus__/4+Ma_1_sq__*gama_plus__+1)+...
            (2*Ma_1_sq__+1));
        c2=(tan_eta_sq__*(Ma_1_qu__*gama_plus__+2*Ma_1_sq__)+(Ma_1_qu__+2*Ma_1_sq__));
        c3=-(tan_eta_sq__*Ma_1_qu__+Ma_1_qu__);
        sin_beta_sq=roots([c3 c2 c1 c0]);
        sin_beta_sq=sin_beta_sq(1:2);
    end
    function tan_eta=functionetaBeta(gamma,Ma_1,beta)
        % function to get tan(eta) by beta
        %
        sin_beta__=sin(beta);
        sin_beta_sq__=sin_beta__*sin_beta__;
        tan_beta__=tan(beta);
        Ma_1_sq__=Ma_1*Ma_1;
        tan_eta=(Ma_1_sq__*sin_beta_sq__-1)/...
            (Ma_1_sq__*((gamma+1)/2-sin_beta_sq__)+1)/tan_beta__;
    end
end

function dX=differConicalFlow...
    (t,X,eta_list,V_r_list,V_eta_list)
% equation of conical flow field
% V=[r;eta;V_r;V_eta]
%
% notice: the characteristic of cone shock flow is the physical quantity is
% the same on the same radial line
%
r=X(1);
eta=X(2);

[y_equal__,~]=interpLinear...
    (eta,eta_list,[V_r_list,V_eta_list]);

V_r=y_equal__(1);
V_eta=y_equal__(2);

dr=V_r;
deta=V_eta/r;
dX=[dr;deta];
end

function dV=differConicalFlowVelocity...
    (eta,V,a_cr_sq,gama_sub,gama_plus)
% equation of conical flow field
% r is the distance to conical tip
% eta is the angle between r and axis of cone
% V=[V_r;V_eta]
%
% notice: the characteristic of cone shock flow is the physical quantity is
% the same on the same radial line
%
% reference: [1] 童秉纲, 孔祥言, 邓国华. 气体动力学[M]. 气体动力学, 2012,
% 158-159.
%
V_r=V(1);
V_eta=V(2);

a_sq=gama_sub/2*(gama_plus/gama_sub*a_cr_sq-V_r*V_r-V_eta*V_eta);

dV_dr=V_eta;
dV_deta=(V_r+V_eta*cot(eta))/(V_eta*V_eta/a_sq-1)-V_r;

dV=[dV_dr;dV_deta];
end

%% auxiliary function

function [Y_pred,X_pred_idx]=interpLinear(X_pred,X,Y)
% a simple linear interpolation to calculate data
% Y can be m x n matrix, n is variable
%
[X,idx]=sort(X);
Y=Y(idx,:);
Y_pred=zeros(length(X_pred),size(Y,2));
X_pred_idx=zeros(length(X_pred),1);
for x_idx=1:length(X_pred)
    x_pred=X_pred(x_idx);
    num=length(X);
    idx=num; % search start from last one, find out X samll than x
    while ((idx > 1) && (X(idx) > x_pred))
        idx=idx-1;
    end

    if (idx == num)
        % out interp
        Y_pred(x_idx,:)=(Y(num,:)-Y(num-1,:))/(X(num)-X(num-1))*(x_pred-X(num))+Y(num,:);
        X_pred_idx(x_idx)=idx;
    elseif (idx == 0)
        Y_pred(x_idx,:)=(Y(2,:)-Y(1,:))/(X(2)-X(1))*(x_pred-X(1))+Y(1,:);
        X_pred_idx(x_idx)=idx;
    else
        % linear interpolation
        Y_pred(x_idx,:)=Y(idx,:)+...
            (Y(idx+1,:)-Y(idx,:))*...
            (x_pred-X(idx,:))/...
            (X(idx+1)-X(idx));
        X_pred_idx(x_idx)=idx;
    end
end

end
