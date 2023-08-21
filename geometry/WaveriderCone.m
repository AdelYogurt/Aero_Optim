classdef WaveriderCone < Body
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
        function self=WaveriderCone(Ma_in,T_in,P_in,beta,...
                lead_edge_fcn,R_0,L_total,W_total)
            % lead_edge_fcn is YOZ plane function
            %
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
            if beta < beta_0,error('getConeGuidedWaverider: beta is too small');end

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

        function surf_total=calSurfaceMatrix(self,U_num,V_par,W_num)
            % calculate surface matrix
            %
            sin_beta=sin(self.beta);
            cos_beta=cos(self.beta);

            if isempty(V_par) || V_par ~= fix(V_par)
                % input is V_torlance
                if isempty(V_par)
                    V_par=1e-4*self.W_total;
                end
                max_level=10;
                lead_edge=@(y) lineLeadEdge(self.lead_edge_fcn,y,sin_beta,cos_beta);

                [Y,F,~]=meshAdapt1D(lead_edge,0,self.W_total/2,V_par,max_level,3);
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

            surf_total={surf_low,surf_up,surf_back};

            function f=lineLeadEdge(lead_edge_fcn,fy,sin_beta,cos_beta)
                fz=lead_edge_fcn(fy);
                fr=sqrt(fz^2+fy^2)/sin_beta;
                fx=fr*cos_beta;
                f=[fx,fz,fr];
            end
        end

        function writeStepOpenShell(self,step_filestr,U_num,V_num,W_num)
            % write surface into step file
            %

            % check file name
            if length(step_filestr) > 4
                if ~strcmp(step_filestr((end-4):end),'.step')
                    step_filestr=[step_filestr,'.step'];
                end
            else
                step_filestr=[step_filestr,'.step'];
            end
            [~,step_filename,~]=fileparts(step_filestr);

            % write head
            step_file=fopen(step_filestr,'w');
            object_index=1;
            object_index=writeStepHead(self,step_file,object_index,step_filename);

            % write surface
            ADVANCED_FACE_index_list=zeros(1,length(self.surface_list));

            surf_total=calSurfaceMatrix(self,U_num,V_num,W_num);

            for surf_idx=1:length(surf_total)
                surf=surf_total{surf_idx};
                surf_name=surf.name;
                if size(surf.X,2) == 2
                    u_list=[0,0,1,1];
                else
                    u_list=[0,0,0,linspace(0,1,size(surf.X,2)),1,1,1];
                end
                if size(surf.X,1) == 2
                    v_list=[0,0,1,1];
                else
                    v_list=[0,0,0,linspace(0,1,size(surf.X,1)),1,1,1];
                end
                surf=SurfaceBSpline(surf_name,[],[],[],surf.X,surf.Y,surf.Z,[],[],u_list,v_list);
                [step_str,object_index,ADVANCED_FACE_index_list(surf_idx)]=surf.getStepNode(object_index);
                fprintf(step_file,step_str);
                fprintf(step_file,'\n');
            end

            % generate OPEN_SHELL
            OPEN_SHELL_index=object_index;
            step_str=[num2str(object_index,'#%d'),' = OPEN_SHELL ',...
                '( ''NONE'', ',...
                '( ',num2str(ADVANCED_FACE_index_list(1:end-1),'#%d, '),' ',num2str(ADVANCED_FACE_index_list(end),'#%d'),' )',...
                ' );\n'];object_index=object_index+1;
            fprintf(step_file,step_str);
            fprintf(step_file,'\n');

            % write model
            SHELL_BASED_SURFACE_MODEL_index=object_index;
            step_str=[num2str(object_index,'#%d'),' = SHELL_BASED_SURFACE_MODEL ',...
                '( ''NONE'', ( ',...
                num2str(OPEN_SHELL_index,'#%d'),' )',...
                ' );\n'];object_index=object_index+1;
            fprintf(step_file,step_str);
            fprintf(step_file,'\n');

            % write end of step file
            writeStepEnd(self,step_file,object_index,step_filename,SHELL_BASED_SURFACE_MODEL_index);

            fclose(step_file);
            clear('step_file');

        end
 
        function drawBody(self,U_num,V_par,W_num)
            % show all surface of body
            %
            surf_total=self.calSurfaceMatrix(U_num,V_par,W_num);

            for surf_index=1:length(surf_total)
                surf=surf_total{surf_index};
                surface(surf.X,surf.Y,surf.Z,LineStyle="none");
            end

            xlabel('x');
            ylabel('y');
            zlabel('z');
            axis equal;
            view(3);
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

% aerodynamci function
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

% auxiliary function
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

function [x_list,fval_list,node_list]=meshAdapt1D(fcn,low_bou,up_bou,torl,max_level,fval_num)
% Binary-tree
% adapt capture function value
% ensure error of linear interplation will less than torl
%
if nargin < 6
    fval_num=1;
end

% node_list which is a matrix store all node
% a node is a array, contain level, index_1, index_2, index_c, node_index_1, node_index_2
% place:
% 1-c-2
% cell:
% 1-2
% if children node is empty, left_index or right_index will be zero
list_add_num=50; % node_list will be extend only when node_list is not enough
node_list=zeros(list_add_num,6,'int64');

% data_list use to sort all float data include coodinate, function value
data_list=zeros(list_add_num,fval_num+1);

% add vertex of cell into data_list first
data_list(1,:)=[low_bou,fcn(low_bou)];
data_list(2,:)=[up_bou,fcn(up_bou)];

% create base root
node_list(1,:)=[0,1,2,0,0,0];

node_num=createNodeTree(1,2); % create node tree from root
node_list=node_list(1:node_num,:);
[x_list,fval_list]=traversalInorder(1); % from small to large get list

% add boundary info
x_list=[data_list(1,1);x_list;data_list(2,1)];
fval_list=[data_list(1,2:end);fval_list;data_list(2,2:end)];

    function node_num=createNodeTree(root_idx,data_num)
        % create node tree from root
        %
        stack=root_idx;
        node_num=root_idx;

        while ~isempty(stack)
            % current node information
            node_idx=stack(end);
            node=node_list(node_idx,:);
            stack=stack(1:end-1);

            if node(1) < max_level
                % judge if linear predict if accptable
                % if not, create children cell
                %
                coord_c=(data_list(node(2),1)+data_list(node(3),1))/2;
                fval_c=fcn(coord_c);
                fval_c_pred=(data_list(node(2),2:end)+data_list(node(3),2:end))/2;
                if any(abs(fval_c-fval_c_pred) > torl)
                    % add 1 new data into data_list
                    data_new_idx=data_num+1;
                    if data_num+1 > size(data_list,1)
                        data_list=[data_list;zeros(list_add_num,fval_num+1)];
                    end
                    data_list(data_new_idx,:)=[coord_c,fval_c];
                    node(4)=data_new_idx;
                    data_num=data_num+1;

                    % add 2 new node to node_list
                    node_new_idx=node_num+[1,2];
                    if node_num+2 > size(node_list,1)
                        node_list=[node_list;zeros(list_add_num,6)];
                    end
                    node_list(node_new_idx,:)=[
                        node(1)+1,node(2),node(4),0,0,0;
                        node(1)+1,node(4),node(3),0,0,0];
                    node([5,6])=node_new_idx;
                    node_num=node_num+2;

                    % add to stack
                    stack=[stack,node_new_idx];
                    node_list(node_idx,:)=node;
                end
            end
        end
    end

    function [x_list,fval_list]=traversalInorder(root_idx)
        % inorder traversal node tree to obtain data
        % inorder will make sure x_list is from small to large
        %
        stack=[];
        x_list=[];
        fval_list=[];
        node_idx=root_idx;

        while ~isempty(stack) || node_idx ~= 0
            while node_idx ~= 0
                stack=[stack,node_idx];

                % node=node.left;
                node_idx=node_list(node_idx,5);
            end

            node_idx=stack(end);
            stack=stack(1:end-1);
            data_idx=node_list(node_idx,4);
            if data_idx ~= 0
                x_list=[x_list;data_list(data_idx,1)];
                fval_list=[fval_list;data_list(data_idx,2:end)];
            end

            % node=node.right;
            node_idx=node_list(node_idx,6);
        end
    end
end
