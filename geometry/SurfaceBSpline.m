classdef SurfaceBSpline < handle
    % B-spline surface
    % define reference to step standard
    %
    properties
        name='';
        u_degree;
        v_degree;

        ctrl_X;
        ctrl_Y;
        ctrl_Z;
        u_ctrl_num;
        v_ctrl_num;

        surface_form='.UNSPECIFIED.'; %
        u_closed='.F.'; % boolean
        v_closed='.F.'; % boolean
        self_intersect='.F.'; % boolean

        u_knot_multi;
        v_knot_multi;
        u_knot_list;
        v_knot_list;
        knot_spec='.UNSPECIFIED.';
    end

    properties
        u_list;
        v_list;

        node_X;
        node_Y;
        node_Z;
        u_node_num;
        v_node_num;
    end

    % define surface
    methods
        function self=SurfaceBSpline(name,point_X,point_Y,point_Z,FLAG_FIT,...
                u_degree,v_degree,u_knot_multi,v_knot_multi,u_knot_list,v_knot_list,...
                u_ctrl_num,v_ctrl_num,U,V)
            % generate BSpline surface by defining control point of fitting point
            %
            % input:
            % name, point_X, point_Y, point_Z(control/fit point), FLAG_FIT,...
            % u_degree(optional), v_degree(optional), ...
            % u_knot_multi(optional), v_knot_multi(optional),...
            % u_knot_list(optional), v_knot_list(optional),...
            % u_ctrl_num(optional), v_ctrl_num(optional),...
            % U(optional), V(optional)
            %
            % output:
            % SurfaceBSpline
            %
            % notice:
            % colume of X, Y, Z is LaWGS format
            % colume of node_X, node_Y, node_Z will by reseve calculate
            % colume direction is u(x), rank direction is v(y)
            %
            if nargin < 15
                V=[];
                if nargin < 14
                    U=[];
                    if nargin < 13
                        v_ctrl_num=[];
                        if nargin < 12
                            u_ctrl_num=[];
                            if nargin < 11
                                v_knot_list=[];
                                if nargin < 10
                                    v_knot_multi = [];
                                    if nargin < 9
                                        u_knot_list=[];
                                        if nargin < 8
                                            u_knot_multi=[];
                                            if nargin < 7
                                                v_degree=[];
                                                if nargin < 6
                                                    u_degree=[];
                                                    if nargin < 5
                                                        FLAG_FIT=[];
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            self.name=name;

            if nargin > 1
                if isempty(FLAG_FIT),FLAG_FIT=false;end

                % check input value size and giving default value
                if FLAG_FIT
                    node_X=point_X;node_Y=point_Y;node_Z=point_Z;
                    [v_node_num,u_node_num]=size(node_X);
                    if any(size(node_Y) ~= [v_node_num,u_node_num]) ||...
                            any(size(node_Z) ~= [v_node_num,u_node_num])
                        error('SurfaceBSpline: size of node_X, node_Y, node_Z do not equal');
                    end

                    if isempty(u_ctrl_num),u_ctrl_num=min(u_node_num,10);end
                    if u_ctrl_num > u_node_num
                        error('SurfaceBSpline: u_control number more than u_node number')
                    end
                    if isempty(v_ctrl_num),v_ctrl_num=min(v_node_num,10);end
                    if v_ctrl_num > v_node_num
                        error('SurfaceBSpline: v_control number more than v_node number')
                    end
                else
                    ctrl_X=point_X;ctrl_Y=point_Y;ctrl_Z=point_Z;
                    [v_ctrl_num,u_ctrl_num]=size(ctrl_X);
                    if isempty(u_degree),u_degree=u_ctrl_num-1;end
                    if isempty(v_degree),v_degree=v_ctrl_num-1;end
                    u_node_num=u_ctrl_num-u_degree+1;
                    v_node_num=v_ctrl_num-v_degree+1;
                    if any(size(ctrl_Y) ~= [v_ctrl_num,u_ctrl_num]) ||...
                            any(size(ctrl_Z) ~= [v_ctrl_num,u_ctrl_num])
                        error('SurfaceBSpline: size of ctrl_X, ctrl_Y, ctrl_Z do not equal');
                    end
                end

                if u_ctrl_num < (u_degree+1) || v_ctrl_num < (v_degree+1)
                    error('SurfaceBSpline: ctrl_num less than degree+1');
                end

                % default value
                if isempty(u_degree),u_degree=u_ctrl_num-1;end
                if isempty(U),U=linspace(0,1,u_node_num);end;U=U(:);
                if isempty(u_knot_multi),u_knot_multi=[u_degree+1,ones(1,u_ctrl_num-u_degree-1),u_degree+1];end
                if isempty(u_knot_list),u_knot_list=interp1(linspace(0,1,u_node_num),U,linspace(0,1,u_ctrl_num-u_degree+1));end
                u_list=baseKnotVec(u_knot_multi,u_knot_list);

                % default value
                if isempty(v_degree),v_degree=v_ctrl_num-1;end
                if isempty(V),V=linspace(0,1,v_node_num);end;V=V(:);
                if isempty(v_knot_multi),v_knot_multi=[v_degree+1,ones(1,v_ctrl_num-v_degree-1),v_degree+1];end
                if isempty(v_knot_list),v_knot_list=interp1(linspace(0,1,v_node_num),V,linspace(0,1,v_ctrl_num-v_degree+1));end
                v_list=baseKnotVec(v_knot_multi,v_knot_list);

                if FLAG_FIT
                    % base on node point list inverse calculate control point list
                    matrix_v=zeros(v_node_num,v_ctrl_num);
                    for node_idx=1:v_node_num
                        v=V(node_idx);
                        for ctrl_idx=1:v_ctrl_num
                            matrix_v(node_idx,ctrl_idx)=baseFcnN(v_list,v,ctrl_idx,v_degree);
                        end
                    end

                    matrix_u=zeros(u_ctrl_num,u_node_num);
                    for node_idx=1:u_node_num
                        u=U(node_idx);
                        for ctrl_idx=1:u_ctrl_num
                            matrix_u(ctrl_idx,node_idx)=baseFcnN(u_list,u,ctrl_idx,u_degree);
                        end
                    end

                    % % old reverse calculate process, maybe no stable...
                    % % while input number of node point is to large.
                    % inv_matrix_hess_v=((matrix_v'*matrix_v)\eye(v_ctrl_num));
                    % inv_matrix_hess_u=(eye(u_ctrl_num)/(matrix_u*matrix_u'));
                    % 
                    % ctrl_X=inv_matrix_hess_v*matrix_v'*node_X*matrix_u'*inv_matrix_hess_u;
                    % ctrl_Y=inv_matrix_hess_v*matrix_v'*node_Y*matrix_u'*inv_matrix_hess_u;
                    % ctrl_Z=inv_matrix_hess_v*matrix_v'*node_Z*matrix_u'*inv_matrix_hess_u;

                    ctrl_X=matrix_v\node_X/matrix_u;
                    ctrl_Y=matrix_v\node_Y/matrix_u;
                    ctrl_Z=matrix_v\node_Z/matrix_u;
                else
                    % generate B spline surface by control point
                end

                % main properties
                self.u_degree=u_degree;
                self.v_degree=v_degree;

                self.ctrl_X=ctrl_X;
                self.ctrl_Y=ctrl_Y;
                self.ctrl_Z=ctrl_Z;
                self.u_ctrl_num=u_ctrl_num;
                self.v_ctrl_num=v_ctrl_num;

                self.u_knot_multi=u_knot_multi;
                self.v_knot_multi=v_knot_multi;
                self.u_knot_list=u_knot_list;
                self.v_knot_list=v_knot_list;

                % calculate node point list
                self.u_list=u_list;
                self.v_list=v_list;

                [U,V]=meshgrid(u_list(u_degree+1:u_ctrl_num+1),v_list(v_degree+1:v_ctrl_num+1));
                [self.node_X,self.node_Y,self.node_Z]=self.calPoint(U,V);
                self.u_node_num=u_node_num;
                self.v_node_num=v_node_num;
            end

        end
    end

    % calculate point
    methods
        function [X,Y,Z,U,V]=calSurface(self,u_param,v_param)
            % generate surface matrix
            %
            % default
            % u_list=linspace(0,1,u_gird_num(default is 20)+1)
            % v_list=linspace(0,1,v_gird_num(default is 20)+1)
            % [U,V]=meshgrid(u_list,v_list), colume is LaWGS format line
            % value_torl=1e-3
            %
            % input:
            % U, V
            % or:
            % u_gird_num, v_gird_num
            % or:
            % value_torl, []
            %
            % output:
            % X,Y,Z (colume is LaWGS format line)
            %
            if nargin < 3
                v_param=[];
                if nargin < 2
                    u_param=[];
                end
            end

            if isempty(u_param),u_param=1e-3;end
            if length(u_param) == 1 && u_param ~= fix(u_param)
                % input is torlance
                value_torl=u_param;
                max_level=50;

                % adapt capture U, V
                low_bou=[0,0];
                up_bou=[1,1];
                [U,V,data_list,~]=meshAdapt2DUV(@(x) coordFcn(self,x),low_bou,up_bou,value_torl,max_level,3);
                X=data_list(:,:,1);
                Y=data_list(:,:,2);
                Z=data_list(:,:,3);
            else
                % input is U, V or u_grid_number, v_grid_number
                if isempty(u_param), u_param=40;end
                if isempty(v_param), v_param=40;end

                if length(u_param) == 1
                    U=[];
                    u_gird_num=u_param;
                else
                    % mean input U matrix
                    U=u_param;
                    u_gird_num=size(U,2)-1;
                end

                if length(v_param) == 1
                    V=[];
                    v_gird_num=v_param;
                else
                    % mean input V matrix
                    V=v_param;
                    v_gird_num=size(V,1)-1;
                end

                % calculate local coordinate matrix
                if isempty(U)
                    U=linspace(0,1,u_gird_num+1);
                    U=repmat(U,v_gird_num+1,1);
                end
                if isempty(V)
                    V=linspace(0,1,v_gird_num+1)';
                    V=repmat(V,1,u_gird_num+1);
                end

                [X,Y,Z]=self.calPoint(U,V);
            end

            function fval=coordFcn(self,x)
                [x,y,z]=self.calPoint(x(1),x(2));
                fval=[x,y,z];
            end
        end

        function [X,Y,Z]=calPoint(self,U,V)
            % according u_x to calculate point
            % u_x_list is u_num x 1 matrix
            % point_list is point_number x dimension matrix
            %
            [rank_num,colume_num]=size(U);
            if any(size(V) ~= [rank_num,colume_num])
                error('SurfaceBSpline.calPoint: size of U_x do not equal to size of V_x');
            end

            [X,Y,Z]=self.calBSpline(U,V);
        end

        function [X,Y,Z]=calBSpline(self,U,V)
            [rank_num,colume_num]=size(U);
            X=zeros(rank_num,colume_num);
            Y=zeros(rank_num,colume_num);
            Z=zeros(rank_num,colume_num);

            N_u_list=zeros(self.u_degree+1,1);
            N_v_list=zeros(1,self.v_degree+1);
            for rank_idx=1:rank_num
                for colume_idx=1:colume_num
                    % local index of u_x in u_list, v_list
                    u_x=U(rank_idx,colume_idx);
                    v_x=V(rank_idx,colume_idx);

                    % [index_end_v,index_end_u]=getIndex(); % y, x

                    idx_end_u=self.u_ctrl_num; % is equal to the section index
                    while idx_end_u > self.u_degree+1 && u_x < self.u_list(idx_end_u)
                        idx_end_u=idx_end_u-1;
                    end
                    idx_end_v=self.v_ctrl_num; % is equal to the section index
                    while idx_end_v > self.v_degree+1 && v_x < self.v_list(idx_end_v)
                        idx_end_v=idx_end_v-1;
                    end

                    idx_start_u=idx_end_u-self.u_degree;
                    idx_start_v=idx_end_v-self.v_degree;

                    % calculate base function
                    for N_idx=1:self.u_degree+1
                        N_u_list(N_idx)=baseFcnN(self.u_list,u_x,N_idx+idx_start_u-1,self.u_degree);
                    end

                    for N_idx=1:self.v_degree+1
                        N_v_list(N_idx)=baseFcnN(self.v_list,v_x,N_idx+idx_start_v-1,self.v_degree);
                    end

                    X(rank_idx,colume_idx)=N_v_list*self.ctrl_X(idx_start_v:idx_end_v,idx_start_u:idx_end_u)*N_u_list;
                    Y(rank_idx,colume_idx)=N_v_list*self.ctrl_Y(idx_start_v:idx_end_v,idx_start_u:idx_end_u)*N_u_list;
                    Z(rank_idx,colume_idx)=N_v_list*self.ctrl_Z(idx_start_v:idx_end_v,idx_start_u:idx_end_u)*N_u_list;
                end
            end
        end

    end

    % calculate coord
    methods
        function [U,V,X,Y,Z]=calCoord(self,X,Y,Z)
            % base on X, Y, Z calculate local coordinate in surface
            %
            XO=X;YO=Y;ZO=Z;geo_torl=100*eps;

            % base on range of node point to preliminary project to U, V
            node_1=[self.node_X(1,1);self.node_Y(1,1);self.node_Z(1,1)];
            node_2=[self.node_X(1,end);self.node_Y(1,end);self.node_Z(1,end)];
            node_3=[self.node_X(end,1);self.node_Y(end,1);self.node_Z(end,1)];
            node_4=[self.node_X(end,end);self.node_Y(end,end);self.node_Z(end,end)];
            node_c=(node_1+node_2+node_3+node_4)/4;
            vector_z=cross(node_4-node_1,node_3-node_2);vector_z=vector_z/norm(vector_z);
            vector_x=(node_2+node_4)/2-node_c;vector_x=vector_x/norm(vector_x);
            vector_y=cross(vector_z,vector_x);
            proj_matrix=[vector_x,vector_y,vector_z]';
            proj_base=node_c;
            proj_node_1=proj_matrix*(node_1-proj_base);
            proj_node_2=proj_matrix*(node_2-proj_base);
            proj_node_3=proj_matrix*(node_3-proj_base);
            proj_node_4=proj_matrix*(node_4-proj_base);
            proj_point=proj_matrix*([X,Y,Z]'-proj_base);
            proj_node=[proj_node_1,proj_node_2,proj_node_3,proj_node_4];

            % project to 2D
            proj_point=proj_point(1:2,:);
            proj_node=proj_node(1:2,:);
            
            % re-deform of uv
            vector_e1=(proj_node(:,2)-proj_node(:,3))/2;
            vector_e1=vector_e1/norm(vector_e1)*max(norm(proj_node(:,2)),norm(proj_node(:,3)));
            vector_e2=(proj_node(:,4)-proj_node(:,1))/2;
            vector_e2=vector_e2/norm(vector_e2)*max(norm(proj_node(:,4)),norm(proj_node(:,1)));
            sqrt2_2=sqrt(2)/2;
            tran_matrix=[sqrt2_2,sqrt2_2;-sqrt2_2,sqrt2_2]/[vector_e1,vector_e2];
            proj_point=tran_matrix*proj_point;

            U=(proj_point(1,:)/2+0.5)';V=(proj_point(2,:)/2+0.5)';

            % use project function to adjust parameter
            [X,Y,Z,U,V]=self.calProject(XO,YO,ZO,U,V,geo_torl);
        end

        function [X,Y,Z,U,V]=calProject(self,XO,YO,ZO,U,V,geo_torl)
            % adjust u, v by Jacobian transformation
            % also can project point to surface
            %
            if nargin < 7,geo_torl=double(eps('single'));end

            iter=0;iter_max=20;

            [X,Y,Z]=self.calPoint(U,V);
            % scatter3(X,Y,Z);
            dX=XO-X;dY=YO-Y;dZ=ZO-Z;
            boolean=(abs(dX)+abs(dY)+abs(dZ)) > geo_torl;
            while any(any(boolean)) && iter < iter_max
                [dX_dU,dY_dU,dZ_dU,dX_dV,dY_dV,dZ_dV]=self.calGradient(U(boolean),V(boolean),X(boolean),Y(boolean),Z(boolean));

                RU_RU=dX_dU.^2+dY_dU.^2+dZ_dU.^2;
                RV_RV=dX_dV.^2+dY_dV.^2+dZ_dV.^2;
                RU_RV=dX_dU.*dX_dV+dY_dU.*dY_dV+dZ_dU.*dZ_dV;
                RU_D=dX_dU.*dX(boolean)+dY_dU.*dY(boolean)+dZ_dU.*dZ(boolean);
                RV_D=dX_dV.*dX(boolean)+dY_dV.*dY(boolean)+dZ_dV.*dZ(boolean);
                RRRR_RR=RU_RU.*RV_RV-(RU_RV).^2;
                dU=(RU_D.*RV_RV-RV_D.*RU_RV)./RRRR_RR;
                dV=(RV_D.*RU_RU-RU_D.*RU_RV)./RRRR_RR;
                dU(isnan(dU) & isinf(dU))=0;
                dV(isnan(dV) & isinf(dV))=0;

                U(boolean)=U(boolean)+dU;
                V(boolean)=V(boolean)+dV;
                U=max(U,0);U=min(U,1);
                V=max(V,0);V=min(V,1);

                [X,Y,Z]=self.calPoint(U,V);
                % scatter3(X,Y,Z);

                dX=XO-X;dY=YO-Y;dZ=ZO-Z;
                iter=iter+1;
                boolean=(abs(dX)+abs(dY)+abs(dZ)) > geo_torl;
            end
        end

        function [dX_dU,dY_dU,dZ_dU,dX_dV,dY_dV,dZ_dV]=calGradient(self,U,V,X,Y,Z)
            % use differ ot calculate gradient
            %
            if nargin < 5
                [X,Y,Z]=self.calPoint(U,V);
            end
            step=100*eps;

            [X_UF,Y_UF,Z_UF]=self.calPoint(U+step,V);
            [X_UB,Y_UB,Z_UB]=self.calPoint(U-step,V);
            Bool_U=(U+step) > 1;
            Bool_D=(U-step) < 0;
            Bool_C=~any([Bool_U,Bool_D],2);
            dX_dU=X;dY_dU=Y;dZ_dU=Z; % allocate memory
            dX_dU(Bool_C)=(X_UF(Bool_C)-X_UB(Bool_C))/2/step;
            dY_dU(Bool_C)=(Y_UF(Bool_C)-Y_UB(Bool_C))/2/step;
            dZ_dU(Bool_C)=(Z_UF(Bool_C)-Z_UB(Bool_C))/2/step;
            dX_dU(Bool_U)=(X(Bool_U)-X_UB(Bool_U))/2/step;
            dY_dU(Bool_U)=(Y(Bool_U)-Y_UB(Bool_U))/2/step;
            dZ_dU(Bool_U)=(Z(Bool_U)-Z_UB(Bool_U))/2/step;
            dX_dU(Bool_D)=(X_UF(Bool_D)-X(Bool_D))/step;
            dY_dU(Bool_D)=(Y_UF(Bool_D)-Y(Bool_D))/step;
            dZ_dU(Bool_D)=(Z_UF(Bool_D)-Z(Bool_D))/step;
            dX_dU=real(dX_dU);dY_dU=real(dY_dU);dZ_dU=real(dZ_dU);

            [X_VF,Y_VF,Z_VF]=self.calPoint(U,V+step);
            [X_VB,Y_VB,Z_VB]=self.calPoint(U,V-step);
            Bool_U=(V+step) > 1;
            Bool_D=(V-step) < 0;
            Bool_C=~any([Bool_U,Bool_D],2);
            dX_dV=X;dY_dV=Y;dZ_dV=Z; % allocate memory
            dX_dV(Bool_C)=(X_VF(Bool_C)-X_VB(Bool_C))/2/step;
            dY_dV(Bool_C)=(Y_VF(Bool_C)-Y_VB(Bool_C))/2/step;
            dZ_dV(Bool_C)=(Z_VF(Bool_C)-Z_VB(Bool_C))/2/step;
            dX_dV(Bool_U)=(X(Bool_U)-X_VB(Bool_U))/2/step;
            dY_dV(Bool_U)=(Y(Bool_U)-Y_VB(Bool_U))/2/step;
            dZ_dV(Bool_U)=(Z(Bool_U)-Z_VB(Bool_U))/2/step;
            dX_dV(Bool_D)=(X_VF(Bool_D)-X(Bool_D))/step;
            dY_dV(Bool_D)=(Y_VF(Bool_D)-Y(Bool_D))/step;
            dZ_dV(Bool_D)=(Z_VF(Bool_D)-Z(Bool_D))/step;
            dX_dV=real(dX_dV);dY_dV=real(dY_dV);dZ_dV=real(dZ_dV);
        end
    end

    % visualizate function
    methods
        function drawSurface(self,axe_hdl,U,V,surface_option,control_option,node_option)
            % draw surface on axes handle
            %
            if nargin < 7
                node_option=[];
                if nargin < 5
                    control_option=[];
                    if nargin < 5
                        surface_option=[];
                        if nargin < 4
                            V=[];
                            if nargin < 3
                                U=[];
                                if nargin < 2
                                    axe_hdl=[];
                                end
                            end
                        end
                    end
                end
            end

            if isempty(axe_hdl),axe_hdl=axes(figure());end

            % default draw option
            if isempty(surface_option)
                surface_option=struct('LineStyle','none');
            end
            if isempty(node_option)
                node_option=struct('Marker','o','MarkerEdgeColor','b','LineStyle','none','FaceAlpha',0);
            end
            if isempty(control_option)
                control_option=struct('Marker','s','MarkerEdgeColor','r','EdgeColor','r','LineStyle','--','FaceAlpha',0);
            end

            % calculate point on surface
            [X,Y,Z]=calSurface(self,U,V);

            % plot surface
            surface(axe_hdl,X,Y,Z,surface_option);
            surface(axe_hdl,self.node_X,self.node_Y,self.node_Z,node_option);
            surface(axe_hdl,self.ctrl_X,self.ctrl_Y,self.ctrl_Z,control_option);
            view(3);
            xlabel('x');
            ylabel('y');
            zlabel('z');

            %             axis equal
            %             x_range=xlim();
            %             y_range=ylim();
            %             z_range=zlim();
            %             center=[mean(x_range),mean(y_range),mean(z_range)];
            %             range=max([x_range(2)-x_range(1),y_range(2)-y_range(1),z_range(2)-z_range(1)])/2;
            %             xlim([center(1)-range,center(1)+range]);
            %             ylim([center(2)-range,center(2)+range]);
            %             zlim([center(3)-range,center(3)+range]);
        end

        function [step_str,object_index,ADVANCED_FACE]=getStep(self,object_index)
            % write BSpline into step file
            %
            if nargin < 2,object_index=1;end
            out_name=self.name;
            if isempty(out_name),out_name='NONE';end

            u_num=self.u_ctrl_num;
            v_num=self.v_ctrl_num;

            % generate CARTESIAN_POINT
            CARTESIAN_POINT=object_index;

            str_control=[];
            for u_idx=1:u_num
                for v_idx=1:v_num
                    str=[num2str(object_index,'#%d ='),' CARTESIAN_POINT ',...
                        '( ',...
                        '''NONE'', ',...
                        '( ',num2str(self.ctrl_X(v_idx,u_idx),'%.16f'),', ',...
                        num2str(self.ctrl_Y(v_idx,u_idx),'%.16f'),', ',...
                        num2str(self.ctrl_Z(v_idx,u_idx),'%.16f'),' )',...
                        ' ) ;\n'];
                    str_control=[str_control,str];
                    object_index=object_index+1;
                end
            end

            % generate B_SPLINE_CURVE_WITH_KNOTS
            B_SPLINE_CURVE_WITH_KNOTS=object_index;
            str_curve=[];
            str_curve=[str_curve,stepCurve((1:v_num) +CARTESIAN_POINT-1,...
                self.v_degree,self.v_knot_multi,self.v_knot_list)];
            str_curve=[str_curve,stepCurve(((1:u_num)*v_num) +CARTESIAN_POINT-1,...
                self.u_degree,self.u_knot_multi,self.u_knot_list)];
            str_curve=[str_curve,stepCurve(((v_num*u_num):-1:(v_num*u_num-v_num+1)) +CARTESIAN_POINT-1,...
                self.v_degree,self.v_knot_multi,self.v_knot_list)];
            str_curve=[str_curve,stepCurve((((u_num-1):-1:0)*v_num+1) +CARTESIAN_POINT-1,...
                self.u_degree,self.u_knot_multi,self.u_knot_list)];

            % generate VERTEX_POINT
            VERTEX_POINT=object_index;
            str_vertex=[];
            str_vertex=[str_vertex,stepVertex(1 +CARTESIAN_POINT-1)];
            str_vertex=[str_vertex,stepVertex(v_num +CARTESIAN_POINT-1)];
            str_vertex=[str_vertex,stepVertex(v_num*u_num +CARTESIAN_POINT-1)];
            str_vertex=[str_vertex,stepVertex(v_num*u_num-v_num+1 +CARTESIAN_POINT-1)];

            % generate EDGE_CURVE
            EDGE_CURVE=object_index;
            str_edge=[];
            str_edge=[str_edge,stepEdge(VERTEX_POINT,VERTEX_POINT+1,B_SPLINE_CURVE_WITH_KNOTS)];
            str_edge=[str_edge,stepEdge(VERTEX_POINT+1,VERTEX_POINT+2,B_SPLINE_CURVE_WITH_KNOTS+1)];
            str_edge=[str_edge,stepEdge(VERTEX_POINT+2,VERTEX_POINT+3,B_SPLINE_CURVE_WITH_KNOTS+2)];
            str_edge=[str_edge,stepEdge(VERTEX_POINT+3,VERTEX_POINT,B_SPLINE_CURVE_WITH_KNOTS+3)];

            % generate ORIENTED_EDGE
            ORIENTED_EDGE=object_index;
            str_oriented=[];
            str_oriented=[str_oriented,stepOriented(EDGE_CURVE)];
            str_oriented=[str_oriented,stepOriented(EDGE_CURVE+1)];
            str_oriented=[str_oriented,stepOriented(EDGE_CURVE+2)];
            str_oriented=[str_oriented,stepOriented(EDGE_CURVE+3)];

            % generate EDGE_LOOP
            EDGE_LOOP=object_index;
            str_loop=[num2str(object_index,'#%d ='),' EDGE_LOOP',...
                ' ( ',...
                '''NONE''',', ',...
                '( ',num2str(ORIENTED_EDGE,'#%d'),', ',...
                num2str(ORIENTED_EDGE+1,'#%d'),', ',...
                num2str(ORIENTED_EDGE+2,'#%d'),', ',...
                num2str(ORIENTED_EDGE+3,'#%d'),') ',...
                ' ) ;\n'];
            object_index=object_index+1;

            % generate FACE_OUTER_BOUND
            FACE_OUTER_BOUND=object_index;
            str_outer=[num2str(object_index,'#%d ='),' FACE_OUTER_BOUND',...
                ' ( ',...
                '''NONE''',', ',...
                num2str(EDGE_LOOP,'#%d'),', ',...
                '.T.',...
                ' ) ;\n'];
            object_index=object_index+1;

            % generate B_SPLINE_SURFACE_WITH_KNOTS
            B_SPLINE_SURFACE_WITH_KNOTS=object_index;
            str_surface=[num2str(object_index,'#%d ='),' B_SPLINE_SURFACE_WITH_KNOTS',...
                ' ( ',...
                '''',out_name,'''',', ',...
                num2str(self.u_degree,'%d'),', ',num2str(self.v_degree,'%d'),', \n'];
            str_surface=[str_surface,'('];
            for u_bias=0:v_num:v_num*(u_num-2)
                str_surface=[str_surface,'( ',...
                    num2str((v_num+u_bias) +CARTESIAN_POINT-1,'#%d'),...
                    num2str((((v_num-1):-1:1)+u_bias) +CARTESIAN_POINT-1,', #%d'),...
                    ' ),\n'];
            end
            u_bias=v_num*(u_num-1);
            str_surface=[str_surface,'( ',...
                num2str((v_num+u_bias) +CARTESIAN_POINT-1,'#%d'),...
                num2str((((v_num-1):-1:1)+u_bias) +CARTESIAN_POINT-1,', #%d'),...
                ' )'];
            str_surface=[str_surface,'),\n'];
            str_surface=[str_surface,'.UNSPECIFIED.',', ','.F.',', ','.F.',', ','.F.',',\n',...
                '( ',num2str(self.u_knot_multi(1),'%d'),num2str(self.u_knot_multi(2:end),', %d'),' ),\n',...
                '( ',num2str(self.v_knot_multi(1),'%d'),num2str(self.v_knot_multi(2:end),', %d'),' ),\n',...
                '( ',num2str(self.u_knot_list(1),'%.16f'),num2str(self.u_knot_list(2:end),', %.16f'),' ),\n',...
                '( ',num2str(self.v_knot_list(1),'%.16f'),num2str(self.v_knot_list(2:end),', %.16f'),' ),\n',...
                '.UNSPECIFIED.',...
                ') ;\n'];
            object_index=object_index+1;

            % generate ADVANCED_FACE
            ADVANCED_FACE=object_index;
            str_face=[num2str(object_index,'#%d ='),' ADVANCED_FACE',...
                ' ( ',...
                '''',out_name,'''',', ',...
                '( ',num2str(FACE_OUTER_BOUND,'#%d'),' )',', '...
                num2str(B_SPLINE_SURFACE_WITH_KNOTS,'#%d'),', ',...
                '.T.',...
                ' ) ;\n'];
            object_index=object_index+1;

            % generate all
            step_str=[str_control,'\n',str_curve,'\n',str_vertex,'\n',...
                str_edge,'\n',str_oriented,'\n',str_loop,'\n',str_outer,'\n',...
                str_surface,'\n',str_face,'\n'];

            function str_curve=stepCurve(point_index,degree,knot_multi,knot_list)
                str_curve=[num2str(object_index,'#%d ='),' B_SPLINE_CURVE_WITH_KNOTS',...
                    ' ( ',...
                    '''NONE''',', ',...
                    num2str(degree,'%d'),', \n',...
                    '( ',num2str(point_index(1),'#%d'),num2str(point_index(2:end),', #%d'),' ),\n',...
                    '.UNSPECIFIED., .F., .F.,\n',...
                    '( ',num2str(knot_multi(1),'%d'),num2str(knot_multi(2:end),', %d'),' ),\n',...
                    '( ',num2str(knot_list(1),'%.16f'),num2str(knot_list(2:end),', %.16f'),' ),\n',...
                    '.UNSPECIFIED.',...
                    ' ) ;\n'];
                object_index=object_index+1;
            end

            function str_vertex=stepVertex(point_index)
                str_vertex=[num2str(object_index,'#%d ='),' VERTEX_POINT',...
                    ' ( ',...
                    '''NONE''',', ',...
                    num2str(point_index,'#%d'),...
                    ' ) ;\n'];
                object_index=object_index+1;
            end

            function str_edge=stepEdge(start_vertex,end_vertex,curve_index)
                str_edge=[num2str(object_index,'#%d ='),' EDGE_CURVE',...
                    ' ( ',...
                    '''NONE''',', ',...
                    num2str(start_vertex,'#%d'),', ',...
                    num2str(end_vertex,'#%d'),', ',...
                    num2str(curve_index,'#%d'),', ',...
                    '.T.',...
                    ' ) ;\n'];
                object_index=object_index+1;
            end

            function str_oriented=stepOriented(edge_index)
                str_oriented=[num2str(object_index,'#%d ='),' ORIENTED_EDGE',...
                    ' ( ',...
                    '''NONE''',', ',...
                    '*',', ','*',', ',...
                    num2str(edge_index,'#%d'),', ',...
                    '.T.',...
                    ' ) ;\n'];
                object_index=object_index+1;
            end
        end
    end
end

%% common function

function [U,V,Fval,node_list]=meshAdapt2DUV(fcn,low_bou,up_bou,torl,max_level,fval_num)
% 2D Omni-Tree
% adapt capture 2 dimemsion function value
% ensure error of linear interplation will less than torl
%
if nargin < 6
    fval_num=1;
end

% node_list which is a matrix store all node
% a node is a array, contain:
% level, idx_1-8(index of data_list), idx_c, node_index_1-4
% place:
% 3-8-4 or 3-4 or 3-8-4
% 1-5-2    6-7    6-c-7
%          1-2    1-5-2
% node:
% 1-2 or 3 or 3-4
%        1    1-2
% if children node is empty, left_index or right_index will be zero
list_add_num=50; % list will be extend only when list is not enough
node_list=zeros(list_add_num,14,'int64');
% data_list use to sort all float data include coodinate, function value
data_list=zeros(list_add_num,fval_num+2);

% add vertex of cell into data_list first
bou_1=[low_bou(1),low_bou(2)];
data_list(1,:)=[bou_1,fcn(bou_1)];
bou_2=[up_bou(1),low_bou(2)];
data_list(2,:)=[bou_2,fcn(bou_2)];
bou_3=[low_bou(1),up_bou(2)];
data_list(3,:)=[bou_3,fcn(bou_3)];
bou_4=[up_bou(1),up_bou(2)];
data_list(4,:)=[bou_4,fcn(bou_4)];

% create base root
node_list(1,:)=[0,1,2,3,4,0,0,0,0,0,0,0,0,0];

% create node tree from root
[node_num,data_num]=createNodeTree(1,4);
node_list=node_list(1:node_num,:);
data_list=data_list(1:data_num,:);

% generate U and V
[u_list,~,u_idx]=unique(data_list(:,1));
[v_list,~,v_idx]=unique(data_list(:,2));
idx_list=[u_idx,v_idx];

[U,V]=meshgrid(u_list,v_list);

% local data to Fval
Fval=nan(length(v_list),length(u_list),fval_num);

for data_index=1:data_num
    if all([idx_list(data_index,2),idx_list(data_index,1)] ~= 0)
        Fval(idx_list(data_index,2),idx_list(data_index,1),:)=data_list(data_index,3:end);
    end
end

% fit nan data
for rank_idx=1:length(v_list)
    for colume_idx=1:length(u_list)
        if isnan(Fval(rank_idx,colume_idx,1))
            Fval(rank_idx,colume_idx,:)=fcn([U(rank_idx,colume_idx),V(rank_idx,colume_idx)]);
        end
    end
end

    function [node_num,data_num]=createNodeTree(root_idx,data_num)
        % create quad tree
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
                [coord_c,coord_5,coord_6,coord_7,coord_8,...
                    fval_c,fval_5,fval_6,fval_7,fval_8]=calCell(node(2),node(3),node(4),node(5));
                [fval_pred_c,fval_pred_5,fval_pred_6,...
                    fval_pred_7,fval_pred_8]=calCellPred(node(2),node(3),node(4),node(5));

                % check u direction
                if any(abs(fval_5-fval_pred_5) > torl) || any(abs(fval_8-fval_pred_8) > torl)
                    add_u_flag=true(1);
                else
                    add_u_flag=false(1);
                end

                % check v direction
                if any(abs(fval_6-fval_pred_6) > torl) || any(abs(fval_7-fval_pred_7) > torl)
                    add_v_flag=true(1);
                else
                    add_v_flag=false(1);
                end

                % check center place
                if ~add_u_flag && ~add_v_flag && any(abs(fval_c-fval_pred_c) > torl)
                    add_u_flag=true(1);
                    add_v_flag=true(1);
                end

                if add_u_flag && add_v_flag
                    % add 5 data into data_list
                    data_new_idx=data_num+(1:5);
                    if data_num+5 > size(data_list,1)
                        data_list=[data_list;zeros(list_add_num,fval_num+2)];
                    end
                    data_list(data_new_idx,:)=[
                        coord_5,fval_5;
                        coord_6,fval_6;
                        coord_7,fval_7;
                        coord_8,fval_8;
                        coord_c,fval_c;];
                    node(6:10)=data_new_idx;
                    data_num=data_num+5;

                    % add 4 new node to node_list
                    node_new_idx=node_num+(1:4);
                    if node_num+4 > size(node_list,1)
                        node_list=[node_list;zeros(list_add_num,14)];
                    end
                    node_list(node_new_idx,:)=[...
                        node(1)+1,node(2),node(6),node(7),node(10),0,0,0,0,0,0,0,0,0;
                        node(1)+1,node(6),node(3),node(10),node(8),0,0,0,0,0,0,0,0,0;
                        node(1)+1,node(7),node(10),node(4),node(9),0,0,0,0,0,0,0,0,0;
                        node(1)+1,node(10),node(8),node(9),node(5),0,0,0,0,0,0,0,0,0;];
                    node(11:14)=node_new_idx;
                    node_num=node_num+4;

                elseif add_u_flag
                    % add 2 data into data_list
                    data_new_idx=data_num+(1:2);
                    if data_num+2 > size(data_list,1)
                        data_list=[data_list;zeros(list_add_num,fval_num+2)];
                    end
                    data_list(data_new_idx,:)=[
                        coord_5,fval_5;
                        coord_8,fval_8;];
                    node([6,9])=data_new_idx;
                    data_num=data_num+2;

                    % add 2 new node to node_list
                    node_new_idx=node_num+(1:2);
                    if node_num+2 > size(node_list,1)
                        node_list=[node_list;zeros(list_add_num,14)];
                    end
                    node_list(node_new_idx,:)=[...
                        node(1)+1,node(2),node(6),node(4),node(9),0,0,0,0,0,0,0,0,0;
                        node(1)+1,node(6),node(3),node(9),node(5),0,0,0,0,0,0,0,0,0;];
                    node([11,12])=node_new_idx;
                    node_num=node_num+2;

                elseif add_v_flag
                    % add 2 data into data_list
                    data_new_idx=data_num+(1:2);
                    if data_num+2 > size(data_list,1)
                        data_list=[data_list;zeros(list_add_num,fval_num+2)];
                    end
                    data_list(data_new_idx,:)=[
                        coord_6,fval_6;
                        coord_7,fval_7;];
                    node([7,8])=data_new_idx;
                    data_num=data_num+2;

                    % add 2 new node to node_list
                    node_new_idx=node_num+(1:2);
                    if node_num+2 > size(node_list,1)
                        node_list=[node_list;zeros(list_add_num,14)];
                    end
                    node_list(node_num+(1:2),:)=[...
                        node(1)+1,node(2),node(3),node(7),node(8),0,0,0,0,0,0,0,0,0;
                        node(1)+1,node(7),node(8),node(4),node(5),0,0,0,0,0,0,0,0,0;];
                    node([11,13])=node_new_idx;
                    node_num=node_num+2;
                else
                    node_new_idx=[];
                end

                % add to stack
                stack=[stack,node_new_idx];
                node_list(node_idx,:)=node;
            end
        end
    end

    function [coord_c,coord_5,coord_6,coord_7,coord_8,...
            fval_c,fval_5,fval_6,fval_7,fval_8]=calCell(vidx_1,vidx_2,vidx_3,vidx_4)
        % calculate index of c, 5, 6, 7, 8 place function value
        % abbreviation:
        % vidx: vertex index
        %
        coord_c=(data_list(vidx_1,[1,2])+data_list(vidx_4,[1,2]))/2;
        fval_c=fcn(coord_c);

        coord_5=(data_list(vidx_1,[1,2])+data_list(vidx_2,[1,2]))/2;
        fval_5=fcn(coord_5);
        coord_6=(data_list(vidx_1,[1,2])+data_list(vidx_3,[1,2]))/2;
        fval_6=fcn(coord_6);
        coord_7=(data_list(vidx_2,[1,2])+data_list(vidx_4,[1,2]))/2;
        fval_7=fcn(coord_7);
        coord_8=(data_list(vidx_3,[1,2])+data_list(vidx_4,[1,2]))/2;
        fval_8=fcn(coord_8);
    end

    function [fval_pred_c,fval_pred_5,fval_pred_6,...
            fval_pred_7,fval_pred_8]=calCellPred(vidx_1,vidx_2,vidx_3,vidx_4)
        % calculate index of c, 5, 6, 7, 8 place linear predict function value
        %
        fval_pred_c=(data_list(vidx_1,3:end)+data_list(vidx_2,3:end)+data_list(vidx_3,3:end)+data_list(vidx_4,3:end))/4;
        fval_pred_5=(data_list(vidx_1,3:end)+data_list(vidx_2,3:end))/2;
        fval_pred_6=(data_list(vidx_1,3:end)+data_list(vidx_3,3:end))/2;
        fval_pred_7=(data_list(vidx_2,3:end)+data_list(vidx_4,3:end))/2;
        fval_pred_8=(data_list(vidx_3,3:end)+data_list(vidx_4,3:end))/2;
    end

end
