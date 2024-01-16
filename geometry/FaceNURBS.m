classdef FaceNURBS < handle
    % Non-Uniform Rational B-Splines Face
    % define reference to step standard
    %
    properties % Explicit Attributes
        name='';
        UDegree; % u_degree
        VDegree; % v_degree
        Poles; % control_points_list
        face_form='.UNSPECIFIED.'; % face_form
        UPeriodic='.F.'; % closed_surface (boolean)
        VPeriodic='.F.'; % closed_surface (boolean)
        self_intersect='.F.'; % self_intersect (boolean)
        Weights; % weights_data

        UMults; % u_knot_multiplicities
        VMults; % v_knot_multiplicities
        UKnots; % UKnots
        VKnots; % VKnots
        knot_spec='.UNSPECIFIED.';

        u_list; % u_knot_vector
        v_list; % v_knot_vector

        reverse=false;
    end

    % define Face
    methods
        function self=FaceNURBS(name,Poles,UDegree,VDegree,...
                UMults,VMults,UKnots,VKnots,Weights)
            % generate Non-Uniform Rational B-Splines Face
            %
            % input:
            % name (str):
            % Poles (matrix): control point, v_pole_num x u_pole_num x dimension matrix
            % UDegree (matrix): optional input
            % VDegree (matrix): optional input
            % UMults (matrix): optional input
            % VMults (matrix): optional input
            % UKnots (matrix): optional input
            % VKnots (matrix): optional input
            % Weights (matrix): optional input
            %
            % output:
            % FaceNURBS
            %
            % notice:
            % colume of X, Y, Z is LaWGS format
            % colume of node_X, node_Y, node_Z will by reseve calculate
            % colume direction is u(x), rank direction is v(y)
            %
            if nargin < 9
                Weights=[];
                if nargin < 8
                    VKnots=[];
                    if nargin < 7
                        VMults=[];
                        if nargin < 6
                            UKnots=[];
                            if nargin < 5
                                UMults=[];
                                if nargin < 4
                                    VDegree=[];
                                    if nargin < 3
                                        UDegree=[];
                                    end
                                end
                            end
                        end
                    end
                end
            end
            self.name=name;

            if nargin > 1
                [v_pole_num,u_pole_num,~]=size(Poles);

                % default value
                if isempty(UDegree),UDegree=u_pole_num-1;end
                if isempty(VDegree),VDegree=v_pole_num-1;end
                if isempty(UMults),UMults=[UDegree+1,ones(1,u_pole_num-UDegree-1),UDegree+1];end
                if isempty(VMults),VMults=[VDegree+1,ones(1,v_pole_num-VDegree-1),VDegree+1];end
                if isempty(UKnots),UKnots=linspace(0,1,u_pole_num-UDegree+1);end
                if isempty(VKnots),VKnots=linspace(0,1,v_pole_num-VDegree+1);end
                if isempty(Weights), Weights=ones(v_pole_num,u_pole_num);end
                u_list=baseKnotVec(UMults,UKnots);
                v_list=baseKnotVec(VMults,VKnots);

                if u_pole_num < (UDegree+1) || v_pole_num < (VDegree+1)
                    error('FaceNURBS: pole_num less than Degree+1');
                end

                if length(u_list) ~= u_pole_num+UDegree+1 || length(v_list) ~= v_pole_num+VDegree+1
                    error('FaceNURBS: knot_num is not equal to pole_num+Degree+1');
                end

                % Explicit Attributes
                self.UDegree=UDegree;
                self.VDegree=VDegree;
                self.Poles=Poles;
                self.UMults=UMults;
                self.VMults=VMults;
                self.UKnots=UKnots;
                self.VKnots=VKnots;
                self.Weights=Weights;

                % Derived Attributes
                self.u_list=u_list;
                self.v_list=v_list;
            end
        end

        function Point=calNURBS(self,U,V)
            % calculate point on Non-Uniform Rational B-Splines Surface
            %
            [rank_num,colume_num]=size(U);
            Point=zeros(rank_num,colume_num,size(self.Poles,3));

            N_u_list=zeros(1,self.UDegree+1);
            N_v_list=zeros(self.VDegree+1,1);
            for rank_idx=1:rank_num
                for colume_idx=1:colume_num
                    % local index of u_x in u_list, v_list
                    u_x=U(rank_idx,colume_idx);
                    v_x=V(rank_idx,colume_idx);

                    % [index_end_v,index_end_u]=getIndex(); % y, x

                    idx_end_u=size(self.Poles,2); % is equal to the section index
                    while idx_end_u > self.UDegree+1 && u_x < self.u_list(idx_end_u)
                        idx_end_u=idx_end_u-1;
                    end
                    idx_end_v=size(self.Poles,1); % is equal to the section index
                    while idx_end_v > self.VDegree+1 && v_x < self.v_list(idx_end_v)
                        idx_end_v=idx_end_v-1;
                    end

                    idx_start_u=idx_end_u-self.UDegree;
                    idx_start_v=idx_end_v-self.VDegree;

                    % calculate base function
                    for N_idx=1:self.UDegree+1
                        N_u_list(N_idx)=baseFcnN(self.u_list,u_x,N_idx+idx_start_u-1,self.UDegree);
                    end

                    for N_idx=1:self.VDegree+1
                        N_v_list(N_idx)=baseFcnN(self.v_list,v_x,N_idx+idx_start_v-1,self.VDegree);
                    end
                    NW_list=N_u_list.*self.Weights(idx_start_v:idx_end_v,idx_start_u:idx_end_u).*N_v_list;

                    Point(rank_idx,colume_idx,:)=sum(NW_list.*self.Poles(idx_start_v:idx_end_v,idx_start_u:idx_end_u,:),[1,2])/sum(NW_list,'all');
                end
            end
        end
    end

    % control Face
    methods
        function addDegree(self,Udeg_tar,Vdeg_tar)
            % increase face degree
            %
            Ctrl=cat(3,self.Poles,self.Weights);
            [v_num,u_num,~]=size(self.Poles);

            % modify V degree
            if self.VDegree < Vdeg_tar
                Ctrl=permute(Ctrl,[1 2 3]);
                Ctrl=reshape(Ctrl,v_num,4*u_num);
                [Ctrl,self.VMults]=GeomApp.addDegree(self.VDegree,Ctrl,self.VMults,self.VKnots,Vdeg_tar);
                % Ctrl=Ctrl./Ctrl(:,end);
                v_num=size(Ctrl,1);
                Ctrl=reshape(Ctrl,[v_num,u_num,4]);
                Ctrl=permute(Ctrl,[1 2 3]);
            end

            % modify U degree
            if self.UDegree < Udeg_tar
                Ctrl=permute(Ctrl,[2 1 3]);
                Ctrl=reshape(Ctrl,u_num,4*v_num);
                [Ctrl,self.UMults]=GeomApp.addDegree(self.UDegree,Ctrl,self.UMults,self.UKnots,Udeg_tar);
                % Ctrl=Ctrl./Ctrl(:,end);
                u_num=size(Ctrl,1);
                Ctrl=reshape(Ctrl,[u_num,v_num,4]);
                Ctrl=permute(Ctrl,[2 1 3]);
            end

            self.Poles=Ctrl(:,:,1:end-1);
            self.Weights=Ctrl(:,:,end);

            self.u_list=baseKnotVec(self.UMults,self.UKnots);
            self.v_list=baseKnotVec(self.VMults,self.VKnots);
        end

        function insertKnot(self,U_ins,V_ins)
            % insert knot to edge
            %
            Ctrl=cat(3,self.Poles,self.Weights);
            U=self.u_list;
            V=self.v_list;
            Vdeg=self.VDegree;
            Udeg=self.UDegree;
            [v_num,u_num,~]=size(self.Poles);

            % Insert knots along the v direction
            if ~isempty(V_ins)
                Ctrl=permute(Ctrl,[1 2 3]);
                Ctrl=reshape(Ctrl,v_num,4*u_num);
                [Ctrl,V]=GeomApp.insertKnot(Vdeg,Ctrl,V,V_ins);
                v_num=size(Ctrl,1);
                Ctrl=reshape(Ctrl,[v_num,u_num,4]);
                Ctrl=permute(Ctrl,[1 2 3]);
            end

            % Insert knots along the u direction
            if ~isempty(U_ins)
                Ctrl=permute(Ctrl,[2 1 3]);
                Ctrl=reshape(Ctrl,u_num,4*v_num);
                [Ctrl,U]=GeomApp.insertKnot(Udeg,Ctrl,U,U_ins);
                u_num=size(Ctrl,1);
                Ctrl=reshape(Ctrl,[u_num,v_num,4]);
                Ctrl=permute(Ctrl,[2 1 3]);
            end

            self.Poles=Ctrl(:,:,1:end-1);
            self.Weights=Ctrl(:,:,end)';

            self.VKnots=unique(V);
            self.VMults=ones(size(self.VKnots));
            for k_idx=1:length(self.VKnots)
                self.VMults(k_idx)=sum(V == self.VKnots(k_idx));
            end

            self.UKnots=unique(U);
            self.UMults=ones(size(self.UKnots));
            for k_idx=1:length(self.UKnots)
                self.UMults(k_idx)=sum(U == self.UKnots(k_idx));
            end

            self.u_list=U;
            self.v_list=V;
        end
    end

    % calculate point
    methods
        function [Point,U,V]=calFace(self,u_param,v_param)
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

            if isempty(u_param)
                pnt=reshape(self.calPoint([0;0;1;1],[0;1;0;1]),4,[]);
                bou_min=min(pnt,[],1);bou_max=max(pnt,[],1);
                u_param=2^-8*mean(bou_max-bou_min);
            end
            if length(u_param) == 1 && u_param ~= fix(u_param)
                % input is torlance
                value_torl=u_param;min_level=2;max_level=8;

                % adapt capture U, V
                low_bou=[0,0];
                up_bou=[1,1];
                dim=numel(self.calPoint(0,0));
                [U,V,data_list,~]=GeomApp.meshAdapt2DUV(@(x) reshape(self.calPoint(x(:,1),x(:,2)),[],dim),low_bou,up_bou,value_torl,min_level,max_level,dim);
                Point=data_list;
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

                Point=self.calPoint(U,V);
            end
        end

        function Point=calPoint(self,U,V)
            % according u_x to calculate point
            % u_x_list is u_num x 1 matrix
            % point_list is point_number x dimension matrix
            %
            [rank_num,colume_num]=size(U);
            if any(size(V) ~= [rank_num,colume_num])
                error('FaceNURBS.calPoint: size of U do not equal to size of V');
            end

            Point=self.calNURBS(U,V);
        end
    end

    % calculate coord
    methods
        function [U,V,Point]=calCoord(self,Point)
            % base on X, Y, Z calculate local coordinate in surface
            %
            Point0=Point;geo_torl=sqrt(eps);
            v_num=size(Point,1);u_num=size(Point,2);

            % generate rough mesh to initialize pre coord
            [Point_base,U_base,V_base]=self.calFace(5,5);
            U=zeros(size(Point,1),size(Point,2));
            V=zeros(size(Point,1),size(Point,2));
            for v_idx=1:v_num
                for u_idx=1:u_num
                    point=Point(v_idx,u_idx,:);
                    dis_abs=sum(abs(Point_base-point),3);
                    [~,idx]=min(dis_abs,[],"all");
                    U(v_idx,u_idx)=U_base(idx(1));
                    V(v_idx,u_idx)=V_base(idx(1));
                end
            end

            % use project function to adjust parameter
            [Point,U,V]=self.calProject(Point0,U,V,geo_torl);
        end

        function [Point,U,V]=calProject(self,Point0,U,V,geo_torl)
            % adjust u, v by Jacobian transformation
            % also can project point to surface
            %
            if nargin < 7,geo_torl=100*eps;end

            iter=0;iter_max=50;

            [v_num,u_num,~]=size(Point0);
            Point0=reshape(Point0,[],1,3);U=U(:);V=V(:);
            
            Point=self.calPoint(U,V);
            dPoint=Point0-Point;
            idx=find(sum(abs(dPoint),3) > geo_torl);
            % scatter3(Point(:,:,1),Point(:,:,2),Point(:,:,3));
            while ~isempty(idx) && iter < iter_max
                [dPoint_dU,dPoint_dV]=self.calGradient(U(idx),V(idx));

                RU_RU=sum(dPoint_dU.^2,3);
                RV_RV=sum(dPoint_dV.^2,3);
                RU_RV=sum(dPoint_dU.*dPoint_dV,3);
                RU_D=0;RV_D=0;
                for d_idx=1:size(self.Poles,3)
                    RU_D=RU_D+dPoint_dU(:,:,d_idx).*dPoint(idx,d_idx);
                    RV_D=RV_D+dPoint_dV(:,:,d_idx).*dPoint(idx,d_idx);
                end
                RRRR_RR=RU_RU.*RV_RV-(RU_RV).^2;
                dU=(RU_D.*RV_RV-RV_D.*RU_RV)./RRRR_RR;
                dV=(RV_D.*RU_RU-RU_D.*RU_RV)./RRRR_RR;
                dU(isnan(dU) | isinf(dU))=0;
                dV(isnan(dV) | isinf(dV))=0;

                U(idx)=U(idx)+dU;
                V(idx)=V(idx)+dV;
                U=max(U,0);U=min(U,1);
                V=max(V,0);V=min(V,1);

                iter=iter+1;

                Point=self.calPoint(U,V);
                dPoint=Point0-Point;
                idx=find(sum(abs(dPoint),3) > geo_torl);
                % scatter3(Point(:,:,1),Point(:,:,2),Point(:,:,3));
            end

            Point=reshape(Point,v_num,u_num,size(self.Poles,3));
            U=reshape(U,v_num,u_num);
            V=reshape(V,v_num,u_num);
        end

        function [dPoint_dU,dPoint_dV]=calGradient(self,U,V,step)
            % use differ to calculate gradient
            %
            if nargin < 4 || isempty(step),step=100*eps;end

            dim=numel(self.calPoint(0,0));
            [v_num,u_num]=size(U);U=U(:);V=V(:);

            Point=self.calPoint(U,V);Point=reshape(Point,[],dim);

            Point_UF=self.calPoint(U+step,V);Point_UF=reshape(Point_UF,[],dim);
            Point_UB=self.calPoint(U-step,V);Point_UB=reshape(Point_UB,[],dim);
            Bool_U=(U+step) > 1;
            Bool_D=(U-step) < 0;
            Bool(:,1)=Bool_U;Bool(:,2)=Bool_D;
            Bool_C=~any(Bool,2);
            dPoint_dU=zeros(v_num*u_num,dim); % allocate memory
            dPoint_dU(Bool_C,:)=(Point_UF(Bool_C,:)-Point_UB(Bool_C,:))/2/step;
            dPoint_dU(Bool_U,:)=(Point(Bool_U,:)-Point_UB(Bool_U,:))/step;
            dPoint_dU(Bool_D,:)=(Point_UF(Bool_D,:)-Point(Bool_D,:))/step;
            dPoint_dU=real(dPoint_dU);
            dPoint_dU=reshape(dPoint_dU,v_num,u_num,dim);

            Point_VF=self.calPoint(U,V+step);Point_VF=reshape(Point_VF,[],dim);
            Point_VB=self.calPoint(U,V-step);Point_VB=reshape(Point_VB,[],dim);
            Bool_U=(V+step) > 1;
            Bool_D=(V-step) < 0;
            Bool(:,1)=Bool_U;Bool(:,2)=Bool_D;
            Bool_C=~any(Bool,2);
            dPoint_dV=zeros(v_num*u_num,dim); % allocate memory
            dPoint_dV(Bool_C,:)=(Point_VF(Bool_C,:)-Point_VB(Bool_C,:))/2/step;
            dPoint_dV(Bool_U,:)=(Point(Bool_U,:)-Point_VB(Bool_U,:))/step;
            dPoint_dV(Bool_D,:)=(Point_VF(Bool_D,:)-Point(Bool_D,:))/step;
            dPoint_dV=real(dPoint_dV);
            dPoint_dV=reshape(dPoint_dV,v_num,u_num,dim);
        end
    end

    % visualizate surface
    methods
        function drawFace(self,axe_hdl,u_param,v_param,srf_option,ctrl_option)
            % draw surface on axes handle
            %
            if nargin < 6
                ctrl_option=[];
                if nargin < 5
                    srf_option=[];
                    if nargin < 4
                        v_param=[];
                        if nargin < 3
                            u_param=[];
                            if nargin < 2
                                axe_hdl=[];
                            end
                        end
                    end
                end
            end

            if isempty(axe_hdl),axe_hdl=axes(figure());end

            % default draw option
            if isempty(srf_option)
                srf_option=struct('LineStyle','none');
            end
            if isempty(ctrl_option)
                ctrl_option=struct('Marker','s','MarkerEdgeColor','r','EdgeColor','r','LineStyle','--','FaceAlpha',0);
            end

            % calculate point on surface
            Point=self.calFace(u_param,v_param);

            % plot surface
            dim=numel(self.calPoint(0,0));
            if dim == 2
                surface(axe_hdl,Point(:,:,1),Point(:,:,2),srf_option);
                if ~isempty(self.Poles)
                    surface(axe_hdl,self.Poles(:,:,1),self.Poles(:,:,2),ctrl_option);
                end
            else
                surface(axe_hdl,Point(:,:,1),Point(:,:,2),Point(:,:,3),srf_option);
                if ~isempty(self.Poles)
                    surface(axe_hdl,self.Poles(:,:,1),self.Poles(:,:,2),self.Poles(:,:,3),ctrl_option);
                end
                zlabel('z');
                view(3);
            end
            xlabel('x');
            ylabel('y');

            % axis equal
            % x_range=xlim();
            % y_range=ylim();
            % z_range=zlim();
            % center=[mean(x_range),mean(y_range),mean(z_range)];
            % range=max([x_range(2)-x_range(1),y_range(2)-y_range(1),z_range(2)-z_range(1)])/2;
            % xlim([center(1)-range,center(1)+range]);
            % ylim([center(2)-range,center(2)+range]);
            % zlim([center(3)-range,center(3)+range]);
        end

        function [step_str,obj_idx,ADVANCED_FACE]=getStep(self,obj_idx,param)
            % write BSpline into step file
            %
            if nargin < 2,obj_idx=1;end
            out_name=self.name;
            if isempty(out_name),out_name='NONE';end

            u_num=size(self.Poles,2);
            v_num=size(self.Poles,1);

            % generate CARTESIAN_POINT
            CARTESIAN_POINT=obj_idx;

            str_control=[];
            for u_idx=1:u_num
                for v_idx=1:v_num
                    str=[num2str(obj_idx,'#%d ='),' CARTESIAN_POINT ',...
                        '( ',...
                        '''NONE'', ',...
                        '( ',num2str(self.Poles(v_idx,u_idx,1),'%.16f'),', ',...
                        num2str(self.Poles(v_idx,u_idx,2),'%.16f'),', ',...
                        num2str(self.Poles(v_idx,u_idx,3),'%.16f'),' )',...
                        ' ) ;\n'];
                    str_control=[str_control,str];
                    obj_idx=obj_idx+1;
                end
            end

            % generate B_SPLINE_CURVE_WITH_KNOTS
            B_SPLINE_CURVE_WITH_KNOTS=obj_idx;
            str_curve=[];
            str_curve=[str_curve,stepCurve((1:v_num) +CARTESIAN_POINT-1,...
                self.VDegree,self.VMults,self.VKnots)];
            str_curve=[str_curve,stepCurve(((1:u_num)*v_num) +CARTESIAN_POINT-1,...
                self.UDegree,self.UMults,self.UKnots)];
            str_curve=[str_curve,stepCurve(((v_num*u_num):-1:(v_num*u_num-v_num+1)) +CARTESIAN_POINT-1,...
                self.VDegree,self.VMults,self.VKnots)];
            str_curve=[str_curve,stepCurve((((u_num-1):-1:0)*v_num+1) +CARTESIAN_POINT-1,...
                self.UDegree,self.UMults,self.UKnots)];

            % generate VERTEX_POINT
            VERTEX_POINT=obj_idx;
            str_vertex=[];
            str_vertex=[str_vertex,stepVertex(1 +CARTESIAN_POINT-1)];
            str_vertex=[str_vertex,stepVertex(v_num +CARTESIAN_POINT-1)];
            str_vertex=[str_vertex,stepVertex(v_num*u_num +CARTESIAN_POINT-1)];
            str_vertex=[str_vertex,stepVertex(v_num*u_num-v_num+1 +CARTESIAN_POINT-1)];

            % generate EDGE_CURVE
            EDGE_CURVE=obj_idx;
            str_edge=[];
            str_edge=[str_edge,stepEdge(VERTEX_POINT,VERTEX_POINT+1,B_SPLINE_CURVE_WITH_KNOTS)];
            str_edge=[str_edge,stepEdge(VERTEX_POINT+1,VERTEX_POINT+2,B_SPLINE_CURVE_WITH_KNOTS+1)];
            str_edge=[str_edge,stepEdge(VERTEX_POINT+2,VERTEX_POINT+3,B_SPLINE_CURVE_WITH_KNOTS+2)];
            str_edge=[str_edge,stepEdge(VERTEX_POINT+3,VERTEX_POINT,B_SPLINE_CURVE_WITH_KNOTS+3)];

            % generate ORIENTED_EDGE
            ORIENTED_EDGE=obj_idx;
            str_oriented=[];
            str_oriented=[str_oriented,stepOriented(EDGE_CURVE)];
            str_oriented=[str_oriented,stepOriented(EDGE_CURVE+1)];
            str_oriented=[str_oriented,stepOriented(EDGE_CURVE+2)];
            str_oriented=[str_oriented,stepOriented(EDGE_CURVE+3)];

            % generate EDGE_LOOP
            EDGE_LOOP=obj_idx;
            str_loop=[num2str(obj_idx,'#%d ='),' EDGE_LOOP',...
                ' ( ',...
                '''NONE''',', ',...
                '( ',num2str(ORIENTED_EDGE,'#%d'),', ',...
                num2str(ORIENTED_EDGE+1,'#%d'),', ',...
                num2str(ORIENTED_EDGE+2,'#%d'),', ',...
                num2str(ORIENTED_EDGE+3,'#%d'),') ',...
                ' ) ;\n'];
            obj_idx=obj_idx+1;

            % generate FACE_OUTER_BOUND
            FACE_OUTER_BOUND=obj_idx;
            str_outer=[num2str(obj_idx,'#%d ='),' FACE_OUTER_BOUND',...
                ' ( ',...
                '''NONE''',', ',...
                num2str(EDGE_LOOP,'#%d'),', ',...
                '.T.',...
                ' ) ;\n'];
            obj_idx=obj_idx+1;

            % generate B_SPLINE_SURFACE_WITH_KNOTS
            B_SPLINE_SURFACE_WITH_KNOTS=obj_idx;
            str_surf=[num2str(obj_idx,'#%d ='),' B_SPLINE_SURFACE_WITH_KNOTS',...
                ' ( ',...
                '''',out_name,'''',', ',...
                num2str(self.UDegree,'%d'),', ',num2str(self.VDegree,'%d'),', \n'];
            str_surf=[str_surf,'('];
            for u_bias=0:v_num:v_num*(u_num-2)
                str_surf=[str_surf,'( ',...
                    num2str((v_num+u_bias) +CARTESIAN_POINT-1,'#%d'),...
                    num2str((((v_num-1):-1:1)+u_bias) +CARTESIAN_POINT-1,', #%d'),...
                    ' ),\n'];
            end
            u_bias=v_num*(u_num-1);
            str_surf=[str_surf,'( ',...
                num2str((v_num+u_bias) +CARTESIAN_POINT-1,'#%d'),...
                num2str((((v_num-1):-1:1)+u_bias) +CARTESIAN_POINT-1,', #%d'),...
                ' )'];
            str_surf=[str_surf,'),\n'];
            str_surf=[str_surf,'.UNSPECIFIED.',', ','.F.',', ','.F.',', ','.F.',',\n',...
                '( ',num2str(self.UMults(1),'%d'),num2str(self.UMults(2:end),', %d'),' ),\n',...
                '( ',num2str(self.VMults(1),'%d'),num2str(self.VMults(2:end),', %d'),' ),\n',...
                '( ',num2str(self.UKnots(1),'%.16f'),num2str(self.UKnots(2:end),', %.16f'),' ),\n',...
                '( ',num2str(self.VKnots(1),'%.16f'),num2str(self.VKnots(2:end),', %.16f'),' ),\n',...
                '.UNSPECIFIED.',...
                ') ;\n'];
            obj_idx=obj_idx+1;

            % generate ADVANCED_FACE
            ADVANCED_FACE=obj_idx;
            str_face=[num2str(obj_idx,'#%d ='),' ADVANCED_FACE',...
                ' ( ',...
                '''',out_name,'''',', ',...
                '( ',num2str(FACE_OUTER_BOUND,'#%d'),' )',', '...
                num2str(B_SPLINE_SURFACE_WITH_KNOTS,'#%d'),', ',...
                '.T.',...
                ' ) ;\n'];
            obj_idx=obj_idx+1;

            % generate all
            step_str=[str_control,'\n',str_curve,'\n',str_vertex,'\n',...
                str_edge,'\n',str_oriented,'\n',str_loop,'\n',str_outer,'\n',...
                str_surf,'\n',str_face,'\n'];

            function str_edge=stepCurve(point_index,Degree,Mults,Knots)
                str_edge=[num2str(obj_idx,'#%d ='),' B_SPLINE_CURVE_WITH_KNOTS',...
                    ' ( ',...
                    '''NONE''',', ',...
                    num2str(Degree,'%d'),', \n',...
                    '( ',num2str(point_index(1),'#%d'),num2str(point_index(2:end),', #%d'),' ),\n',...
                    '.UNSPECIFIED., .F., .F.,\n',...
                    '( ',num2str(Mults(1),'%d'),num2str(Mults(2:end),', %d'),' ),\n',...
                    '( ',num2str(Knots(1),'%.16f'),num2str(Knots(2:end),', %.16f'),' ),\n',...
                    '.UNSPECIFIED.',...
                    ' ) ;\n'];
                obj_idx=obj_idx+1;
            end

            function str_vertex=stepVertex(point_index)
                str_vertex=[num2str(obj_idx,'#%d ='),' VERTEX_POINT',...
                    ' ( ',...
                    '''NONE''',', ',...
                    num2str(point_index,'#%d'),...
                    ' ) ;\n'];
                obj_idx=obj_idx+1;
            end

            function str_edge=stepEdge(start_vertex,end_vertex,edge_index)
                str_edge=[num2str(obj_idx,'#%d ='),' EDGE_CURVE',...
                    ' ( ',...
                    '''NONE''',', ',...
                    num2str(start_vertex,'#%d'),', ',...
                    num2str(end_vertex,'#%d'),', ',...
                    num2str(edge_index,'#%d'),', ',...
                    '.T.',...
                    ' ) ;\n'];
                obj_idx=obj_idx+1;
            end

            function str_oriented=stepOriented(edge_index)
                str_oriented=[num2str(obj_idx,'#%d ='),' ORIENTED_EDGE',...
                    ' ( ',...
                    '''NONE''',', ',...
                    '*',', ','*',', ',...
                    num2str(edge_index,'#%d'),', ',...
                    '.T.',...
                    ' ) ;\n'];
                obj_idx=obj_idx+1;
            end
        end
    end
end
