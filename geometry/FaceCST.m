classdef FaceCST < FaceNURBS
    % CST surface
    %
    % notice:
    % shape function in origin reference is local coordinate,
    % point calculate by class and shape will be translate to global
    %
    % reference: [1] Kulfan B. A Universal Parametric Geometry
    % Representation Method - "CST" [C]. 45th AIAA Aerospace Sciences
    % Meeting and Exhibit, Reno, Nevada, U.S.A.
    %
    properties
        % base parameter
        sym_x; % if true, class fcn will start from 0.5 to 1
        sym_y; % if true, class fcn will start from 0.5 to 1

        % origin parameter
        C_par_X;
        C_par_Y;
        C_par_ZV;
        C_par_ZU;

        class_fcn_X; % (V)
        class_fcn_Y; % (U)
        class_fcn_ZV; % (U,V), direction is V
        class_fcn_ZU; % (U,V), direction is U

        shape_fcn; % (U, V), default is LX, LY, LZ, output SX, SY, SZ
        class_fcn; % (U, V), default is CX, CY, CZV.*CZU, output CX, CY, CZ
    end

    properties
        % deform parameter
        deform_fcn_X=[]; % (U,V)
        deform_fcn_Y=[]; % (U,V)
        deform_fcn_Z=[]; % (U,V)

        % rotate parameter
        rotate_matrix=[];

        % translate parameter
        translate=[];
    end

    properties
        fit_data;
    end

    % define surface
    methods
        function self=FaceCST(name,C_par_X,C_par_Y,C_par_ZV,C_par_ZU,sym_x,sym_y,LX,LY,LZ)
            % generate 3D CST surface by LX, LY, LZ, C_par_X, C_par_Y, C_par_ZU, C_par_ZV
            %
            % u, x=LX*C(v)*S(v)
            % v, y=LY*C(u)*S(u)
            % w, z=LZ*C(u,v)*S(u,v)
            %
            % input:
            % name, LX, LY, LZ, C_par_X, C_par_Y, C_par_ZU, C_par_ZV, sym_y
            %
            % notice:
            % if input N1, N2 == 0, or C_par is empty,...
            % class_fcn will equal to 1
            % class_fcn_X(V)
            % class_fcn_Y(U)
            % class_fcn_ZU(U,V)(mainly control U direction)
            % class_fcn_ZV(U,V)(mainly control V direction)
            %
            self=self@FaceNURBS(name);
            if nargin < 10
                LZ=[];
                if nargin < 9
                    LY=[];
                    if nargin < 8
                        LX=[];
                        if nargin < 7
                            sym_y=[];
                            if nargin < 6
                                sym_x=[];
                                if nargin < 5
                                    C_par_ZU=[];
                                    if nargin < 4
                                        C_par_ZV=[];
                                        if nargin < 3
                                            C_par_Y=[];
                                            if nargin < 2
                                                C_par_X=[];
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end

            if isempty(sym_x),self.sym_x=false;
            else,self.sym_x=sym_x; end
            if isempty(sym_y),self.sym_y=false;
            else,self.sym_y=sym_y; end

            if isempty(C_par_X),C_par_X=[0,0];end
            if isempty(C_par_Y),C_par_Y=[0,0];end
            if isempty(C_par_ZV),C_par_ZV=[0,0];end
            if isempty(C_par_ZU),C_par_ZU=[0,0];end

            if isempty(LX),LX=1;end
            if isempty(LY),LY=1;end
            if isempty(LZ),LZ=1;end

            % class function
            if isnumeric(C_par_X)
                self.C_par_X=C_par_X;
                self.class_fcn_X=@(V) baseFcnClass(V,C_par_X(1),C_par_X(2));
            else
                self.C_par_X=[];
                self.class_fcn_X=C_par_X;
            end
            if isnumeric(C_par_Y)
                self.C_par_Y=C_par_Y;
                self.class_fcn_Y=@(U) baseFcnClass(U,C_par_Y(1),C_par_Y(2));
            else
                self.C_par_Y=[];
                self.class_fcn_Y=C_par_Y;
            end
            if isnumeric(C_par_ZV)
                self.C_par_ZV=C_par_ZV;
                self.class_fcn_ZV=@(U,V) baseFcnClass(V,C_par_ZV(1),C_par_ZV(2));
            else
                self.C_par_ZV=[];
                self.class_fcn_ZV=C_par_ZV;
            end
            if isnumeric(C_par_ZU)
                self.C_par_ZU=C_par_X;
                self.class_fcn_ZU=@(U,V) baseFcnClass(U,C_par_ZU(1),C_par_ZU(2));
            else
                self.C_par_ZU=[];
                self.class_fcn_ZU=C_par_ZU;
            end
            self.class_fcn=@(U,V) cat(3,self.class_fcn_X(V),self.class_fcn_Y(U),self.class_fcn_ZV(U,V).*self.class_fcn_ZU(U,V));
            self.shape_fcn=@(U,V) cat(3,U*LX,V*LY,ones(size(U))*LZ);
        end
    
        function Point=calCST(self,U,V)
            % calculate class
            % calculate class
            U_class=U;V_class=V;
            if self.sym_y,V_class=(V_class/2)+0.5;end
            if self.sym_x,U_class=(U_class/2)+0.5;end
            CPoint=self.class_fcn(U_class,V_class);

            % calculate shape
            SPoint=self.shape_fcn(U,V);

            Point=SPoint.*CPoint;
        end
    end

    % fit shape
    methods
        function addNURBS(self,Poles,UDegree,VDegree,...
                UMults,VMults,UKnots,VKnots,Weights)
            % add NURBS as shape function
            %
            if nargin < 9
                Weights=[];
                if nargin < 8
                    VKnots=[];
                    if nargin < 7
                        VMults = [];
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
            [v_pole_num,u_pole_num,dimension]=size(Poles);

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
                error('FaceCST.addNURBS: pole_num less than Degree+1');
            end

            if length(u_list) ~= u_pole_num+UDegree+1 || length(v_list) ~= v_pole_num+VDegree+1
                error('FaceCST.addNURBS: knot_num is not equal to pole_num+Degree+1');
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
            self.u_list=u_list;
            self.v_list=v_list;

            % use upper class BSpline calculate
            self.shape_fcn=@(U,V) self.calNURBS(U,V);
        end
    
        function fitNURBS(self,Nodes,UDegree,VDegree,u_pole_num,v_pole_num,U_node,V_node)
            % fit NURBS base on CST class function
            %
            if nargin < 8
                V_node=[];
                if nargin < 7
                    U_node=[];
                    if nargin < 6
                        v_pole_num=[];
                        if nargin < 5
                            u_pole_num=[];
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

            [v_node_num,u_node_num,dimension]=size(Nodes);
            if isempty(u_pole_num),u_pole_num=u_node_num;end
            if isempty(v_pole_num),v_pole_num=v_node_num;end
            if u_pole_num > u_node_num || v_pole_num > v_node_num
                error('FaceCST.fitNURBS: pole_num more than node_num')
            end

            % default value
            if isempty(U_node)
                U_node=vecnorm(Nodes(:,2:end,:)-Nodes(:,1:end-1,:),2,3);
                U_node=mean(U_node,1);
                U_node=[0;cumsum(U_node')];U_node=U_node/U_node(end);
            end
            U_node=U_node(:)';
            if isempty(V_node)
                V_node=vecnorm(Nodes(2:end,:,:)-Nodes(1:end-1,:,:),2,3);
                V_node=mean(V_node,2);
                V_node=[0;cumsum(V_node)];V_node=V_node/V_node(end);
            end
            V_node=V_node(:);

            UMults=[UDegree+1,ones(1,u_pole_num-UDegree-1),UDegree+1];
            UKnots=linspace(0,1,u_node_num-UDegree+1);
            for j=2:u_node_num-UDegree
                UKnots(j)=mean(U_node(j:j+UDegree-1));
            end
            % modify
            UKnots=interp1(linspace(0,1,u_node_num-UDegree+1),UKnots,linspace(0,1,u_pole_num-UDegree+1));
            u_list=baseKnotVec(UMults,UKnots);

            VMults=[VDegree+1,ones(1,v_pole_num-VDegree-1),VDegree+1];
            VKnots=linspace(0,1,v_node_num-VDegree+1);
            for j=2:v_node_num-VDegree
                VKnots(j)=mean(V_node(j:j+VDegree-1));
            end
            % modify
            VKnots=interp1(linspace(0,1,v_node_num-VDegree+1),VKnots,linspace(0,1,v_pole_num-VDegree+1));
            v_list=baseKnotVec(VMults,VKnots);

            % translate node to local coordinate
            [U_matrix,V_matrix]=meshgrid(U_node,V_node);
            Nodes=self.axisGlobalToLocal(Nodes,U_matrix,V_matrix);

            % base on node point list inverse calculate control point list
            matrix_v=zeros(v_node_num,v_pole_num);
            for node_idx=1:v_node_num
                v=V_node(node_idx);
                for ctrl_idx=1:v_pole_num
                    matrix_v(node_idx,ctrl_idx)=baseFcnN(v_list,v,ctrl_idx,VDegree);
                end
            end

            matrix_u=zeros(u_pole_num,u_node_num);
            for node_idx=1:u_node_num
                u=U_node(node_idx);
                for ctrl_idx=1:u_pole_num
                    matrix_u(ctrl_idx,node_idx)=baseFcnN(u_list,u,ctrl_idx,UDegree);
                end
            end

            % add class coefficient
            if self.sym_y,V_node=(V_node/2)+0.5;end
            if self.sym_x,U_node=(U_node/2)+0.5;end

            matrix_v_class_X=matrix_v.*self.class_fcn_X(V_node);
            matrix_v_class_Y=matrix_v;
            matrix_v_class_ZV=matrix_v.*self.class_fcn_ZV(0.5*ones(size(V_node)),V_node);

            matrix_u_class_X=matrix_u;
            matrix_u_class_Y=matrix_u.*self.class_fcn_Y(U_node);
            matrix_u_class_ZU=matrix_u.*self.class_fcn_ZU(U_node,0.5*ones(size(U_node)));

            % reverse calculate control point
            Poles(:,:,1)=matrix_v_class_X\Nodes(:,:,1)/matrix_u_class_X;
            Poles(:,:,2)=matrix_v_class_Y\Nodes(:,:,2)/matrix_u_class_Y;
            Poles(:,:,3)=matrix_v_class_ZV\Nodes(:,:,3)/matrix_u_class_ZU;
            Weights=ones(v_pole_num,u_pole_num);

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

            % fit data
            self.fit_data.Nodes=Nodes;
            self.fit_data.U_node=U_node;
            self.fit_data.V_node=V_node;
            self.fit_data.matrix_u=matrix_u;
            self.fit_data.matrix_v=matrix_v;

            % use upper class BSpline calculate
            self.shape_fcn=@(U,V) self.calNURBS(U,V);
        end
    end

    % deform, rotate, translate
    methods
        function addDeform(self,deform_fcn_X,deform_fcn_Y,deform_fcn_Z)
            % base on local coordinate deform surface
            %
            % input:
            % tran_fcn_X(V), tran_fcn_Y(U), tran_fcn_Z(U,V)
            %
            self.deform_fcn_X=deform_fcn_X;
            self.deform_fcn_Y=deform_fcn_Y;
            self.deform_fcn_Z=deform_fcn_Z;
        end

        function addRotate(self,ang_x,ang_y,ang_z)
            % base on angle to rotate surface
            % rotate order:
            % y, z, x
            %
            % input:
            % ang_x(deg), ang_y(deg), ang_z(deg)
            %
            matrix=eye(3);

            % process rotate
            if ang_y ~= 0
                cRY=cos(ang_y/180*pi);sRY=sin(ang_y/180*pi);
                matrix=[
                    cRY 0 sRY;
                    0 1 0
                    -sRY 0 cRY]*matrix;
            end

            if ang_z ~= 0
                cRZ=cos(ang_z/180*pi);sRZ=sin(ang_z/180*pi);
                matrix=[
                    cRZ -sRZ 0
                    sRZ cRZ 0
                    0 0 1]*matrix;
            end

            if ang_x ~= 0
                cRX=cos(ang_x/180*pi);sRX=sin(ang_x/180*pi);
                matrix=[
                    1 0 0;
                    0 cRX -sRX
                    0 sRX cRX]*matrix;
            end

            self.rotate_matrix=matrix;
        end

        function addTranslate(self,tran_x,tran_y,tran_z)
            % base on angle to rotate surface
            %
            self.translate=[tran_x,tran_y,tran_z];
        end

        function Point=axisLocalToGlobal(self,Point,U,V)
            % deform surface
            if ~isempty(self.deform_fcn_X)
                Point(:,:,1)=Point(:,:,1)+self.deform_fcn_X(U,V);
            end
            if ~isempty(self.deform_fcn_Y)
                Point(:,:,2)=Point(:,:,2)+self.deform_fcn_Y(U,V);
            end
            if ~isempty(self.deform_fcn_Z)
                Point(:,:,3)=Point(:,:,3)+self.deform_fcn_Z(U,V);
            end

            % rotate surface
            if ~isempty(self.rotate_matrix)
                matrix=self.rotate_matrix;
                X_old=Point(:,:,1);Y_old=Point(:,:,2);Z_old=Point(:,:,3);
                Point(:,:,1)=matrix(1,1)*X_old+matrix(1,2)*Y_old+matrix(1,3)*Z_old;
                Point(:,:,2)=matrix(2,1)*X_old+matrix(2,2)*Y_old+matrix(2,3)*Z_old;
                Point(:,:,3)=matrix(3,1)*X_old+matrix(3,2)*Y_old+matrix(3,3)*Z_old;
            end

            % translate surface
            if ~isempty(self.translate)
                Point(:,:,1)=Point(:,:,1)+self.translate(1);
                Point(:,:,2)=Point(:,:,2)+self.translate(2);
                Point(:,:,3)=Point(:,:,3)+self.translate(3);
            end
        end

        function Point=axisGlobalToLocal(self,Point,U,V)
            % re-translate surface
            if ~isempty(self.translate)
                Point(:,:,1)=Point(:,:,1)-self.translate(1);
                Point(:,:,2)=Point(:,:,2)-self.translate(2);
                Point(:,:,3)=Point(:,:,3)-self.translate(3);
            end

            % re-rotate surface
            if ~isempty(self.rotate_matrix)
                matrix=self.rotate_matrix';
                X_old=Point(:,:,1);Y_old=Point(:,:,2);Z_old=Point(:,:,3);
                Point(:,:,1)=matrix(1,1)*X_old+matrix(1,2)*Y_old+matrix(1,3)*Z_old;
                Point(:,:,2)=matrix(2,1)*X_old+matrix(2,2)*Y_old+matrix(2,3)*Z_old;
                Point(:,:,3)=matrix(3,1)*X_old+matrix(3,2)*Y_old+matrix(3,3)*Z_old;
            end

            % re-deform surface
            if ~isempty(self.deform_fcn_X)
                Point(:,:,1)=Point(:,:,1)-self.deform_fcn_X(U,V);
            end
            if ~isempty(self.deform_fcn_Y)
                Point(:,:,2)=Point(:,:,2)-self.deform_fcn_Y(U,V);
            end
            if ~isempty(self.deform_fcn_Z)
                Point(:,:,3)=Point(:,:,3)-self.deform_fcn_Z(U,V);
            end
        end

    end

    % calculate point
    methods
        function Point=calPoint(self,U,V)
            % calculate point on surface
            %
            % u, x=LX*C(v)*S(v)
            % v, y=LY*C(u)*S(u)
            % w, z=LZ*C(u,v)*S(u,v)
            %
            Point=calCST(self,U,V);
            Point=self.axisLocalToGlobal(Point,U,V);
        end
    end

    % calculate coord
    methods
        function [U,V,Point]=calCoord(self,Point)
            % base on X, Y, Z calculate local coordinate in surface
            %
            PointO=Point;geo_torl=100*eps;

            % undone translate surface
            if ~isempty(self.translate)
                Point(:,:,1)=Point(:,:,1)-self.translate(1);
                Point(:,:,2)=Point(:,:,2)-self.translate(2);
                Point(:,:,3)=Point(:,:,3)-self.translate(3);
            end

            % undone rotate surface
            if ~isempty(self.rotate_matrix)
                matrix=self.rotate_matrix';
                X_old=Point(:,:,1);Y_old=Point(:,:,2);Z_old=Point(:,:,3);
                Point(:,:,1)=matrix(1,1)*X_old+matrix(1,2)*Y_old+matrix(1,3)*Z_old;
                Point(:,:,2)=matrix(2,1)*X_old+matrix(2,2)*Y_old+matrix(2,3)*Z_old;
                Point(:,:,3)=matrix(3,1)*X_old+matrix(3,2)*Y_old+matrix(3,3)*Z_old;
            end

            Point_S=self.shape_fcn(0,0);
            Point_E=self.shape_fcn(1,1);

            V=Point(:,:,2)./abs(Point_E(:,:,2)-Point_S(:,:,2));
            U=Point(:,:,1)./abs(Point_E(:,:,1)-Point_S(:,:,1));

            U=max(U,0);U=min(U,1);
            V=max(V,0);V=min(V,1);

            % use project function to adjust parameter
            [Point,U,V]=self.calProject(PointO,U,V,geo_torl);
        end
    end

    % visualizate function
    methods
        function gplot(self,axe_hdl,u_param,v_param,srf_option,ctrl_option)
            % draw surface on axes handle
            %
            if nargin < 5
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
            surface(axe_hdl,Point(:,:,1),Point(:,:,2),Point(:,:,3),srf_option);
            if ~isempty(self.Poles)
                u_pole_num=size(self.Poles,2);v_pole_num=size(self.Poles,1);
                U_pole=interp1(linspace(0,1,u_pole_num-self.UDegree+1),self.u_list(self.UDegree+1:u_pole_num+1),linspace(0,1,u_pole_num));
                V_pole=interp1(linspace(0,1,v_pole_num-self.VDegree+1),self.v_list(self.VDegree+1:v_pole_num+1),linspace(0,1,v_pole_num));
                [U_pole,V_pole]=meshgrid(U_pole,V_pole);
                if self.sym_y,V_pole=(V_pole/2)+0.5;end
                if self.sym_x,U_pole=(U_pole/2)+0.5;end
                Poles=self.Poles.*self.class_fcn(U_pole,V_pole);
                Poles=self.axisLocalToGlobal(Poles,U_pole,V_pole);
                surface(axe_hdl,Poles(:,:,1),Poles(:,:,2),Poles(:,:,3),ctrl_option);
            end
            view(3);
            xlabel('x');
            ylabel('y');
            zlabel('z');

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

        function fce=getNURBS(self,u_param,v_param,u_pole_num,v_pole_num)
            % convert CST Face into NURBS Face
            %
            % input:
            % u_param, v_param
            % or:
            % value_torl, []
            %
            if nargin < 5
                v_pole_num=[];
                if nargin < 4
                    u_pole_num=[];
                    if nargin < 3
                        v_param=[];
                        if nargin < 2
                            u_param=[];
                        end
                    end
                end
            end

            if ~isempty(u_param) && length(u_param) > 1 && u_param(1,1) > u_param(1,2)
                u_param=fliplr(u_param);
            end
            if ~isempty(v_param) && length(v_param) > 1 && v_param(1,1) > v_param(2,1)
                v_param=flipud(v_param);
            end
            [Nodes,U_node,V_node]=calFace(self,u_param,v_param);

            UDegree=min([size(U_node,2)-1,3,u_pole_num-1]);VDegree=min([size(V_node,1)-1,3,v_pole_num-1]);
            fce=GeomApp.VertexToFace(self.name,Nodes,UDegree,VDegree,u_pole_num,v_pole_num,U_node(1,:),V_node(:,1));
        end

        function [step_str,obj_idx,ADVANCED_FACE]=getStep(self,obj_idx,param)
            % interface of BSpline surface getStep function
            %
            fce=self.getNURBS(param);
            [step_str,obj_idx,ADVANCED_FACE]=fce.getStep(obj_idx,param);
        end
    end
end
