classdef SurfaceCST < SurfaceBSpline
    % CST surface
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

        shape_fcn; % (U, V), default is LX, LY, LZ
    end

    properties
        % deform parameter
        deform_fcn_X=[]; % (V)
        deform_fcn_Y=[]; % (U)
        deform_fcn_Z=[]; % (U,V)

        % rotation parameter
        rotate_matrix=[];

        % translation parameter
        translation=[];
    end

    properties
        fit_data;
    end

    % define surface
    methods
        function self=SurfaceCST(name,C_par_X,C_par_Y,C_par_ZV,C_par_ZU,sym_x,sym_y,LX,LY,LZ)
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
            self=self@SurfaceBSpline(name);
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
                self.class_fcn_X=@(V) defcnClass(V,C_par_X(1),C_par_X(2));
            else
                self.C_par_X=[];
                self.class_fcn_X=C_par_X;
            end
            if isnumeric(C_par_Y)
                self.C_par_Y=C_par_Y;
                self.class_fcn_Y=@(U) defcnClass(U,C_par_Y(1),C_par_Y(2));
            else
                self.C_par_Y=[];
                self.class_fcn_Y=C_par_Y;
            end
            if isnumeric(C_par_ZV)
                self.C_par_X=C_par_ZV;
                self.class_fcn_ZV=@(U,V) defcnClass(V,C_par_ZV(1),C_par_ZV(2));
            else
                self.C_par_ZV=[];
                self.class_fcn_ZV=C_par_ZV;
            end
            if isnumeric(C_par_ZU)
                self.C_par_ZU=C_par_X;
                self.class_fcn_ZU=@(U,V) defcnClass(U,C_par_ZU(1),C_par_ZU(2));
            else
                self.C_par_ZU=[];
                self.class_fcn_ZU=C_par_ZU;
            end

            self.shape_fcn=@(U,V) defcnShape(U,V,LX,LY,LZ);

            function [X,Y,Z]=defcnShape(U,V,LX,LY,LZ)
                X=U*LX;Y=V*LY;Z=ones(size(U))*LZ;
            end
        end
    end

    % fit shape
    methods
        function addShapeBSpline(self,point_X,point_Y,point_Z,FLAG_FIT,...
                u_degree,v_degree,u_knot_multi,v_knot_multi,u_knot_list,v_knot_list,...
                u_ctrl_num,v_ctrl_num,U,V)
            % fit node point to generate shape function
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

            if isempty(FLAG_FIT),FLAG_FIT=false;end

            % check input value size and giving default value
            if FLAG_FIT
                node_X=point_X;node_Y=point_Y;node_Z=point_Z;
                [v_node_num,u_node_num]=size(node_X);
                if any(size(node_Y) ~= [v_node_num,u_node_num]) ||...
                        any(size(node_Z) ~= [v_node_num,u_node_num])
                    error('SurfaceCST.addShapeBSpline: size of node_X, node_Y, node_Z do not equal');
                end

                if isempty(u_ctrl_num),u_ctrl_num=u_node_num;end
                if u_ctrl_num > u_node_num
                    error('SurfaceCST.addShapeBSpline: u_control number more than u_node number')
                end
                if isempty(v_ctrl_num),v_ctrl_num=v_node_num;end
                if v_ctrl_num > v_node_num
                    error('SurfaceCST.addShapeBSpline: v_control number more than v_node number')
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
                    error('SurfaceCST.addShapeBSpline: size of ctrl_X, ctrl_Y, ctrl_Z do not equal');
                end
            end

            if u_ctrl_num < (u_degree+1) || v_ctrl_num < (v_degree+1)
                error('SurfaceCST.addShapeBSpline: ctrl_num less than degree+1');
            end

            % default value
            if isempty(u_degree),u_degree=u_ctrl_num-1;end
            if isempty(U),U=linspace(0,1,u_node_num);end;U=U(:);
            if isempty(u_knot_multi),u_knot_multi=[u_degree+1,ones(1,u_ctrl_num-u_degree-1),u_degree+1];end
            if isempty(u_knot_list),u_knot_list=interp1(linspace(0,1,u_node_num),U,linspace(0,1,u_ctrl_num-u_degree+1));end
            u_list=getKnotVec(u_knot_multi,u_knot_list);

            % default value
            if isempty(v_degree),v_degree=v_ctrl_num-1;end
            if isempty(V),V=linspace(0,1,v_node_num);end;V=V(:);
            if isempty(v_knot_multi),v_knot_multi=[v_degree+1,ones(1,v_ctrl_num-v_degree-1),v_degree+1];end
            if isempty(v_knot_list),v_knot_list=interp1(linspace(0,1,v_node_num),V,linspace(0,1,v_ctrl_num-v_degree+1));end
            v_list=getKnotVec(v_knot_multi,v_knot_list);

            if FLAG_FIT
                % process symmetry
                U_class=U;V_class=V;
                if self.sym_y,V_class=(V_class/2)+0.5;end
                if self.sym_x,U_class=(U_class/2)+0.5;end
    
                % base on node point list inverse calculate control point list
                matrix_v=zeros(v_node_num,v_ctrl_num);
                for node_idx=1:v_node_num
                    v=V(node_idx);
                    for ctrl_idx=1:v_ctrl_num
                        matrix_v(node_idx,ctrl_idx)=baseFcnN(v_list,v,ctrl_idx,v_degree);
                    end
                end
    
                matrix_v_class_X=matrix_v.*self.class_fcn_X(V_class);
                matrix_v_class_Y=matrix_v;
                matrix_v_class_ZV=matrix_v.*self.class_fcn_ZV(U_class,V_class);
    
                matrix_u=zeros(u_ctrl_num,u_node_num);
                for node_idx=1:u_node_num
                    u=U(node_idx);
                    for ctrl_idx=1:u_ctrl_num
                        matrix_u(ctrl_idx,node_idx)=baseFcnN(u_list,u,ctrl_idx,u_degree);
                    end
                end
    
                matrix_u_class_X=matrix_u;
                matrix_u_class_Y=matrix_u.*self.class_fcn_Y(U_class');
                matrix_u_class_ZU=matrix_u.*self.class_fcn_ZU(U_class',V_class');
    
                % undone translation surface
                if ~isempty(self.translation)
                    node_X=node_X-self.translation(1);
                    node_Y=node_Y-self.translation(2);
                    node_Z=node_Z-self.translation(3);
                end
    
                % undone rotation surface
                if ~isempty(self.rotate_matrix)
                    matrix=self.rotate_matrix';
                    X_old=node_X;Y_old=node_Y;Z_old=Z;
                    node_X=matrix(1,1)*X_old+matrix(1,2)*Y_old+matrix(1,3)*Z_old;
                    node_Y=matrix(2,1)*X_old+matrix(2,2)*Y_old+matrix(2,3)*Z_old;
                    node_Z=matrix(3,1)*X_old+matrix(3,2)*Y_old+matrix(3,3)*Z_old;
                end
    
                % undone deform surface
                [U_Mat,V_Mat]=meshgrid(U,V);
                if ~isempty(self.deform_fcn_Y),node_Y=node_Y-self.deform_fcn_Y(U_Mat);end
                if ~isempty(self.deform_fcn_X),node_X=node_X-self.deform_fcn_X(V_Mat);end
                if ~isempty(self.deform_fcn_Z),node_Z=node_Z-self.deform_fcn_Z(U_Mat,V_Mat);end
 
                ctrl_X=matrix_v_class_X\node_X/matrix_u_class_X;
                ctrl_Y=matrix_v_class_Y\node_Y/matrix_u_class_Y;
                ctrl_Z=matrix_v_class_ZV\node_Z/matrix_u_class_ZU;

                self.fit_data.node_X=node_X;
                self.fit_data.node_Y=node_Y;
                self.fit_data.node_Z=node_Z;
                self.fit_data.U=U;
                self.fit_data.V=V;
                self.fit_data.matrix_u=matrix_u;
                self.fit_data.matrix_v=matrix_v;
            elseif ~isempty(ctrl_X) && ~isempty(ctrl_Y) && ~isempty(ctrl_Z)
                % generate B spline surface by control point
                self.fit_data=[];
            else
                error('SurfaceCST.addShapeBSpline: error input, lack control point or node point');
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

            % shape function
            self.u_list=u_list;
            self.v_list=v_list;
            self.shape_fcn=@(U,V) self.calPointShape(U,V);
        end

        function [X,Y,Z]=calPointShape(self,U,V)
            % according u_x to calculate point
            % u_x_list is u_x_number x 1 matrix
            % point_list is point_number x dimension matrix
            %
            [rank_num,colume_num]=size(U);
            if any(size(V) ~= [rank_num,colume_num])
                error('SurfaceBSpline.calPoint: size of U_x do not equal to size of V_x');
            end

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
            % base on angle to rotation surface
            % rotation order:
            % y, z, x
            %
            % input:
            % ang_x(deg), ang_y(deg), ang_z(deg)
            %
            matrix=eye(3);

            % process rotation
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
            % base on angle to rotation surface
            %
            self.translation=[tran_x,tran_y,tran_z];
        end

    end

    % calculate point
    methods
        function [X,Y,Z]=calPoint(self,U,V)
            % calculate point on surface
            %
            % u, x=LX*C(v)*S(v)
            % v, y=LY*C(u)*S(u)
            % w, z=LZ*C(u,v)*S(u,v)
            %

            % calculate origin surface
            [X,Y,Z]=self.shape_fcn(U,V);

            % calculate class
            U_class=U;V_class=V;
            if self.sym_y,V_class=(V_class/2)+0.5;end
            if self.sym_x,U_class=(U_class/2)+0.5;end
            X=X.*self.class_fcn_X(V_class);
            Y=Y.*self.class_fcn_Y(U_class);
            Z=Z.*self.class_fcn_ZV(U_class,V_class).*self.class_fcn_ZU(U_class,V_class);

            % deform surface
            if ~isempty(self.deform_fcn_X)
                X=X+self.deform_fcn_X(V);
            end
            if ~isempty(self.deform_fcn_Y)
                Y=Y+self.deform_fcn_Y(U);
            end
            if ~isempty(self.deform_fcn_Z)
                Z=Z+self.deform_fcn_Z(U,V);
            end

            % rotation surface
            if ~isempty(self.rotate_matrix)
                matrix=self.rotate_matrix;
                X_old=X;Y_old=Y;Z_old=Z;
                X=matrix(1,1)*X_old+matrix(1,2)*Y_old+matrix(1,3)*Z_old;
                Y=matrix(2,1)*X_old+matrix(2,2)*Y_old+matrix(2,3)*Z_old;
                Z=matrix(3,1)*X_old+matrix(3,2)*Y_old+matrix(3,3)*Z_old;
            end

            % translation surface
            if ~isempty(self.translation)
                X=X+self.translation(1);
                Y=Y+self.translation(2);
                Z=Z+self.translation(3);
            end
        end
    end

    % calculate coord
    methods
        function [U,V,X,Y,Z]=calCoord(self,X,Y,Z)
            % base on X, Y, Z calculate local coordinate in surface
            %
            XO=X;YO=Y;ZO=Z;geo_torl=100*eps;

            % undone translation surface
            if ~isempty(self.translation)
                X=X-self.translation(1);
                Y=Y-self.translation(2);
                Z=Z-self.translation(3);
            end

            % undone rotation surface
            if ~isempty(self.rotate_matrix)
                matrix=self.rotate_matrix';
                X_old=X;Y_old=Y;Z_old=Z;
                X=matrix(1,1)*X_old+matrix(1,2)*Y_old+matrix(1,3)*Z_old;
                Y=matrix(2,1)*X_old+matrix(2,2)*Y_old+matrix(2,3)*Z_old;
                Z=matrix(3,1)*X_old+matrix(3,2)*Y_old+matrix(3,3)*Z_old;
            end

            [X_S,Y_S]=self.shape_fcn(0,0);
            [X_E,Y_E]=self.shape_fcn(1,1);

            V=Y./abs(Y_E-Y_S);
            U=X./abs(X_E-X_S);

            U=max(U,0);U=min(U,1);
            V=max(V,0);V=min(V,1);

            % % re deform surface
            % if ~isempty(self.deform_fcn_Y),Y=Y-self.deform_fcn_Y(U);end
            % if ~isempty(self.deform_fcn_X),X=X-self.deform_fcn_X(V);end

            % use project function to adjust parameter
            [X,Y,Z,U,V]=self.calProject(XO,YO,ZO,U,V,geo_torl);
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
            if ~isempty(self.node_X) && ~isempty(self.node_Y) && ~isempty(self.node_Z)
                [U,V]=meshgrid(self.u_list(self.u_degree+1:self.u_ctrl_num+1),self.v_list(self.v_degree+1:self.v_ctrl_num+1));
                [X,Y,Z]=calPointPure(self.node_X,self.node_Y,self.node_Z,U,V);
                surface(axe_hdl,X,Y,Z,node_option);
            end
            if ~isempty(self.ctrl_X) && ~isempty(self.ctrl_Y) && ~isempty(self.ctrl_Z)
                u_list=interp1(linspace(0,1,self.u_ctrl_num-self.u_degree+1),self.u_list(self.u_degree+1:self.u_ctrl_num+1),linspace(0,1,self.u_ctrl_num));
                v_list=interp1(linspace(0,1,self.v_ctrl_num-self.v_degree+1),self.v_list(self.v_degree+1:self.v_ctrl_num+1),linspace(0,1,self.v_ctrl_num));
                [U,V]=meshgrid(u_list,v_list);
                [X,Y,Z]=calPointPure(self.ctrl_X,self.ctrl_Y,self.ctrl_Z,U,V);
                surface(axe_hdl,X,Y,Z,control_option);
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

            function [X,Y,Z]=calPointPure(PX,PY,PZ,U,V)
                % calculate class
                U_class=U;V_class=V;
                if self.sym_y,V_class=(V_class/2)+0.5;end
                if self.sym_x,U_class=(U_class/2)+0.5;end
                
                X=PX.*self.class_fcn_X(V_class);
                Y=PY.*self.class_fcn_Y(U_class);
                Z=PZ.*self.class_fcn_ZV(U_class,V_class).*self.class_fcn_ZU(U_class,V_class);

                % deform surface
                if ~isempty(self.deform_fcn_X)
                    X=X+self.deform_fcn_X(V);
                end
                if ~isempty(self.deform_fcn_Y)
                    Y=Y+self.deform_fcn_Y(U);
                end
                if ~isempty(self.deform_fcn_Z)
                    Z=Z+self.deform_fcn_Z(U,V);
                end

                % rotation surface
                if ~isempty(self.rotate_matrix)
                    matrix=self.rotate_matrix;
                    X_old=X;Y_old=Y;Z_old=Z;
                    X=matrix(1,1)*X_old+matrix(1,2)*Y_old+matrix(1,3)*Z_old;
                    Y=matrix(2,1)*X_old+matrix(2,2)*Y_old+matrix(2,3)*Z_old;
                    Z=matrix(3,1)*X_old+matrix(3,2)*Y_old+matrix(3,3)*Z_old;
                end

                % translation surface
                if ~isempty(self.translation)
                    X=X+self.translation(1);
                    Y=Y+self.translation(2);
                    Z=Z+self.translation(3);
                end
            end
        end

        function surface_BSpline=getSurfaceBSpline(self,u_param,v_param)
            % convert CST surface into BSpline surface
            %
            % input:
            % U,V
            % u_grid_number, v_grid_number
            % value torlance
            %
            if nargin < 3
                v_param=[];
                if nargin < 2
                    u_param=[];
                end
            end

            if ~isempty(u_param) && length(u_param) > 1 && u_param(1,1) > u_param(1,2)
                u_param=fliplr(u_param);
            end
            if ~isempty(v_param) && length(v_param) > 1 && v_param(1,1) > v_param(2,1)
                v_param=flipud(v_param);
            end
            [X,Y,Z,u_param,v_param]=calSurface(self,u_param,v_param);

            u_degree=min(size(u_param,2)-1,3);v_degree=min(size(v_param,1)-1,3);
            surface_BSpline=SurfaceBSpline(self.name,X,Y,Z,true,u_degree,v_degree);
        end

        function [step_str,object_index,ADVANCED_FACE]=getStep(self,object_index)
            % interface of BSpline surface getStep function
            %
            surf_BSpline=self.getSurfaceBSpline(1e-2);
            [step_str,object_index,ADVANCED_FACE]=surf_BSpline.getStep(object_index);
        end
    end
end

%% common function

function C=defcnClass(U,N1,N2)
% default class function of X, Y
%
NP=calNormPar(N1,N2);
C=U.^N1.*(1-U).^N2/NP;
end

function nomlz_par=calNormPar(N1,N2)
% calculate normailize class function parameter by N1, N2
%
nomlz_par=(N1./(N1+N2)).^N1.*(N2./(N1+N2)).^N2;
nomlz_par((N1 == 0) & (N2 == 0))=1;
end

function knot_vec=getKnotVec(knot_multi,knot_list)
% base on knot_multi and knot_list to create knot vector
%
knot_vec=zeros(1,sum(knot_multi));
start_idx=1;end_idx=knot_multi(1);
for n_idx=1:length(knot_list)
    knot_vec(start_idx:end_idx)=knot_list(n_idx);
    start_idx=start_idx+knot_multi(n_idx);
    end_idx=end_idx+knot_multi(n_idx);
end
end

function N=baseFcnN(u_list,u_x,i,k)
% base function of BSpline curve
%
if k == 0
    if ((u_list(i) <= u_x) && (u_x <= u_list(i+1)))
        if any(u_list == u_x) && u_x ~= u_list(1) && u_x ~= u_list(end)
            N=0.5;
        else
            N=1;
        end
    else
        N=0;
    end
else
    if u_list(i+k) == u_list(i)
        A=0;
    else
        A=(u_x-u_list(i))/(u_list(i+k)-u_list(i));
    end

    if u_list(i+k+1) == u_list(i+1)
        B=0;
    else
        B=(u_list(i+k+1)-u_x)/(u_list(i+k+1)-u_list(i+1));
    end

    N=A*baseFcnN(u_list,u_x,i,k-1)+B*baseFcnN(u_list,u_x,i+1,k-1);
end
end
