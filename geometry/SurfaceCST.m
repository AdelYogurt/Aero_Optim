classdef SurfaceCST < handle
    % Class-Shape Transformation (CST) Surface
    %
    % notice:
    % origin reference of shape function is local coordinate
    % point calculate by class and shape will be translate to global
    %
    % reference: [1] Kulfan B. A Universal Parametric Geometry
    % Representation Method - "CST" [C]. 45th AIAA Aerospace Sciences
    % Meeting and Exhibit, Reno, Nevada, U.S.A.
    %
    properties % CST parameter
        sym_x; % if true, U_class fcn will start from 0.5 to 1
        sym_y; % if true, V_class fcn will start from 0.5 to 1

        C_par_X; % origin parameter, default is [0.0, 0.0]
        C_par_Y; % origin parameter, default is [0.0, 0.0]
        C_par_ZV; % origin parameter, default is [0.0, 0.0]
        C_par_ZU; % origin parameter, default is [0.0, 0.0]

        spline; % 3D BSpline surface

        class_fcn_X; % (V)
        class_fcn_Y; % (U)
        class_fcn_ZV; % (U, V), default is baseFcnClass(V, C_par_ZV(1), C_par_ZV(2))
        class_fcn_ZU; % (U, V), default is baseFcnClass(U, C_par_ZU(1), C_par_ZU(2))

        class_fcn; % (U, V), output are CX, CY, CZ(default is CZV.*CZU)
        shape_fcn; % (U, V), default are LX, LY, LZ, output are SX, SY, SZ
    end

    properties % basical surface parameter
        deform_fcn_X=[]; % deform parameter: deform_fcn_X(U,V)
        deform_fcn_Y=[]; % deform parameter: deform_fcn_Y(U,V)
        deform_fcn_Z=[]; % deform parameter: deform_fcn_Z(U,V)
        rotate_matrix=[]; % rotate parameter
        translate=[]; % translate parameter
    end

    properties % fit point set parameter
        fit_data;
    end

    methods % define surface
        function self=SurfaceCST(C_par_X,C_par_Y,C_par_ZV,C_par_ZU,sym_x,sym_y,LX,LY,LZ)
            % generate 3D CST surface by LX, LY, LZ, C_par_X, C_par_Y, C_par_ZU, C_par_ZV
            % class_fcn(u,v)=[class_fcn_X(v),class_fcn_Y(u),class_fcn_ZV(u,v).*class_fcn_ZU(u,v)]
            % X=LX*class_fcn(u,v)*shape_fcn(u,v)
            % Y=LY*class_fcn(u,v)*shape_fcn(u,v)
            %
            % input:
            % LX, LY, LZ, C_par_X, C_par_Y, C_par_ZU, C_par_ZV, sym_y
            %
            % notice:
            % if C_par is [0.0, 0.0] or empty, class_fcn will equal to 1.0
            % class_fcn_X(V)
            % class_fcn_Y(U)
            % class_fcn_ZV(U,V)(mainly control V direction)
            % class_fcn_ZU(U,V)(mainly control U direction)
            %
            if nargin < 9
                LZ=[];
                if nargin < 8
                    LY=[];
                    if nargin < 7
                        LX=[];
                        if nargin < 6
                            sym_y=[];
                            if nargin < 5
                                sym_x=[];
                                if nargin < 4
                                    C_par_ZU=[];
                                    if nargin < 3
                                        C_par_ZV=[];
                                        if nargin < 2
                                            C_par_Y=[];
                                            if nargin < 1
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

            % shape function
            self.spline=Surface(cat(3,[0,LX;0,LX],[0,0;LY,LY],[LZ,LZ;LZ,LZ]));
            self.shape_fcn=@(U,V) self.spline.calPoint(U,V);
        end

        function Points=calPoint(self,U_x,V_x)
            % calculate point on curve
            %
            U_class=U_x;V_class=V_x;
            if self.sym_y,V_class=(V_class/2)+0.5;end
            if self.sym_x,U_class=(U_class/2)+0.5;end

            % calculate class
            CPoint=self.class_fcn(U_class,V_class);

            % calculate shape
            SPoint=self.shape_fcn(U_x,V_x);

            Points=SPoint.*CPoint;
            Points=self.axisLocalToGlobal(Points,U_x,V_x);
        end
        
        function [Points,dPoints_dU,dPoints_dV]=calGradient(self,U_x,V_x,step)
            % use differ to calculate gradient
            %
            if nargin < 4 || isempty(step),step=100*eps;end

            dim=numel(self.calPoint(0,0));
            [v_num,u_num]=size(U_x);U_x=U_x(:);V_x=V_x(:);

            Points=self.calPoint(U_x,V_x);Points=reshape(Points,[],dim);

            Point_UF=self.calPoint(min(U_x+step,1),V_x);Point_UF=reshape(Point_UF,[],dim);
            Point_UB=self.calPoint(max(U_x-step,0),V_x);Point_UB=reshape(Point_UB,[],dim);
            Bool_U=(U_x+step) >= 1;
            Bool_D=(U_x-step) <= 0;
            Bool(:,1)=Bool_U;Bool(:,2)=Bool_D;
            Bool_C=~any(Bool,2);
            dPoints_dU=zeros(v_num*u_num,dim); % allocate memory
            dPoints_dU(Bool_C,:)=(Point_UF(Bool_C,:)-Point_UB(Bool_C,:))/2/step;
            dPoints_dU(Bool_U,:)=(Points(Bool_U,:)-Point_UB(Bool_U,:))/step;
            dPoints_dU(Bool_D,:)=(Point_UF(Bool_D,:)-Points(Bool_D,:))/step;
            dPoints_dU=real(dPoints_dU);
            dPoints_dU=reshape(dPoints_dU,v_num,u_num,dim);

            Point_VF=self.calPoint(U_x,min(V_x+step,1));Point_VF=reshape(Point_VF,[],dim);
            Point_VB=self.calPoint(U_x,max(V_x-step,0));Point_VB=reshape(Point_VB,[],dim);
            Bool_U=(V_x+step) >= 1;
            Bool_D=(V_x-step) <= 0;
            Bool(:,1)=Bool_U;Bool(:,2)=Bool_D;
            Bool_C=~any(Bool,2);
            dPoints_dV=zeros(v_num*u_num,dim); % allocate memory
            dPoints_dV(Bool_C,:)=(Point_VF(Bool_C,:)-Point_VB(Bool_C,:))/2/step;
            dPoints_dV(Bool_U,:)=(Points(Bool_U,:)-Point_VB(Bool_U,:))/step;
            dPoints_dV(Bool_D,:)=(Point_VF(Bool_D,:)-Points(Bool_D,:))/step;
            dPoints_dV=real(dPoints_dV);
            dPoints_dV=reshape(dPoints_dV,v_num,u_num,dim);

            Points=reshape(Points,v_num,u_num,dim);
        end

        function srf=convertSpline(self)
            % convert CST Surface into Spline Surface
            %
            warning('SurfaceCST.convertSpline: can not convert to spline surface exactly');
            [U_node,V_node]=meshgrid(linspace(0,1,21),linspace(0,1,21));
            Nodes=self.calPoint(U_node,V_node);UDegree=3;VDegree=3;
            srf=GeomApp.interpPointToSurface(Nodes,UDegree,VDegree);
        end

        function srf_hdl=displayGeom(self,axe_hdl,srf_option,u_param,v_param)
            % draw surface on axes handle
            %
            if nargin < 5
                v_param=[];
                if nargin < 4
                    u_param=[];
                    if nargin < 3
                        srf_option=[];
                        if nargin < 2
                            axe_hdl=[];
                        end
                    end
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end

            % default draw option
            if isempty(srf_option),srf_option=struct('LineStyle','none');end

            % calculate point on surface
            if isempty(u_param), u_param=101;end
            if length(u_param) == 1, u_param=linspace(0,1,u_param);end
            if isempty(v_param), v_param=101;end
            if length(v_param) == 1, v_param=linspace(0,1,v_param);end
            [u_param,v_param]=meshgrid(u_param,v_param);
            Points=self.calPoint(u_param,v_param);

            % plot surface
            srf_hdl=surface(axe_hdl,Points(:,:,1),Points(:,:,2),Points(:,:,3),srf_option);
            view(3);
            xlabel('x');
            ylabel('y');
            zlabel('z');
        end

        function srf_hdl=displayPoles(self,axe_hdl,pole_option)
            % draw surface on axes handle
            %
            if nargin < 3
                pole_option=[];
                if nargin < 2
                    axe_hdl=[];
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end

            % default draw option
            if isempty(pole_option)
                pole_option=struct('Marker','s','MarkerEdgeColor','r','EdgeColor','r','LineStyle','--','FaceAlpha',0);
            end

            % plot surface
            if ~isempty(self.spline)
                % load properties
                Poles=self.spline.Poles;
                [v_pole_num,u_pole_num,~]=size(Poles);
                UDegree=self.spline.UDegree;
                VDegree=self.spline.VDegree;
                u_list=self.spline.u_list;
                v_list=self.spline.v_list;

                % calculate U_pole, V_pole
                U_pole=interp1(linspace(0,1,u_pole_num-UDegree+1),u_list(UDegree+1:u_pole_num+1),linspace(0,1,u_pole_num));
                V_pole=interp1(linspace(0,1,v_pole_num-VDegree+1),v_list(VDegree+1:v_pole_num+1),linspace(0,1,v_pole_num));
                [U_pole,V_pole]=meshgrid(U_pole,V_pole);
                if self.sym_y,V_pole=(V_pole/2)+0.5;end
                if self.sym_x,U_pole=(U_pole/2)+0.5;end
                Poles=self.axisLocalToGlobal(Poles.*self.class_fcn(U_pole,V_pole),U_pole,V_pole);

                srf_hdl=surface(axe_hdl,Poles(:,:,1),Poles(:,:,2),Poles(:,:,3),pole_option);
            end
            view(3);
            xlabel('x');
            ylabel('y');
            zlabel('z');
        end

        function crv=getBoundCurve(self,idx)
            switch idx
                case 1
                    crv=Curve(reshape(self.Poles(1,:,:),size(self.Poles,2),[]),self.UDegree,self.UMults,self.UKnots,self.Weights(1,:));
                case 2
                    crv=Curve(reshape(self.Poles(end,:,:),size(self.Poles,2),[]),self.UDegree,self.UMults,self.UKnots,self.Weights(end,:));
                case 3
                    crv=Curve(reshape(self.Poles(:,1,:),size(self.Poles,1),[]),self.VDegree,self.VMults,self.VKnots,self.Weights(:,1));
                case 4
                    crv=Curve(reshape(self.Poles(:,end,:),size(self.Poles,1),[]),self.VDegree,self.VMults,self.VKnots,self.Weights(:,end));
            end
        end
    end

    methods % add spline shape function
        function addSpline(self,spline)
            % add Spline as shape function
            %
            if size(spline.Poles,3) > 3
                warning('SurfaceCST.addSpline: dimension of spline is large than 3, spline will be project to 3D')
                spline.Poles=spline.Poles(:,:,1:3);
            elseif size(spline.Poles,3) == 2
                [v_pole_num,u_pole_num,~]=size(spline.Poles);
                spline.Poles=cat(3,repmat(linspace(0,1,u_pole_num),v_pole_num,1),spline.Poles);
            elseif size(spline.Poles,3) == 1
                [v_pole_num,u_pole_num,~]=size(spline.Poles);
                [U_pole,V_pole]=meshgrid(linspace(0,1,u_pole_num),linspace(0,1,v_pole_num));
                spline.Poles=cat(3,U_pole,V_pole,spline.Poles);
            end
            if  min(spline.UKnots) ~= 0.0 || max(spline.UKnots) ~= 1.0
                warning('SurfaceCST.addSpline: UKnots of spline will be rescale to [0.0, 1.0]')
                spline.UKnots=(spline.UKnots-min(spline.UKnots))/(max(spline.UKnots)-min(spline.UKnots));
            end
            if  min(spline.VKnots) ~= 0.0 || max(spline.VKnots) ~= 1.0
                warning('SurfaceCST.addSpline: VKnots of spline will be rescale to [0.0, 1.0]')
                spline.VKnots=(spline.VKnots-min(spline.VKnots))/(max(spline.VKnots)-min(spline.VKnots));
            end
            self.spline=spline;

            % add shape function
            self.shape_fcn=@(U_x) self.spline.calPoint(U_x);
        end

        function fitSpline(self,Nodes,UDegree,VDegree,u_pole_num,v_pole_num,U_node,V_node)
            % fit Spline base on CST class function
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
                error('SurfaceCST.fitSpline: pole_num more than node_num')
            end

            % default value
            if isempty(U_node)
                U_node=vecnorm(Nodes(:,2:end,:)-Nodes(:,1:end-1,:),2,3);
                U_node=mean(U_node,1);U_node=[0;cumsum(U_node')];
            end
            U_node=(U_node(:)-min(U_node))/(max(U_node)-min(U_node));
            if isempty(V_node)
                V_node=vecnorm(Nodes(2:end,:,:)-Nodes(1:end-1,:,:),2,3);
                V_node=mean(V_node,2);V_node=[0;cumsum(V_node)];
            end
            V_node=(V_node(:)-min(V_node))/(max(V_node)-min(V_node));

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
            v_fit_matrix=zeros(v_node_num,v_pole_num);
            [N_list,idx_srt,idx_end]=baseFcnN(V_node,VDegree,v_list);
            for deg_idx=1:VDegree+1
                idx=sub2ind([v_node_num,v_pole_num],(1:(v_node_num))',idx_srt+(deg_idx-1));
                v_fit_matrix(idx)=N_list(:,deg_idx);
            end
            u_fit_matrix=zeros(u_node_num,u_pole_num);
            [N_list,idx_srt,idx_end]=baseFcnN(U_node,UDegree,u_list);
            for deg_idx=1:UDegree+1
                idx=sub2ind([u_node_num,u_pole_num],(1:(u_node_num))',idx_srt+(deg_idx-1));
                u_fit_matrix(idx)=N_list(:,deg_idx);
            end
            u_fit_matrix=u_fit_matrix';

            % add class coefficient
            V_class=V_node;
            if self.sym_y,V_class=(V_node/2)+0.5;end
            U_class=U_node;
            if self.sym_x,U_class=(U_node/2)+0.5;end

            matrix_v_class_X=v_fit_matrix.*self.class_fcn_X(V_class);
            matrix_v_class_Y=v_fit_matrix;
            matrix_v_class_ZV=v_fit_matrix.*self.class_fcn_ZV(0.5*ones(size(V_class)),V_class);

            matrix_u_class_X=u_fit_matrix;
            matrix_u_class_Y=u_fit_matrix.*self.class_fcn_Y(U_class');
            matrix_u_class_ZU=u_fit_matrix.*self.class_fcn_ZU(U_class',0.5*ones(size(U_class')));

            % reverse calculate control point
            Poles(:,:,1)=matrix_v_class_X\Nodes(:,:,1)/matrix_u_class_X;
            Poles(:,:,2)=matrix_v_class_Y\Nodes(:,:,2)/matrix_u_class_Y;
            Poles(:,:,3)=matrix_v_class_ZV\Nodes(:,:,3)/matrix_u_class_ZU;
            Weights=ones(v_pole_num,u_pole_num);

            % add spline and shape function
            self.spline=Surface(Poles,UDegree,VDegree,UMults,VMults,UKnots,VKnots,Weights);
            self.shape_fcn=@(U_x,V_x) self.spline.calPoint(U_x,V_x);

            % fit data
            self.fit_data.Nodes=Nodes;
            self.fit_data.U_node=U_node;
            self.fit_data.V_node=V_node;
            self.fit_data.u_fit_matrix=u_fit_matrix;
            self.fit_data.v_fit_matrix=v_fit_matrix;
        end
    end

    methods % deform, rotate, translate
        function self=addDeform(self,deform_fcn_X,deform_fcn_Y,deform_fcn_Z)
            % base on local coordinate deform surface
            %
            % input:
            % tran_fcn_X(V), tran_fcn_Y(U), tran_fcn_Z(U,V)
            %
            self.deform_fcn_X=deform_fcn_X;
            self.deform_fcn_Y=deform_fcn_Y;
            self.deform_fcn_Z=deform_fcn_Z;
        end

        function self=addRotate(self,ang_x,ang_y,ang_z)
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

        function self=addTranslate(self,tran_x,tran_y,tran_z)
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

    methods % calculate coord
        function [U,V]=calCoordinate(self,Points,geom_torl)
            % base on X, Y, Z calculate local coordinate in surface
            %
            if nargin < 3, geom_torl=[];end
            if isempty(geom_torl), geom_torl=sqrt(eps);end

            % find point to start
            [U,V]=findNearest(self,Points,20);

            % use project function to adjust parameter
            [U,V]=self.projectPoint(Points,geom_torl,U,V);
        end

        function [U,V]=projectPoint(self,Points_init,geom_torl,U,V)
            % adjust U, V by Jacobian transformation
            % also can project point to surface
            %
            if nargin < 5
                U=[];V=[];
                if nargin < 3
                    geom_torl=[];
                end
            end

            % find point to start
            if isempty(U) && isempty(V)
                [U,V]=findNearest(self,Points_init,20);
            end

            [v_num,u_num,~]=size(Points_init);
            num=u_num*v_num;
            Points_init=reshape(Points_init,num,[]);U=U(:);V=V(:);

            % iteration
            iter=0;iter_max=50;
            done=false;idx=1:v_num*u_num;
            while ~done
                [Points,dPoints_dU,dPoints_dV]=self.calGradient(U(idx),V(idx));
                dPoint=reshape(Points_init(idx,:),length(idx),1,[])-Points;

                RU_RU=sum(dPoints_dU.^2,3);
                RV_RV=sum(dPoints_dV.^2,3);
                RU_RV=sum(dPoints_dU.*dPoints_dV,3);
                RU_D=sum(dPoints_dU.*dPoint,3);
                RV_D=sum(dPoints_dV.*dPoint,3);
                RRRR_RR=RU_RU.*RV_RV-(RU_RV).^2;
                dU=(RU_D.*RV_RV-RV_D.*RU_RV)./RRRR_RR;
                dV=(RV_D.*RU_RU-RU_D.*RU_RV)./RRRR_RR;
                dU(isnan(dU) | isinf(dU))=0;
                dV(isnan(dV) | isinf(dV))=0;

                U(idx)=U(idx)+dU;
                V(idx)=V(idx)+dV;
                U=max(U,0);U=min(U,1);
                V=max(V,0);V=min(V,1);

                idx=find(sum(abs(dPoint),3) > geom_torl);

                % Points_inv=self.calPoint(U,V);
                % scatter3(Points_inv(:,:,1),Points_inv(:,:,2),Points_inv(:,:,3));

                iter=iter+1;
                if isempty(idx) || iter >= iter_max
                    done=true;
                end
            end

            U=reshape(U,v_num,u_num);
            V=reshape(V,v_num,u_num);
        end

        function [U,V]=findNearest(self,Points,param)
            % find nearest U, V in grid
            %
            if nargin < 3, param=[];end
            if isempty(param), param=20;end
            [v_num,u_num,~]=size(Points);

            % generate rough mesh to initialize pre coord
            sample_grid=((0:(param-1))+(1:param))/2/param;
            [U_base,V_base]=meshgrid(sample_grid,sample_grid);
            Point_base=self.calPoint(U_base,V_base);
            U=zeros(v_num,u_num);
            V=zeros(v_num,u_num);
            for v_idx=1:v_num
                for u_idx=1:u_num
                    point=Points(v_idx,u_idx,:);
                    dis_abs=sum(abs(Point_base-point),3);
                    [~,idx]=min(dis_abs,[],"all");
                    U(v_idx,u_idx)=U_base(idx(1));
                    V(v_idx,u_idx)=V_base(idx(1));
                end
            end
        end

    end
end
