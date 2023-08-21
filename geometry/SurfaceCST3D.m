classdef SurfaceCST3D < handle
    % generate surface by CST parameter
    % process order:
    % generate origin surface, deform, rotation(global y z x), translation
    %
    properties
        % base parameter
        name='';
        LX;
        LY;
        LZ;
        deform_ID; % 1 is x, 2 is y
        symmetry_y;

        % origin parameter
        shape_fcn_X=[]; % (V)
        shape_fcn_Y=[]; % (U)
        shape_fcn_Z=[]; % (U, V)
        N1_fcn_Z=[]; % (U)
        N2_fcn_Z=[]; % (U)
        class_fcn_Z=[]; % (N1_fcn_Z, N2_fcn_Z, V)

        % origin gradient parameter
        dshape_fcn_X=[]; % (V)
        dshape_fcn_Y=[]; % (U)
        dshape_fcn_Z=[]; % (U, V)
        dclass_fcn_Z=[]; % (N1_fcn_Z, N2_fcn_Z, V)

        % deform parameter
        tran_fcn_X=[]; % (V)
        tran_fcn_Y=[]; % (U)
        tran_fcn_Z=[]; % (U,V)

        % gradient deform parameter
        dtran_fcn_X=[]; % (V)
        dtran_fcn_Y=[]; % (U)
        dtran_fcn_Z=[]; % (U,V)

        % rotation parameter
        rotation_matrix=[];

        % translation parameter
        translation=[];
    end

    % define surface function
    methods
        function self=SurfaceCST3D(name,LX,LY,LZ,shape_par_X,shape_par_Y,shape_par_Z,class_par_Z,symmetry_y)
            % generate 3D CST surface by LX,LY,LZ,shape_par_Y,shape_par_Z,class_par_Z
            %
            % notice:
            % if input N1, N2 == 0, shape_fun will equal to 1
            % shape_fcn_X or shape_par_Y should be empty
            % default function is (u/v)^N1*(1-u/v)^N2
            % u, x=X(v), shape_fcn_X(V)
            % v, y=Y(u), shape_fcn_Y(U)
            % w, z=Z(u,v), shape_fcn_Z(U,V)(mainly control U direction), class_fcn_Z{N1_fcn_Z(U), N2_fcn_Z(U)}(mainly control V direction)
            %
            % input:
            % name, LX, LY, LZ, shape_par_X, shape_par_Y, shape_par_Z, class_par_Z, symmetry_y
            %
            if nargin < 9 || isempty(symmetry_y)
                self.symmetry_y=false(1);
                if nargin < 8
                    class_par_Z=[];
                    if nargin < 7
                        shape_par_Z=[];
                        if nargin < 6
                            shape_par_Y=[];
                            if nargin < 5
                                shape_par_X=[];
                            end
                        end
                    end
                end
            else
                self.symmetry_y=symmetry_y;
            end
            self.name=name;
            self.LX=LX;
            self.LY=LY;
            self.LZ=LZ;

            if ~isempty(shape_par_Y) && ~isempty(shape_par_X)
                error('SurfaceCST3D: donot support both X and Y deform')
            elseif ~isempty(shape_par_X)
                % deform direction x
                self.deform_ID=1;
            elseif ~isempty(shape_par_Y)
                % deform direction y
                self.deform_ID=2;
            else
                self.deform_ID=0;
            end

            if ~isempty(shape_par_X) && isnumeric(shape_par_X)
                % shape_par_Y only one cross-section parameter
                self.shape_fcn_X=@(V) self.defunShape(V,shape_par_X(1),shape_par_X(2));
                % self.dshape_fcn_X=@(V) self.defundShape(V,shape_par_X(1),shape_par_X(2));
                self.dshape_fcn_X=@(V) self.differFcn(self.shape_fcn_X,V);
            elseif isa(shape_par_X, 'function_handle')
                % function handle
                self.shape_fcn_X=shape_par_X;
                self.dshape_fcn_X=@(V) self.differFcn(shape_par_X,V);
            end

            if ~isempty(shape_par_Y) && isnumeric(shape_par_Y)
                % shape_par_Y only one cross-section parameter
                self.shape_fcn_Y=@(U) self.defunShape(U,shape_par_Y(1),shape_par_Y(2));
                % self.dshape_fcn_Y=@(U) self.defundShape(U,shape_par_Y(1),shape_par_Y(2));
                self.dshape_fcn_Y=@(U) self.differFcn(self.shape_fcn_Y,U);
            else
                % function handle
                self.shape_fcn_Y=shape_par_Y;
                self.dshape_fcn_Y=@(U) self.differFcn(shape_par_Y,U);
            end

            if isempty(shape_par_Z) && isempty(class_par_Z)
                self.shape_fcn_Z=[];
                self.N1_fcn_Z=[];
                self.N2_fcn_Z=[];
                self.class_fcn_Z=[];
            else
                if ~isempty(shape_par_Z) && (isnumeric(shape_par_Z) || iscell(shape_par_Z))
                    % shape_par_Z can have two parameter V direction cross-section parameter
                    % {V section 1 parameter, V section 2 parameter}
                    if isnumeric(shape_par_Z)
                        shape_par_Z={shape_par_Z};
                    end
                    if length(shape_par_Z) == 1
                        shape_par_Z=[shape_par_Z,shape_par_Z];
                    end
                    self.shape_fcn_Z=@(U,V) self.defunShapeZ(U,V,shape_par_Z{1},shape_par_Z{2});
                    % self.dshape_fcn_Z=@(U,V) self.defundShapeZ(U,V,shape_par_Z{1},shape_par_Z{2});
                    self.dshape_fcn_Z=@(U,V) self.differFcn2(self.shape_fcn_Z,U,V);
                elseif isa(shape_par_Z, 'function_handle')
                    % function handle
                    self.shape_fcn_Z=shape_par_Z;
                    self.dshape_fcn_Z=@(U,V) self.differFcn2(shape_par_Z,U,V);
                end

                if isa(class_par_Z, 'function_handle')
                    self.class_fcn_Z=class_par_Z;
                elseif ~isempty(class_par_Z)
                    % class_par_Z can have two parameter U direction cross-section parameter
                    % {U section 1 parameter, U section 2 parameter}
                    if isnumeric(class_par_Z)
                        class_par_Z={class_par_Z,class_par_Z};
                    elseif iscell(class_par_Z) && length(class_par_Z) == 1
                        class_par_Z=[class_par_Z,class_par_Z];
                    end

                    if isnumeric(class_par_Z{1})
                        self.N1_fcn_Z=@(U) class_par_Z{1}(1).*(1-U)+class_par_Z{2}(1).*U;
                    else
                        % function handle
                        self.N1_fcn_Z=class_par_Z{1};
                    end

                    if isnumeric(class_par_Z{2})
                        self.N2_fcn_Z=@(U) class_par_Z{2}(2).*(1-U)+class_par_Z{1}(2).*U;
                    else
                        % function handle
                        self.N2_fcn_Z=class_par_Z{2};
                    end
                    self.class_fcn_Z=@(U,V) self.defunClassZ(U,V);
                end
                self.dclass_fcn_Z=@(U,V) self.differFcn2(self.class_fcn_Z,U,V);
            end

        end

        function addDeform(self,tran_fcn_X,tran_fcn_Y,tran_fcn_Z)
            % base on local coordinate deform or translation surface
            %
            % input:
            % tran_fcn_X(V), tran_fcn_Y(U), tran_fcn_Z(U,V)
            %
            if self.deform_ID == 1
                if ~isempty(tran_fcn_Y)
                    error('SurfaceCST3D: donot support both X and Y deform')
                end
            elseif self.deform_ID == 2
                if ~isempty(tran_fcn_X)
                    error('SurfaceCST3D: donot support both X and Y deform')
                end
            end

            self.tran_fcn_X=tran_fcn_X;
            self.tran_fcn_Y=tran_fcn_Y;
            self.tran_fcn_Z=tran_fcn_Z;
            self.dtran_fcn_X=@(V) self.differFcn(tran_fcn_X,V);
            self.dtran_fcn_Y=@(U) self.differFcn(tran_fcn_Y,U);
            self.dtran_fcn_Z=@(U,V) self.differFcn2(tran_fcn_Z,U,V);
        end

        function addRotation(self,ang_x,ang_y,ang_z)
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

            self.rotation_matrix=matrix;
        end

        function addTranslation(self,tran_x,tran_y,tran_z)
            % base on angle to rotation surface
            %
            self.translation=[tran_x,tran_y,tran_z];
        end

    end

    % common function
    methods
        function S=defunShape(self,U,N1,N2)
            % default shape function of X and Y
            %
            NP=self.calNormPar(N1,N2);
            S=U.^N1.*(1-U).^N2/NP;
        end

        function G=defundShape(self,U,N1,N2)
            % default shape function of X and Y
            %
            NP=self.calNormPar(N1,N2);
            G=(N1*U.^(N1-1).*(1-U).^N2-N2*U.^N1.*(1-U).^(N2-1))/NP;

            % numerical stable
            G(G < -1e4)=-1e4;
            G(G > 1e4)=1e4;
        end

        function S=defunShapeZ(self,U,V,NS1,NS2)
            % default shape function of Z
            % default N1, N2 is relation with V
            %
            N1=NS1(1).*(1-V)+NS2(1).*V;
            N2=NS1(2).*(1-V)+NS2(2).*V;
            NP=self.calNormPar(N1,N2);
            S=U.^N1.*(1-U).^N2./NP;
        end

        function [G_U,G_V]=defundShapeZ(self,U,V,NS1,NS2)
            % default shape function of Z
            % default N1, N2 is relation with V
            %
            N1=NS1(1).*(1-V)+NS2(1).*V;
            N2=NS1(2).*(1-V)+NS2(2).*V;
            NP=self.calNormPar(N1,N2);
            G_U=(N1.*U.^(N1-1).*(1-U).^N2-N2.*U.^N1.*(1-U).^(N2-1))./NP;

            dN1_dV=-NS1(1)+NS2(1);
            dN2_dV=-NS1(2)+NS2(2);
            G_V=(dN1_dV.*log(U).*U.^N1.*(1-U).^N2+dN2_dV.*log(1-U).*U.^N1.*(1-U).^N2)./NP;
            G_V(U == 0 | U == 1)=-1e4;

            % numerical stable
            G_U(G_U < -1e4)=-1e4;
            G_U(G_U > 1e4)=1e4;
            G_V(G_V < -1e4)=-1e4;
            G_V(G_V > 1e4)=1e4;
        end

        function C=defunClassZ(self,U,V)
            % default class function of Z
            % default N1, N2 is relation with U
            %
            N1=self.N1_fcn_Z(U);
            N2=self.N2_fcn_Z(U);
            NP=self.calNormPar(N1,N2);
            C=V.^N1.*(1-V).^N2./NP;
        end

        function nomlz_par=calNormPar(self,N1,N2)
            % calculate normailize class function parameter by N1, N2
            %
            nomlz_par=(N1./(N1+N2)).^N1.*(N2./(N1+N2)).^N2;
            nomlz_par((N1 == 0) & (N2 == 0))=1;
        end

        function G=differFcn(self,fcn,U)
            % differ to get gradient
            %
            step=1e-5;
            S=fcn(U);
            S_F=fcn(U+step);
            G=(S_F-S)/step;
            bool=~isreal(G) | isnan(G);
            if any(bool) %% try back walk differ
                S_B=fcn(U(bool)-step);
                G(bool)=(S(bool)-S_B)/step;
            end

            % numerical stable
            G(G < -1e4)=-1e4;
            G(G > 1e4)=1e4;
        end

        function [G_U,G_V]=differFcn2(self,fcn,U,V)
            % differ to get gradient
            %
            step=1e-5;
            S=fcn(U,V);

            S_F=fcn(U+step,V);
            G_U=(S_F-S)/step;
            bool=~isreal(G_U) | isnan(G_U);
            if any(bool) %% try back walk differ
                S_B=fcn(U(bool)-step,V(bool));
                G_U(bool)=(S(bool)-S_B)/step;
            end

            S_F=fcn(U,V+step);
            G_V=(S_F-S)/step;
            bool=~isreal(G_V) | isnan(G_V);
            if any(bool) %% try back walk differ
                S_B=fcn(U(bool),V(bool)-step);
                G_V(bool)=(S(bool)-S_B)/step;
            end

            % numerical stable
            G_U(G_U < -1e4)=-1e4;
            G_U(G_U > 1e4)=1e4;
            G_V(G_V < -1e4)=-1e4;
            G_V(G_V > 1e4)=1e4;
        end

    end

    % calculate grid coordinate function
    methods
        function [X,Y,Z,U,V]=calSurface(self,varargin)
            % generate 3D CST matrix by U, V or u_gird_number, v_gird_number
            % X=CX*U
            % Y=CY*shape_fcn_Y.*V, if symmetry_y, Y=CY*shape_fcn_Y.*(V-0.5)
            % Z=CZ*shape_fcn_Z.*class_fcn_Z
            %
            % default
            % u_list=linspace(0,1,u_gird_num(default is 20)+1)
            % v_list=linspace(0,1,v_gird_num(default is 20)+1), if symmetry_y v_list=linspace(0.5,1,v_gird_num+1);
            % [U,V]=meshgrid(u_list,v_list), colume is LaWGS format line
            % value_torl=1e-3
            %
            % u, x=X
            % v, y=Y(u)
            % w, z=Z(u,v)
            %
            % input:
            % U, V
            % u_grid_number, v_grid_number
            % value_torl
            %
            % notice:
            % shape_fcn_Y(U), shape_fcn_Z(U,V), class_fcn_Z{N1_fcn_Z(U), N2_fcn_Z(U)}
            %
            % output:
            % X,Y,Z (colume is LaWGS format line)
            %

            if nargin <= 2
                % input is torlance
                if nargin < 2
                    value_torl=1e-3;
                else
                    value_torl=varargin{1};
                end
                max_level=50;

                % adapt capture U, V
                if self.symmetry_y
                    low_bou=[0,0.5];
                else
                    low_bou=[0,0];
                end
                up_bou=[1,1];
                [U,V,Fval,~]=meshAdapt2DM(@(x) coordFcn(self,x),low_bou,up_bou,value_torl,max_level,3);
                X=Fval(:,:,1);
                Y=Fval(:,:,2);
                Z=Fval(:,:,3);
            elseif nargin == 3
                % input is U, V or u_grid_number, v_grid_number
                if length(varargin{1}) == 1
                    U=[];
                    u_grid_numebr=varargin{1};
                else
                    % mean input U matrix
                    U=varargin{1};
                    u_grid_numebr=size(U,2)-1;
                end

                if length(varargin{2}) == 1
                    V=[];
                    v_grid_numebr=varargin{2};
                else
                    % mean input V matrix
                    V=varargin{2};
                    v_grid_numebr=size(V,1)-1;
                end

                % calculate local coordinate matrix
                if isempty(U)
                    u_list=linspace(0,1,u_grid_numebr+1);
                    U=repmat(u_list,v_grid_numebr+1,1);
                end
                if isempty(V)
                    v_list=linspace(0,1,v_grid_numebr+1)';
                    if self.symmetry_y
                        v_list=v_list/2+0.5;
                    end
                    V=repmat(v_list,1,u_grid_numebr+1);
                end
            
                [X,Y,Z]=self.calPoint(U,V);
            end
            
            function fval=coordFcn(self,x)
                [x,y,z]=self.calPoint(x(1),x(2));
                fval=[x,y,z];
            end
        end

        function [X,Y,Z]=calPoint(self,U,V)
            % calculate point on surface
            %

            % calculate origin surface matrix
            X=self.LX*U;
            if ~isempty(self.shape_fcn_X)
                X=X.*self.shape_fcn_X(V);
            end

            if self.symmetry_y
                Y=self.LY*(V-0.5);
            else
                Y=self.LY*V;
            end
            if ~isempty(self.shape_fcn_Y)
                Y=Y.*self.shape_fcn_Y(U);
            end

            if isempty(self.shape_fcn_Z) && isempty(self.class_fcn_Z)
                Z=zeros(size(X));
            elseif isempty(self.shape_fcn_Z)
                Z=self.LZ*self.class_fcn_Z(U,V);
            elseif isempty(self.class_fcn_Z)
                Z=self.LZ*self.shape_fcn_Z(U,V);
            else
                Z=self.LZ*self.shape_fcn_Z(U,V).*self.class_fcn_Z(U,V);
            end

            % deform surface
            if ~isempty(self.tran_fcn_X)
                X=X+self.tran_fcn_X(V);
            end
            if ~isempty(self.tran_fcn_Y)
                Y=Y+self.tran_fcn_Y(U);
            end
            if ~isempty(self.tran_fcn_Z)
                Z=Z+self.tran_fcn_Z(U,V);
            end

            % rotation surface
            if ~isempty(self.rotation_matrix)
                matrix=self.rotation_matrix;
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

        function [U,V,X,Y,Z]=calCoordinate(self,X,Y,Z)
            % base on X, Y, Z calculate local coordinate in surface
            %
            XO=X;YO=Y;ZO=Z;geo_torl=1e-6;

            % undone translation surface
            if ~isempty(self.translation)
                X=X-self.translation(1);
                Y=Y-self.translation(2);
                Z=Z-self.translation(3);
            end

            % undone rotation surface
            if ~isempty(self.rotation_matrix)
                matrix=self.rotation_matrix';
                X_old=X;Y_old=Y;Z_old=Z;
                X=matrix(1,1)*X_old+matrix(1,2)*Y_old+matrix(1,3)*Z_old;
                Y=matrix(2,1)*X_old+matrix(2,2)*Y_old+matrix(2,3)*Z_old;
            end

            % identify base local coordinate
            if self.deform_ID == 1 || self.deform_ID == 0
                if self.LY > 0
                    Y=max(Y,0);Y=min(Y,self.LY);
                else
                    Y=max(Y,self.LY);Y=min(Y,0);
                end

                V=Y./self.LY;
                if self.symmetry_y
                    V=V+0.5;
                end

                % re deform surface
                if ~isempty(self.tran_fcn_X)
                    X=X-self.tran_fcn_X(V);
                end

                if self.LX > 0
                    X=max(X,0);
                else
                    X=min(X,0);
                end
                U=X./self.LX;
                if ~isempty(self.shape_fcn_X)
                    shape=self.shape_fcn_X(V);
                    shape(shape==0)=1;
                    U=U./shape;
                end

                % judge local coordinate
                U=min(U,1);
            elseif self.deform_ID == 2
                if self.LX > 0
                    X=max(X,0);X=min(X,self.LX);
                else
                    X=max(X,self.LX);X=min(X,0);
                end

                U=X./self.LX;

                % re deform surface
                if ~isempty(self.tran_fcn_Y)
                    Y=Y-self.tran_fcn_Y(U);
                end

                V=Y./self.LY;
                if ~isempty(self.shape_fcn_Y)
                    shape=self.shape_fcn_Y(U);
                    shape(shape==0)=1;
                    V=V./shape;
                end

                % judge local coordinate
                if self.symmetry_y
                    V=V+0.5;
                end
                V=min(V,1);
            end

            U=max(U,0);U=min(U,1);
            V=max(V,0);V=min(V,1);

            % use project function to adjust parameter
            [X,Y,Z,U,V]=self.calProject(XO,YO,ZO,U,V,geo_torl);
        end
    end

    % calculate project point coordinate function
    methods
        function [dX_dU,dY_dU,dZ_dU,dX_dV,dY_dV,dZ_dV]=calGradient(self,U,V)
            % calculate x, y, z gradient to u, v
            %
            
            % calculate origin surface matrix gradient
            dX_dU=self.LX*ones(size(U));
            if ~isempty(self.shape_fcn_X)
                dX_dU=dX_dU.*(self.shape_fcn_X(V));
                dX_dV=dX_dU.*U.*self.dshape_fcn_X(V);
            else
                dX_dV=zeros(size(U));
            end

            dY_dV=self.LY*ones(size(V));
            if ~isempty(self.shape_fcn_Y)
                dY_dV=dY_dV.*(self.shape_fcn_Y(U));
                dY_dU=dY_dV.*V.*self.dshape_fcn_Y(U);
            else
                dY_dU=zeros(size(V));
            end

            if isempty(self.shape_fcn_Z) && isempty(self.class_fcn_Z)
                dZ_dU=zeros(size(U));dZ_dV=zeros(size(V));
            elseif isempty(self.shape_fcn_Z)
                [dZ_dU,dZ_dV]=self.dclass_fcn_Z(U,V);
                dZ_dU=dZ_dU*self.LZ;dZ_dV=dZ_dV*self.LZ;
            elseif isempty(self.class_fcn_Z)
                [dZ_dU,dZ_dV]=self.dshape_fcn_Z(U,V);
                dZ_dU=dZ_dU*self.LZ;dZ_dV=dZ_dV*self.LZ;
            else
                [dSZ_dU,dSZ_dV]=self.dshape_fcn_Z(U,V);
                [dCZ_dU,dCZ_dV]=self.dclass_fcn_Z(U,V);
                dZ_dU=(dSZ_dU.*self.class_fcn_Z(U,V)+self.shape_fcn_Z(U,V).*dCZ_dU)*self.LZ;
                dZ_dV=(dSZ_dV.*self.class_fcn_Z(U,V)+self.shape_fcn_Z(U,V).*dCZ_dV)*self.LZ;
            end

            % calculate deform surface gradient
            if ~isempty(self.tran_fcn_X)
                dX_dV=dX_dV+self.dtran_fcn_X(V);
            end
            if ~isempty(self.tran_fcn_Y)
                dY_dU=dY_dU+self.dtran_fcn_Y(U);
            end
            if ~isempty(self.tran_fcn_Z)
                [dTZ_dU,dTZ_dV]=self.dtran_fcn_Z(U,V);
                dZ_dU=dZ_dU+dTZ_dU;
                dZ_dV=dZ_dU+dTZ_dV;
            end

            % calculate rotation surface gradient
            if ~isempty(self.rotation_matrix)
                matrix=self.rotation_matrix;
                dX_dU_old=dX_dU;dY_dU_old=dY_dU;dZ_dU_old=dZ_dU;
                dX_dV_old=dX_dV;dY_dV_old=dY_dV;dZ_dV_old=dZ_dV;
                dX_dU=matrix(1,1)*dX_dU_old+matrix(1,2)*dY_dU_old+matrix(1,3)*dZ_dU_old;
                dX_dV=matrix(1,1)*dX_dV_old+matrix(1,2)*dY_dV_old+matrix(1,3)*dZ_dV_old;
                dY_dU=matrix(2,1)*dX_dU_old+matrix(2,2)*dY_dU_old+matrix(2,3)*dZ_dU_old;
                dY_dV=matrix(2,1)*dX_dV_old+matrix(2,2)*dY_dV_old+matrix(2,3)*dZ_dV_old;
                dZ_dU=matrix(3,1)*dX_dU_old+matrix(3,2)*dY_dU_old+matrix(3,3)*dZ_dU_old;
                dZ_dV=matrix(3,1)*dX_dV_old+matrix(3,2)*dY_dV_old+matrix(3,3)*dZ_dV_old;
            end

        end

        function [dX_dU,dY_dU,dZ_dU,dX_dV,dY_dV,dZ_dV]=calGradientDiffer(self,U,V,X,Y,Z)
            % use differ ot calculate gradient
            %
            if nargin < 5
                [X,Y,Z]=self.calPoint(U,V);
            end
            step=1e-5;

            [X_UF,Y_UF,Z_UF]=self.calPoint(U+step,V);
            dX_dU=(X_UF-X)/step;dY_dU=(Y_UF-Y)/step;dZ_dU=(Z_UF-Z)/step;
            bool=~isreal(dX_dU) | isnan(dX_dU) |...
                ~isreal(dY_dU) | isnan(dY_dU) |...
                ~isreal(dZ_dU) | isnan(dZ_dU);
            if any(bool) %% try back walk differ
                [X_UB,Y_UB,Z_UB]=self.calPoint(U(bool)-step,V(bool));
                dX_dU(bool)=(X(bool)-X_UB)/step;dY_dU(bool)=(Y(bool)-Y_UB)/step;dZ_dU(bool)=(Z(bool)-Z_UB)/step;
            end

            [X_VF,Y_VF,Z_VF]=self.calPoint(U,V+step);
            dX_dV=(X_VF-X)/step;dY_dV=(Y_VF-Y)/step;dZ_dV=(Z_VF-Z)/step;
            bool=~isreal(dX_dV) | isnan(dX_dV) |...
                ~isreal(dY_dV) | isnan(dY_dV) |...
                ~isreal(dZ_dV) | isnan(dZ_dV);
            if any(bool) %% try back walk differ
                [X_VB,Y_VB,Z_VB]=self.calPoint(U(bool),V(bool)-step);
                dX_dV(bool)=(X(bool)-X_VB)/step;dY_dV(bool)=(Y(bool)-Y_VB)/step;dZ_dV(bool)=(Z(bool)-Z_VB)/step;
            end
        end

        function [X,Y,Z,U,V]=calProject(self,XO,YO,ZO,U,V,geo_torl)
            % adjust u, v by Jacobian transformation
            % also can project point to surface
            %
            if nargin < 7
                geo_torl=1e-6;
            end

            iter=0;
            iter_max=10;

            [X,Y,Z]=self.calPoint(U,V);
%             scatter3(X,Y,Z);
            dX=XO-X;dY=YO-Y;dZ=ZO-Z;
            boolean=(abs(dX)+abs(dY)+abs(dZ)) > geo_torl;
            while any(any(boolean)) && iter < iter_max
                [dX_dU,dY_dU,dZ_dU,dX_dV,dY_dV,dZ_dV]=self.calGradientDiffer(U(boolean),V(boolean),X(boolean),Y(boolean),Z(boolean));

                RU_RU=dX_dU.^2+dY_dU.^2+dZ_dU.^2;
                RV_RV=dX_dV.^2+dY_dV.^2+dZ_dV.^2;
                RU_RV=dX_dU.*dX_dV+dY_dU.*dY_dV+dZ_dU.*dZ_dV;
                RU_D=dX_dU.*dX(boolean)+dY_dU.*dY(boolean)+dZ_dU.*dZ(boolean);
                RV_D=dX_dV.*dX(boolean)+dY_dV.*dY(boolean)+dZ_dV.*dZ(boolean);
                RRRR_RR=RU_RU.*RV_RV-(RU_RV).^2;
                dU=(RU_D.*RV_RV-RV_D.*RU_RV)./RRRR_RR;
                dV=(RV_D.*RU_RU-RU_D.*RU_RV)./RRRR_RR;

                U(boolean)=U(boolean)+dU;
                V(boolean)=V(boolean)+dV;
                U=max(U,0);U=min(U,1);
                V=max(V,0);V=min(V,1);

                [X,Y,Z]=self.calPoint(U,V);
%                 scatter3(X,Y,Z);

                dX=XO-X;dY=YO-Y;dZ=ZO-Z;
                iter=iter+1;
                boolean=(abs(dU)+abs(dV)) > geo_torl;
            end
        end

    end

    % external function
    methods
        function drawSurface()

        end
        
        function surface_BSpline=getSurfaceBSpline(self,varargin)
            % convert CST surface into BSpline surface
            %
            % input:
            % U,V
            % u_grid_number, v_grid_number
            % value torlance
            %

            if nargin == 1
                [X,Y,Z,U,V]=calSurface(self);
            elseif nargin == 2
                [X,Y,Z,U,V]=calSurface(self,varargin{1});
            elseif nargin == 3
                U=varargin{1};
                V=varargin{2};
                if length(U) > 1 && U(1,1) > U(1,2)
                    U=fliplr(U);
                end
                if length(V) > 1 && V(1,1) > V(2,1)
                    V=flipud(V);
                end
                [X,Y,Z,U,V]=calSurface(self,U,V);
            else
                error('SurfaceCST3D.getBSpline: error input');
            end

            if size(X,2) == 2
                u_list=[0,0,1,1];
            else
                u_list=[0,0,0,U(1,:),1,1,1];
            end
            if size(X,1) == 2
                v_list=[0,0,1,1];
            else
                v_list=[0,0,0,V(:,1)',1,1,1];
            end
            surface_BSpline=SurfaceBSpline(self.name,[],[],[],X,Y,Z,[],[],u_list,v_list);
        end

        function [step_str,object_index,OPEN_SHELL_index]=getStepNode(self,object_index,varargin)
            % interface of BSpline surface getStepNode function
            %
            if nargin == 2
                surf_BSpline=self.getSurfaceBSpline();
            elseif nargin == 3
                surf_BSpline=self.getSurfaceBSpline(varargin{1});
            else
                surf_BSpline=self.getSurfaceBSpline(varargin{1},varargin{2});
            end

            [step_str,object_index,OPEN_SHELL_index]=surf_BSpline.getStepNode(object_index);
        end
    end
end

function Y=calDicrete(fcn,U,U_bou)
Y_bou=fcn(U_bou);
u_number=length(U);
Y=U;
for u_index=1:u_number
    Y(u_index)=interpLinear(U(u_index),U_bou,Y_bou);
end
end

function y_pred=interpLinear(x_pred,X,Y)
num=length(X);
idx=num; % search start from last one, find out X samll than x
while ((idx > 1) && (X(idx) > x_pred))
    idx=idx-1;
end

if (idx == num)
    % out interp
    y_pred=(Y(end)-Y(end-1))/(X(end)-X(end-1))*(x_pred-X(end))+Y(end);
elseif (idx == 0)
    y_pred=(Y(2)-Y(1))/(X(2)-X(1))*(x_pred-X(1))+Y(1);
else
    % linear interpolation
    y_pred=Y(idx)+...
        (Y(idx+1)-Y(idx))*...
        (x_pred-X(idx))/...
        (X(idx+1)-X(idx));
end
end

function [U,V,Fval,node_list]=meshAdapt2DM(fcn,low_bou,up_bou,torl,max_level,fval_num)
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

[node_num,data_num]=createNodeTree(1,4); % create node tree from root
node_list=node_list(1:node_num,:);
data_list=data_list(1:data_num,:);

% % inorder traversal to obtain u and v aus coordinate
% idx_list=zeros(data_num,2); % data idx in u_list
% u_list=traversalInorder(1,1,[11,13],[12,14],[6,9,10]);
% v_list=traversalInorder(1,2,[11,12],[13,14],[7,8,10]);
% 
% % add boundary
% u_list=[low_bou(1),u_list,up_bou(1)];
% v_list=[low_bou(2),v_list,up_bou(2)];

[u_list,~,u_idx]=unique(data_list(:,1));
[v_list,~,v_idx]=unique(data_list(:,2));
idx_list=[u_idx,v_idx];

[U,V]=meshgrid(u_list,v_list);

% local data to Fval
Fval=nan(length(v_list),length(u_list),fval_num);

for data_idx=1:data_num
    if all([idx_list(data_idx,2),idx_list(data_idx,1)] ~= 0)
        Fval(idx_list(data_idx,2),idx_list(data_idx,1),:)=data_list(data_idx,3:end);
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

    function x_list=traversalInorder(root_idx,coord_idx,idx_1,idx_2,place_c)
        % inorder traversal node tree to get u_list or v_list
        % if input coord_idx=1, idx_1=[cell_1,cell_3], idx_2=[cell_2,cell_4], traversal is u_list
        % if input coord_idx=2, idx_1=[cell_1,cell_2], idx_2=[cell_3,cell_4], traversal is v_list
        %
        stack={};
        x_list=[];
        node_idx=root_idx;

        while length(stack)~=0 || ~isempty(node_idx)
            while ~isempty(node_idx)
                % only add one node_idx is enough
                % because we only need center coordinate of one cell
                stack=[stack,{node_idx}];

                % left node
                node_idx=node_list(node_idx,idx_1);
                node_idx=node_idx(node_idx~=0);
            end

            node_idx=stack{end};
            stack=stack(1:end-1);
            % obtain center place coordinate of node
            data_idx=node_list(node_idx,place_c);
            data_idx=data_idx(data_idx~=0);

            if ~isempty(data_idx)
                x_list=[x_list,data_list(data_idx(1),coord_idx)];
                idx_list(data_idx,coord_idx)=length(x_list)+1;
            end

            % right node
            node_idx=node_list(node_idx,idx_2);
            node_idx=node_idx(node_idx~=0);
        end
    end

end
