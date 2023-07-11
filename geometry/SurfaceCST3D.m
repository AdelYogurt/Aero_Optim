classdef SurfaceCST3D < handle
    % generate surface by CST parameter
    % process order:
    % generate origin surface, deform, rotation(global y z x), translation
    properties
        LX;
        LY;
        LZ;

        % origin parameter
        deform_ID; % 1 is x, 2 is y
        shape_fcn_X; % (PSI)
        shape_fcn_Y; % (XI)
        shape_fcn_Z; % (XI,PSI)
        N1_fcn_Z; % (XI)
        N2_fcn_Z; % (XI)
        class_fcn_Z;
        symmetry_y;

        % deform parameter
        tran_fcn_X=[]; % (PSI)
        tran_fcn_Y=[]; % (XI)
        tran_fcn_Z=[]; % (XI,PSI)

        % rotation parameter
        rotation_matrix=[];

        % translation parameter
        translation=[];
    end

    % define surface function
    methods
        function self=SurfaceCST3D(LX,LY,LZ,shape_par_X,shape_par_Y,shape_par_Z,class_par_Z,symmetry_y)
            % generate 3D CST surface by LX,LY,LZ,shape_par_Y,shape_par_Z,class_par_Z
            %
            % xi, x=X
            % psi, y=Y(xi)
            % zeta, z=Z(xi,psi)
            %
            % input:
            % LX,LY,LZ,shape_par_Y,shape_par_Z,class_par_Z,symmetry_y
            %
            % notice:
            % shape_fcn_X(PSI), shape_fcn_Y(XI), shape_fcn_Z(XI,PSI), class_fcn_Z{N1_fcn_Z(XI), N2_fcn_Z(XI)}
            % if input N1, N2 == 0, shape_fun will equal to 1
            % one of shape_fcn_X and shape_par_Y should be empty
            %
            if nargin < 8 || isempty(symmetry_y)
                self.symmetry_y=false(1);
                if nargin < 7
                    class_par_Z=[];
                    if nargin < 6
                        shape_par_Z=[];
                        if nargin < 5
                            shape_par_Y=[];
                            if nargin < 4
                                shape_par_X=[];
                            end
                        end
                    end
                end
            else
                self.symmetry_y=symmetry_y;
            end
            self.LX=LX;
            self.LY=LY;
            self.LZ=LZ;

            if ~isempty(shape_par_Y) && ~isempty(shape_par_X)
                error('CST3DSurface: donot support both X and Y deform')
            elseif ~isempty(shape_par_X)
                % deform direction x
                self.deform_ID=1;
            elseif ~isempty(shape_par_Y)
                % deform direction y
                self.deform_ID=2;
            else
                self.deform_ID=0;
            end

            if isempty(shape_par_X)
                self.shape_fcn_X=[];
            elseif isnumeric(shape_par_X)
                % shape_par_Y only one cross-section parameter
                self.shape_fcn_X=@(XI) self.defunShape...
                    (XI,shape_par_X(1),shape_par_X(2));
            else
                % function handle
                self.shape_fcn_X=shape_par_X;
            end

            if isempty(shape_par_Y)
                self.shape_fcn_Y=[];
            elseif isnumeric(shape_par_Y)
                % shape_par_Y only one cross-section parameter
                self.shape_fcn_Y=@(XI) self.defunShape...
                    (XI,shape_par_Y(1),shape_par_Y(2));
            else
                % function handle
                self.shape_fcn_Y=shape_par_Y;
            end

            if isempty(shape_par_Z) && isempty(class_par_Z)
                self.shape_fcn_Z=[];
                self.N1_fcn_Z=[];
                self.N2_fcn_Z=[];
                self.class_fcn_Z=[];
            else
                if isempty(shape_par_Z)
                    self.shape_fcn_Z=[];
                elseif isnumeric(shape_par_Z)
                    % shape_par_Z can have two parameter cross-section parameter
                    % default shape function donot
                    if size(shape_par_Z,1) == 1
                        shape_par_Z=[shape_par_Z;shape_par_Z];
                    end
                    self.shape_fcn_Z=@(XI,PSI) self.defunShapeZ...
                        (XI,PSI,shape_par_Z(1,1),shape_par_Z(1,2),shape_par_Z(2,1),shape_par_Z(2,2));
                else
                    % function handle
                    self.shape_fcn_Z=shape_par_Z;
                end

                if isa(class_par_Z, 'function_handle')
                    self.class_fcn_Z=class_par_Z;
                else
                    % class_par_Z is {N1,N2}, N1/N2 hace two parameter(means two cross-section)
                    if isnumeric(class_par_Z)
                        class_par_Z={class_par_Z,class_par_Z};
                    elseif iscell(class_par_Z) && length(class_par_Z) == 1
                        class_par_Z={class_par_Z{1},class_par_Z{1}};
                    elseif ~iscell(class_par_Z)
                        error('calCST: z class function parameter error input')
                    end

                    if isnumeric(class_par_Z{1})
                        self.N1_fcn_Z=@(XI) class_par_Z{1}(2).*XI+class_par_Z{1}(1).*(1-XI);
                    else
                        % function handle
                        self.N1_fcn_Z=class_par_Z{1};
                    end

                    if isnumeric(class_par_Z{2})
                        self.N2_fcn_Z=@(XI) class_par_Z{2}(2).*XI+class_par_Z{2}(1).*(1-XI);
                    else
                        % function handle
                        self.N2_fcn_Z=class_par_Z{2};
                    end
                    self.class_fcn_Z=@(XI,PSI) self.defunClassZ(XI,PSI);
                end
            end

        end

        function addDeform(self,tran_fcn_X,tran_fcn_Y,tran_fcn_Z)
            % base on local coordinate deform or translation surface
            %
            % input:
            % tran_fcn_X(PSI), tran_fcn_Y(XI), tran_fcn_Z(XI,PSI)
            %
            if self.deform_ID == 1
                if ~isempty(tran_fcn_Y)
                    error('CST3DSurface: donot support both X and Y deform')
                end
            elseif self.deform_ID == 2
                if ~isempty(tran_fcn_X)
                    error('CST3DSurface: donot support both X and Y deform')
                end
            end

            self.tran_fcn_X=tran_fcn_X;
            self.tran_fcn_Y=tran_fcn_Y;
            self.tran_fcn_Z=tran_fcn_Z;
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
        function S=defunShape(self,XI,N1,N2)
            % default shape function of X and Y
            %
            NP=self.calNormPar(N1,N2);
            S=XI.^N1.*(1-XI).^N2/NP;
        end

        function S=defunShapeZ(self,XI,PSI,N11,N12,N21,N22)
            % default shape function of Z
            % default N1, N2 is relation with PSI
            %
            N1=N21.*PSI+N11.*(1-PSI);
            N2=N22.*PSI+N12.*(1-PSI);
            NP=self.calNormPar(N1,N2);
            S=XI.^N1.*(1-XI).^N2./NP;
        end

        function C=defunClassZ(self,XI,PSI)
            % default class function of Z
            % default N1, N2 is relation with XI
            %
            N1=self.N1_fcn_Z(XI);
            N2=self.N2_fcn_Z(XI);
            NP=self.calNormPar(N1,N2);
            C=PSI.^N1.*(1-PSI).^N2./NP;
        end

        function nomlz_par=calNormPar(self,N1,N2)
            % calculate normailize class function parameter by N1, N2
            %
            nomlz_par=(N1./(N1+N2)).^N1.*(N2./(N1+N2)).^N2;
            nomlz_par((N1 == 0) & (N2 == 0))=1;
        end
    end

    % calculate grid coordinate function
    methods
        function [X,Y,Z,XI,PSI]=calSurface(self,varargin)
            % generate 3D CST matrix by XI, PSI or xi_gird_number, psi_gird_number
            % X=CX*XI
            % Y=CY*shape_fcn_Y.*PSI, if symmetry_y Y=CY*shape_fcn_Y.*(PSI-0.5)
            % Z=CZ*shape_fcn_Z.*class_fcn_Z
            %
            % default
            % xi_list=linspace(0,1,xi_gird_num(default is 20)+1)
            % psi_list=linspace(0,1,psi_gird_num(default is 20)+1), if symmetry_y psi_list=linspace(0.5,1,psi_gird_num+1);
            % [XI,PSI]=meshgrid(xi_list,psi_list), colume is LaWGS format line
            % value_torl=1e-3
            %
            % xi, x=X
            % psi, y=Y(xi)
            % zeta, z=Z(xi,psi)
            %
            % input:
            % XI, PSI
            % xi_grid_number, psi_grid_number
            % value_torl
            %
            % notice:
            % shape_fcn_Y(XI), shape_fcn_Z(XI,PSI), class_fcn_Z{N1_fcn_Z(XI), N2_fcn_Z(XI)}
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

                % adapt capture y
                if self.symmetry_y
                    low_bou_y=0.5;
                else
                    low_bou_y=0;
                end
                if ~isempty(self.shape_fcn_X)
                    [psi_list,~,~]=girdAdapt1D(self.shape_fcn_X,low_bou_y,1,value_torl,max_level);
                else
                    psi_list=[];
                end
                if length(psi_list) < 21
                    psi_list=linspace(low_bou_y,1,21);
                end

                % adapt capture x
                if ~isempty(self.shape_fcn_Y)
                    [xi_list,~,~]=girdAdapt1D(self.shape_fcn_Y,0,1,value_torl,max_level);
                else
                    xi_list=[];
                end
                if length(xi_list) < 21
                    xi_list=linspace(0,1,21);
                end

                [XI,PSI]=meshgrid(xi_list,psi_list);
            elseif nargin == 3
                % input is XI, PSI or xi_grid_number, psi_grid_number
                if length(varargin{1}) == 1
                    XI=[];
                    xi_grid_numebr=varargin{1};
                else
                    % mean input XI matrix
                    XI=varargin{1};
                    xi_grid_numebr=size(XI,2)-1;
                end

                if length(varargin{2}) == 1
                    PSI=[];
                    psi_grid_numebr=varargin{2};
                else
                    % mean input PSI matrix
                    PSI=varargin{2};
                    psi_grid_numebr=size(PSI,1)-1;
                end

                % calculate local coordinate matrix
                if isempty(XI)
                    xi_list=linspace(0,1,xi_grid_numebr+1);
                    XI=repmat(xi_list,psi_grid_numebr+1,1);
                end
                if isempty(PSI)
                    psi_list=linspace(0,1,psi_grid_numebr+1)';
                    if self.symmetry_y
                        psi_list=psi_list/2+0.5;
                    end
                    PSI=repmat(psi_list,1,xi_grid_numebr+1);
                end
            end

            [X,Y,Z]=self.calPoint(XI,PSI);
        end

        function [X,Y,Z]=calPoint(self,XI,PSI)
            % calculate point on surface
            %

            % calculate origin surface matrix
            X=self.LX*XI;
            if ~isempty(self.shape_fcn_X)
                X=X.*self.shape_fcn_X(PSI);
            end

            if self.symmetry_y
                Y=self.LY*(PSI-0.5);
            else
                Y=self.LY*PSI;
            end
            if ~isempty(self.shape_fcn_Y)
                Y=Y.*self.shape_fcn_Y(XI);
            end

            if isempty(self.shape_fcn_Z) && isempty(self.class_fcn_Z)
                Z=zeros(size(X));
            elseif isempty(self.shape_fcn_Z)
                Z=self.LZ*self.class_fcn_Z(XI,PSI);
            elseif isempty(self.class_fcn_Z)
                Z=self.LZ*self.shape_fcn_Z(XI,PSI);
            else
                Z=self.LZ*self.shape_fcn_Z(XI,PSI).*self.class_fcn_Z(XI,PSI);
            end

            % deform surface
            if ~isempty(self.tran_fcn_X)
                X=X+self.tran_fcn_X(PSI);
            end

            if ~isempty(self.tran_fcn_Y)
                Y=Y+self.tran_fcn_Y(XI);
            end

            if ~isempty(self.tran_fcn_Z)
                Z=Z+self.tran_fcn_Z(XI,PSI);
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

        function [XI,PSI]=calCoordinate(self,X,Y,Z)
            % base on X, Y, Z calculate local coordinate in surface
            %

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

                PSI=Y./self.LY;
                if self.symmetry_y
                    PSI=PSI+0.5;
                end

                % re deform surface
                if ~isempty(self.tran_fcn_X)
                    X=X-self.tran_fcn_X(PSI);
                end

                if self.LX > 0
                    X=max(X,0);
                else
                    X=min(X,0);
                end
                XI=X./self.LX;
                if ~isempty(self.shape_fcn_X)
                    shape=self.shape_fcn_X(PSI);
                    shape(shape==0)=1;
                    XI=XI./shape;
                end

                % judge local coordinate
                XI=min(XI,1);
            elseif self.deform_ID == 2
                if self.LX > 0
                    X=max(X,0);X=min(X,self.LX);
                else
                    X=max(X,self.LX);Y=min(X,0);
                end

                XI=X./self.LX;

                % re deform surface
                if ~isempty(self.tran_fcn_Y)
                    Y=Y-self.tran_fcn_Y(XI);
                end

                PSI=Y./self.LY;
                if ~isempty(self.shape_fcn_Y)
                    shape=self.shape_fcn_Y(XI);
                    shape(shape==0)=1;
                    PSI=PSI./shape;
                end

                % judge local coordinate
                if self.symmetry_y
                    PSI=PSI+0.5;
                end
                PSI=min(PSI,1);
            end

        end

    end

    % discrete surface function
    methods
        function [X,Y,Z]=calPointDiscrete(self,XI,PSI,XI_surf,PSI_surf)
            % calculate point on surface
            %
            [X_surf,Y_surf,Z_surf]=calSurface(self,XI_surf,PSI_surf);

            xi_gird_num=length(XI_surf)-1;
            psi_gird_num=length(PSI_surf)-1;

            point_num=length(XI);
            X=zeros(point_num,1);Y=zeros(point_num,1);Z=zeros(point_num,1);

            for point_idx=1:point_num
                xi=XI(point_idx);
                psi=PSI(point_idx);

                % locate point in which grid
                xi_idx=1; % search start from last one, find out X samll than x
                while ((xi_idx < xi_gird_num) && (XI_surf(xi_idx+1) < xi))
                    xi_idx=xi_idx+1;
                end
                psi_idx=1; % search start from last one, find out X samll than x
                while ((psi_idx < psi_gird_num) && (PSI_surf(psi_idx+1) < psi))
                    psi_idx=psi_idx+1;
                end

                grid_xi=(xi-XI_surf(xi_idx))/(XI_surf(xi_idx+1)-XI_surf(xi_idx));
                grid_psi=(psi-PSI_surf(psi_idx))/(PSI_surf(psi_idx+1)-PSI_surf(psi_idx));
                N1=(1-grid_xi)*(1-grid_psi);
                N2=grid_xi*(1-grid_psi);
                N3=grid_xi*grid_psi;
                N4=(1-grid_xi)*grid_psi;

                idx_1=(psi_gird_num+1)*(xi_idx-1)+psi_idx;
                idx_2=(psi_gird_num+1)*(xi_idx)+psi_idx;
                idx_3=(psi_gird_num+1)*(xi_idx)+psi_idx+1;
                idx_4=(psi_gird_num+1)*(xi_idx-1)+psi_idx+1;

                X(point_idx)=X_surf(idx_1)*N1+X_surf(idx_2)*N2+X_surf(idx_3)*N3+X_surf(idx_4)*N4;
                Y(point_idx)=Y_surf(idx_1)*N1+Y_surf(idx_2)*N2+Y_surf(idx_3)*N3+Y_surf(idx_4)*N4;
                Z(point_idx)=Z_surf(idx_1)*N1+Z_surf(idx_2)*N2+Z_surf(idx_3)*N3+Z_surf(idx_4)*N4;
            end
        end

        function [XI,PSI]=calCoordinateDiscrete(self,X,Y,Z,XI_bou,PSI_bou,check_flag)
            % base on X, Y, Z calculate local coordinate in surface
            % calculate local coordinate of point which surface is discrete
            %
            if nargin < 7
                check_flag=false(1);
            end
            geometry_torlance=1e-3;
            X_origin=X;Y_origin=Y;Z_origin=Z;

            if length(XI_bou) == 1
                XI_bou=linspace(0,1,XI_bou+1);
            end

            if length(PSI_bou) == 1
                PSI_bou=linspace(0,1,PSI_bou+1)';
            end

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

                PSI=Y./self.LY;
                if self.symmetry_y
                    PSI=PSI+0.5;
                end

                % re deform surface
                if ~isempty(self.tran_fcn_X)
                    % how to solve?
                    X=X-calDicrete(self.tran_fcn_X,PSI,PSI_bou);
                end

                if self.LX > 0
                    X=max(X,0);
                else
                    X=min(X,0);
                end
                % base on up bou to re calculate PSI local coordinate
                XI=X./self.LX;
                if ~isempty(self.shape_fcn_X)
                    shape=calDicrete(self.shape_fcn_X,PSI,PSI_bou);
                    shape(shape==0)=1;
                    XI=XI./shape;
                end

                % judge local coordinate
                XI=min(XI,1);

            elseif self.deform_ID == 2
                if self.LX > 0
                    X=max(X,0);X=min(X,self.LX);
                else
                    X=max(X,self.LX);Y=min(X,0);
                end

                XI=X./self.LX;

                % re deform surface
                if ~isempty(self.tran_fcn_Y)
                    Y=Y-calDicrete(self.tran_fcn_Y,XI,XI_bou);
                end

                if self.LY > 0
                    Y=max(Y,0);
                else
                    Y=min(Y,0);
                end
                % base on up bou to re calculate PSI local coordinate
                PSI=Y./self.LY;
                if ~isempty(self.shape_fcn_Y)
                    shape=calDicrete(self.shape_fcn_Y,XI,XI_bou);
                    shape(shape==0)=1;
                    PSI=PSI./shape;
                end

                % judge local coordinate
                if self.symmetry_y
                    PSI=PSI+0.5;
                end
                PSI=min(PSI,1);

            end

            if check_flag
                [X_new,Y_new,Z_new]=self.calPointDiscrete(XI,PSI,XI_bou,PSI_bou);
                boolean=(abs(X_origin-X_new)+abs(Y_origin-Y_new)+abs(Z_origin-Z_new)) > geometry_torlance;
                recal_index=find(boolean);

                % re project point to surface and calculate coordinate
                for idx=1:length(recal_index)
                    obj_fcn=@(x) objFcn(x,X_origin(recal_index(idx)),Y_origin(recal_index(idx)),Z_origin(recal_index(idx)));
                    x=fmincon(obj_fcn,[XI(recal_index(idx)),PSI(recal_index(idx))],[],[],[],[],[0,0],[1,1],[],optimoptions('fmincon','Display','none'));
                    XI(recal_index(idx))=x(1);PSI(recal_index(idx))=x(2);
                end
            end

            function obj=objFcn(x,x_origin,y_origin,z_origin)
                [x_new,y_new,z_new]=self.calPointDiscrete(x(1),x(2),XI_bou,PSI_bou);
                obj=(abs(x_origin-x_new)+abs(y_origin-y_new)+abs(z_origin-z_new));
            end
        end

    end

    % external function
    methods
        function surface_BSpline=getBSpline(self,varargin)
            % convert CST surface into BSpline surface
            %
            % input:
            % XI,PSI
            % xi_grid_number, psi_grid_number
            % value torlance
            %

            if nargin == 1
                [X,Y,Z,XI,PSI]=calSurface(self);
            elseif nargin == 2
                [X,Y,Z,XI,PSI]=calSurface(self,varargin{1});
            elseif nargin == 3
                XI=varargin{1};
                PSI=varargin{2};
                if length(XI) > 1 && XI(1,1) > XI(1,2)
                    XI=fliplr(XI);
                end
                if length(PSI) > 1 && PSI(1,1) > PSI(2,1)
                    PSI=flipud(PSI);
                end
                [X,Y,Z,XI,PSI]=calSurface(self,XI,PSI);
            else
                error('SurfaceCST3D.getBSpline: error input');
            end

            u_list=[0,0,0,XI(1,:),0,0,0];
            v_list=[0,0,0,PSI(:,1)',0,0,0];
            surface_BSpline=SurfaceBSpline([],[],[],X,Y,Z,u_list,v_list);
        end
    end
end

function Y=calDicrete(fcn,XI,XI_bou)
Y_bou=fcn(XI_bou);
xi_number=length(XI);
Y=XI;
for xi_index=1:xi_number
    Y(xi_index)=interpLinear(XI(xi_index),XI_bou,Y_bou);
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

function [x_list,fval_list,node_list]=girdAdapt1D(fcn,low_bou,up_bou,torl,max_level)
% adapt capture function value
% ensure error of linear interplation will less than torl
%

% node_list which is a matrix store all node
% a node is a array, contain level, x, fval, node_low_bou_index, node_up_bou_index, left_index, right_index
% if children node is empty, left_index or right_index will be zero
node_add_num=50; % node_list will be extend only when node_list is not enough
node_list=zeros(node_add_num,7);

% add low_bou data and up_bou data
% root will start from index 3
node_list(1,:)=[0,low_bou,fcn(low_bou),0,0,0,0];
node_list(2,:)=[0,up_bou,fcn(up_bou),0,0,0,0];

x=(low_bou+up_bou)/2;
fval=fcn(x);
level=0;
node_list(3,:)=[level,x,fval,1,2,0,0]; % root

node_num=createNodeTree(3); % create node tree from root
node_list=node_list(1:node_num,:);
[x_list,fval_list]=traversalInorder(3); % from small to large get list

% add boundary info
x_list=[node_list(1,2),x_list,node_list(2,2)];
fval_list=[node_list(1,3),fval_list,node_list(2,3)];

    function node_num=createNodeTree(root_idx)
        %
        %
        stack=root_idx;
        node_num=root_idx;

        while ~isempty(stack)
            % current node information
            node_idx=stack(end);
            node=node_list(node_idx,:);
            stack=stack(1:end-1);

            % check left
            x_left=(node_list(node(4),2)+node(2))/2;
            fval_left=fcn(x_left);
            fval_left_pred=(node_list(node(4),3)+node(3))/2;
            if abs(fval_left-fval_left_pred) > torl && node(1) < max_level
                % create new node
                node_left_index=node_num+1;
                node_left=[node(1)+1,x_left,fval_left,node(4),node_idx,0,0];

                if node_left_index > size(node_list,1)
                    node_list=[node_list;zeros(node_add_num,7)];
                end
                node_list(node_left_index,:)=node_left;
                node_num=node_num+1;

                stack=[stack,node_left_index];
            else
                node_left_index=0;
            end

            % check right
            x_right=(node_list(node(5),2)+node(2))/2;
            fval_right=fcn(x_right);
            fval_right_pred=(node_list(node(5),3)+node(3))/2;
            if abs(fval_right-fval_right_pred) > torl && node(1) < max_level
                % create new node
                node_right_index=node_num+1;
                node_right=[node(1)+1,x_right,fval_right,node_idx,node(5),0,0];

                if node_right_index > size(node_list,1)
                    node_list=[node_list;zeros(node_add_num,7)];
                end
                node_list(node_right_index,:)=node_right;
                node_num=node_num+1;

                stack=[stack,node_right_index];
            else
                node_right_index=0;
            end

            % add children index into current node data
            node_list(node_idx,6:7)=[node_left_index,node_right_index];
        end
    end

    function [x_list,fval_list]=traversalInorder(root_idx)
        stack=[];
        x_list=[];
        fval_list=[];
        node_idx=root_idx;

        while ~isempty(stack) || node_idx ~= 0
            while node_idx ~= 0
                stack=[stack,node_idx];

                % node=node.left;
                node_idx=node_list(node_idx,6);
            end

            node_idx=stack(end);
            stack=stack(1:end-1);
            x_list=[x_list,node_list(node_idx,2)];
            fval_list=[fval_list,node_list(node_idx,3)];

            % node=node.right;
            node_idx=node_list(node_idx,7);
        end
    end
end

