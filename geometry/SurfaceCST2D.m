classdef SurfaceCST2D
    properties
        LX;
        LY;
        class_fun_X;
        shape_fcn_X;
        deform_ID; % 1 or 2, 1 is x, 2 is y
        class_fun_Y;
        shape_fcn_Y;
        symmetry_x;
        symmetry_y;
    end
    methods
        function self=SurfaceCST2D(LX,LY,shape_par_X,class_par_X,shape_par_Y,class_par_Y,symmetry_x,symmetry_y)
            % generate 2D surface by LX, LY, shape_par_X,class_par_X,shape_par_Y,class_par_Y
            %
            % xi, x=X(psi)
            % psi, y=Y(xi)
            %
            % input:
            % LX,LY,[],[],shape_par_Y,class_par_Y,symmetry_x,symmetry_y
            % LX,LY,shape_par_X,class_par_X,[],[],symmetry_x,symmetry_y
            %
            % notice:
            % shape_fcn_X(XI), class_fun_X[N1, N2], shape_fcn_Y(XI), class_fun_Y[N1, N2]
            %
            if nargin < 7 || isempty(symmetry_x)
                self.symmetry_x=false(1);
            else
                self.symmetry_x=symmetry_x;
            end
            if nargin < 8 || isempty(symmetry_y)
                self.symmetry_y=false(1);
            else
                self.symmetry_y=symmetry_y;
            end
            self.LX=LX;
            self.LY=LY;

            if isempty(shape_par_X) && isempty(class_par_X)
                % deform direction
                self.deform_ID=2;
            elseif isempty(shape_par_Y) && isempty(class_par_Y)
                % deform direction
                self.deform_ID=1;
            else
                error('CST2DSurface: donot support both X and Y deform')
            end

            % generate shape function and class function
            if isempty(shape_par_X)
                self.shape_fcn_X=@(PSI) 1;
            elseif isnumeric(shape_par_X)
                self.shape_fcn_X=@(PSI) self.defunClass...
                    (PSI,shape_par_X(1),shape_par_X(2));
            else
                % function handle
                self.shape_fcn_X=shape_par_X;
            end

            if isempty(class_par_X)
                self.class_fun_X=@(PSI) 1;
            elseif isnumeric(class_par_X)
                self.class_fun_X=@(PSI) self.defunClass...
                    (PSI,class_par_X(1),class_par_X(2));
            else
                % function handle
                self.class_fun_X=class_par_X;
            end

            if isempty(shape_par_Y)
                self.shape_fcn_Y=@(XI) 1;
            elseif isnumeric(shape_par_Y)
                self.shape_fcn_Y=@(XI) self.defunClass...
                    (XI,shape_par_Y(1),shape_par_Y(2));
            else
                % function handle
                self.shape_fcn_Y=shape_par_Y;
            end

            if isempty(class_par_Y)
                self.class_fun_Y=@(XI) 1;
            elseif isnumeric(class_par_Y)
                self.class_fun_Y=@(XI) self.defunClass...
                    (XI,class_par_Y(1),class_par_Y(2));
            else
                % function handle
                self.class_fun_Y=class_par_Y;
            end
        end

        function [X,Y]=calSurface(self,varargin)
            % generate 2D surface by XI,PSI
            % X=CX*shape_fcn_X.*XI, if symmetry_x X=CX*shape_fcn_X.*(XI-0.5)
            % Y=CY*shape_fcn_Y.*PSI, if symmetry_y Y=CY*shape_fcn_Y.*(PSI-0.5)
            %
            % default
            % xi_list=linspace(0,1,xi_gird_num+1) if symmetry_x xi_list=linspace(0.5,1,xi_gird_num+1)
            % psi_list=linspace(0,1,psi_gird_num+1) if symmetry_y psi_list=linspace(0.5,1,psi_gird_num+1)
            %
            % xi, x=X(psi)
            % psi, y=Y(xi)
            %
            % input:
            % LX,LY,shape_par_X,class_par_X,shape_par_Y,class_par_Y,symmetry_x,symmetry_y,XI,PSI
            % LX,LY,shape_par_X,class_par_X,shape_par_Y,class_par_Y,symmetry_x,symmetry_y,xi_grid_number,psi_gird_num
            %
            % notice:
            % shape_fcn_X(XI), class_fun_X[N1, N2], shape_fcn_Y(XI), class_fun_Y[N1, N2]
            %
            % output:
            % X,Y
            %

            % xi
            if nargin < 2
                % mean input xi_grid_numebr
                XI=[];
                xi_grid_numebr=20;
            else
                if length(varargin{1}) == 1
                    XI=[];
                    xi_grid_numebr=varargin{1};
                else
                    % mean input XI matrix
                    XI=varargin{1};
                    xi_grid_numebr=size(XI,2)-1;
                end
            end

            % psi
            if nargin < 3
                % mean input xi_grid_numebr
                PSI=[];
                psi_grid_numebr=20;
            else
                if length(varargin{2}) == 1
                    PSI=[];
                    psi_grid_numebr=varargin{2};
                else
                    % mean input PSI matrix
                    PSI=varargin{2};
                    psi_grid_numebr=size(PSI,1)-1;
                end
            end

            % calculate local coordinate matrix
            if isempty(XI)
                if self.symmetry_x
                    XI=repmat(linspace(0.5,1,xi_grid_numebr+1),psi_grid_numebr+1,1);
                else
                    XI=repmat(linspace(0,1,xi_grid_numebr+1),psi_grid_numebr+1,1);
                end
            end
            if isempty(PSI)
                if self.symmetry_y
                    PSI=repmat(linspace(0.5,1,psi_grid_numebr+1)',1,xi_grid_numebr+1);
                else
                    PSI=repmat(linspace(0,1,psi_grid_numebr+1)',1,xi_grid_numebr+1);
                end
            end

            % calculate coordination
            if self.symmetry_x
                X=self.LX*self.shape_fcn_X(PSI).*self.class_fun_X(PSI).*(XI-0.5);
            else
                X=self.LX*self.shape_fcn_X(PSI).*self.class_fun_X(PSI).*XI;
            end
            if self.symmetry_y
                Y=self.LY*self.shape_fcn_Y(XI).*self.class_fun_Y(XI).*(PSI-0.5);
            else
                Y=self.LY*self.shape_fcn_Y(XI).*self.class_fun_Y(XI).*PSI;
            end

        end

        function [XI,PSI]=calCoordinate(self,X,Y)
            % base on X, Y calculate local coordinate in surface
            %
            if self.deform_ID == 1
                if self.symmetry_y
                    PSI=Y./self.LY+0.5;
                else
                    PSI=Y./self.LY;
                end

                deform_matrix=self.shape_fcn_X(PSI).*self.class_fun_X(PSI);
                deform_matrix(deform_matrix == 0)=1;

                XI=X./self.LX./deform_matrix;
                if self.symmetry_x
                    XI=XI+0.5;
                end
            elseif self.deform_ID == 2
                if self.symmetry_x
                    XI=X./self.LX+0.5;
                else
                    XI=X./self.LX;
                end

                deform_matrix=self.shape_fcn_Y(XI).*self.class_fun_Y(XI);
                deform_matrix(deform_matrix == 0)=1;

                PSI=Y./self.LY./deform_matrix;
                if self.symmetry_y
                    PSI=PSI+0.5;
                end
            end
        end

        function C=defunClass(self,X,N1,N2)
            % default class function
            %
            NP=self.calNormPar(N1,N2);
            C=X.^N1.*(1-X).^N2/NP;
        end

        function nomlz_par=calNormPar(self,N1,N2)
            % calculate normailize class function parameter by N1, N2
            %
            nomlz_par=(N1./(N1+N2)).^N1.*(N2./(N1+N2)).^N2;
            nomlz_par((N1 == 0) & (N2 == 0))=1;
        end
    end
end