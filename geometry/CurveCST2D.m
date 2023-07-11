classdef CurveCST2D < handle
    % generate curve by CST parameter
    properties
        LX;
        LY;

        % origin parameter
        shape_fcn_Y;
        class_fun_Y;
        symmetry_x;
        
        % deform parameter
        tran_fun_Y=[]; % (XI)

        % rotation parameter
        rotation_matrix=[];

        % translation parameter
        translation=[];
    end

    methods
        function self=CurveCST2D...
                (LX,LY,shape_par_Y,class_par_Y,symmetry_x)
            % generate 2D CST line by LX, LY, shape_par_Y, class_par_Y
            %
            % xi, x=X
            % psi, y=Y(xi)
            %
            % input:
            % LX,LY,shape_par_Y,class_par_Y,symmetry_x
            %
            % notice:
            % shape_fcn_Y(XI), class_fun_Y[N1, N2]
            %
            if nargin < 5 || isempty(symmetry_x)
                self.symmetry_x=false(1);
            else
                self.symmetry_x=symmetry_x;
            end
            self.LX=LX;
            self.LY=LY;

            if isempty(shape_par_Y)
                self.shape_fcn_Y=@(XI) 1;
            elseif isnumeric(shape_par_Y)
                % shape_par_Y only one cross-section parameter
                self.shape_fcn_Y=@(XI) self.defunClass...
                    (XI,shape_par_Y(1),shape_par_Y(2));
            else
                % function handle
                self.shape_fcn_Y=shape_par_Y;
            end

            self.class_fun_Y=@(XI) self.defunClass...
                (XI,class_par_Y(1),class_par_Y(2));

        end
    end

    % common function
    methods
        function C=defunClass(self,XI,N1,N2)
            % default of Z class function
            % default N1, N2 is relation with XI
            %
            NP=self.calNormPar(N1,N2);
            C=XI.^N1.*(1-XI).^N2./NP;
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
        function [X,Y]=calCurve(self,varargin)
            % calculate 2D CST vector by XI
            % X=CX*XI, if symmetry_x X=CX*(XI-0.5)
            % Y=CY*shape_fcn_Y.*class_fun_Y
            %
            % default
            % xi_list=linspace(0,1,xi_gird_num+1) if symmetry_x xi_list=linspace(0.5,1,xi_gird_num+1)
            %
            % xi, x=X
            % psi, y=Y(xi)
            %
            % input:
            % XI
            % xi_grid_number
            %
            % notice:
            % shape_fcn_Y(XI), class_fun_Y[N1, N2]
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
                    xi_grid_numebr=length(XI)-1;
                end
            end

            if isempty(XI)
                if self.symmetry_x
                    XI=linspace(0.5,1,xi_grid_numebr+1);
                else
                    XI=linspace(0,1,xi_grid_numebr+1);
                end
            end

            [X,Y]=self.calPoint(XI);
        end

        function [X,Y]=calPoint(self,XI)
            % calculate point on curve
            %

            % calculate origin surface matrix
            if self.symmetry_x
                X=self.LX.*(XI-0.5);
            else
                X=self.LX.*XI;
            end
            
            Y=self.LY*self.shape_fcn_Y(XI).*self.class_fun_Y(XI);
        end

        function XI=calCoordinate(self,X)
            % base on X, Y calculate local coordinate in surface
            %
            XI=X./self.LX;
        end
    end

end
