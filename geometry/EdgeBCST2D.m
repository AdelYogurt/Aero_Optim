classdef EdgeBCST2D < EdgeCST2D
    %
    %
    properties
        blunt_vec;
        radius;
        blunt_fcn;
    end

    methods
        function self=EdgeBCST2D(name,C_par,sym,LX,LY)
            % generate 2D CST line with blunt edge
            %
            % u, x=LX*S(u)
            % v, y=LY*C(u)*S(u)
            %
            % input:
            % name, C_par, sym, radius, LX, LY
            %
            % notice:
            % if input N1, N2 == 0, or C_par is empty,...
            % class_fcn will equal to 1
            % class_fcn(U)
            %
            if nargin < 5
                LY=[];
                if nargin < 4
                    LX=[];
                    if nargin < 3
                        sym=[];
                        if nargin < 2
                            C_par=[];
                        end
                    end
                end
            end
            self=self@EdgeCST2D(name,C_par,sym,LX,LY)

            self.blunt_fcn=[];
        end

        function addBlunt(self,blunt_vec,radius)
            % calculate blunt parameter of CST curve
            % blunt_vec: -1, blunting to -y; 1, blunting to +y
            % radius: radius of edge blunting
            % joint_k: slope of tangent line
            %
            if isempty(blunt_vec),blunt_vec=1;end
            self.blunt_vec=blunt_vec;
            self.radius=radius;

            % calculate cross point
            [LX,LY]=self.shape_fcn(1);
            fminbnd_option=optimset('Display','none');
            if LX-2*radius < 0
                error('EdgeBCST2D.addBlunt: blunt radius is too large');
            end

            % find out grad equal u in curve
            u_min=fminbnd(@(u) (self.calPoint(u)-(LX-2*radius))^2,0,1,fminbnd_option);
            solve_fcn=@(u) solveFcn(self,u,blunt_vec,radius,LX,u_min,1);
            equal_U=fminbnd(solve_fcn,u_min,1,fminbnd_option);
            [SX,SY]=self.shape_fcn(equal_U);
            [CX,CY]=self.class_fcn(equal_U);
            X=CX.*SX;Y=CY.*SY;
            RX=X-(LX-radius);DX=X;
            DY=blunt_vec*sqrt(radius*radius-RX*RX+eps)-Y;
            
            if blunt_vec*DY > 0
                self.blunt_fcn=@(X,Y) defcnBlunt(X,Y,blunt_vec,radius,DX,DY,LX);
            else
                self.blunt_fcn=[];
            end

            function obj=solveFcn(self,u,blunt_vec,radius,LX,u_min,u_max)
                u=max(u,u_min);u=min(u,u_max);
                [sx,sy]=self.shape_fcn(u);
                [cx,cy]=self.class_fcn(u);
                x=cx.*sx;
                rx=x-(LX-radius);
                ry=blunt_vec*sqrt(radius*radius-rx*rx+eps);
                grad_radius=-rx/ry;
                [dx_du,dy_du]=self.calGradient(u);
                grad_curve=dy_du./dx_du;
                obj=(grad_radius-grad_curve)^2;
            end

            function [X,Y]=defcnBlunt(X,Y,blunt_vec,radius,DX,DY,LX)
                % blunt bias of curve
                %
                radius_sq=radius*radius;
                Bool_radius=X(:) >= DX;BY=zeros(size(Y));
                BY(Bool_radius)=blunt_vec*sqrt(radius_sq-(X(Bool_radius)-(LX-radius)).^2+eps)-Y(Bool_radius);
                BY(~Bool_radius)=DY*sin(pi/2/DX*X(~Bool_radius));
                Y=Y+BY;
            end
        end
    end

    % calculate point
    methods
        function [X,Y,Z]=calPoint(self,U)
            % calculate point on curve
            %
            
            % calculate origin curve
            [SX,SY]=self.shape_fcn(U);
            [CX,CY]=self.class_fcn(U);
            X=CX.*SX;Y=CY.*SY;
            if ~isempty(self.blunt_fcn)
                [X,Y]=self.blunt_fcn(X,Y);
            end
            [X,Y]=self.axisLocalToGlobal(X,Y,U);
            Z=zeros(size(X));
        end
    end
end
