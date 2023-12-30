classdef EdgeBCST2D < EdgeCST2D
    % blunt edge version of CST
    %
    properties
        blunt_vec; % blunt direction
        radius; % blunt radius
        blunt_fcn; % [X,Y]=blunt_fcn(X,Y), blunt point deform function
    end

    methods
        function self=EdgeBCST2D(name,C_par,sym,LX,LY)
            % generate 2D BCST line with blunt edge
            % X=LX*shape_fcn(u)
            % Y=LY*class_fcn(u)*shape_fcn(u)
            %
            % input:
            % name:
            % C_par:
            % sym;
            % LX:
            % LY:
            %
            % output:
            % EdgeCST2D
            %
            % notice:
            % if input N1, N2 == 0, or C_par is empty,...
            % class_fcn will equal to 1
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
            self.blunt_fcn=[];

            % calculate cross point
            point=self.shape_fcn(1);LX=point(1);
            fminbnd_option=optimset('Display','none','TolX',100*eps);
            if LX-2*radius < 0
                error('EdgeBCST2D.addBlunt: blunt radius is too large');
            end

            % find out grad equal u in curve
            function obj=solveCross(self,u,LX,radius)
                p=self.calCST(u);
                obj=(p(1)-(LX-2*radius))^2;
            end
            cross_fcn=@(u) solveCross(self,u,LX,radius);
            u_min=fminbnd(cross_fcn,0,1,fminbnd_option)+1e-6;
            
            function obj=solveGrad(self,u,blunt_vec,radius,LX)
                p=self.calCST(u);
                rx=p(1)-(LX-radius);
                ry=blunt_vec*sqrt(radius*radius-rx*rx+eps);
                grad_radius=-rx/ry;
                dp_du=self.calGradient(u,sqrt(eps));
                grad_curve=dp_du(2)/dp_du(1);
                obj=(grad_radius-grad_curve)^2;
            end
            grad_fcn=@(u) solveGrad(self,u,blunt_vec,radius,LX);
            equal_u=fminbnd(grad_fcn,u_min,1,fminbnd_option);

            % calculate class
            point=self.calCST(equal_u);
            RX=point(1)-(LX-radius);DX=point(1);
            DY=blunt_vec*sqrt(radius*radius-RX*RX+eps)-point(2);
            
            if blunt_vec*DY > 0
                self.blunt_fcn=@(Point) defcnBlunt(Point,blunt_vec,radius,DX,DY,LX);
            else
                self.blunt_fcn=[];
            end

            function Point=defcnBlunt(Point,blunt_vec,radius,DX,DY,LX)
                % blunt bias of curve
                %
                X=Point(:,1);Y=Point(:,2);
                radius_sq=radius*radius;
                Bool_radius=X(:) >= DX;BY=zeros(size(Y));
                BY(Bool_radius)=blunt_vec*sqrt(radius_sq-(X(Bool_radius)-(LX-radius)).^2+eps)-Y(Bool_radius);
                BY(~Bool_radius)=DY*sin(pi/2/DX*X(~Bool_radius));
                Y=Y+BY;
                Point=[X,Y];
            end
        end
    end

    % calculate point
    methods
        function Point=calPoint(self,U)
            % calculate point on curve
            %
            Point=calCST(self,U);
            if ~isempty(self.blunt_fcn)
                Point=self.blunt_fcn(Point);
            end
            Point=self.axisLocalToGlobal(Point,U);
        end
    end
end
