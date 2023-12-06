classdef CurveBCST2D < CurveCST2D
    %
    %
    properties
        blunt_vec;
        radius;
        joint_k;

        blunt_fcn;
    end

    methods
        function self=CurveBCST2D(name,C_par,sym,LX,LY)
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
            self=self@CurveCST2D(name,C_par,sym,LX,LY)

            self.blunt_fcn=@(X,Y) defcnBlunt(X,Y);

            function [X,Y]=defcnBlunt(X,Y),end
        end

        function addBlunt(self,blunt_vec,radius,joint_k)
            % calculate blunt parameter of CST curve
            % blunt_vec: -1, blunting to -y; 1, blunting to +y
            % radius: radius of edge blunting
            % joint_k: slope of tangent line
            %
            if isempty(blunt_vec),blunt_vec=1;end
            self.blunt_vec=blunt_vec;
            self.radius=radius;
            self.joint_k=joint_k;

            % calculate cross point
            [LX,LY]=self.shape_fcn(1);
            fsolve_option=optimoptions('fsolve','Display','none');
            if LX-2*radius < 0
                error('CurveBCST2D.addBlunt: blunt radius is too large');
            end

            if isinf(joint_k)
                tang_x=LX-2*radius;
                tang_y=0;
            else
                tang_x=LX-(joint_k/sqrt(1+joint_k*joint_k)+1)*radius;
                tang_y=blunt_vec*radius/sqrt(1+joint_k*joint_k);
            end
            solve_fcn=@(U) self.shape_fcn(U).*self.class_fcn(U)-tang_x;
            cross_u=fsolve(solve_fcn,tang_x/LX,fsolve_option);
            [SX,SY]=self.shape_fcn(cross_u);
            [CX,CY]=self.class_fcn(cross_u);
            cross_x=SX.*CX;cross_y=SY.*CY;
            dtang_y=tang_y-cross_y;
            [dX_dU,dY_dU]=self.calGradient(cross_u);
            grad_tang=dY_dU/dX_dU;

            % if ~isinf(joint_k)
            %     k=joint_k;
            %     b=-joint_k*tang_x+tang_y;
            %     solve_fcn=@(U) solveFcn(U,k,b,self);
            %     [cross_u,fval,exitflag,output]=fsolve(solve_fcn,tang_x/LX,fsolve_option);
            %     if cross_u > 1 || cross_u < 0 || ~isreal(cross_u) || abs(fval) > 0.1
            %         error('CurveBCST2D.addBlunt: blunt radius is too large');
            %     end
            %     [SX,SY]=self.shape_fcn(cross_u);
            %     [CX,CY]=self.class_fcn(cross_u);
            %     cross_x=SX.*CX;cross_y=SY.*CY;
            % end
            % 
            % self.blunt_fcn=@(X,Y) defcnBlunt(X,Y,blunt_vec,radius,tang_x,cross_x,LX,dtang_y);

            % calculate cross point of x axis
            djoint_k=joint_k-grad_tang;
            if djoint_k == 0
                error('CurveBCST2D.addBlunt: blunt radius is too large');
            end
            cross_x=(djoint_k*tang_x-dtang_y)/djoint_k;
            if cross_x < 0 || cross_x > LX
                error('CurveBCST2D.addBlunt: blunt radius is too large');
            end

            self.blunt_fcn=@(X,Y) defcnBlunt(X,Y,blunt_vec,radius,tang_x,cross_x,LX,dtang_y);

            function obj=solveFcn(u,k,b,self)
                [SX__,SY__]=self.shape_fcn(u);
                [CX__,CY__]=self.class_fcn(u);
                obj=SX__.*CX__*k+b-SY__.*CY__;
            end

            function [X,Y]=defcnBlunt(X,Y,blunt_vec,radius,tang_x,cross_x,LX,dtang_y)
                % blunt bias of curve
                %
                radius_sq=radius*radius;
                Bool_radius=X >= tang_x;
                BY(Bool_radius)=blunt_vec*sqrt(radius_sq-(X(Bool_radius)-(LX-radius)).^2+eps)-Y(Bool_radius);
                A=tang_x-2*cross_x+0;B=-2*tang_x+2*cross_x;C=tang_x-X(~Bool_radius);
                Delta=sqrt(B.*B-4*A.*C);
                U_nega=(-B-Delta)/2./A;U_posi=(-B+Delta)/2./A;
                U=U_nega;U(U < 0 | U > 1)=U_posi(U < 0 | U > 1);
                BY(~Bool_radius)=blunt_vec*(1-U).^2*dtang_y;
                if blunt_vec == 1,BY=max(BY,0);
                else,BY=min(BY,0);end
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
            [X,Y]=self.blunt_fcn(X,Y);
            [X,Y]=self.axisLocalToGlobal(X,Y,U);
            Z=[];
        end
    end
end