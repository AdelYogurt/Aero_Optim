classdef CurveCST < handle & matlab.mixin.Copyable
    % Class-Shape Transformation (CST) Curve
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
        sym; % if true, U_class will start from 0.5 to 1

        C_par; % origin parameter, default is [0.0, 0.0]
        spline; % 2D BSpline curve

        class_fcn; % (U), output are CX, CY
        shape_fcn; % (U), default are LX, LY, output are SX, SY
    end

    properties % basical curve parameter
        deform_fcn=[]; % deform parameter: deform_fcn(U), only y direction
        rotate_matrix=[]; % rotate parameter
        translate=[]; % translate parameter
    end

    properties % fit point set parameter
        fit_data;
    end

    methods % define curve
        function self=CurveCST(C_par,sym,LX,LY)
            % generate 2D CST curve by LX, LY, C_par
            % X=LX*class_fcn(u)*shape_fcn(u)
            % Y=LY*class_fcn(u)*shape_fcn(u)
            %
            % input:
            % C_par (array): N1 and N2 array([N1, N2]), default is [0.0, 0.0]
            % sym (boolean): whether using half of U_x while calculate class_fcn
            % LX (double): total length of curve
            % LY (double): total height of curve
            %
            % output:
            % CurveCST
            %
            % notice:
            % if C_par is [0.0, 0.0] or empty, class_fcn will equal to 1.0
            %
            if nargin < 4
                LY=[];
                if nargin < 3
                    LX=[];
                    if nargin < 2
                        sym=[];
                        if nargin < 1
                            C_par=[];
                        end
                    end
                end
            end

            if isempty(sym),self.sym=false;
            else,self.sym=sym; end
            if isempty(C_par),C_par=[0,0];end
            if isempty(LX),LX=1;end
            if isempty(LY),LY=1;end

            % class function
            if isnumeric(C_par)
                self.C_par=C_par;
                self.class_fcn=@(U) [ones(size(U)),baseFcnClass(U,C_par(1),C_par(2))];
            else
                self.C_par=[];
                self.class_fcn=C_par;
            end

            % shape function
            self.spline=Curve([0,LY;LX,LY]);
            self.shape_fcn=@(U) self.spline.calPoint(U);
        end

        function Points=calPoint(self,U_x)
            % calculate point on curve
            %
            U_x=U_x(:);U_class=U_x;
            if self.sym,U_class=(U_class/2)+0.5;end

            % calculate class
            CPoint=self.class_fcn(U_class);

            % calculate shape
            SPoint=self.shape_fcn(U_x);

            Points=SPoint.*CPoint;
            Points=self.axisLocalToGlobal(Points,U_x);
        end
        
        function [Points,dPoints_dU]=calGradient(self,U_x,step)
            % use differ to calculate gradient
            %
            if nargin < 3 || isempty(step),step=100*eps;end

            dim=numel(self.calPoint(0));
            U_x=U_x(:);

            Points=self.calPoint(U_x);

            Points_UF=self.calPoint(min(U_x+step,1));
            Points_UB=self.calPoint(max(U_x-step,0));
            Bool_U=(U_x+step) > 1;
            Bool_D=(U_x-step) < 0;
            Bool(:,1)=Bool_U;Bool(:,2)=Bool_D;
            Bool_C=~any(Bool,2);
            dPoints_dU=zeros(size(U_x,1),dim); % allocate memory
            dPoints_dU(Bool_C,:)=(Points_UF(Bool_C,:)-Points_UB(Bool_C,:))/2/step;
            dPoints_dU(Bool_U,:)=(Points(Bool_U,:)-Points_UB(Bool_U,:))/step;
            dPoints_dU(Bool_D,:)=(Points_UF(Bool_D,:)-Points(Bool_D,:))/step;
            dPoints_dU=real(dPoints_dU);
        end

        function [k1,k2,x1,x2]=calTangTorl(self,torl)
            % calculate differ tangent by control discrete error in torl
            %
            if ~isempty(self.C_par)
                N1=self.C_par(1);
                N2=self.C_par(2);
                nomlz_par=(N1./(N1+N2)).^N1.*(N2./(N1+N2)).^N2;
                K1=self.shape_fcn(0);KY1=K1(2)/nomlz_par;
                K2=self.shape_fcn(1);KY2=K2(2)/nomlz_par;
                KX=abs(K2(1)-K1(1));
                if self.sym,KX=KX*2;end

                if N1 == 1
                    k1=N1*KY1/KX;
                    x1=0;
                else
                    k1=N1*KY1/KX*(torl/abs(KY1)/abs(1-N1))^((N1-1)/N1);
                    x1=KX*(k1*KX/N1/abs(KY1))^(1/(N1-1));
                end

                if N2 == 1
                    k2=N2*KY2/KX;
                    x2=0;
                else
                    k2=N2*KY2/KX*(torl/abs(KY2)/abs(1-N2))^((N2-1)/N2);
                    x2=KX*(k1*KX/N2/abs(KY2))^(1/(N2-1));
                end
            else
                k1=[];k2=[];
                x1=[];x2=[];
            end
        end

        function crv=convertSpline(self)
            % convert CST Curve to Spline Curve
            %
            % reference: [1] Marshall D D. Creating Exact Bezier
            % Representations of CST Shapes [C]. 21st AIAA Computational
            % Fluid Dynamics Conference, San Diego, CA.
            %
            if ~isempty(self.C_par) && ~isempty(self.spline)
                if self.C_par(1) == 0.5 && self.C_par(2) == 1.0
                    if ~isempty(self.deform_fcn), KY=self.deform_fcn(1);
                    else, KY=0;end

                    SP=self.shape_fcn(1);
                    LX=SP(1);PY=SP(2);
                    if ~isempty(self.spline.Poles), PY=self.spline.Poles(:,2);end

                    [Poles,Degree,~]=getBezier(self.C_par(1),self.C_par(2),LX,PY,KY);
                    crv=Curve(Poles,Degree);
                else
                    warning('CurveCST.convertSpline: can not convert to spline curve exactly');
                    Nodes=self.calPoint(linspace(0,1,21));Degree=3;
                    crv=GeomApp.interpPointToCurve(Nodes,Degree);
                end
            end

            function [Poles,Degree,B_mat]=getBezier(N1,N2,LX,PY,KY)
                % convert CST into Bezier excatly
                %
                c_nomlz=(N1./(N1+N2)).^N1.*(N2./(N1+N2)).^N2;
                c_nomlz((N1 == 0) & (N2 == 0))=1;
                PY=PY(:)';
                Degree=length(PY)-1;
                fac_list_temp=[1,cumprod(1:2*Degree+3)];

                if Degree == 0
                    B_mat=zeros(2,2*Degree+3+1);

                    % X
                    B_mat(1,2+1)=LX;
                    B_mat(2,1+1)=PY/c_nomlz;
                    B_mat(2,2+1)=KY;
                    B_mat(2,3+1)=-PY/c_nomlz;
                else
                    % calculate a list
                    fac_list=fac_list_temp(1:Degree+1);
                    ir_mat=zeros(Degree+1);
                    sign_mat=ones(Degree+1);
                    i=0:Degree;ni_list=fac_list(Degree+1)./fac_list(Degree-i+1)./fac_list(i+1);
                    for i=0:Degree
                        r=0:i;
                        ir_mat(1:(i+1),i+1)=fac_list(i+1)./fac_list(i-r+1)./fac_list(r+1);
                        sign_mat(1:(i+1),i+1)=(-1).^(i-r);
                    end
                    a_list=sum(ni_list.*ir_mat.*sign_mat.*PY',1);

                    % calculate B matrix
                    B_mat=zeros(2,2*Degree+3+1);

                    % X
                    B_mat(1,2+1)=LX;

                    % Y
                    B_mat(2,1+1)=a_list(1)/c_nomlz;
                    B_mat(2,2+1)=KY;
                    i=1:Degree;
                    B_mat(2,i*2+1+1)=(a_list(i+1)-a_list(i))/c_nomlz;
                    B_mat(2,2*Degree+3+1)=-a_list(Degree+1)/c_nomlz;
                end
                % calculate new Poles
                fac_list=fac_list_temp;
                Degree=2*Degree+3;
                i=0:Degree;mi_list=fac_list(Degree+1)./fac_list(Degree-i+1)./fac_list(i+1);
                miij_mat=zeros(Degree+1);
                for i=0:Degree
                    j=0:i;
                    miij_mat(i+1,1:(i+1))=fac_list(Degree-j+1)./fac_list((Degree-j)-(i-j)+1)./fac_list((i-j)+1);
                end

                Poles=zeros(size(B_mat))';
                Poles(:,1)=sum(miij_mat.*B_mat(1,:),2)./mi_list';
                Poles(:,2)=sum(miij_mat.*B_mat(2,:),2)./mi_list';
            end
        end

        function ln_hdl=displayGeom(self,axe_hdl,crv_option,u_param)
            % display curve on axes
            %
            if nargin < 4
                u_param=[];
                if nargin < 3
                    crv_option=[];
                    if nargin < 2
                        axe_hdl=[];
                    end
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end
                
            if isempty(crv_option), crv_option=struct();end

            % calculate Points on curve
            if isempty(u_param),u_param=101;end
            if length(u_param) == 1,u_param=linspace(0,1,u_param);end
            Points=self.calPoint(u_param);

            % draw Points on axe_hdl
            ln_hdl=line(axe_hdl,Points(:,1),Points(:,2),crv_option);
            xlabel('x');
            ylabel('y');
        end

        function ln_hdl=displayPoles(self,axe_hdl,pole_option)
            % draw curve on figure handle
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
                pole_option=struct('Marker','s','LineStyle','--','Color','r');
            end

            % draw Poles on axe_hdl
            if ~isempty(self.spline)
                % load properties
                Poles=self.spline.Poles;
                pole_num=size(Poles,1);
                Degree=self.spline.Degree;
                u_list=self.spline.u_list;

                % calculate U_pole
                U_pole=interp1(linspace(0,1,pole_num-Degree+1),u_list(Degree+1:pole_num+1),linspace(0,1,pole_num))';
                if self.sym;U_pole=U_pole/2+0.5;end
                Poles=self.axisLocalToGlobal(Poles.*self.class_fcn(U_pole),U_pole);

                ln_hdl=line(axe_hdl,Poles(:,1),Poles(:,2),pole_option);
                xlabel('x');
                ylabel('y');
            end
        end
    end

    methods % add spline shape function
        function addSpline(self,spline)
            % add spline as shape function
            %
            if size(spline.Poles,2) > 2
                warning('CurveCST.addSpline: dimension of spline is large than 2, spline will be project to 2D')
                spline.Poles=spline.Poles(:,1:2);
            elseif size(spline.Poles,2) == 1
                spline.Poles=[linspace(0,1,size(spline.Poles,1))',spline.Poles];
            end
            if  min(spline.Knots) ~= 0.0 || max(spline.Knots) ~= 1.0
                warning('SurfaceCST.addSpline: Knots of spline will be rescale to [0.0, 1.0]')
                spline.Knots=(spline.Knots-min(spline.Knots))/(max(spline.Knots)-min(spline.Knots));
            end
            self.spline=spline;

            % add shape function
            self.shape_fcn=@(U_x) self.spline.calPoint(U_x);
        end

        function fitSpline(self,Nodes,Degree,pole_num,U_node)
            % fit Spline base on CST class function
            %
            if nargin < 5
                U_node=[];
                if nargin < 4
                    pole_num = [];
                    if nargin < 3
                        Degree=[];
                    end
                end
            end

            node_num=size(Nodes,1);
            if size(Nodes,2) == 1
                Nodes=[linspace(0,1,node_num)',Nodes];
            end
            if isempty(pole_num),pole_num=node_num;end
            if pole_num > node_num
                error('CurveCST.fitSpline: pole_num can not more than node_num')
            end

            % default value
            if isempty(U_node)
                U_node=vecnorm(Nodes(2:end,:)-Nodes(1:end-1,:),2,2);
                U_node=[0;cumsum(U_node)];
            end
            U_node=(U_node(:)-min(U_node))/(max(U_node)-min(U_node));
            U_node_knots=U_node;
            U_node_knots=interp1(linspace(0,1,length(U_node_knots)),U_node_knots,linspace(0,1,pole_num));

            Mults=[Degree+1,ones(1,pole_num-Degree-1),Degree+1];
            Knots=linspace(0,1,node_num-Degree+1);
            for j=2:node_num-Degree
                Knots(j)=mean(U_node_knots(j:j+Degree-1));
            end
            Knots=interp1(linspace(0,1,node_num-Degree+1),Knots,linspace(0,1,pole_num-Degree+1));
            u_list=baseKnotVec(Mults,Knots);

            % translate node to local coordinate
            Nodes=self.axisGlobalToLocal(Nodes,U_node);

            % base on node point list inverse calculate control point list
            fit_matrix=zeros(node_num,pole_num);
            [N_list,idx_srt,idx_end]=baseFcnN(U_node,Degree,u_list);
            for deg_idx=1:Degree+1
                idx=sub2ind([node_num,pole_num],(1:(node_num))',idx_srt+(deg_idx-1));
                fit_matrix(idx)=N_list(:,deg_idx);
            end

            % add class coefficient
            U_class=U_node;
            if self.sym,U_class=(U_node/2)+0.5;end
            CPoints=self.class_fcn(U_class);
            matrix_class=fit_matrix.*CPoints(:,2);

            % reverse calculate control point
            Poles=[fit_matrix\Nodes(:,1),matrix_class\Nodes(:,2)];
            Weights=ones(pole_num,1);

            % add spline and shape function
            self.spline=Curve(Poles,Degree,Mults,Knots,Weights);
            self.shape_fcn=@(U_x) self.spline.calPoint(U_x);

            % fit data
            self.fit_data.Nodes=Nodes;
            self.fit_data.U_node=U_node;
            self.fit_data.fit_matrix=fit_matrix;
        end

        function fit_err=optimClass(self,optim_option)
            % optimization coefficient of class function
            %
            if nargin < 2
                optim_option=[];
            end

            if isempty(optim_option)
                optim_option=optimoptions('fminunc','Display','none');
            end

            if isempty(self.C_par)
                error('CurveCST.optimClass: input C_par is handle, cannot optimize coefficient of class function');
            end

            C_par_low_bou=[0,0];
            C_par_up_bou=[1e3,1e3];
            C_par_init=self.C_par;
            obj_fit=@(C_par) self.fitError(C_par,C_par_low_bou,C_par_up_bou);
            [C_par_optim,fit_err]=fminunc(obj_fit,C_par_init,optim_option);

            C_par_optim=max(C_par_optim,C_par_low_bou);
            C_par_optim=min(C_par_optim,C_par_up_bou);
            self.C_par=C_par_optim;
            self.class_fcn=@(U) [ones(size(U)),baseFcnClass(U,C_par_optim(1),C_par_optim(2))];
        end

        function RMSE=fitError(self,C_par,C_par_low_bou,C_par_up_bou)
            % fit error of C_par
            %
            if nargin < 2
                C_par=self.C_par;
            else
                C_par=max(C_par,C_par_low_bou);
                C_par=min(C_par,C_par_up_bou);
            end

            % updata C_par
            self.class_fcn=@(U) [ones(size(U)),baseFcnClass(U,C_par(1),C_par(2))];

            % load fit data
            U_node=self.fit_data.U_node;
            Nodes=self.fit_data.Nodes;
            fit_matrix=self.fit_data.fit_matrix;

            % add class coefficient
            U_class=U_node;
            if self.sym,U_class=(U_node/2)+0.5;end
            Point=self.class_fcn(U_class);
            matrix_class=fit_matrix.*Point(:,2);

            % reverse calculate control point
            Poles=[fit_matrix\Nodes(:,1),matrix_class\Nodes(:,2)];
            self.spline.Poles=Poles;

            Point_cal=self.calPoint(U_node);
            RMSE=sqrt(mean((Nodes-Point_cal).^2,"all"));
        end
    end

    methods % deform, rotate, translate
        function self=addDeform(self,deform_fcn_Y)
            % base on local coordinate deform curve
            %
            % input:
            % deform_fcn_Y(U)
            %
            self.deform_fcn=deform_fcn_Y;
        end

        function self=addRotate(self,ang)
            % base on angle to rotate surface
            % rotate is anti-clock
            %
            % input:
            % ang(deg)
            %
            cR=cos(ang/180*pi);sR=sin(ang/180*pi);
            matrix=[
                cR,-sR
                sR,cR];

            self.rotate_matrix=matrix;
        end

        function self=addTranslate(self,tran_x,tran_y)
            % base on angle to rotate surface
            %
            self.translate=[tran_x,tran_y];
        end

        function Point=axisLocalToGlobal(self,Point,U)
            % deform curve
            if ~isempty(self.deform_fcn)
                Point(:,2)=Point(:,2)+self.deform_fcn(U);
            end

            % rotate curve
            if ~isempty(self.rotate_matrix)
                matrix=self.rotate_matrix;
                X_old=Point(:,1);Y_old=Point(:,2);
                Point(:,1)=matrix(1,1)*X_old+matrix(1,2)*Y_old;
                Point(:,2)=matrix(2,1)*X_old+matrix(2,2)*Y_old;
            end

            % translate curve
            if ~isempty(self.translate)
                Point=Point+self.translate;
            end
        end

        function Point=axisGlobalToLocal(self,Point,U)
            % re-translate curve
            if ~isempty(self.translate)
                Point=Point-self.translate;
            end

            % re-rotate curve
            if ~isempty(self.rotate_matrix)
                matrix=self.rotate_matrix';
                X_old=Point(:,1);Y_old=Point(:,2);
                Point(:,1)=matrix(1,1)*X_old+matrix(1,2)*Y_old;
                Point(:,2)=matrix(2,1)*X_old+matrix(2,2)*Y_old;
            end

            % re-deform curve
            if ~isempty(self.deform_fcn)
                Point(:,2)=Point(:,2)-self.deform_fcn(U);
            end
        end

    end
    
    methods % calculate coord
        function U=calCoordinate(self,Points,geom_torl)
            % base on X, Y, Z calculate local coordinate in curve
            %
            if nargin < 3, geom_torl=[];end
            if isempty(geom_torl), geom_torl=sqrt(eps);end

            % find point to start
            U=findNearest(self,Points,20);

            % use project function to adjust parameter
            U=self.projectPoint(Points,geom_torl,U);
        end

        function U=projectPoint(self,Points_init,geom_torl,U)
            % adjust U by Jacobian transformation
            % also can project point to curve
            %
            if nargin < 4
                U=[];
                if nargin < 3
                    geom_torl=[];
                end
            end

            % find point to start
            if isempty(U)
                U=findNearest(self,Points_init,20);
            end
            [num,~]=size(Points_init);U=U(:);

            iter=0;iter_max=50;
            done=false;idx=1:num;
            while ~done
                [Points,dPoints_dU]=self.calGradient(U(idx));
                dPoint=Points_init(idx,:)-Points;

                RU_RU=sum(dPoints_dU.^2,2);
                RU_D=sum(dPoints_dU.*dPoint,2);

                dU=RU_D./RU_RU;
                dU(isnan(dU) | isinf(dU))=0;

                U(idx)=U(idx)+dU;
                U=max(U,0);U=min(U,1);

                idx=find(abs(RU_D) > geom_torl);

                % Points_inv=self.calPoint(U);
                % scatter3(Points_inv(:,1),Points_inv(:,2),Points_inv(:,3));

                iter=iter+1;
                if isempty(idx) || iter >= iter_max
                    done=true;
                end
            end
        end

        function U=findNearest(self,Points,param)
            % find nearest U in grid
            %
            if nargin < 3, param=[];end
            if isempty(param), param=20;end
            num=size(Points,1);

            % generate rough mesh to initialize pre coord
            U_base=((0:(param-1))+(1:param))/2/param;
            Point_base=self.calPoint(U_base);
            U=zeros(num,1);
            for p_idx=1:num
                point=Points(p_idx,:);
                dis_abs=sum(abs(Point_base-point),2);
                [~,idx]=min(dis_abs,[],"all");
                U(p_idx)=U_base(idx(1));
            end
        end
    
    end
end