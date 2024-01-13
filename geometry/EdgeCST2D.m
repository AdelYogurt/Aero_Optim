classdef EdgeCST2D < EdgeNURBS
    % CST curve
    %
    properties
        sym; % if true, class U will start from 0.5 to 1
        C_par; % origin parameter
        shape_fcn; % (U), default is LX, LY, output is SX, SY
        class_fcn; % (U), default is 1, CY, output is CX, CY
    end

    properties
        deform_fcn=[]; % deform parameter: deform_fcn(U), only y direction
        rotate_matrix=[]; % rotate parameter
        translate=[]; % translate parameter
    end

    properties
        fit_data;
    end

    % define curve
    methods
        function self=EdgeCST2D(name,C_par,sym,LX,LY)
            % generate 2D CST line by LX, LY, C_par
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
            self=self@EdgeNURBS(name)

            if isempty(sym),self.sym=false;
            else,self.sym=sym; end
            if isempty(C_par),C_par=[0,0];end
            if isempty(LX),LX=1;end
            if isempty(LY),LY=1;end

            self.shape_fcn=@(U) [U*LX,ones(size(U))*LY];

            % class function
            if isnumeric(C_par)
                self.C_par=C_par;
                self.class_fcn=@(U) [ones(size(U)),baseFcnClass(U,C_par(1),C_par(2))];
            else
                self.C_par=[];
                self.class_fcn=C_par;
            end

            self.dimension=2;
        end
    
        function Point=calCST(self,U)
            % calculate class
            U_class=U;
            if self.sym,U_class=(U_class/2)+0.5;end
            CPoint=self.class_fcn(U_class);

            % calculate shape
            SPoint=self.shape_fcn(U);

            Point=SPoint.*CPoint;
        end
    end

    % add NURBS shape function
    methods
        function addNURBS(self,Poles,Degree,Mults,Knots,Weights)
            % add NURBS as shape function
            %
            if nargin < 6
                Weights=[];
                if nargin < 5
                    Knots = [];
                    if nargin < 4
                        Mults=[];
                        if nargin < 3
                            Degree=[];
                        end
                    end
                end
            end
            Poles=Poles(:,1:2);
            pole_num=size(Poles,1);

            % default value
            if isempty(Degree),Degree=pole_num-1;end
            if isempty(Mults),Mults=[Degree+1,ones(1,pole_num-Degree-1),Degree+1];end
            if isempty(Knots),Knots=linspace(0,1,pole_num-Degree+1);end
            if isempty(Weights), Weights=ones(1,pole_num);end;Weights=Weights(:)';
            u_list=baseKnotVec(Mults,Knots);

            if pole_num < (Degree+1)
                error('EdgeCST2D.addNURBS: pole_num less than Degree+1');
            end

            if length(u_list) ~= pole_num+Degree+1
                error('EdgeCST2D.addNURBS: knot_num is not equal to pole_num+Degree+1');
            end

            % Explicit Attributes
            self.Degree=Degree;
            self.Poles=Poles;
            self.Mults=Mults;
            self.Knots=Knots;
            self.Weights=Weights;

            % Derived Attributes
            self.pole_num=pole_num;
            self.u_list=u_list;

            % add shape function
            self.shape_fcn=@(U) self.calNURBS(U);
        end

        function fitNURBS(self,Nodes,Degree,pole_num,U_node)
            % fit NURBS base on CST class function
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
            if isempty(pole_num),pole_num=node_num;end
            if pole_num > node_num
                error('EdgeCST2D.fitNURBS: pole_num can not more than node_num')
            end

            % default value
            if isempty(U_node)
                U_node=vecnorm(Nodes(2:end,:)-Nodes(1:end-1,:),2,2);
                U_node=[0;cumsum(U_node)];U_node=U_node/U_node(end);
            end
            U_node=U_node(:);

            Mults=[Degree+1,ones(1,pole_num-Degree-1),Degree+1];
            Knots=linspace(0,1,node_num-Degree+1);
            for j=2:node_num-Degree
                Knots(j)=mean(U_node(j:j+Degree-1));
            end
            Knots=interp1(linspace(0,1,node_num-Degree+1),Knots,linspace(0,1,pole_num-Degree+1));
            u_list=baseKnotVec(Mults,Knots);

            % translate node to local coordinate
            Nodes=self.axisGlobalToLocal(Nodes,U_node);

            % base on node point list inverse calculate control point list
            fit_matrix=zeros(node_num,pole_num);
            for node_idx=1:node_num
                u=U_node(node_idx);
                for ctrl_idx=1:pole_num
                    fit_matrix(node_idx,ctrl_idx)=baseFcnN(u_list,u,ctrl_idx,Degree);
                end
            end

            % add class coefficient
            if self.sym,U_node=(U_node/2)+0.5;end
            Point=self.class_fcn(U_node);
            matrix_class=fit_matrix.*Point(:,2);

            % reverse calculate control point
            Poles=[fit_matrix\Nodes(:,1),matrix_class\Nodes(:,2)];
            Weights=ones(1,pole_num);

            % Explicit Attributes
            self.Degree=Degree;
            self.Poles=Poles;
            self.Mults=Mults;
            self.Knots=Knots;
            self.Weights=Weights;

            % Derived Attributes
            self.pole_num=pole_num;
            self.u_list=u_list;

            % add shape function
            self.shape_fcn=@(U) self.calNURBS(U);

            % fit data
            self.fit_data.Nodes=Nodes;
            self.fit_data.U_node=U_node;
            self.fit_data.fit_matrix=fit_matrix;
        end

        function optimClass(self,optim_option)
            % optimization coefficient of class function
            %
            if nargin < 2
                optim_option=[];
            end

            if isempty(optim_option)
                optim_option=optimoptions('fminunc','Display','none');
            end

            if isempty(self.C_par)
                error('EdgeCST2D.optimClass: input C_par is handle, cannot optimize coefficient of class function');
            end

            C_par_low_bou=[0,0];
            C_par_up_bou=[1e3,1e3];
            C_par_init=self.C_par;
            obj_fit=@(C_par) self.fitError(C_par,C_par_low_bou,C_par_up_bou);
            C_par_optim=fminunc(obj_fit,C_par_init,optim_option);

            C_par_optim=max(C_par_optim,C_par_low_bou);
            C_par_optim=min(C_par_optim,C_par_up_bou);
            self.C_par=C_par_optim;
            self.class_fcn=@(U) [ones(size(U)),baseFcnClass(U,C_par_optim(1),C_par_optim(2))];
        end

        function RMSE=fitError(self,C_par,C_par_low_bou,C_par_up_bou)
            % fit error of C_par
            %
            C_par=max(C_par,C_par_low_bou);
            C_par=min(C_par,C_par_up_bou);

            % updata C_par
            self.class_fcn=@(U) [ones(size(U)),baseFcnClass(U,C_par(1),C_par(2))];

            % load fit data
            U_node=self.fit_data.U_node;
            Nodes=self.fit_data.Nodes;
            fit_matrix=self.fit_data.fit_matrix;
            
            % add class coefficient
            Point=self.class_fcn(U_node);
            matrix_class=fit_matrix.*Point(:,2);

            % reverse calculate control point
            Poles=[fit_matrix\Nodes(:,1),matrix_class\Nodes(:,2)];
            self.Poles=Poles;

            Point_cal=self.calCST(U_node);
            RMSE=sqrt(mean((Nodes-Point_cal).^2,"all"));
        end
    end

    % deform, rotate, translate
    methods
        function addDeform(self,deform_fcn_Y)
            % base on local coordinate deform curve
            %
            % input:
            % deform_fcn_Y(U)
            %
            self.deform_fcn=deform_fcn_Y;
        end

        function addRotate(self,ang)
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

        function addTranslate(self,tran_x,tran_y)
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

    % calculate point
    methods
        function Point=calPoint(self,U)
            % calculate point on curve
            %
            U=U(:);
            Point=calCST(self,U);
            Point=self.axisLocalToGlobal(Point,U);
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
    end

    % calculate coord
    methods
        function U=calCoord(self,Point)
            % base on X, Y calculate local coordinate in surface
            %

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

            LP=self.shape_fcn(0)-self.shape_fcn(1);
            U=Point(:,1)./abs(LP(1));
            U=max(U,0);U=min(U,1);
        end
    end

    % visualizate function
    methods
        function drawEdge(self,axe_hdl,u_param,crv_option,ctrl_option)
            % draw curve on figure handle
            %
            if nargin < 5
                ctrl_option=[];
                if nargin < 4
                    crv_option=[];
                    if nargin < 3
                        u_param=[];
                        if nargin < 2
                            axe_hdl=[];
                        end
                    end
                end
            end

            if isempty(axe_hdl),axe_hdl=axes(figure());end

            % default draw option
            if isempty(crv_option)
                crv_option=struct();
            end
            if isempty(ctrl_option)
                ctrl_option=struct('Marker','s','LineStyle','--','Color','r');
            end

            % calculate point on curve
            Points=self.calEdge(u_param);

            % plot line
            if self.dimension == 2
                line(axe_hdl,Points(:,1),Points(:,2),crv_option);
                if ~isempty(self.Poles)
                    u_num=self.pole_num-self.Degree+1;
                    U_pole=interp1(linspace(0,1,u_num),self.u_list(self.Degree+1:self.pole_num+1),linspace(0,1,self.pole_num))';
                    if self.sym;U_pole=U_pole/2+0.5;end
                    Poles=self.Poles.*self.class_fcn(U_pole);
                    Poles=self.axisLocalToGlobal(Poles);
                    line(axe_hdl,Poles(:,1),Poles(:,2),ctrl_option);
                end
            end
            xlabel('x');
            ylabel('y');

            % axis equal
            % x_range=xlim();
            % y_range=ylim();
            % center=[mean(x_range),mean(y_range)];
            % range=max([x_range(2)-x_range(1),y_range(2)-y_range(1)])/2;
            % xlim([center(1)-range,center(1)+range]);
            % ylim([center(2)-range,center(2)+range]);
        end

        function edg=getNURBS(self,u_param)
            % convert CST Edge into NURBS Edge
            %
            % input:
            % u_param
            % or:
            % value_torl
            %
            if nargin < 2
                u_param=[];
            end

            if ~isempty(u_param) && length(u_param) > 1 && u_param(1,1) > u_param(1,2)
                u_param=fliplr(u_param);
            end
            [Nodes,U_node]=self.calEdge(u_param);

            Degree=min(length(U_node)-1,3);
            edg=GeomApp.VertexToEdge(self.name,Nodes,Degree,[],U_node);
        end
    end
end
