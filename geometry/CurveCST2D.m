classdef CurveCST2D < CurveBSpline
    % CST curve
    %
    properties
        sym; % if true, class U will start from 0.5 to 1

        % origin parameter
        C_par;

        class_fcn_2D; % (U)

        shape_fcn; % (U), default is LX, LY, output is SX, SY
        class_fcn; % (U), default is 1, CY, output is CX, CY
    end

    properties
        % deform parameter
        deform_fcn=[]; % (U), only y direction

        % rotate parameter
        rotate_matrix=[];

        % translate parameter
        translate=[];
    end

    properties
        fit_data;
    end

    % define curve
    methods
        function self=CurveCST2D(name,C_par,sym,LX,LY)
            % generate 2D CST line by LX, LY, C_par
            %
            % u, x=LX*S(u)
            % v, y=LY*C(u)*S(u)
            %
            % input:
            % name, C_par, sym, LX, LY
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
            self=self@CurveBSpline(name)

            if isempty(sym),self.sym=false;
            else,self.sym=sym; end
            if isempty(C_par),C_par=[0,0];end
            if isempty(LX),LX=1;end
            if isempty(LY),LY=1;end

            % class function
            if isnumeric(C_par)
                self.C_par=C_par;
                self.class_fcn_2D=@(U) baseFcnClass(U,C_par(1),C_par(2));
            else
                self.C_par=[];
                self.class_fcn_2D=C_par;
            end

            self.shape_fcn=@(U) defcnShape(U,LX,LY);
            self.class_fcn=@(U) defcnClass(U,self);
            self.dimension=2;

            function [X,Y]=defcnShape(U,LX,LY)
                X=U*LX;Y=ones(size(X))*LY;
            end

            function [X,Y]=defcnClass(U,self)
                [Y]=self.calClass(U);
                X=ones(size(U));
            end
        end
    end

    % add BSpline shape function
    methods
        function addShapeBSpline(self,point_X,point_Y,FLAG_FIT,...
                degree,knot_multi,knot_list,ctrl_num,U)
            % fit node point to generate shape function
            if nargin < 9
                U=[];
                if nargin < 8
                    ctrl_num=[];
                    if nargin < 7
                        knot_list=[];
                        if nargin < 6
                            knot_multi = [];
                            if nargin < 5
                                degree=[];
                                if nargin < 4
                                    FLAG_FIT=[];
                                end
                            end
                        end
                    end
                end
            end

            if isempty(FLAG_FIT),FLAG_FIT=false;end

            % check input value size and giving default value
            if FLAG_FIT
                node_X=point_X(:);node_Y=point_Y(:);
                node_num=length(node_X);
                if length(node_Y) ~= node_num
                    error('CurveCST2D.addShapeBSpline: size of node_X, node_Y do not equal');
                end
                if isempty(ctrl_num),ctrl_num=node_num;end
                if ctrl_num > node_num
                    error('CurveCST2D.addShapeBSpline: ctrl_num can not more than node_num')
                end
            else
                ctrl_X=point_X(:);ctrl_Y=point_Y(:);
                ctrl_num=length(ctrl_X);
                if isempty(degree),degree=ctrl_num-1;end
                node_num=ctrl_num-degree+1;
                if length(ctrl_Y) ~= ctrl_num
                    error('CurveCST2D.addShapeBSpline: size of ctrl_X, ctrl_Y do not equal');
                end
            end

            % default value
            if isempty(degree),degree=ctrl_num-1;end
            if isempty(U),U=linspace(0,1,node_num);end;U=U(:)';
            if isempty(knot_multi),knot_multi=[degree+1,ones(1,ctrl_num-degree-1),degree+1];end
            if isempty(knot_list),knot_list=interp1(linspace(0,1,node_num),U,linspace(0,1,ctrl_num-degree+1));end
            u_list=baseKnotVec(knot_multi,knot_list);
            if node_num > 5, du_coord=1/(node_num-3);u_coord=[0,du_coord/2,linspace(du_coord,1-du_coord,node_num-4),1-du_coord/2,1];
            else, u_coord=linspace(0,1,node_num);end
            U=interp1(linspace(0,1,length(knot_list)),knot_list,u_coord);

            if ctrl_num < (degree+1)
                error('CurveCST2D.addShapeBSpline: ctrl_num less than degree+1');
            end

            if FLAG_FIT
                % process symmetry
                [node_X,node_Y]=self.axisGlobalToLocal(node_X,node_Y);

                [CY]=self.calClass(U);

                % base on node point list inverse calculate control point list
                fit_matrix=zeros(node_num,ctrl_num);
                for node_idx=1:node_num
                    u=U(node_idx);
                    for ctrl_idx=1:ctrl_num
                        fit_matrix(node_idx,ctrl_idx)=baseFcnN(u_list,u,ctrl_idx,degree);
                    end
                end

                matrix_class=fit_matrix.*CY;

                % undone translate surface
                if ~isempty(self.translate)
                    node_X=node_X-self.translate(1);
                    node_Y=node_Y-self.translate(2);
                end

                % undone rotate surface
                if ~isempty(self.rotate_matrix)
                    matrix=self.rotate_matrix';
                    X_old=node_X;Y_old=node_Y;
                    node_X=matrix(1,1)*X_old+matrix(1,2)*Y_old;
                    node_Y=matrix(2,1)*X_old+matrix(2,2)*Y_old;
                end

                if ~isempty(self.deform_fcn),node_Y=node_Y-self.deform_fcn(U);end

                ctrl_X=fit_matrix\node_X;
                ctrl_Y=matrix_class\node_Y;

                self.fit_data.node_X=node_X;
                self.fit_data.node_Y=node_Y;
                self.fit_data.U=U;
                self.fit_data.matrix=fit_matrix;
            else
                % generate B spline curve by control point
                self.fit_data=[];
            end

            % main properties
            self.degree=degree;
            self.ctrl_X=ctrl_X;
            self.ctrl_Y=ctrl_Y;
            self.ctrl_Z=zeros(size(ctrl_Y));
            self.ctrl_num=ctrl_num;
            self.knot_multi=knot_multi;
            self.knot_list=knot_list;

            % shape function
            self.u_list=u_list;

            self.shape_fcn=@(U) self.calBSpline(U);
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
                error('CurveCST2D.optimClass: input C_par is handle, cannot optimize coefficient of class function');
            end
            
            C_par_low_bou=[0,0];
            C_par_up_bou=[1e3,1e3];
            C_par_init=self.C_par;
            obj_fit=@(C_par) self.fitError(C_par,C_par_low_bou,C_par_up_bou);
            C_par_optim=fminunc(obj_fit,C_par_init,optim_option);

            self.C_par=C_par_optim;
            self.class_fcn_2D=@(U) baseFcnClass(U,C_par_optim(1),C_par_optim(2));
        end

        function RMSE=fitError(self,C_par,C_par_low_bou,C_par_up_bou)
            % fit error of C_par
            %
            C_par=max(C_par,C_par_low_bou);
            C_par=min(C_par,C_par_up_bou);
            self.class_fcn_2D=@(U) baseFcnClass(U,C_par(1),C_par(2));
            U=self.fit_data.U;

            node_X=self.fit_data.node_X;
            node_Y=self.fit_data.node_Y;
            matrix=self.fit_data.matrix;
            matrix_class=matrix.*self.class_fcn_2D(U);

            % reverse calculate control point
            ctrl_X=matrix\node_X;
            ctrl_Y=matrix_class\self.fit_data.node_Y;
            self.ctrl_X=ctrl_X;
            self.ctrl_Y=ctrl_Y;

            [X_cal,Y_cal]=self.calPoint(U);
            RMSE=sqrt(mean([(node_X-X_cal).^2;(node_Y-Y_cal).^2]));
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
                cR -sR
                sR cR];

            self.rotate_matrix=matrix;
        end

        function addTranslate(self,tran_x,tran_y)
            % base on angle to rotate surface
            %
            self.translate=[tran_x,tran_y];
        end

        function [X,Y]=axisLocalToGlobal(self,X,Y,U)
            % deform curve
            if ~isempty(self.deform_fcn)
                Y=Y+self.deform_fcn(U);
            end

            % rotate curve
            if ~isempty(self.rotate_matrix)
                matrix=self.rotate_matrix;
                X_old=X;Y_old=Y;
                X=matrix(1,1)*X_old+matrix(1,2)*Y_old;
                Y=matrix(2,1)*X_old+matrix(2,2)*Y_old;
            end

            % translate curve
            if ~isempty(self.translate)
                X=X+self.translate(1);
                Y=Y+self.translate(2);
            end
        end

        function [X,Y]=axisGlobalToLocal(self,X,Y,U)
            % re-rotate curve
            if ~isempty(self.rotate_matrix)
                matrix=self.rotate_matrix';
                X_old=X;Y_old=Y;
                X=matrix(1,1)*X_old+matrix(1,2)*Y_old;
                Y=matrix(2,1)*X_old+matrix(2,2)*Y_old;
            end

            % translate curve
            if ~isempty(self.translate)
                X=X-self.translate(1);
                Y=Y-self.translate(2);
            end

            % re-deform curve
            if ~isempty(self.deform_fcn)
                Y=Y-self.deform_fcn(U);
            end
        end

    end

    % calculate point
    methods
        function [X,Y]=calPoint(self,U)
            % calculate point on curve
            %
            
            % calculate origin curve
            [SX,SY]=self.shape_fcn(U);
            [CX,CY]=self.class_fcn(U);
            X=CX.*SX;Y=CY.*SY;
            [X,Y]=self.axisLocalToGlobal(X,Y,U);
            Z=[];
        end

        function [Y]=calClass(self,U)
            % calculate class
            U_class=U;
            if self.sym,U_class=(U_class/2)+0.5;end
            Y=self.class_fcn_2D(U_class);
        end

    end

    % calculate coord
    methods
        function U=calCoord(self,X,Y)
            % base on X, Y calculate local coordinate in surface
            %
            
            % undone translate surface
            if ~isempty(self.translate)
                X=X-self.translate(1);
                Y=Y-self.translate(2);
            end

            % undone rotate surface
            if ~isempty(self.rotate_matrix)
                matrix=self.rotate_matrix';
                X_old=X;Y_old=Y;
                X=matrix(1,1)*X_old+matrix(1,2)*Y_old;
                Y=matrix(2,1)*X_old+matrix(2,2)*Y_old;
            end
            
            U=X./abs(self.shape_fcn(0)-self.shape_fcn(1));
            U=max(U,0);U=min(U,1);
        end
    end

    % visualizate function
    methods
        function drawCurve(self,axe_hdl,u_param,line_option,ctrl_option,node_option)
            % draw curve on figure handle
            %
            if nargin < 6
                node_option=[];
                if nargin < 5
                    ctrl_option=[];
                    if nargin < 4
                        line_option=[];
                        if nargin < 3
                            u_param=[];
                            if nargin < 2
                                axe_hdl=[];
                            end
                        end
                    end
                end
            end

            if isempty(axe_hdl),axe_hdl=axes(figure());end

            % default draw option
            if isempty(line_option)
                line_option=struct();
            end
            if isempty(node_option)
                node_option=struct('Marker','o','LineStyle','none');
            end
            if isempty(ctrl_option)
                ctrl_option=struct('Marker','s','LineStyle','--','Color','r');
            end

            % calculate point on curve
            if self.dimension == 2
                [X,Y]=self.calCurve(u_param);

                % plot line
                line(axe_hdl,X,Y,line_option);
                if ~isempty(self.ctrl_X) && ~isempty(self.ctrl_Y)
                    U=interp1(linspace(0,1,self.ctrl_num-self.degree+1),self.u_list(self.degree+1:self.ctrl_num+1),linspace(0,1,self.ctrl_num))';
                    [CY]=self.calClass(U);
                    X=self.ctrl_X;Y=CY.*self.ctrl_Y;
                    [X,Y]=self.axisLocalToGlobal(X,Y);
                    line(axe_hdl,X,Y,ctrl_option);
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

    end
end
