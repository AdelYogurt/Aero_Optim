classdef CurveCST2D < CurveBSpline
    % CST curve
    %
    properties
        LX;
        LY;
        sym; % if true, class U will start from 0.5 to 1

        % origin parameter
        C_par;

        class_fcn; % (U)

        shape_fcn; % (U), default is 1
    end

    properties
        % deform parameter
        deform_fcn=[]; % (U), only y direction

        % rotation parameter
        rotation_matrix=[];

        % translation parameter
        translation=[];
    end

    properties
        fit_data;
    end

    % define curve
    methods
        function self=CurveCST2D(name,LX,LY,C_par,sym)
            % generate 2D CST line by LX, LY, C_par
            %
            % u, x=LX*S(u)
            % v, y=LY*C(u)*S(u)
            %
            % input:
            % name, LX, LY, C_par_Y, symmetry
            %
            % notice:
            % if input N1, N2 == 0, or C_par is empty,...
            % class_fcn will equal to 1
            % class_fcn(U)
            %
            self=self@CurveBSpline(name)
            if nargin < 5
                sym=[];
                if nargin < 4
                    C_par=[];
                end
            end
            self.LX=LX;
            self.LY=LY;

            if isempty(sym),self.sym=false;
            else,self.sym=sym; end

            if isempty(C_par),C_par=[0,0];end

            % class function
            if isnumeric(C_par)
                self.C_par=C_par;
                self.class_fcn=@(U) defcnClass(U,C_par(1),C_par(2));
            else
                self.C_par=[];
                self.class_fcn=C_par;
            end

            self.shape_fcn=@(U) defcnShape(U);
            self.dimension=2;

            function [X,Y]=defcnShape(U)
                X=U;Y=ones(size(X));
            end
        end
    end

    % add BSpline shape function
    methods
        function addShapeBSpline(self,ctrl_X,ctrl_Y,...
                node_X,node_Y,degree,...
                knot_multi,knot_list,ctrl_num,U)
            % fit node point to generate shape function
            if nargin < 10
                U=[];
                if nargin < 9
                    ctrl_num=[];
                    if nargin < 8
                        knot_list=[];
                        if nargin < 7
                            knot_multi = [];
                            if nargin < 6
                                degree=[];
                            end
                        end
                    end
                end
            end

            % check input value size and giving default value
            if isempty(ctrl_X) && isempty(ctrl_Y)
                node_X=node_X(:);node_Y=node_Y(:);
                node_num=length(node_X);
                if length(node_Y) ~= node_num
                    error('CurveCST2D.addShapeBSpline: size of node_X, node_Y do not equal');
                end

                if isempty(ctrl_num),ctrl_num=node_num;end
                if ctrl_num > node_num
                    error('CurveCST2D.addShapeBSpline: ctrl_num can not more than node_num')
                end
            else
                ctrl_X=ctrl_X(:);ctrl_Y=ctrl_Y(:);
                ctrl_num=length(ctrl_X);
                if isempty(degree),degree=ctrl_num-1;end
                node_num=ctrl_num-degree+1;
                if length(ctrl_Y) ~= ctrl_num
                    error('CurveCST2D.addShapeBSpline: size of ctrl_X, ctrl_Y do not equal');
                end
            end

            % default value u
            if isempty(degree),degree=ctrl_num-1;end
            if isempty(U),U=linspace(0,1,node_num);end;U=U(:);
            if isempty(knot_multi),knot_multi=[degree+1,ones(1,ctrl_num-degree-1),degree+1];end
            if isempty(knot_list),knot_list=interp1(linspace(0,1,node_num),U,linspace(0,1,ctrl_num-degree+1));end
            u_list=getKnotVec(knot_multi,knot_list);

            if ctrl_num < (degree+1)
                error('CurveCST2D.addShapeBSpline: ctrl_num less than degree+1');
            end

            if isempty(ctrl_X) && isempty(ctrl_Y)
                % process symmetry
                U_class=U;
                if self.sym,U_class=U_class+0.5;end

                % base on node point list inverse calculate control point list
                matrix=zeros(node_num,ctrl_num);
                for node_idx=1:node_num
                    u=U(node_idx);
                    for ctrl_idx=1:ctrl_num
                        matrix(node_idx,ctrl_idx)=baseFcnN(u_list,u,ctrl_idx,degree);
                    end
                end

                matrix_class=matrix.*self.class_fcn(U_class);

                % undone translation surface
                if ~isempty(self.translation)
                    node_X=node_X-self.translation(1);
                    node_Y=node_Y-self.translation(2);
                end

                % undone rotation surface
                if ~isempty(self.rotation_matrix)
                    matrix=self.rotation_matrix';
                    X_old=node_X;Y_old=node_Y;
                    node_X=matrix(1,1)*X_old+matrix(1,2)*Y_old;
                    node_Y=matrix(2,1)*X_old+matrix(2,2)*Y_old;
                end

                if ~isempty(self.deform_fcn),node_Y=node_Y-self.deform_fcn(U);end

                node_Y=node_Y./self.LY;
                node_X=node_X./self.LX;

                ctrl_X=matrix\node_X;
                ctrl_Y=matrix_class\node_Y;

                self.fit_data.node_X=node_X;
                self.fit_data.node_Y=node_Y;
                self.fit_data.U=U;
                self.fit_data.matrix=matrix;
            elseif ~isempty(ctrl_X) && ~isempty(ctrl_Y)
                % generate B spline curve by control point
                self.fit_data=[];
            else
                error('CurveBSpline: error input, lack control point or node point');
            end

            % main properties
            self.degree=degree;
            self.ctrl_X=ctrl_X;
            self.ctrl_Y=ctrl_Y;
            self.ctrl_num=ctrl_num;
            self.knot_multi=knot_multi;
            self.knot_list=knot_list;

            % shape function
            self.u_list=u_list;
            self.shape_fcn=@(U) self.calPointShape(U);
        end

        function [X,Y]=calPointShape(self,U)
            % according U to calculate point
            %
            point_num=length(U);

            X=zeros(point_num,1);
            Y=zeros(point_num,1);

            N_list=zeros(1,self.degree+1);
            for point_idx=1:point_num
                % local u_x in u_list index
                u=U(point_idx);
                idx_end=self.ctrl_num; % is equal to the section index
                while idx_end > self.degree+1 && u < self.u_list(idx_end)
                    idx_end=idx_end-1;
                end
                idx_start=idx_end-self.degree;
                for N_idx=1:self.degree+1
                    N_list(N_idx)=baseFcnN(self.u_list,u,N_idx+idx_start-1,self.degree);
                end

                X(point_idx)=N_list*self.ctrl_X(idx_start:idx_end);
                Y(point_idx)=N_list*self.ctrl_Y(idx_start:idx_end);
            end
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
            self.class_fcn=@(U) defcnClass(U,C_par_optim(1),C_par_optim(2));
            self.shape_fcn=@(U) self.calPointShape(U);
        end

        function RMSE=fitError(self,C_par,C_par_low_bou,C_par_up_bou)
            % fit error of C_par
            %
            C_par=max(C_par,C_par_low_bou);
            C_par=min(C_par,C_par_up_bou);
            self.class_fcn=@(U) defcnClass(U,C_par(1),C_par(2));
            U=self.fit_data.U;

            node_X=self.fit_data.node_X;
            node_Y=self.fit_data.node_Y;
            matrix=self.fit_data.matrix;
            matrix_class=matrix.*self.class_fcn(U);

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
            % base on angle to rotation surface
            % rotation is anti-clock
            %
            % input:
            % ang(deg)
            %
            cR=cos(ang/180*pi);sR=sin(ang/180*pi);
            matrix=[
                cR -sR
                sR cR];

            self.rotation_matrix=matrix;
        end

        function addTranslate(self,tran_x,tran_y)
            % base on angle to rotation surface
            %
            self.translation=[tran_x,tran_y];
        end
    end

    % calculate point
    methods
        function [X,Y,U]=calCurve(self,u_param)
            % calculate 2D CST vector by U
            % if sys, u will start from 0.5 to 1
            %
            % default
            % u_list=linspace(0,1,u_gird_num+1) if sym, u_list=linspace(0.5,1,u_gird_num+1)
            % value_torl=1e-3
            %
            % u, x=LX*u
            % v, y=LY*C(u)*S(u)
            %
            % input:
            % U
            % u_grid_number
            % value_torl
            %
            % output:
            % X, Y, U
            %
            if nargin < 2 || isempty(u_param)
                u_param=1e-3;
            end

            if  u_param ~= fix(u_param) && length(u_param) == 1
                value_torl=u_param;
                max_level=50;

                % adapt capture U
                low_bou=0;
                up_bou=1;
                if self.sym
                    low_bou(1)=-0.5;up_bou(1)=0.5;
                end
                [U,data_list,~]=meshAdapt1D(@(x) coordFcn(self,x),low_bou,up_bou,value_torl,max_level,2);
                X=data_list(:,1);
                Y=data_list(:,2);
            else
                if length(u_param) == 1
                    U=[];
                    u_grid_num=u_param;
                else
                    % mean input U matrix
                    U=u_param;
                    u_grid_num=length(U)-1;
                end

                if isempty(U)
                    U=linspace(0,1,u_grid_num+1)';
                    if self.sym,U=U-0.5;end
                end

                [X,Y]=self.calPoint(U);
            end

            function fval=coordFcn(self,x)
                [X_cal,Y_cal]=self.calPoint(x);
                fval=[X_cal,Y_cal];
            end
        end

        function [X,Y]=calPoint(self,U)
            % calculate point on curve
            %
            U=U(:);
            if self.sym,U=U-0.5;end

            % calculate origin surface matrix
            U_class=U;
            if self.sym,U_class=U_class+0.5;end
            X=self.LX;
            Y=self.LY*self.class_fcn(U_class);

            [X_cal,Y_cal]=self.shape_fcn(U);
            X=X.*X_cal;
            Y=Y.*Y_cal;

            % deform curve
            if ~isempty(self.deform_fcn)
                Y=Y+self.deform_fcn(U);
            end

            % rotation curve
            if ~isempty(self.rotation_matrix)
                matrix=self.rotation_matrix;
                X_old=X;Y_old=Y;
                X=matrix(1,1)*X_old+matrix(1,2)*Y_old;
                Y=matrix(2,1)*X_old+matrix(2,2)*Y_old;
            end

            % translation curve
            if ~isempty(self.translation)
                X=X+self.translation(1);
                Y=Y+self.translation(2);
            end
        end
    end

    % calculate coord
    methods
        function U=calCoord(self,X,Y)
            % base on X, Y calculate local coordinate in surface
            %
            
            % undone translation surface
            if ~isempty(self.translation)
                X=X-self.translation(1);
                Y=Y-self.translation(2);
            end

            % undone rotation surface
            if ~isempty(self.rotation_matrix)
                matrix=self.rotation_matrix';
                X_old=X;Y_old=Y;
                X=matrix(1,1)*X_old+matrix(1,2)*Y_old;
                Y=matrix(2,1)*X_old+matrix(2,2)*Y_old;
            end
            
            U=X./self.LX;
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
                if ~isempty(self.node_X) && ~isempty(self.node_Y)
                    line(axe_hdl,self.node_X,self.node_Y,node_option);
                end
                if ~isempty(self.ctrl_X) && ~isempty(self.ctrl_Y)
                    line(axe_hdl,self.ctrl_X*self.LX,self.ctrl_Y*self.LY,ctrl_option);
                end
            end
            xlabel('x');
            ylabel('y');

%             axis equal
%             x_range=xlim();
%             y_range=ylim();
%             center=[mean(x_range),mean(y_range)];
%             range=max([x_range(2)-x_range(1),y_range(2)-y_range(1)])/2;
%             xlim([center(1)-range,center(1)+range]);
%             ylim([center(2)-range,center(2)+range]);
        end
    end
end

%% common function

function C=defcnClass(U,N1,N2)
% default of class function
%
NP=calNormPar(N1,N2);
C=U.^N1.*(1-U).^N2./NP;
end

function nomlz_par=calNormPar(N1,N2)
% calculate normailize class function parameter by N1, N2
%
nomlz_par=(N1./(N1+N2)).^N1.*(N2./(N1+N2)).^N2;
nomlz_par((N1 == 0) & (N2 == 0))=1;
end

function knot_vec=getKnotVec(knot_multi,knot_list)
% base on knot_multi and knot_list to create knot vector
%
knot_vec=zeros(1,sum(knot_multi));
start_idx=1;end_idx=knot_multi(1);
for n_idx=1:length(knot_list)
    knot_vec(start_idx:end_idx)=knot_list(n_idx);
    start_idx=start_idx+knot_multi(n_idx);
    end_idx=end_idx+knot_multi(n_idx);
end
end

function N=baseFcnN(u_list,u_x,i,k)
% base function of BSpline curve
%
if k == 0
    if ((u_list(i) <= u_x) && (u_x <= u_list(i+1)))
        if any(u_list == u_x) && u_x ~= u_list(1) && u_x ~= u_list(end)
            N=0.5;
        else
            N=1;
        end
    else
        N=0;
    end
else
    if u_list(i+k) == u_list(i)
        A=0;
    else
        A=(u_x-u_list(i))/(u_list(i+k)-u_list(i));
    end

    if u_list(i+k+1) == u_list(i+1)
        B=0;
    else
        B=(u_list(i+k+1)-u_x)/(u_list(i+k+1)-u_list(i+1));
    end

    N=A*baseFcnN(u_list,u_x,i,k-1)+B*baseFcnN(u_list,u_x,i+1,k-1);
end
end

function [x_list,fval_list,node_list]=meshAdapt1D(fcn,low_bou,up_bou,torl,max_level,fval_num)
% Binary-tree
% adapt capture function value
% ensure error of linear interplation will less than torl
%
if nargin < 6
    fval_num=1;
end

% node_list which is a matrix store all node
% a node is a array, contain level, index_1, index_2, index_c, node_index_1, node_index_2
% place:
% 1-c-2
% cell:
% 1-2
% if children node is empty, left_index or right_index will be zero
list_add_num=50; % node_list will be extend only when node_list is not enough
node_list=zeros(list_add_num,6,'int64');

% data_list use to sort all float data include coodinate, function value
data_list=zeros(list_add_num,fval_num+1);

% add vertex of cell into data_list first
data_list(1,:)=[low_bou,fcn(low_bou)];
data_list(2,:)=[up_bou,fcn(up_bou)];

% create base root
node_list(1,:)=[0,1,2,0,0,0];

node_num=createNodeTree(1,2); % create node tree from root
node_list=node_list(1:node_num,:);
[x_list,fval_list]=traversalInorder(1); % from small to large get list

% add boundary info
x_list=[data_list(1,1);x_list;data_list(2,1)];
fval_list=[data_list(1,2:end);fval_list;data_list(2,2:end)];

    function node_num=createNodeTree(root_idx,data_num)
        % create node tree from root
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
                coord_c=(data_list(node(2),1)+data_list(node(3),1))/2;
                fval_c=fcn(coord_c);
                fval_c_pred=(data_list(node(2),2:end)+data_list(node(3),2:end))/2;
                if any(abs(fval_c-fval_c_pred) > torl)
                    % add 1 new data into data_list
                    data_new_idx=data_num+1;
                    if data_num+1 > size(data_list,1)
                        data_list=[data_list;zeros(list_add_num,fval_num+1)];
                    end
                    data_list(data_new_idx,:)=[coord_c,fval_c];
                    node(4)=data_new_idx;
                    data_num=data_num+1;

                    % add 2 new node to node_list
                    node_new_idx=node_num+[1,2];
                    if node_num+2 > size(node_list,1)
                        node_list=[node_list;zeros(list_add_num,6)];
                    end
                    node_list(node_new_idx,:)=[
                        node(1)+1,node(2),node(4),0,0,0;
                        node(1)+1,node(4),node(3),0,0,0];
                    node([5,6])=node_new_idx;
                    node_num=node_num+2;

                    % add to stack
                    stack=[stack,node_new_idx];
                    node_list(node_idx,:)=node;
                end
            end
        end
    end

    function [x_list,fval_list]=traversalInorder(root_idx)
        % inorder traversal node tree to obtain data
        % inorder will make sure x_list is from small to large
        %
        stack=[];
        x_list=[];
        fval_list=[];
        node_idx=root_idx;

        while ~isempty(stack) || node_idx ~= 0
            while node_idx ~= 0
                stack=[stack,node_idx];

                % node=node.left;
                node_idx=node_list(node_idx,5);
            end

            node_idx=stack(end);
            stack=stack(1:end-1);
            data_idx=node_list(node_idx,4);
            if data_idx ~= 0
                x_list=[x_list;data_list(data_idx,1)];
                fval_list=[fval_list;data_list(data_idx,2:end)];
            end

            % node=node.right;
            node_idx=node_list(node_idx,6);
        end
    end
end
