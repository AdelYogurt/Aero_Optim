classdef EdgeBSpline < handle
    % B-spline curve
    % define reference to step standard
    %
    properties
        name='';
        degree;

        ctrl_X=[]; % X(U)
        ctrl_Y=[]; % Y(U)
        ctrl_Z=[]; % Z(U)
        ctrl_num;

        curve_form='.UNSPECIFIED.'; %
        closed_curve='.F.'; % boolean
        self_intersect='.F.'; % boolean

        knot_multi;
        knot_list;
        knot_spec='.UNSPECIFIED.';
    end

    properties
        dimension;
        u_list;

        node_X=[]; % U
        node_Y=[]; % U
        node_Z=[]; % U
        node_num;
    end

    % define curve
    methods
        function self=EdgeBSpline(name,point_X,point_Y,point_Z,FLAG_FIT,...
                degree,knot_multi,knot_list,ctrl_num,U)
            % generate BSpline curve by defining control point of fitting point
            %
            % input:
            % name, point_X, point_Y, point_Z(control/fit point), FLAG_FIT,...
            % degree(optional), knot_multi(optional), knot_list(optional),...
            % ctrl_num(optional), U(optional)
            %
            % output:
            % EdgeBSpline
            %
            % notice:
            % if input degree is empty, degree default is ctrl_num-1,...
            % which is Bezier curve
            %
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
                                if nargin < 5
                                    FLAG_FIT=[];
                                end
                            end
                        end
                    end
                end
            end
            self.name=name;

            if nargin > 1
                if isempty(FLAG_FIT),FLAG_FIT=false;end
                
                % check input value size and giving default value
                if FLAG_FIT
                    node_X=point_X(:);node_Y=point_Y(:);node_Z=point_Z(:);
                    node_num=length(node_X);
                    if length(node_Y) ~= node_num || length(node_Z) ~= node_num
                        error('EdgeBSpline: size of node_X, node_Y, node_Z do not equal');
                    end
                    if isempty(ctrl_num),ctrl_num=node_num;end
                    if ctrl_num > node_num
                        error('EdgeBSpline: ctrl_num more than node_num')
                    end
                else
                    ctrl_X=point_X(:);ctrl_Y=point_Y(:);ctrl_Z=point_Z(:);
                    ctrl_num=length(ctrl_X);
                    if isempty(degree),degree=ctrl_num-1;end
                    node_num=ctrl_num-degree+1;
                    if length(ctrl_Y) ~= ctrl_num || length(ctrl_Z) ~= ctrl_num
                        error('EdgeBSpline: size of ctrl_X, ctrl_Y, ctrl_Z do not equal');
                    end
                end

                % default value
                if isempty(degree),degree=ctrl_num-1;end
                if isempty(U),U=linspace(0,1,node_num);end;U=U(:)';
                if isempty(knot_multi),knot_multi=[degree+1,ones(1,ctrl_num-degree-1),degree+1];end
                if isempty(knot_list),knot_list=interp1(linspace(0,1,length(U)),U,linspace(0,1,ctrl_num-degree+1));end
                u_list=baseKnotVec(knot_multi,knot_list);
                if node_num > 5, du_coord=1/(node_num-3);u_coord=[0,du_coord/2,linspace(du_coord,1-du_coord,node_num-4),1-du_coord/2,1];
                else, u_coord=linspace(0,1,node_num);end
                U=interp1(linspace(0,1,length(knot_list)),knot_list,u_coord);

                if ctrl_num < (degree+1)
                    error('EdgeBSpline: ctrl_num less than degree+1');
                end

                if length(u_list) ~= ctrl_num+degree+1
                    error('EdgeBSpline: number of knot num do not equal to ctrl_num+degree');
                end

                if FLAG_FIT
                    % base on node point list inverse calculate control point list
                    fit_matrix=zeros(node_num,ctrl_num);
                    for node_idx=1:node_num
                        u=U(node_idx);
                        for ctrl_idx=1:ctrl_num
                            fit_matrix(node_idx,ctrl_idx)=baseFcnN(u_list,u,ctrl_idx,degree);
                        end
                    end
                    ctrl_X=fit_matrix\node_X;
                    ctrl_Y=fit_matrix\node_Y;
                    ctrl_Z=fit_matrix\node_Z;
                else
                    % generate B spline curve by control point
                end

                % main properties
                self.degree=degree;
                self.ctrl_X=ctrl_X;
                self.ctrl_Y=ctrl_Y;
                self.ctrl_Z=ctrl_Z;
                self.ctrl_num=ctrl_num;
                self.knot_multi=knot_multi;
                self.knot_list=knot_list;

                % calculate node point list
                self.dimension=3;
                self.u_list=u_list;
                
                [self.node_X,self.node_Y,self.node_Z]=self.calPoint(u_list(degree+1:ctrl_num+1));
                self.node_num=node_num;
            end
            
        end
    end

    % control curve
    methods
        function reverse(self)
            % revese direction of curve
            %
            self.ctrl_X=flipud(self.ctrl_X);
            self.ctrl_Y=flipud(self.ctrl_Y);
            self.ctrl_Z=flipud(self.ctrl_Z);
            self.knot_multi=fliplr(self.knot_multi);
            self.knot_list=1-fliplr(self.knot_list);
            self.u_list=1-fliplr(self.u_list);
            self.node_X=flipud(self.node_X);
            self.node_Y=flipud(self.node_Y);
            self.node_Z=flipud(self.node_Z);
        end

        function changeDegree(self,degree_target)
            if degree_target > self.degree
                degree_old=self.degree;
                ctrl_old=[self.ctrl_X,self.ctrl_Y,self.ctrl_Z];
                ctrl_num_old=self.ctrl_num;
                u_list_old=self.u_list;
                knot_multi_old=self.knot_multi;

                for degree_temp=self.degree+1:degree_target
                    % add repeat node
                    degree_new=degree_old+1;
                    ctrl_num_new=ctrl_num_old+1;
                    knot_multi_new=knot_multi_old+1;
                    u_list_new=baseKnotVec(knot_multi_new,self.knot_list);

                    % calculate new ctrl
                    matrix=zeros(ctrl_num_new,ctrl_num_old);
                    for j=1:ctrl_num_new
                        for i=1:ctrl_num_old
                            matrix(j,i)=baseFcnL(u_list_old,u_list_new,i,j,degree_old);
                        end
                    end
                    matrix=1/degree_new*matrix;
                    ctrl_new=matrix*ctrl_old;

                    % sort old curve data
                    degree_old=degree_new;
                    ctrl_num_old=ctrl_num_new;
                    knot_multi_old=knot_multi_new;
                    u_list_old=u_list_new;
                    ctrl_old=ctrl_new;
                end

                self.degree=degree_new;

                self.ctrl_X=ctrl_new(:,1);
                self.ctrl_Y=ctrl_new(:,2);
                self.ctrl_Z=ctrl_new(:,3);
                self.ctrl_num=ctrl_num_new;

                self.knot_multi=knot_multi_new;

                self.u_list=u_list_new;
                self.node_num=self.ctrl_num-self.degree+1;
                [self.node_X,self.node_Y,self.node_Z]=self.calPoint(u_list_new(degree_new+1:ctrl_num_new+1));
            elseif degree_target < self.degree

            end
        end
    end

    % calculate point
    methods
        function [X,Y,Z,U]=calEdge(self,u_param)
            % generate curve matrix by u_x_list or point_number
            %
            if nargin < 2 || isempty(u_param)
                u_param=1e-3;
            end

            if length(u_param) == 1 && u_param ~= fix(u_param) 
                value_torl=u_param;min_level=5;max_level=50;
                low_bou=0;up_bou=1;
                [U,data_list,~]=meshAdapt1D(@(x) coordFcn(self,x),low_bou,up_bou,value_torl,min_level,max_level,3);
                X=data_list(:,1);
                Y=data_list(:,2);
                Z=data_list(:,3);
            else
                if length(u_param) == 1
                    U=[];
                    u_grid_numebr=u_param;
                else
                    % mean input U matrix
                    U=u_param;
                    u_grid_numebr=length(U)-1;
                end

                if isempty(U)
                    U=linspace(0,1,u_grid_numebr+1);
                end

                if self.dimension == 2
                    [X,Y]=self.calPoint(U);Z=[];
                else
                    [X,Y,Z]=self.calPoint(U);
                end
            end

            function fval=coordFcn(self,x)
                [X_cal,Y_cal,Z_cal]=self.calPoint(x);
                fval=[X_cal,Y_cal,Z_cal];
            end

            function fval=coordFcn2D(self,x)
                [X_cal,Y_cal]=self.calPoint(x);
                fval=[X_cal,Y_cal];
            end
        end

        function [X,Y,Z]=calPoint(self,U)
            % according U to calculate point
            %
            [X,Y,Z]=calBSpline(self,U);
        end

        function [X,Y,Z]=calBSpline(self,U)
            point_num=numel(U);

            X=zeros(size(U));
            Y=zeros(size(U));
            Z=zeros(size(U));

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
                Z(point_idx)=N_list*self.ctrl_Z(idx_start:idx_end);
            end
        end
    end

    % calculate coord
    methods
        function [U,X,Y,Z]=calCoord(self,X,Y,Z)

        end

        function [X,Y,Z,U]=calPorject(self,XO,YO,ZO,U,geo_torl)
            % adjust u by Jacobian transformation
            % also can project point to surface
            %
            if nargin < 6,geo_torl=double(eps('single'));end
        end

        function [dX_dU,dY_dU,dZ_dU]=calGradient(self,U)
            % use differ to calculate gradient
            %
            [X,Y,Z]=self.calPoint(U);
            step=100*eps;

            [X_UF,Y_UF,Z_UF]=self.calPoint(U+step);
            [X_UB,Y_UB,Z_UB]=self.calPoint(U-step);
            Bool_U=(U+step) > 1;
            Bool_D=(U-step) < 0;
            Bool(:,:,1)=Bool_U;Bool(:,:,2)=Bool_D;
            Bool_C=~any(Bool,3);
            dX_dU=X;dY_dU=Y;dZ_dU=Z; % allocate memory
            dX_dU(Bool_C)=(X_UF(Bool_C)-X_UB(Bool_C))/2/step;
            dX_dU(Bool_U)=(X(Bool_U)-X_UB(Bool_U))/2/step;
            dX_dU(Bool_D)=(X_UF(Bool_D)-X(Bool_D))/step;
            dX_dU=real(dX_dU);
            dY_dU(Bool_C)=(Y_UF(Bool_C)-Y_UB(Bool_C))/2/step;
            dY_dU(Bool_U)=(Y(Bool_U)-Y_UB(Bool_U))/2/step;
            dY_dU(Bool_D)=(Y_UF(Bool_D)-Y(Bool_D))/step;
            dY_dU=real(dY_dU);
            if ~isempty(Z)
                dZ_dU(Bool_C)=(Z_UF(Bool_C)-Z_UB(Bool_C))/2/step;
                dZ_dU(Bool_U)=(Z(Bool_U)-Z_UB(Bool_U))/2/step;
                dZ_dU(Bool_D)=(Z_UF(Bool_D)-Z(Bool_D))/step;
                dZ_dU=real(dZ_dU);
            end
        end
    end

    % visualizate curve
    methods
        function drawEdge(self,axe_hdl,u_param,line_option,ctrl_option,node_option)
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
                [X,Y]=self.calEdge(u_param);

                % plot line
                line(axe_hdl,X,Y,line_option);
                if ~isempty(self.node_X) && ~isempty(self.node_Y)
                    line(axe_hdl,self.node_X,self.node_Y,node_option);
                end
                if ~isempty(self.ctrl_X) && ~isempty(self.ctrl_Y)
                    line(axe_hdl,self.ctrl_X,self.ctrl_Y,ctrl_option);
                end
            else
                [X,Y,Z]=self.calEdge(u_param);

                % plot line
                line(axe_hdl,X,Y,Z,line_option);
                if ~isempty(self.node_X) && ~isempty(self.node_Y) && ~isempty(self.node_Z)
                    line(axe_hdl,self.node_X,self.node_Y,self.node_Z,node_option);
                end
                if ~isempty(self.ctrl_X) && ~isempty(self.ctrl_Y) && ~isempty(self.ctrl_Z)
                    line(axe_hdl,self.ctrl_X,self.ctrl_Y,self.ctrl_Z,ctrl_option);
                end
                zlabel('z');
                view(3);
            end
            xlabel('x');
            ylabel('y');

            % axis equal
            % x_range=xlim();
            % y_range=ylim();
            % z_range=zlim();
            % center=[mean(x_range),mean(y_range),mean(z_range)];
            % range=max([x_range(2)-x_range(1),y_range(2)-y_range(1),z_range(2)-z_range(1)])/2;
            % xlim([center(1)-range,center(1)+range]);
            % ylim([center(2)-range,center(2)+range]);
            % zlim([center(3)-range,center(3)+range]);
        end

        function [step_str,object_index]=getStep(self,object_index)
            % write BSpline into step file
            %
            if nargin < 2,object_index=1;end
            out_name=self.name;
            if isempty(out_name),out_name='NONE';end

            control_num=self.ctrl_num;

            % generate CARTESIAN_POINT
            CARTESIAN_POINT=object_index;

            str_control=[];control=[self.ctrl_X,zeros(self.ctrl_num,3-self.dimension)];
            for control_idx=1:control_num
                str=[num2str(object_index,'#%d ='),' CARTESIAN_POINT',...
                    ' ( ',...
                    '''NONE''',', ',...
                    '( ',num2str(control(control_idx,1),'%.16f'),', ',...
                    num2str(control(control_idx,2),'%.16f'),', ',...
                    num2str(control(control_idx,2),'%.16f'),' )',...
                    ' ) ;\n'];
                str_control=[str_control,str];
                object_index=object_index+1;
            end

            % generate curve
            point_index=((1:control_num)-1)+CARTESIAN_POINT;
            str_curve=[num2str(object_index,'#%d ='),' B_SPLINE_CURVE_WITH_KNOTS',...
                ' ( ',...
                '''',out_name,'''',', ',...
                num2str(self.degree,'%d'),', ',...
                '( ',num2str(point_index(1),'#%d'),num2str(point_index(2:end),', #%d'),' ),\n',...
                self.curve_form,', ',self.closed_curve,', ',self.self_intersect,',\n',...
                '( ',num2str(self.knot_multi(1),'%d'),num2str(self.knot_multi(2:end),', %d'),' ),\n',...
                '( ',num2str(self.knot_list(1),'%.16f'),num2str(self.knot_list(2:end),', %.16f'),' ),\n',...
                self.knot_spec,...
                ' ) ;\n'];

            step_str=[str_control,'\n',str_curve,'\n'];

        end

    end
end

%% common function

function [x_list,fval_list,node_list]=meshAdapt1D(fcn,low_bou,up_bou,torl,min_level,max_level,fval_num)
% Binary-tree
% adapt capture function value
% ensure error of linear interplation will less than torl
%
if nargin < 7,fval_num=1;end

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

                if any(abs(fval_c-fval_c_pred) > torl) || node(1) <= min_level
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

function L=baseFcnL(u_list,u_new_list,i,j,k)
% base function of increase degree
%
if k == 0
    if ((u_list(i) <= u_new_list(j)) && (u_new_list(j) <= u_list(i+1)))
        if any(u_list == u_new_list(j)) && u_new_list(j) ~= u_list(1) && u_new_list(j) ~= u_list(end)
            L=0.5;
        else
            L=1;
        end
    else
        L=0;
    end
else
    if u_list(i+k) == u_list(i)
        A=0;
    else
        A=(u_new_list(j+k+1)-u_list(i))/(u_list(i+k)-u_list(i));
    end

    if u_list(i+k+1) == u_list(i+1)
        B=0;
    else
        B=(u_list(i+k+1)-u_new_list(j+k+1))/(u_list(i+k+1)-u_list(i+1));
    end

    L=A*baseFcnL(u_list,u_new_list,i,j,k-1)+B*baseFcnL(u_list,u_new_list,i+1,j,k-1)+baseFcna(u_list,u_new_list,i,j,k);
end
end

function a=baseFcna(u_list,u_new_list,i,j,k)
% discrete base function of BSpline curve
%
if k == 0
    if ((u_list(i) <= u_new_list(j)) && (u_new_list(j) < u_list(i+1)))
        if any(u_list == u_new_list(j)) && u_new_list(j) ~= u_list(1) && u_new_list(j) ~= u_list(end)
            a=0.5;
        else
            a=1;
        end
    else
        a=0;
    end
else
    if u_list(i+k) == u_list(i)
        A=0;
    else
        A=(u_new_list(j+k)-u_list(i))/(u_list(i+k)-u_list(i));
    end

    if u_list(i+k+1) == u_list(i+1)
        B=0;
    else
        B=(u_list(i+k+1)-u_new_list(j+k))/(u_list(i+k+1)-u_list(i+1));
    end

    a=A*baseFcna(u_list,u_new_list,i,j,k-1)+B*baseFcna(u_list,u_new_list,i+1,j,k-1);
end
end