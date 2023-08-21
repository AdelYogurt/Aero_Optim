classdef CurveCST2D < handle
    % generate curve by CST parameter
    %
    properties
        name='';
        LX;
        LY;

        % origin parameter
        shape_fcn_Y;
        class_fun_Y;
        symmetry_x;
        
        % deform parameter
        tran_fun_Y=[]; % (U)

        % rotation parameter
        rotation_matrix=[];

        % translation parameter
        translation=[];
    end

    % define function
    methods
        function self=CurveCST2D(name,LX,LY,shape_par_Y,class_par_Y,symmetry_x)
            % generate 2D CST line by LX, LY, shape_par_Y, class_par_Y
            %
            % u, x=X
            % v, y=Y(u)
            %
            % input:
            % LX,LY,shape_par_Y,class_par_Y,symmetry_x
            %
            % notice:
            % shape_fcn_Y(U), class_fun_Y[N1, N2]
            %
            if nargin < 6 || isempty(symmetry_x)
                self.symmetry_x=false(1);
            else
                self.symmetry_x=symmetry_x;
            end
            self.name=name;
            self.LX=LX;
            self.LY=LY;

            if isempty(shape_par_Y)
                self.shape_fcn_Y=@(U) 1;
            elseif isnumeric(shape_par_Y)
                % shape_par_Y only one cross-section parameter
                self.shape_fcn_Y=@(U) self.defunClass...
                    (U,shape_par_Y(1),shape_par_Y(2));
            else
                % function handle
                self.shape_fcn_Y=shape_par_Y;
            end

            self.class_fun_Y=@(U) self.defunClass...
                (U,class_par_Y(1),class_par_Y(2));

        end
    end

    % common function
    methods
        function C=defunClass(self,U,N1,N2)
            % default of Z class function
            % default N1, N2 is relation with U
            %
            NP=self.calNormPar(N1,N2);
            C=U.^N1.*(1-U).^N2./NP;
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
        function [X,Y,U]=calCurve(self,varargin)
            % calculate 2D CST vector by U
            % X=CX*U, if symmetry_x, X=CX*(U-0.5)
            % Y=CY*shape_fcn_Y.*class_fun_Y
            %
            % default
            % u_list=linspace(0,1,u_gird_num+1) if symmetry_x u_list=linspace(0.5,1,u_gird_num+1)
            % value_torl=1e-3
            %
            % u, x=X
            % v, y=Y(u)
            %
            % input:
            % U
            % u_grid_number
            % value_torl
            %
            % notice:
            % shape_fcn_Y(U), class_fun_Y[N1, N2]
            %
            % output:
            % X,Y
            %

            if nargin < 2 || isempty(varargin{1})
                varargin{1}=1e-3;
            end

            if  varargin{1} ~= fix(varargin{1}) && length(varargin{1}) == 1
                value_torl=varargin{1};
                max_level=50;

                % adapt capture U
                if self.symmetry_x
                    low_bou=0.5;
                else
                    low_bou=0;
                end
                up_bou=1;
                [U,fval_list,~]=meshAdapt1D(@(x) coordFcn(self,x),low_bou,up_bou,value_torl,max_level,2);
                X=fval_list(1,:);
                Y=fval_list(2,:);
            else
                if length(varargin{1}) == 1
                    U=[];
                    u_grid_numebr=varargin{1};
                else
                    % mean input U matrix
                    U=varargin{1};
                    u_grid_numebr=length(U)-1;
                end

                if isempty(U)
                    if self.symmetry_x
                        U=linspace(0.5,1,u_grid_numebr+1);
                    else
                        U=linspace(0,1,u_grid_numebr+1);
                    end
                end

                [X,Y]=self.calPoint(U);
            end

            function fval=coordFcn(self,x)
                [x,y]=self.calPoint(x);
                fval=[x,y];
            end
        end

        function [X,Y]=calPoint(self,U)
            % calculate point on curve
            %

            % calculate origin surface matrix
            if self.symmetry_x
                X=self.LX.*(U-0.5);
            else
                X=self.LX.*U;
            end
            
            Y=self.LY*self.shape_fcn_Y(U).*self.class_fun_Y(U);
        end

        function U=calCoordinate(self,X)
            % base on X, Y calculate local coordinate in surface
            %
            U=X./self.LX;
        end
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
x_list=[data_list(1,1),x_list,data_list(2,1)];
fval_list=[data_list(1,2:end)',fval_list,data_list(2,2:end)'];

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
                x_list=[x_list,data_list(data_idx,1)];
                fval_list=[fval_list,data_list(data_idx,2:end)'];
            end

            % node=node.right;
            node_idx=node_list(node_idx,6);
        end
    end
end
