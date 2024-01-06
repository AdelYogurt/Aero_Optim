classdef GeomApp
    methods(Static)
        function edg=VertexToEdge(name,Nodes,Degree,pole_num,U_node)
            % generate BSpline curve by defined fitting points
            %
            % input:
            % name:
            % Nodes (matrix): fit point, node_num x dimension matrix
            % Degree (optional):
            % Mults (optional):
            % Knots (optional):
            % pole_num (optional):
            % U_node (optional):
            %
            % output:
            % EdgeNURBS
            %
            % notice:
            % if input Degree is empty, Degree default is pole_num-1, which is Bezier curve
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
                error('GeomApp.VertexToEdge: pole_num more than node_num')
            end
            % default value
            if isempty(U_node)
                U_node=vecnorm(Nodes(2:end,:)-Nodes(1:end-1,:),2,2);
                U_node=[0;cumsum(U_node)];U_node=U_node/U_node(end);
            end
            U_node=U_node(:);

            % reference:
            % [1] 施法中. 计算机辅助几何设计与非均匀有理B样条[M]. 208-281.
            % [2] https://blog.csdn.net/he_nan/article/details/107134271
            Mults=[Degree+1,ones(1,pole_num-Degree-1),Degree+1];
            Knots=linspace(0,1,node_num-Degree+1);
            for j=2:node_num-Degree
                Knots(j)=mean(U_node(j:j+Degree-1));
            end
            % modify
            Knots=interp1(linspace(0,1,node_num-Degree+1),Knots,linspace(0,1,pole_num-Degree+1));
            % c=(node_num)/(pole_num-Degree);
            % for j=1:pole_num-Degree-1
            %     i=floor(j*c);
            %     a=(j*c-i);
            %     Knots(j+1)=(1-a)*u_node(i)+a*u_node(i+1);
            % end
            u_list=baseKnotVec(Mults,Knots);

            % base on node point list inverse calculate control point list
            fit_matrix=zeros(node_num,pole_num);
            for node_idx=1:node_num
                u=U_node(node_idx);
                for ctrl_idx=1:pole_num
                    fit_matrix(node_idx,ctrl_idx)=baseFcnN(u_list,u,ctrl_idx,Degree);
                end
            end
            Poles=fit_matrix\Nodes;

            edg=EdgeNURBS(name,Poles,Degree,Mults,Knots);
        end

        function fce=VertexToFace(name,Nodes,UDegree,VDegree,u_pole_num,v_pole_num,U_node,V_node)
            % generate BSpline curve by defined fitting points
            %
            % input:
            % name:
            % Nodes (matrix): fit point, u_node_num x v_node_num x dimension matrix
            % Degree (optional):
            % Mults (optional):
            % Knots (optional):
            % u_pole_num (optional):
            % v_pole_num (optional):
            % U_node (optional):
            % V_node (optional):
            %
            % output:
            % FaceNURBS
            %
            % notice:
            % if input Degree is empty, Degree default is pole_num-1, which is Bezier curve
            %
            if nargin < 8
                V_node=[];
                if nargin < 7
                    U_node=[];
                    if nargin < 6
                        v_pole_num=[];
                        if nargin < 5
                            u_pole_num=[];
                            if nargin < 4
                                VDegree=[];
                                if nargin < 3
                                    UDegree=[];
                                end
                            end
                        end
                    end
                end
            end

            [v_node_num,u_node_num,dimension]=size(Nodes);
            if isempty(u_pole_num),u_pole_num=u_node_num;end
            if isempty(v_pole_num),v_pole_num=v_node_num;end
            if u_pole_num > u_node_num || v_pole_num > v_node_num
                error('GeomApp.VertexToEdge: pole_num more than node_num')
            end

            % default value
            if isempty(U_node)
                U_node=vecnorm(Nodes(:,2:end,:)-Nodes(:,1:end-1,:),2,3);
                U_node=mean(U_node,1);
                U_node=[0;cumsum(U_node')];U_node=U_node/U_node(end);
            end
            U_node=U_node(:);
            if isempty(V_node)
                V_node=vecnorm(Nodes(2:end,:,:)-Nodes(1:end-1,:,:),2,3);
                V_node=mean(V_node,2);
                V_node=[0;cumsum(V_node)];V_node=V_node/V_node(end);
            end
            V_node=V_node(:);
            
            UMults=[UDegree+1,ones(1,u_pole_num-UDegree-1),UDegree+1];
            UKnots=linspace(0,1,u_node_num-UDegree+1);
            for j=2:u_node_num-UDegree
                UKnots(j)=mean(U_node(j:j+UDegree-1));
            end
            % modify
            UKnots=interp1(linspace(0,1,u_node_num-UDegree+1),UKnots,linspace(0,1,u_pole_num-UDegree+1));
            u_list=baseKnotVec(UMults,UKnots);

            VMults=[VDegree+1,ones(1,v_pole_num-VDegree-1),VDegree+1];
            VKnots=linspace(0,1,v_node_num-VDegree+1);
            for j=2:v_node_num-VDegree
                VKnots(j)=mean(V_node(j:j+VDegree-1));
            end
            % modify
            VKnots=interp1(linspace(0,1,v_node_num-VDegree+1),VKnots,linspace(0,1,v_pole_num-VDegree+1));
            v_list=baseKnotVec(VMults,VKnots);

            % base on node point list inverse calculate control point list
            matrix_v=zeros(v_node_num,v_pole_num);
            for node_idx=1:v_node_num
                v=V_node(node_idx);
                for ctrl_idx=1:v_pole_num
                    matrix_v(node_idx,ctrl_idx)=baseFcnN(v_list,v,ctrl_idx,VDegree);
                end
            end
            matrix_u=zeros(u_pole_num,u_node_num);
            for node_idx=1:u_node_num
                u=U_node(node_idx);
                for ctrl_idx=1:u_pole_num
                    matrix_u(ctrl_idx,node_idx)=baseFcnN(u_list,u,ctrl_idx,UDegree);
                end
            end
            Poles=pagemrdivide(pagemldivide(matrix_v,Nodes),matrix_u);

            fce=FaceNURBS(name,Poles,UDegree,VDegree,UMults,VMults,UKnots,VKnots);
        end

    end

    methods(Static)
        function edg=JointEdge(name,edg_list,geom_torl)
            % connect BSpline curve into one curve
            %
            if nargin < 3
                geom_torl=[];
            end

            if isempty(geom_torl),geom_torl=100*eps;end

            % search max Degree
            Degree=0;
            edg_idx=1;
            while edg_idx <= length(edg_list)
                edg_curr=edg_list{edg_idx};
                if isempty(edg_curr)
                    edg_list(edg_idx)=[];
                else
                    if Degree < edg_curr.Degree
                        Degree=edg_curr.Degree;
                    end
                end
                edg_idx=edg_idx+1;
            end

            % increase Degree of edg and load data
            edg_num=length(edg_list);
            Poles=[];
            Mults=[];
            Knots=[];
            Weights=[];
            total_length=0;
            for edg_idx=1:edg_num
                edg=edg_list{edg_idx};
                edg.addDegree(Degree);
                Mults_new=edg.Mults;
                Mults_new(1)=Mults_new(1)-1;
                Weights_new=edg.Weights;
                edg_length=sum(vecnorm(diff(edg.Poles),2,2));

                Poles=[Poles(1:end-1,:);edg.Poles];
                Mults=[Mults(1:end-1),Mults_new];
                Knots=[Knots(1:end-1),edg.Knots*edg_length+total_length];
                Weights=[Weights(1:end-1),Weights_new/Weights_new(1)];

                Weights=Weights/Weights(end);
                total_length=total_length+edg_length;
            end
            Mults(1)=Mults(1)+1;

            % connect
            Knots=Knots/total_length;
            Knots=Knots/max(Knots);

            edg=EdgeNURBS(name,Poles,Degree,Mults,Knots,Weights);
        end

        function edg_list=correctEdge(edg_list,geom_torl)
            % correct line point order
            % order should be anti clockwise and start from first line
            %
            edg_num=length(edg_list);

            % first edg check if connext new edg
            edg_curr=edg_list{1};
            ctrl_curr=[edg_curr.ctrl_X,edg_curr.ctrl_Y,edg_curr.ctrl_Z];
            edg_next=edg_list{2};
            ctrl_next=[edg_next.ctrl_X,edg_next.ctrl_Y,edg_next.ctrl_Z];
            if norm(ctrl_curr(end,:)-ctrl_next(1,:)) > geom_torl
                edg_list=fliplr(edg_list);
            end

            % start from first edg
            for edg_idx=1:edg_num-1
                edg_curr=edg_list{edg_idx};
                ctrl_curr=[edg_curr.ctrl_X,edg_curr.ctrl_Y,edg_curr.ctrl_Z];
                edg_next=edg_list{edg_idx+1};
                ctrl_next=[edg_next.ctrl_X,edg_next.ctrl_Y,edg_next.ctrl_Z];
                if norm(ctrl_curr(end,:)-ctrl_next(1,:)) > geom_torl
                    % load remain line all vertex
                    vertex_list=zeros((edg_num-edg_idx)*2,3);
                    for remain_idx=edg_idx+1:edg_num
                        edg_rema=edg_list{remain_idx};
                        ctrl_rema=[edg_rema.ctrl_X,edg_rema.ctrl_Y,edg_rema.ctrl_Z];
                        vertex_list(2*(remain_idx-edg_idx)-1,:)=ctrl_rema(1,:);
                        vertex_list(2*(remain_idx-edg_idx),:)=ctrl_rema(end,:);
                    end

                    % search next connect point
                    dist=vecnorm((ctrl_curr(end,:)-vertex_list),2,2);
                    overlap_idx=find(dist < geom_torl);
                    if ~any(overlap_idx)
                        error('geomMapGrid: correctLine: line not connect');
                    end
                    exchange_idx=ceil(overlap_idx/2)+edg_idx;

                    % exchange line
                    line_temp=edg_list{exchange_idx};
                    edg_list{exchange_idx}=edg_next;
                    edg_list{edg_idx+1}=line_temp;

                    % reorder point in line order
                    if mod(overlap_idx,2) == 0
                        edg_list{edg_idx+1}.reverse;
                    end
                end
            end
        end

        function Point=MapGrid(line_u0,line_u1,line_v0,line_v1,U,V)
            % generate discrete area by mapping
            % using Lagrange polynomial mapping method
            % local parameter u and v is equispaced
            %
            % input:
            % line_u0 (matrix): u_num x dimension
            % line_u1 (matrix): u_num x dimension
            % line_v0 (matrix): v_num x dimension
            % line_v1 (matrix): v_num x dimension
            %
            % output:
            % Point (matrix): v_num x u_num x dimension
            %
            if nargin < 6
                V=[];
                if nargin < 5
                    U=[];
                end
            end

            geom_torl=100*eps;
            line_list={line_u0,line_v1,flipud(line_u1),flipud(line_v0)};
            line_list=GeomApp.correctLine(line_list,geom_torl);
            line_u0=line_list{1};
            line_v1=line_list{2};
            line_u1=flipud(line_list{3});
            line_v0=flipud(line_list{4});

            u_num=size(line_u0,1);v_num=size(line_v0,1);
            if u_num ~= size(line_u1,1) || v_num ~= size(line_v1,1)
                error('GeomApp.MapGrid: opposite line of boundary discrete number no equal')
            end
            dimension=size(line_u0,2);
            if isempty(U),U=linspace(0,1,u_num);end
            if isempty(V),V=linspace(0,1,v_num);end

            Point=zeros(v_num,u_num,dimension);

            % preproces boundary
            Point(1,:,:)=line_u0(:,:);
            Point(end,:,:)=line_u1(:,:);
            Point(:,1,:)=line_v0(:,:);
            Point(:,end,:)=line_v1(:,:);

            if u_num > 2 && v_num > 2
                H=diff(U);G=diff(V);

                H_i_j=H(1:end-1)+H(2:end);
                H_ij=H(1:end-1).*H(2:end).*H_i_j/2;
                G_i_j=G(1:end-1)+G(2:end);
                G_ij=G(1:end-1).*G(2:end).*G_i_j/2;

                % solve prepare
                % D_xx
                d_xx=spdiags([H(2:end)./H_ij;-H_i_j./H_ij;H(1:end-1)./H_ij]',-1:1,v_num-2,v_num-2)';
                I_M=speye(u_num-2,u_num-2);
                D_xx=kron(I_M,d_xx);

                % D_yy
                d_yy=spdiags([G(2:end)./G_ij;-G_i_j./G_ij;G(1:end-1)./G_ij]',-1:1,u_num-2,u_num-2)';
                I_K=speye(v_num-2,v_num-2);
                D_yy=kron(d_yy,I_K);

                % Laplace matrix
                Lap=(D_xx+D_yy);

                Point_B=zeros(v_num-2,u_num-2,dimension);
                Point_B(1,:,:)=Point_B(1,:,:)-Point(1,2:u_num-1,:)/G(2)^2;
                Point_B(end,:,:)=Point_B(end,:,:)-Point(end,2:u_num-1,:)/G(end-1)^2;
                Point_B(:,1,:)=Point_B(:,1,:)-Point(2:v_num-1,1,:)/H(2)^2;
                Point_B(:,end,:)=Point_B(:,end,:)-Point(2:v_num-1,end,:)/H(end-1)^2;
                for dim_idx=1:dimension
                    Point_B_page=Point_B(:,:,dim_idx);
                    Point(2:v_num-1,2:u_num-1,dim_idx)=reshape(Lap\Point_B_page(:),v_num-2,u_num-2);
                end
            end
        end

        function line_list=correctLine(line_list,geom_torl)
            % correct line point order
            % order should be anti clockwise and start from first line
            %
            line_num=length(line_list);

            % start from first line
            for line_idx=1:line_num-1
                line_curr=line_list{line_idx};
                line_next=line_list{line_idx+1};
                if norm(line_curr(end,:)-line_next(1,:)) > geom_torl
                    % load remain line all vertex
                    vertex_list=zeros((line_num-line_idx)*2,size(line_curr,2));
                    for remain_idx=line_idx+1:line_num
                        line_rema=line_list{remain_idx};
                        vertex_list(2*(remain_idx-line_idx)-1,:)=line_rema(1,:);
                        vertex_list(2*(remain_idx-line_idx),:)=line_rema(end,:);
                    end

                    % search next connect point
                    dist=vecnorm((line_curr(end,:)-vertex_list),2,2);
                    overlap_idx=find(dist < geom_torl);
                    if ~any(overlap_idx)
                        save();
                        error('geomMapGrid: correctLine: line not connect');
                    end
                    exchange_idx=ceil(overlap_idx/2)+line_idx;

                    % exchange line
                    line_temp=line_list{exchange_idx};
                    line_list{exchange_idx}=line_next;
                    line_list{line_idx+1}=line_temp;

                    % reorder point in line order
                    if mod(overlap_idx,2) == 0
                        line_list{line_idx+1}=flipud(line_list{line_idx+1});
                    end
                end
            end
        end
    end

    methods(Static)
        function [U,Fval,node_list]=meshAdapt1D(fcn,low_bou,up_bou,torl,min_level,max_level,fval_num)
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
            [U,Fval]=traversalInorder(1); % from small to large get list

            % add boundary info
            U=[data_list(1,1);U;data_list(2,1)];
            Fval=[data_list(1,2:end);Fval;data_list(2,2:end)];

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

        function [U,V,Fval,node_list]=meshAdapt2DUV(fcn,low_bou,up_bou,torl,min_level,max_level,fval_num)
            % 2D Omni-Tree
            % adapt capture 2 dimemsion function value
            % ensure error of linear interplation will less than torl
            %
            if nargin < 7,fval_num=1;end

            % node_list which is a matrix store all node
            % a node is a array, contain:
            % level, idx_1-8(index of data_list), idx_c, node_index_1-4
            % place:
            % 3-8-4 or 3-4 or 3-8-4
            % 1-5-2    6-7    6-c-7
            %          1-2    1-5-2
            % node:
            % 1-2 or 3 or 3-4
            %        1    1-2
            % if children node is empty, left_index or right_index will be zero
            list_add_num=50; % list will be extend only when list is not enough
            node_list=zeros(list_add_num,14,'int64');
            % data_list use to sort all float data include coodinate, function value
            data_list=zeros(list_add_num,fval_num+2);

            % add vertex of cell into data_list first
            bou_1=[low_bou(1),low_bou(2)];
            data_list(1,:)=[bou_1,fcn(bou_1)];
            bou_2=[up_bou(1),low_bou(2)];
            data_list(2,:)=[bou_2,fcn(bou_2)];
            bou_3=[low_bou(1),up_bou(2)];
            data_list(3,:)=[bou_3,fcn(bou_3)];
            bou_4=[up_bou(1),up_bou(2)];
            data_list(4,:)=[bou_4,fcn(bou_4)];

            % create base root
            node_list(1,:)=[0,1,2,3,4,0,0,0,0,0,0,0,0,0];

            % create node tree from root
            [node_num,data_num]=createNodeTree(1,4);
            node_list=node_list(1:node_num,:);
            data_list=data_list(1:data_num,:);

            % generate U and V
            [u_list,~,u_idx]=unique(data_list(:,1));
            [v_list,~,v_idx]=unique(data_list(:,2));
            idx_list=[u_idx,v_idx];

            [U,V]=meshgrid(u_list,v_list);

            % local data to Fval
            Fval=nan(length(v_list),length(u_list),fval_num);

            for data_index=1:data_num
                if all([idx_list(data_index,2),idx_list(data_index,1)] ~= 0)
                    Fval(idx_list(data_index,2),idx_list(data_index,1),:)=data_list(data_index,3:end);
                end
            end

            % fit nan data
            for rank_idx=1:length(v_list)
                for colume_idx=1:length(u_list)
                    if isnan(Fval(rank_idx,colume_idx,1))
                        Fval(rank_idx,colume_idx,:)=fcn([U(rank_idx,colume_idx),V(rank_idx,colume_idx)]);
                    end
                end
            end

            function [node_num,data_num]=createNodeTree(root_idx,data_num)
                % create quad tree
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
                        [coord_c,coord_5,coord_6,coord_7,coord_8,...
                            fval_c,fval_5,fval_6,fval_7,fval_8]=calCell(node(2),node(3),node(4),node(5));
                        [fval_pred_c,fval_pred_5,fval_pred_6,...
                            fval_pred_7,fval_pred_8]=calCellPred(node(2),node(3),node(4),node(5));

                        % check u direction
                        if any(abs(fval_5-fval_pred_5) > torl) || any(abs(fval_8-fval_pred_8) > torl)
                            add_u_flag=true(1);
                        else
                            add_u_flag=false(1);
                        end

                        % check v direction
                        if any(abs(fval_6-fval_pred_6) > torl) || any(abs(fval_7-fval_pred_7) > torl)
                            add_v_flag=true(1);
                        else
                            add_v_flag=false(1);
                        end

                        % check center place
                        if any(abs(fval_c-fval_pred_c) > torl)
                            add_u_flag=true(1);
                            add_v_flag=true(1);
                        end

                        if (add_u_flag && add_v_flag) || node(1) < min_level
                            % add 5 data into data_list
                            data_new_idx=data_num+(1:5);
                            if data_num+5 > size(data_list,1)
                                data_list=[data_list;zeros(list_add_num,fval_num+2)];
                            end
                            data_list(data_new_idx,:)=[
                                coord_5,fval_5;
                                coord_6,fval_6;
                                coord_7,fval_7;
                                coord_8,fval_8;
                                coord_c,fval_c;];
                            node(6:10)=data_new_idx;
                            data_num=data_num+5;

                            % add 4 new node to node_list
                            node_new_idx=node_num+(1:4);
                            if node_num+4 > size(node_list,1)
                                node_list=[node_list;zeros(list_add_num,14)];
                            end
                            node_list(node_new_idx,:)=[...
                                node(1)+1,node(2),node(6),node(7),node(10),0,0,0,0,0,0,0,0,0;
                                node(1)+1,node(6),node(3),node(10),node(8),0,0,0,0,0,0,0,0,0;
                                node(1)+1,node(7),node(10),node(4),node(9),0,0,0,0,0,0,0,0,0;
                                node(1)+1,node(10),node(8),node(9),node(5),0,0,0,0,0,0,0,0,0;];
                            node(11:14)=node_new_idx;
                            node_num=node_num+4;

                        elseif add_u_flag
                            % add 2 data into data_list
                            data_new_idx=data_num+(1:2);
                            if data_num+2 > size(data_list,1)
                                data_list=[data_list;zeros(list_add_num,fval_num+2)];
                            end
                            data_list(data_new_idx,:)=[
                                coord_5,fval_5;
                                coord_8,fval_8;];
                            node([6,9])=data_new_idx;
                            data_num=data_num+2;

                            % add 2 new node to node_list
                            node_new_idx=node_num+(1:2);
                            if node_num+2 > size(node_list,1)
                                node_list=[node_list;zeros(list_add_num,14)];
                            end
                            node_list(node_new_idx,:)=[...
                                node(1)+1,node(2),node(6),node(4),node(9),0,0,0,0,0,0,0,0,0;
                                node(1)+1,node(6),node(3),node(9),node(5),0,0,0,0,0,0,0,0,0;];
                            node([11,12])=node_new_idx;
                            node_num=node_num+2;

                        elseif add_v_flag
                            % add 2 data into data_list
                            data_new_idx=data_num+(1:2);
                            if data_num+2 > size(data_list,1)
                                data_list=[data_list;zeros(list_add_num,fval_num+2)];
                            end
                            data_list(data_new_idx,:)=[
                                coord_6,fval_6;
                                coord_7,fval_7;];
                            node([7,8])=data_new_idx;
                            data_num=data_num+2;

                            % add 2 new node to node_list
                            node_new_idx=node_num+(1:2);
                            if node_num+2 > size(node_list,1)
                                node_list=[node_list;zeros(list_add_num,14)];
                            end
                            node_list(node_num+(1:2),:)=[...
                                node(1)+1,node(2),node(3),node(7),node(8),0,0,0,0,0,0,0,0,0;
                                node(1)+1,node(7),node(8),node(4),node(5),0,0,0,0,0,0,0,0,0;];
                            node([11,13])=node_new_idx;
                            node_num=node_num+2;
                        else
                            node_new_idx=[];
                        end

                        % add to stack
                        stack=[stack,node_new_idx];
                        node_list(node_idx,:)=node;
                    end
                end
            end

            function [coord_c,coord_5,coord_6,coord_7,coord_8,...
                    fval_c,fval_5,fval_6,fval_7,fval_8]=calCell(vidx_1,vidx_2,vidx_3,vidx_4)
                % calculate index of c, 5, 6, 7, 8 place function value
                % abbreviation:
                % vidx: vertex index
                %
                coord_c=(data_list(vidx_1,[1,2])+data_list(vidx_4,[1,2]))/2;
                fval_c=fcn(coord_c);

                coord_5=(data_list(vidx_1,[1,2])+data_list(vidx_2,[1,2]))/2;
                fval_5=fcn(coord_5);
                coord_6=(data_list(vidx_1,[1,2])+data_list(vidx_3,[1,2]))/2;
                fval_6=fcn(coord_6);
                coord_7=(data_list(vidx_2,[1,2])+data_list(vidx_4,[1,2]))/2;
                fval_7=fcn(coord_7);
                coord_8=(data_list(vidx_3,[1,2])+data_list(vidx_4,[1,2]))/2;
                fval_8=fcn(coord_8);
            end

            function [fval_pred_c,fval_pred_5,fval_pred_6,...
                    fval_pred_7,fval_pred_8]=calCellPred(vidx_1,vidx_2,vidx_3,vidx_4)
                % calculate index of c, 5, 6, 7, 8 place linear predict function value
                %
                fval_pred_c=(data_list(vidx_1,3:end)+data_list(vidx_2,3:end)+data_list(vidx_3,3:end)+data_list(vidx_4,3:end))/4;
                fval_pred_5=(data_list(vidx_1,3:end)+data_list(vidx_2,3:end))/2;
                fval_pred_6=(data_list(vidx_1,3:end)+data_list(vidx_3,3:end))/2;
                fval_pred_7=(data_list(vidx_2,3:end)+data_list(vidx_4,3:end))/2;
                fval_pred_8=(data_list(vidx_3,3:end)+data_list(vidx_4,3:end))/2;
            end

        end

    end
end