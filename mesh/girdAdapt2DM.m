function [XI,PSI,Fval,node_list]=girdAdapt2DM(fcn,low_bou,up_bou,torl,max_level)
% adapt capture 2 dimemsion function value
% ensure error of linear interplation will less than torl
%

% node_list which is a matrix store all node
% a node is a array, contain level, x, fval_c, fval_l, fval_u, node_low_bou_index, node_up_bou_index, left_index, right_index
% if children node is empty, left_index or right_index will be zero
node_add_num=50; % node_list will be extend only when node_list is not enough
node_list=zeros(node_add_num,9,2); % x tree, y_tree

% add low_bou data and up_bou data
% root will start from index 3
fval_bou_1=fcn([low_bou(1),low_bou(2)]);
fval_bou_2=fcn([low_bou(1),up_bou(2)]);
fval_bou_3=fcn([up_bou(1),low_bou(2)]);
fval_bou_4=fcn([up_bou(1),up_bou(2)]);

node_list(1,:,1)=[0,low_bou(1),fval_bou_1,0,fval_bou_3,0,0,0,0];
node_list(2,:,1)=[0,up_bou(1),fval_bou_2,0,fval_bou_4,0,0,0,0];
node_list(1,:,2)=[0,low_bou(2),fcn(up_bou),0,0,0,0];
node_list(2,:,2)=[0,up_bou(2),fcn(up_bou),0,0,0,0];

x=(low_bou+up_bou)/2;
fval=fcn(x);
level=0;
node_list(3,:)=[level,x,fval,1,2,0,0]; % root

node_num=createNodeTree(3); % create node tree from root
node_list=node_list(1:node_num,:);
[X,Fval]=traversalInorder(3); % from small to large get list

% add boundary info
X=[node_list(1,2),X,node_list(2,2)];
Fval=[node_list(1,3),Fval,node_list(2,3)];

    function node_num=createNodeTree(root_idx)
        %
        %
        stack=root_idx;
        node_num=root_idx;

        while ~isempty(stack)
            % current node information
            node_idx=stack(end);
            node=node_list(node_idx,:);
            stack=stack(1:end-1);

            % check left
            x_left=(node_list(node(4),2)+node(2))/2;
            fval_left=fcn(x_left);
            fval_left_pred=(node_list(node(4),3)+node(3))/2;
            if abs(fval_left-fval_left_pred) > torl && node(1) < max_level
                % create new node
                node_left_index=node_num+1;
                node_left=[node(1)+1,x_left,fval_left,node(4),node_idx,0,0];

                if node_left_index > size(node_list,1)
                    node_list=[node_list;zeros(node_add_num,7)];
                end
                node_list(node_left_index,:)=node_left;
                node_num=node_num+1;

                stack=[stack,node_left_index];
            else
                node_left_index=0;
            end

            % check right
            x_right=(node_list(node(5),2)+node(2))/2;
            fval_right=fcn(x_right);
            fval_right_pred=(node_list(node(5),3)+node(3))/2;
            if abs(fval_right-fval_right_pred) > torl && node(1) < max_level
                % create new node
                node_right_index=node_num+1;
                node_right=[node(1)+1,x_right,fval_right,node_idx,node(5),0,0];

                if node_right_index > size(node_list,1) 
                    node_list=[node_list;zeros(node_add_num,7)];
                end
                node_list(node_right_index,:)=node_right;
                node_num=node_num+1;

                stack=[stack,node_right_index];
            else
                node_right_index=0;
            end

            % add children index into current node data
            node_list(node_idx,6:7)=[node_left_index,node_right_index];
        end
    end

    function [x_list,fval_list]=traversalInorder(root_idx)
        stack=[];
        x_list=[];
        fval_list=[];
        node_idx=root_idx;

        while ~isempty(stack) || node_idx ~= 0
            while node_idx ~= 0
                stack=[stack,node_idx];

                % node=node.left;
                node_idx=node_list(node_idx,6);
            end

            node_idx=stack(end);
            stack=stack(1:end-1);
            x_list=[x_list,node_list(node_idx,2)];
            fval_list=[fval_list,node_list(node_idx,3)];

            % node=node.right;
            node_idx=node_list(node_idx,7);
        end
    end
end
