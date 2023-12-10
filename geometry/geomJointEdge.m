function edge=geomJointEdge(name,edge_list,geom_torl)
% connect BSpline curve into one curve
%
if nargin < 3
    geom_torl=[];
end

if isempty(geom_torl),geom_torl=100*eps;end

% check if control point is connect
edge_list=correctEdge(edge_list,geom_torl);

% search max degree
edge_num=length(edge_list);
degree=0;
for edge_idx=1:edge_num
    edge_curr=edge_list{edge_idx};
    if degree < edge_curr.degree
        degree=edge_curr.degree;
    end
end

% increase degree of edge and load data
ctrl_X=[];
ctrl_Y=[];
ctrl_Z=[];
knot_multi=[];
knot_list=[];
for edge_idx=1:edge_num
    edge=edge_list{edge_idx};
    edge.changeDegree(degree);
    ctrl_X=[ctrl_X(1:end-1);edge.ctrl_X];
    ctrl_Y=[ctrl_Y(1:end-1);edge.ctrl_Y];
    ctrl_Z=[ctrl_Z(1:end-1);edge.ctrl_Z];
    knot_multi_new=edge.knot_multi;
    knot_multi_new(1)=knot_multi_new(1)-1;
    knot_multi=[knot_multi(1:end-1),knot_multi_new];
    knot_list=[knot_list(1:end-1),edge.knot_list+edge_idx-1];
end
knot_multi(1)=knot_multi(1)+1;

% connect
knot_list=knot_list/max(knot_list);

edge=EdgeBSpline(name,ctrl_X,ctrl_Y,ctrl_Z,[],degree,knot_multi,knot_list);
end

function edge_list=correctEdge(edge_list,geom_torl)
% correct line point order
% order should be anti clockwise and start from first line
%
edge_num=length(edge_list);

% start from first edge
for edge_idx=1:edge_num-1
    edge_curr=edge_list{edge_idx};
    ctrl_curr=[edge_curr.ctrl_X,edge_curr.ctrl_Y,edge_curr.ctrl_Z];
    edge_next=edge_list{edge_idx+1};
    ctrl_next=[edge_next.ctrl_X,edge_next.ctrl_Y,edge_next.ctrl_Z];
    if norm(ctrl_curr(end,:)-ctrl_next(1,:)) > geom_torl
        % load remain line all vertex
        vertex_list=zeros((edge_num-edge_idx)*2,3);
        for remain_idx=edge_idx+1:edge_num
            edge_rema=edge_list{remain_idx};
            ctrl_rema=[edge_rema.ctrl_X,edge_rema.ctrl_Y,edge_rema.ctrl_Z];
            vertex_list(2*(remain_idx-edge_idx)-1,:)=ctrl_rema(1,:);
            vertex_list(2*(remain_idx-edge_idx),:)=ctrl_rema(end,:);
        end

        % search next connect point
        dist=vecnorm((ctrl_curr(end,:)-vertex_list),2,2);
        overlap_idx=find(dist < geom_torl);
        if ~any(overlap_idx)
            save();
            error('geomMapGrid: correctLine: line not connect');
        end
        exchange_idx=ceil(overlap_idx/2)+edge_idx;

        % exchange line
        line_temp=edge_list{exchange_idx};
        edge_list{exchange_idx}=edge_next;
        edge_list{edge_idx+1}=line_temp;

        % reorder point in line order
        if mod(overlap_idx,2) == 0
            edge_list{edge_idx+1}=edge_list{edge_idx+1}.reverse;
        end
    end
end
end