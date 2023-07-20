function [index_front,index_back,point_front,point_back]=divideFB(index,point,threshold)
% divide point front and back by threshold
%
if nargin < 3
    threshold=0;
end
Bool=point(:,1) >= threshold;

index_front=index(Bool,:);
index_back=index(~Bool,:);

point_front=point(Bool,:);
point_back=point(~Bool,:);

end