function [index_low,index_up,point_low,point_up]=dividePoint(index,point,dimension,threshold)
% divide point up and down by threshold
%
if nargin < 4
    threshold=0;
end
Bool=point(:,dimension) >= threshold;

index_up=index(Bool,:);
index_low=index(~Bool,:);

point_up=point(Bool,:);
point_low=point(~Bool,:);

end