function [index_up,index_low,point_up,point_low]=divideUL(index,point,threshold)
% divide point up and down by threshold
%
if nargin < 3
    threshold=0;
end
Bool=point(:,3) >= threshold;

index_up=index(Bool,:);
index_low=index(~Bool,:);

point_up=point(Bool,:);
point_low=point(~Bool,:);

end