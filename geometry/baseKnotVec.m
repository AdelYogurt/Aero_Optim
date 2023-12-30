function knot_vec=baseKnotVec(Mults,Knots)
% base on Mults and Knots to create knot vector
%
if length(Mults) ~= length(Knots)
    error('baseKnotVec: length of Mults and Knots is not equal')
end
knot_vec=zeros(1,sum(Mults));
start_idx=1;end_idx=Mults(1);
for n_idx=1:length(Knots)
    knot_vec(start_idx:end_idx)=Knots(n_idx);
    start_idx=start_idx+Mults(n_idx);
    end_idx=end_idx+Mults(n_idx);
end
end