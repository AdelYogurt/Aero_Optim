function knot_vec=baseKnotVec(knot_multi,knot_list)
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