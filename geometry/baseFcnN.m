function [N_list,idx_srt,idx_end]=baseFcnN(U_x,k,u_list)
% Basis function for B-Spline
%
% input:
% U_x: parametric points
% k: spline degree
% u_list: knot sequence
%
% output:
% N-Basis functions vector(numel(U_x)*(k+1))
%

idx_end=zeros(size(U_x));
for j=1:numel(U_x)
    if (U_x(j)==u_list(k+1)),idx_end(j)=k+1; continue,end
    idx_end(j)=find(U_x(j) <= u_list,1,'first')-1;
end
idx_srt=idx_end-k;

u_num=length(U_x);
N_list=zeros(u_num,k+1);
N=zeros(1,k+1);
for jj=1:u_num
    i=idx_end(jj); %% findspan uses 0-based numbering
    u=U_x(jj);
    left=zeros(k+1,1);
    right=zeros(k+1,1);
    N(1)=1;
    for j=1:k
        left(j+1)=u-u_list(i+1-j);
        right(j+1)=u_list(i+j)-u;
        saved=0;
        for r=0:j-1
            temp=N(r+1)/(right(r+2)+left(j-r+1));
            N(r+1)=saved+right(r+2)*temp;
            saved=left(j-r+1)*temp;
        end
        N(j+1)=saved;
    end
    N_list(jj,:)=N;
end
end

%% origin BSpline base function

% function N=baseFcnN(u_list,u_x,i,k)
% % base function of BSpline curve
% %
% if k == 0
%     if ((u_list(i) <= u_x) && (u_x <= u_list(i+1)))
%         if any(u_list == u_x) && u_x ~= u_list(1) && u_x ~= u_list(end) % avoid start and end knot repeat
%             N=0.5;
%         else
%             N=1;
%         end
%     else
%         N=0;
%     end
% else
%     if u_list(i+k) == u_list(i)
%         A=0;
%     else
%         A=(u_x-u_list(i))/(u_list(i+k)-u_list(i));
%     end
% 
%     if u_list(i+k+1) == u_list(i+1)
%         B=0;
%     else
%         B=(u_list(i+k+1)-u_x)/(u_list(i+k+1)-u_list(i+1));
%     end
% 
%     N=A*baseFcnN(u_list,u_x,i,k-1)+B*baseFcnN(u_list,u_x,i+1,k-1);
% end
% end
