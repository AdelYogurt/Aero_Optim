classdef GeomBSpline
    methods(Static)
        function Pnts=bspeval(Ctrls,Degree,u_list,U_x)
            % evaluate B-Spline at control points.
            %
            % input:
            % d-Degree of the B-Spline.
            % c-Control Points, matrix of size (dim,nc).
            % k-Knot sequence, row vector of size nk.
            % u-Parametric evaluation points, row vector of size nu.
            %
            % output:
            % p-Evaluated points, matrix of size (dim,nu)
            %
            U_x=U_x(:);
            u_num=length(U_x);
            [ctrl_num,dim]=size(Ctrls);
            idx_end=findspan(ctrl_num-1,Degree,U_x,u_list);
            N_list=baseFcnN(U_x,Degree,u_list,idx_end);
            idx_start=idx_end-Degree;
            Pnts=zeros(u_num,dim);
            for deg_idx=1:Degree+1
                Pnts=Pnts+N_list(:,deg_idx).*Ctrls(idx_start+deg_idx,:);
            end
        end
    end
end

function s=findspan(n,p,U_x,u_list)
s=zeros(size(U_x));
for j=1:numel(U_x)
    if (U_x(j)==u_list(n+2)), s(j)=n; continue, end
    s(j)=find(U_x(j) >= u_list,1,'last')-1;
end
end

function N_list=baseFcnN(U_x,Degree,u_list,u_idx)
% Basis function for B-Spline
%
% input:
% iv: knot span
% uv: parametric points
% p: spline degree
% U: knot sequence
%
% output:
% N-Basis functions vector(numel(uv)*(p+1))
%

u_num=length(U_x);
N_list=zeros(u_num, Degree+1);
for jj=1:u_num
    i=u_idx(jj) + 1; %% findspan uses 0-based numbering
    u=U_x(jj);
    left=zeros(Degree+1,1);
    right=zeros(Degree+1,1);
    N(1)=1;
    for j=1:Degree
        left(j+1)=u-u_list(i+1-j);
        right(j+1)=u_list(i+j)-u;
        saved=0;
        for r=0:j-1
            temp=N(r+1)/(right(r+2) + left(j-r+1));
            N(r+1)=saved + right(r+2)*temp;
            saved=left(j-r+1)*temp;
        end
        N(j+1)=saved;
    end
    N_list (jj, :)=N;
end
end
