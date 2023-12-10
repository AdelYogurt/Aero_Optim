function [X,Y,Z]=geomMapGrid(line_u0,line_u1,line_v0,line_v1,U,V)
% generate discrete area by mapping
% using Lagrange polynomial mapping method
% local parameter u and v is equispaced
%
% input:
% line_u0, line_u1, line_v0, line_v1
%
% output:
% X,Y,Z
%
if nargin < 6
    V=[];
    if nargin < 5
        U=[];
    end
end

geom_torl=100*eps;
line_list={line_u0,line_v1,flipud(line_u1),flipud(line_v0)};
line_list=correctLine(line_list,geom_torl);
line_u0=line_list{1};
line_v1=line_list{2};
line_u1=flipud(line_list{3});
line_v0=flipud(line_list{4});

u_num=size(line_u0,1);v_num=size(line_u1,1);
if u_num ~= size(line_v0,1) || v_num ~= size(line_v1,1)
    error('geomMapGrid: opposite line of boundary discrete number no equal')
end
dimension=size(line_u0,2);
if isempty(U),U=linspace(0,1,u_num);end
if isempty(V),V=linspace(0,1,v_num);end

X=zeros(v_num,u_num);Y=zeros(v_num,u_num);Z=zeros(v_num,u_num);

% preproces boundary
X(1,:)=line_u0(:,1);X(end,:)=line_u1(:,1);
X(:,1)=line_v0(:,1);X(:,end)=line_v1(:,1);
Y(1,:)=line_u0(:,2);Y(end,:)=line_u1(:,2);
Y(:,1)=line_v0(:,2);Y(:,end)=line_v1(:,2);
if dimension == 3
    Z(1,:)=line_u0(:,3);Z(end,:)=line_u1(:,3);
    Z(:,1)=line_v0(:,3);Z(:,end)=line_v1(:,3);
end

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

    X_B=zeros(v_num-2,u_num-2);
    X_B(1,:)=X_B(1,:)-X(1,2:u_num-1)/G(2)^2;
    X_B(end,:)=X_B(end,:)-X(end,2:u_num-1)/G(end-1)^2;
    X_B(:,1)=X_B(:,1)-X(2:v_num-1,1)/H(2)^2;
    X_B(:,end)=X_B(:,end)-X(2:v_num-1,end)/H(end-1)^2;
    X(2:v_num-1,2:u_num-1)=reshape(Lap\X_B(:),v_num-2,u_num-2);

    Y_B=zeros(v_num-2,u_num-2);
    Y_B(1,:)=Y_B(1,:)-Y(1,2:u_num-1)/G(2)^2;
    Y_B(end,:)=Y_B(end,:)-Y(end,2:u_num-1)/G(end-1)^2;
    Y_B(:,1)=Y_B(:,1)-Y(2:v_num-1,1)/H(2)^2;
    Y_B(:,end)=Y_B(:,end)-Y(2:v_num-1,end)/H(end-1)^2;
    Y(2:v_num-1,2:u_num-1)=reshape(Lap\Y_B(:),v_num-2,u_num-2);

    if dimension == 3
        Z_B=zeros(v_num-2,u_num-2);
        Z_B(1,:)=Z_B(1,:)-Z(1,2:u_num-1)/G(2)^2;
        Z_B(end,:)=Z_B(end,:)-Z(end,2:u_num-1)/G(end-1)^2;
        Z_B(:,1)=Z_B(:,1)-Z(2:v_num-1,1)/H(2)^2;
        Z_B(:,end)=Z_B(:,end)-Z(2:v_num-1,end)/H(end-1)^2;
        Z(2:v_num-1,2:u_num-1)=reshape(Lap\Z_B(:),v_num-2,u_num-2);
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
