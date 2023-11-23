function [X,Y,Z]=geomMapGrid(line_1,line_2,line_3,line_4)
% generate discrete area by mapping
% using Lagrange polynomial mapping method
% local parameter u and v is equispaced
%
% input:
% line_1, line_2, line_3, line_4
%
% output:
% X,Y,Z
%
% notice:
% input line should connect each other,...
% order is end of line_1 equal to start of line_2...
%
geom_torl=100*eps;
[line_1,line_2,line_3,line_4]=correctLine(line_1,line_2,line_3,line_4,geom_torl);

u_num=size(line_1,1);
v_num=size(line_2,1);
if u_num ~= size(line_3,1) || v_num ~= size(line_4,1)
    error('geomMapGrid: opposite line of boundary discrete number no equal')
end
dimension=size(line_1,2);

X=zeros(v_num,u_num);Y=zeros(v_num,u_num);Z=zeros(v_num,u_num);

% preproces boundary
X(1,:)=line_1(:,1);X(v_num,:)=line_3(end:-1:1,1);
X(:,u_num)=line_2(:,1);X(:,1)=line_4(end:-1:1,1);
Y(1,:)=line_1(:,2);Y(v_num,:)=line_3(end:-1:1,2);
Y(:,u_num)=line_2(:,2);Y(:,1)=line_4(end:-1:1,2);
if dimension == 3
    Z(1,:)=line_1(:,3);Z(v_num,:)=line_3(end:-1:1,3);
    Z(:,u_num)=line_2(:,3);Z(:,1)=line_4(end:-1:1,3);
end

if u_num > 2 && v_num > 2

    h=1/(u_num-1);
    g=1/(v_num-1);

    % solve prepare
    % D_xx
    d_xx=spdiags([ones(v_num-2,1) -2*ones(v_num-2,1) ones(v_num-2,1)],-1:1,v_num-2,v_num-2);
    d_xx=d_xx/h^2;
    I_M=speye(u_num-2,u_num-2);
    D_xx=kron(I_M,d_xx);

    % D_yy
    d_yy=spdiags([ones(u_num-2,1) -2*ones(u_num-2,1) ones(u_num-2,1)],-1:1,u_num-2,u_num-2);
    d_yy=d_yy/g^2;
    I_K=speye(v_num-2,v_num-2);
    D_yy=kron(d_yy,I_K);

    % Laplace matrix
    Lap=(D_xx+D_yy);

    X_B=zeros(v_num-2,u_num-2);
    X_B(1,:)=X_B(1,:)-X(1,2:u_num-1)/g^2;
    X_B(end,:)=X_B(end,:)-X(end,2:u_num-1)/g^2;
    X_B(:,1)=X_B(:,1)-X(2:v_num-1,1)/h^2;
    X_B(:,end)=X_B(:,end)-X(2:v_num-1,end)/h^2;
    X(2:v_num-1,2:u_num-1)=reshape(Lap\X_B(:),v_num-2,u_num-2);

    Y_B=zeros(v_num-2,u_num-2);
    Y_B(1,:)=Y_B(1,:)-Y(1,2:u_num-1)/g^2;
    Y_B(end,:)=Y_B(end,:)-Y(end,2:u_num-1)/g^2;
    Y_B(:,1)=Y_B(:,1)-Y(2:v_num-1,1)/h^2;
    Y_B(:,end)=Y_B(:,end)-Y(2:v_num-1,end)/h^2;
    Y(2:v_num-1,2:u_num-1)=reshape(Lap\Y_B(:),v_num-2,u_num-2);

    if dimension == 3
        Z_B=zeros(v_num-2,u_num-2);
        Z_B(1,:)=Z_B(1,:)-Z(1,2:u_num-1)/g^2;
        Z_B(end,:)=Z_B(end,:)-Z(end,2:u_num-1)/g^2;
        Z_B(:,1)=Z_B(:,1)-Z(2:v_num-1,1)/h^2;
        Z_B(:,end)=Z_B(:,end)-Z(2:v_num-1,end)/h^2;
        Z(2:v_num-1,2:u_num-1)=reshape(Lap\Z_B(:),v_num-2,u_num-2);
    end

end

end

function [line_1,line_2,line_3,line_4]=correctLine(line_1,line_2,line_3,line_4,geom_torl)
% correct line point order
% order should be anti clockwise and start from line_1
%

% start from line_1
if norm(line_1(end,:)-line_2(1,:)) > geom_torl
    % search next connect point
    dist=vecnorm((line_1(end,:)-[...
        line_2(1,:);line_2(end,:);
        line_3(1,:);line_3(end,:);
        line_4(1,:);line_4(end,:)]),[],2);
    overlap_idx=find(dist < geom_torl);
    if ~any(overlap_idx)
        save();
        error('geomMapGrid: correctLine: line not connect');
    end

    % exchange line
    switch ceil(overlap_idx/2)+1
        case 3
            line_temp=line_2;
            line_2=line_3;
            line_3=line_temp;
        case 4
            line_temp=line_2;
            line_2=line_4;
            line_4=line_temp;
    end

    % reorder point in line order
    if mod(overlap_idx,2) == 0
        line_2=flipud(line_2);
    end
end

if norm(line_2(end,:)-line_3(1,:)) > geom_torl
    % search next connect point
    dist=vecnorm((line_2(end,:)-[...
        line_3(1,:);line_3(end,:);
        line_4(1,:);line_4(end,:)]),[],2);
    overlap_idx=find(dist < geom_torl);
    if length(overlap_idx) ~= 1
        save();
        error('geomMapGrid: correctLine: line not connect');
    end

    % exchange line
    switch ceil(overlap_idx/2)+1
        case 4
            line_temp=line_3;
            line_3=line_4;
            line_4=line_temp;
    end

    % reorder point in line order
    if mod(overlap_idx,2) == 0
        line_3=flipud(line_3);
    end
end

if norm(line_3(end,:)-line_4(1,:)) > geom_torl
    line_4=flipud(line_4);
end

if norm(line_4(end,:)-line_1(1,:)) > geom_torl
    save();
    error('geomMapGrid: correctLine: line not connect');
end
end
