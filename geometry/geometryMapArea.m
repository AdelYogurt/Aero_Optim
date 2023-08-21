function [X,Y]=geometryMapArea(line_list)
% generate discrete area(include all discrete point on area) by mapping
% line1, line2 in line_list... is discrete segment point list[x1,y1;x2,y2...]
% four point of area is base, form the xi, eta coordinate
% and disturbance xi, eta from boundary(line1, line2)
%

geometry_torlance=1e-4;

correctLine();

line_list{3}=line_list{3}(end:-1:1,:);
line_list{4}=line_list{4}(end:-1:1,:);

xi_num=size(line_list{1},1);
eta_num=size(line_list{2},1);
if xi_num ~= size(line_list{3},1) || eta_num ~= size(line_list{4},1)
    error('geometryMapArea: opposite line of boundary discrete number no equal')
end

point_base=[line_list{1}(1,:);line_list{2}(1,:);line_list{3}(end,:);line_list{4}(end,:)];
disturbance=cell(1,4);
% calclulate each line disturbance
point_index_list=[1,2;2,3;4,3;1,4];
for line_index=1:4
    % a=y2*x1-x2*y1, b=y1-y2, c=x2-x1
    bc=[
        point_base(point_index_list(line_index,1),2)-point_base(point_index_list(line_index,2),2);
        point_base(point_index_list(line_index,2),1)-point_base(point_index_list(line_index,1),1)];
    a=point_base(point_index_list(line_index,2),2)*point_base(point_index_list(line_index,1),1)-...
        point_base(point_index_list(line_index,2),1)*point_base(point_index_list(line_index,1),2);
    disturbance{line_index}=(a+line_list{line_index}*bc)/...
        sqrt(bc(1)^2+bc(2)^2);
end

% xi, eta to global coordinate is
% x=(1-xi)*(1-eta)*x1+xi*(1-eta)*x2+xi*eta*x3+(1-xi)*eta*x4;
% base on Jacobian determinant
% xi, eta point xi, eta base vector is
% dx_dxi=-(1-eta)*x1+(1-eta)*x2+eta*x3-eta*x4; x
% dy_dxi=-(1-eta)*y1+(1-eta)*y2+eta*y3-eta*y4; y
% dx_deta=-(1-xi)*x1-xi*x2+xi*x3+(1-xi)*x4; x
% dy_deta=-(1-xi)*y1-xi*y2+xi*y3+(1-xi)*y4; y

disturbance{1}=disturbance{1}';
disturbance{3}=disturbance{3}';

xi_list=linspace(0,1,xi_num);
eta_list=linspace(0,1,eta_num);

[XI,ETA]=meshgrid(xi_list,eta_list);

% base location
X=...
    (1-XI).*(1-ETA)*point_base(1,1)+...
    XI.*(1-ETA)*point_base(2,1)+...
    XI.*ETA*point_base(3,1)+...
    (1-XI).*ETA*point_base(4,1);
Y=...
    (1-XI).*(1-ETA)*point_base(1,2)+...
    XI.*(1-ETA)*point_base(2,2)+...
    XI.*ETA*point_base(3,2)+...
    (1-XI).*ETA*point_base(4,2);

dX_dXI=-(1-ETA)*point_base(1,1)+(1-ETA)*point_base(2,1)+...
    ETA*point_base(3,1)-ETA*point_base(4,1);
dY_dXI=-(1-ETA)*point_base(1,2)+(1-ETA)*point_base(2,2)+...
    ETA*point_base(3,2)-ETA*point_base(4,2);
dX_dETA=-(1-XI)*point_base(1,1)-XI*point_base(2,1)+...
    XI*point_base(3,1)+(1-XI)*point_base(4,1);
dY_dETA=-(1-XI)*point_base(1,2)-XI*point_base(2,2)+...
    XI*point_base(3,2)+(1-XI)*point_base(4,2);

norm_dXI=sqrt(dX_dXI.^2+dY_dXI.^2);
norm_dETA=sqrt(dX_dETA.^2+dY_dETA.^2);

dX_dXI=dX_dXI./norm_dXI;
dY_dXI=dY_dXI./norm_dXI;
dX_dETA=dX_dETA./norm_dETA;
dY_dETA=dY_dETA./norm_dETA;

% disturbance by boundary
X=X+...
    -dX_dXI.*(disturbance{4}.*(1-XI)+disturbance{2}.*XI)+...
    dX_dETA.*(disturbance{1}.*(1-ETA)+disturbance{3}.*ETA);

X(1,:)=line_list{1}(:,1)';
X(:,end)=line_list{2}(:,1);
X(end,:)=line_list{3}(:,1)';
X(:,1)=line_list{4}(:,1);

Y=Y+...
    -dY_dXI.*(disturbance{4}.*(1-XI)+disturbance{2}.*XI)+...
    dY_dETA.*(disturbance{1}.*(1-ETA)+disturbance{3}.*ETA);

Y(1,:)=line_list{1}(:,2)';
Y(:,end)=line_list{2}(:,2);
Y(end,:)=line_list{3}(:,2)';
Y(:,1)=line_list{4}(:,2);

    function correctLine()
       % correct line point order(anti clockwise from line1 order)
       %
       
       % start from line1
       if norm(line_list{1}(end,:)-line_list{2}(1,:)) > geometry_torlance
           % search next connect point
           distance__=sum(abs(line_list{1}(end,:)-[...
               line_list{2}(1,:);line_list{2}(end,:);
               line_list{3}(1,:);line_list{3}(end,:);
               line_list{4}(1,:);line_list{4}(end,:)]),2);
           line_index__=find(distance__ < geometry_torlance);
           if length(line_index__) ~= 1
               save();
               error('geometryMapArea: correctLine: line not connect');
           end
           
           % exchange line
           temp__=line_list{1+ceil(line_index__/2)};
           line_list{1+ceil(line_index__/2)}=line_list{2};
           line_list{2}=temp__;
           
           % reorder point in line order
           if mod(index,2) == 0
               line_list{2}=flipud(line_list{2});
           end
       end
       
       if norm(line_list{2}(end,:)-line_list{3}(1,:)) > geometry_torlance
           % search next connect point
           distance__=sum(abs(line_list{2}(end,:)-[...
               line_list{3}(1,:);line_list{3}(end,:);
               line_list{4}(1,:);line_list{4}(end,:)]),2);
           line_index__=find(distance__ < geometry_torlance);
           if length(line_index__) ~= 1
               save();
               error('geometryMapArea: correctLine: line not connect');
           end
           
           % exchange line
           temp__=line_list{2+ceil(line_index__/2)};
           line_list{2+ceil(line_index__/2)}=line_list{3};
           line_list{3}=temp__;
           
           % reorder point in line order
           if mod(line_index__,2) == 0
               line_list{3}=flipud(line_list{3});
           end
       end
       
       if norm(line_list{3}(end,:)-line_list{4}(1,:)) > geometry_torlance
           line_list{4}=flipud(line_list{4});
       end
       
       if norm(line_list{4}(end,:)-line_list{1}(1,:)) > geometry_torlance
           save();
           error('geometryMapArea: correctLine: line not connect');
       end
    end
end
