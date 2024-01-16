classdef GeomApp
    methods(Static)
        function edg=VertexToEdge(name,Nodes,Degree,pole_num,U_node)
            % generate BSpline curve by defined fitting points
            %
            % input:
            % name:
            % Nodes (matrix): fit point,node_num x dimension matrix
            % Degree (optional):
            % Mults (optional):
            % Knots (optional):
            % pole_num (optional):
            % U_node (optional):
            %
            % output:
            % EdgeNURBS
            %
            % notice:
            % if input Degree is empty,Degree default is pole_num-1,which is Bezier curve
            %
            if nargin < 5
                U_node=[];
                if nargin < 4
                    pole_num=[];
                    if nargin < 3
                        Degree=[];
                    end
                end
            end

            node_num=size(Nodes,1);
            if isempty(pole_num),pole_num=node_num;end
            if pole_num > node_num
                error('GeomApp.VertexToEdge: pole_num more than node_num')
            end
            % default value
            if isempty(U_node)
                U_node=vecnorm(Nodes(2:end,:)-Nodes(1:end-1,:),2,2);
                U_node=[0;cumsum(U_node)];U_node=U_node/U_node(end);
            end
            U_node=U_node(:);

            % reference:
            % [1] 施法中. 计算机辅助几何设计与非均匀有理B样条[M]. 208-281.
            % [2] https://blog.csdn.net/he_nan/article/details/107134271
            Mults=[Degree+1,ones(1,pole_num-Degree-1),Degree+1];
            Knots=linspace(0,1,node_num-Degree+1);
            for j=2:node_num-Degree
                Knots(j)=mean(U_node(j:j+Degree-1));
            end
            % modify
            Knots=interp1(linspace(0,1,node_num-Degree+1),Knots,linspace(0,1,pole_num-Degree+1));
            % c=(node_num)/(pole_num-Degree);
            % for j=1:pole_num-Degree-1
            %     i=floor(j*c);
            %     a=(j*c-i);
            %     Knots(j+1)=(1-a)*u_node(i)+a*u_node(i+1);
            % end
            u_list=baseKnotVec(Mults,Knots);

            % base on node point list inverse calculate control point list
            fit_matrix=zeros(node_num,pole_num);
            for node_idx=1:node_num
                u=U_node(node_idx);
                for ctrl_idx=1:pole_num
                    fit_matrix(node_idx,ctrl_idx)=baseFcnN(u_list,u,ctrl_idx,Degree);
                end
            end
            Poles=fit_matrix\Nodes;

            edg=EdgeNURBS(name,Poles,Degree,Mults,Knots);
        end

        function fce=VertexToFace(name,Nodes,UDegree,VDegree,u_pole_num,v_pole_num,U_node,V_node)
            % generate BSpline curve by defined fitting points
            %
            % input:
            % name:
            % Nodes (matrix): fit point,u_node_num x v_node_num x dimension matrix
            % Degree (optional):
            % Mults (optional):
            % Knots (optional):
            % u_pole_num (optional):
            % v_pole_num (optional):
            % U_node (optional):
            % V_node (optional):
            %
            % output:
            % FaceNURBS
            %
            % notice:
            % if input Degree is empty,Degree default is pole_num-1,which is Bezier curve
            %
            if nargin < 8
                V_node=[];
                if nargin < 7
                    U_node=[];
                    if nargin < 6
                        v_pole_num=[];
                        if nargin < 5
                            u_pole_num=[];
                            if nargin < 4
                                VDegree=[];
                                if nargin < 3
                                    UDegree=[];
                                end
                            end
                        end
                    end
                end
            end

            [v_node_num,u_node_num,dimension]=size(Nodes);
            if isempty(u_pole_num),u_pole_num=u_node_num;end
            if isempty(v_pole_num),v_pole_num=v_node_num;end
            if u_pole_num > u_node_num || v_pole_num > v_node_num
                error('GeomApp.VertexToEdge: pole_num more than node_num')
            end

            % default value
            if isempty(U_node)
                U_node=vecnorm(Nodes(:,2:end,:)-Nodes(:,1:end-1,:),2,3);
                U_node=mean(U_node,1);
                U_node=[0;cumsum(U_node')];U_node=U_node/U_node(end);
            end
            U_node=U_node(:);
            if isempty(V_node)
                V_node=vecnorm(Nodes(2:end,:,:)-Nodes(1:end-1,:,:),2,3);
                V_node=mean(V_node,2);
                V_node=[0;cumsum(V_node)];V_node=V_node/V_node(end);
            end
            V_node=V_node(:);
            
            UMults=[UDegree+1,ones(1,u_pole_num-UDegree-1),UDegree+1];
            UKnots=linspace(0,1,u_node_num-UDegree+1);
            for j=2:u_node_num-UDegree
                UKnots(j)=mean(U_node(j:j+UDegree-1));
            end
            % modify
            UKnots=interp1(linspace(0,1,u_node_num-UDegree+1),UKnots,linspace(0,1,u_pole_num-UDegree+1));
            u_list=baseKnotVec(UMults,UKnots);

            VMults=[VDegree+1,ones(1,v_pole_num-VDegree-1),VDegree+1];
            VKnots=linspace(0,1,v_node_num-VDegree+1);
            for j=2:v_node_num-VDegree
                VKnots(j)=mean(V_node(j:j+VDegree-1));
            end
            % modify
            VKnots=interp1(linspace(0,1,v_node_num-VDegree+1),VKnots,linspace(0,1,v_pole_num-VDegree+1));
            v_list=baseKnotVec(VMults,VKnots);

            % base on node point list inverse calculate control point list
            matrix_v=zeros(v_node_num,v_pole_num);
            for node_idx=1:v_node_num
                v=V_node(node_idx);
                for ctrl_idx=1:v_pole_num
                    matrix_v(node_idx,ctrl_idx)=baseFcnN(v_list,v,ctrl_idx,VDegree);
                end
            end
            matrix_u=zeros(u_pole_num,u_node_num);
            for node_idx=1:u_node_num
                u=U_node(node_idx);
                for ctrl_idx=1:u_pole_num
                    matrix_u(ctrl_idx,node_idx)=baseFcnN(u_list,u,ctrl_idx,UDegree);
                end
            end
            Poles=pagemrdivide(pagemldivide(matrix_v,Nodes),matrix_u);

            fce=FaceNURBS(name,Poles,UDegree,VDegree,UMults,VMults,UKnots,VKnots);
        end

        function fce=EdgeToFace(name,edg_u0,edg_u1,edg_0v,edg_1v)
            % mapping edge to generate face
            %
            if nargin < 5
                edg_1v=[];
                if nargin < 4
                    edg_0v=[];
                end
            end

            if (~isempty(edg_u0) && ~isempty(edg_u1)) || (~isempty(edg_0v) && ~isempty(edg_1v))
                if ~isempty(edg_u0) && ~isempty(edg_u1)
                    [edg_u0,edg_u1,UDegree,UMults,UKnots]=GeomApp.matchEdge(edg_u0,edg_u1);
                    Ctrl_u0=[edg_u0.Poles,edg_u0.Weights];
                    Ctrl_u1=[edg_u1.Poles,edg_u1.Weights];
                end

                if ~isempty(edg_0v) && ~isempty(edg_1v)
                    [edg_0v,edg_1v,VDegree,VMults,VKnots]=GeomApp.matchEdge(edg_0v,edg_1v);
                    Ctrl_0v=[edg_0v.Poles,edg_0v.Weights];
                    Ctrl_1v=[edg_1v.Poles,edg_1v.Weights];
                end

                if isempty(edg_0v) && isempty(edg_1v)
                    Ctrl_0v=[Ctrl_u0(1,:);Ctrl_u1(1,:)];
                    Ctrl_1v=[Ctrl_u0(end,:);Ctrl_u1(end,:)];
                    VDegree=1;VMults=[2,2];VKnots=[0,1];
                elseif isempty(edg_u0) && isempty(edg_u1)
                    Ctrl_u0=[Ctrl_0v(1,:);Ctrl_1v(1,:)];
                    Ctrl_u1=[Ctrl_0v(end,:);Ctrl_1v(end,:)];
                    UDegree=1;UMults=[2,2];UKnots=[0,1];
                elseif isempty(edg_0v) || isempty(edg_1v)
                    edg_v=[edg_0v,edg_1v];
                    VDegree=edg_v.Degree;
                    VMults=edg_v.Mults;
                    VKnots=edg_v.Knots;

                    if ~isempty(edg_0v)
                        Ctrl_0v=[edg_0v.Poles,edg_0v.Weights];
                        % generate new ctrl
                        Ctrl_1v=getMatchCtrl(Ctrl_u0(end,:),Ctrl_u1(end,:),VDegree,VMults,VKnots);
                    elseif ~isempty(edg_1v)
                        Ctrl_1v=[edg_1v.Poles,edg_1v.Weights];
                        % generate new ctrl
                        Ctrl_0v=getMatchCtrl(Ctrl_u0(1,:),Ctrl_u1(1,:),VDegree,VMults,VKnots);
                    end
                elseif isempty(edg_u0) || isempty(edg_u1)
                    edg_u=[edg_u0,edg_u1];
                    UDegree=edg_u.Degree;
                    UMults=edg_u.Mults;
                    UKnots=edg_u.Knots;

                    if ~isempty(edg_u0)
                        Ctrl_u0=[edg_u0.Poles,edg_u0.Weights];
                        % generate new ctrl
                        Ctrl_u1=getMatchCtrl(Ctrl_0v(end,:),Ctrl_1v(end,:),UDegree,UMults,UKnots);
                    elseif ~isempty(edg_u1)
                        Ctrl_u1=[edg_u1.Poles,edg_u1.Weights];
                        % generate new ctrl
                        Ctrl_u0=getMatchCtrl(Ctrl_0v(1,:),Ctrl_1v(1,:),UDegree,UMults,UKnots);
                    end
                end
            else
                error('GeomApp.EdgeToFace: error input edge');
            end
            u_list=baseKnotVec(UMults,UKnots);
            v_list=baseKnotVec(VMults,VKnots);
            u_pole_num=size(Ctrl_u0,1);v_pole_num=size(Ctrl_0v,1);
            U_pole=interp1(linspace(0,1,u_pole_num-UDegree+1),u_list(UDegree+1:u_pole_num+1),linspace(0,1,u_pole_num));
            V_pole=interp1(linspace(0,1,v_pole_num-VDegree+1),v_list(VDegree+1:v_pole_num+1),linspace(0,1,v_pole_num));
            [U,V]=meshgrid(U_pole,V_pole);
            % Ctrl=GeomApp.MapGrid(Ctrl_u0,Ctrl_u1,Ctrl_0v,Ctrl_1v,U_pole,V_pole);

%             UPoles=cat(3,[Ctrl_u0(:,1)';Ctrl_u1(:,1)'],[Ctrl_u0(:,2)';Ctrl_u1(:,2)'],[Ctrl_u0(:,3)';Ctrl_u1(:,3)']);
%             UWeights=[Ctrl_u0(:,4)';Ctrl_u1(:,4)'];
%             Ufce=FaceNURBS('',UPoles,UDegree,[],UMults,[],UKnots,[],UWeights);
%             
%             VPoles=cat(3,[Ctrl_0v(:,1);Ctrl_1v(:,1)],[Ctrl_0v(:,2);Ctrl_1v(:,2)],[Ctrl_0v(:,3);Ctrl_1v(:,3)]);
%             VWeights=[Ctrl_0v(:,4),Ctrl_1v(:,4)];
%             Vfce=FaceNURBS('',VPoles,[],VDegree,[],VMults,[],VKnots,VWeights);
% 
%             UVPoles=cat(3,[Ctrl_0v(:,1);Ctrl_1v(:,1)],[Ctrl_0v(:,2);Ctrl_1v(:,2)],[Ctrl_0v(:,3);Ctrl_1v(:,3)]);
%             UVWeights=[Ctrl_0v(:,4),Ctrl_1v(:,4)];
%             UVfce=FaceNURBS('',VPoles,[],VDegree,[],VMults,[],VKnots,VWeights);

            % generata U, V and UV direction Poles
            Ctrl_u0=reshape(Ctrl_u0,1,size(Ctrl_u0,1),size(Ctrl_u0,2));
            Ctrl_u1=reshape(Ctrl_u1,1,size(Ctrl_u1,1),size(Ctrl_u1,2));
            Ctrl_0v=reshape(Ctrl_0v,size(Ctrl_0v,1),1,size(Ctrl_0v,2));
            Ctrl_1v=reshape(Ctrl_1v,size(Ctrl_1v,1),1,size(Ctrl_1v,2));
            ctrl_00=Ctrl_u0(1,1,:);ctrl_10=Ctrl_u0(1,end,:);
            ctrl_01=Ctrl_u1(1,1,:);ctrl_11=Ctrl_u1(1,end,:);
            U_Ctrl=Ctrl_u0.*(1-V)+Ctrl_u1.*V;
            V_Ctrl=Ctrl_0v.*(1-U)+Ctrl_1v.*U;
            UV_Ctrl=ctrl_00.*(1-V).*(1-U)+ctrl_10.*(1-V).*U+...
                ctrl_01.*V.*(1-U)+ctrl_11.*V.*U;

            % combine coefficient to construct Coons surface
            Ctrl=U_Ctrl+V_Ctrl-UV_Ctrl;

            Poles=Ctrl(:,:,1:end-1);
            Weights=Ctrl(:,:,end);

            fce=FaceNURBS(name,Poles,UDegree,VDegree,UMults,VMults,UKnots,VKnots,Weights);
            
            function Ctrl=getMatchCtrl(Ctrl_0,Ctrl_1,Degree,Mults,Knots)
                u_num=length(baseKnotVec(Mults,Knots));
                pole_num=u_num-Degree-1;
                UC=linspace(0,1,pole_num);
                Ctrl=Ctrl_0*(1-UC)+Ctrl_1*UC;
            end
        end
    end

    methods(Static)
        function [Ctrl_new,Mults_new]=addDegree(deg,Ctrl,Mults,Knots,deg_tar)
            % increase curve degree
            %
            if deg_tar <= deg
                Ctrl_new=Ctrl;
                Mults_new=Mults;
                return;
            end

            U=baseKnotVec(Mults,Knots);
            for deg_temp=deg+1:deg_tar
                % add repeat node
                Deg_new=deg+1;
                Mults_new=Mults+1;
                U_new=baseKnotVec(Mults_new,Knots);

                % calculate new ctrl
                matrix=zeros(size(Ctrl,1)+length(Knots)-1,size(Ctrl,1));
                for j=1:size(Ctrl,1)+length(Knots)-1
                    for i=1:size(Ctrl,1)
                        matrix(j,i)=baseFcnL(U,U_new,i,j,deg);
                    end
                end
                matrix=1/Deg_new*matrix;
                Ctrl_new=(matrix*Ctrl);

                % sort old curve data
                Ctrl=Ctrl_new;
                deg=Deg_new;
                Mults=Mults_new;
                U=U_new;
            end
        end

        function [Ctrl_new,U_new]=insertKnot(deg,Ctrl,U,U_ins)
            % insert knot to NURBS
            %

            % allocate memory
            [ctrl_num,dim]=size(Ctrl);
            U_ins =sort(U_ins);
            nu=length(U_ins);
            u_num=length(U);
            ctrl_num_new=ctrl_num+nu;
            Ctrl_new=zeros(ctrl_num_new,dim);
            U_new=zeros(1,u_num+nu);

            n=ctrl_num-1;
            r=nu-1;
            m=n+deg+1;

            % locate add knots position
            if (U_ins(1) == U(ctrl_num+1)),a=ctrl_num-1;
            else,a=find(U_ins(1) >= U,1,'last')-1;end
            if (U_ins(r+1) == U(ctrl_num+1)),b=ctrl_num-1;
            else,b=find(U_ins(r+1) >= U,1,'last')-1;end
            b=b+1;

            Ctrl_new(1:a-deg+1,:)=Ctrl(1:a-deg+1,:);
            Ctrl_new(b+nu:ctrl_num+nu,:)=Ctrl(b:ctrl_num,:);

            U_new(1:a+1)=U(1:a+1);
            U_new(b+deg+nu+1:m+nu+1)=U(b+deg+1:m+1);

            ii=b+deg-1;
            ss=ii+nu;

            % calculate new ctrl point
            for jj=r:-1:0
                ind=(a+1):ii;
                ind=ind(U_ins(jj+1)<=U(ind+1));
                Ctrl_new(ind+ss-ii-deg,:)=Ctrl(ind-deg,:);
                U_new(ind+ss-ii+1)=U(ind+1);
                ii=ii-length(ind);
                ss=ss-length(ind);
                Ctrl_new(ss-deg,:)=Ctrl_new(ss-deg+1,:);
                for l=1:deg
                    ind=ss-deg+l;
                    alfa=U_new(ss+l+1)-U_ins(jj+1);
                    if abs(alfa) == 0
                        Ctrl_new(ind,:)=Ctrl_new(ind+1,:);
                    else
                        alfa=alfa/(U_new(ss+l+1)-U(ii-deg+l+1));
                        tmp=(1-alfa)*Ctrl_new(ind+1,:);
                        Ctrl_new(ind,:)=alfa*Ctrl_new(ind,:)+tmp;
                    end
                end
                U_new(ss+1)=U_ins(jj+1);
                ss=ss-1;
            end
        end
    end

    methods(Static)
        function edg=JointEdge(name,edg_list,geom_torl)
            % connect BSpline curve into one curve
            %
            if nargin < 3,geom_torl=[];end

            if isempty(geom_torl),geom_torl=100*eps;end

            % search max Degree
            Degree=0;
            edg_idx=1;
            while edg_idx <= length(edg_list)
                edg_curr=edg_list{edg_idx};
                if isempty(edg_curr)
                    edg_list(edg_idx)=[];
                else
                    if Degree < edg_curr.Degree
                        Degree=edg_curr.Degree;
                    end
                end
                edg_idx=edg_idx+1;
            end
            edg_num=length(edg_list);

            % load all Poles to correct line order
            line_list={};
            for edg_idx=1:edg_num
                line_list=[line_list,{[edg_list{edg_idx}.Poles,edg_list{edg_idx}.Weights]}];
            end
            [~,map_list,order_list]=GeomApp.correctLine(line_list,geom_torl);
            edg_list=edg_list(map_list);

            % increase Degree of edg and load data
            Poles=[];
            Mults=[];
            Knots=[];
            Weights=[];
            total_length=0;
            for edg_idx=1:edg_num
                % load edge
                edg=edg_list{edg_idx};
                if order_list(edg_idx),edg.reverse();end

                edg.addDegree(Degree);
                Mults_new=edg.Mults;
                Mults_new(1)=Mults_new(1)-1;
                Weights_new=edg.Weights;
                edg_length=sum(vecnorm(diff(edg.Poles),2,2));

                Poles=[Poles(1:end-1,:);edg.Poles];
                Mults=[Mults(1:end-1),Mults_new];
                Knots=[Knots(1:end-1),edg.Knots*edg_length+total_length];
                Weights=[Weights(1:end-1);Weights_new/Weights_new(1)];

                Weights=Weights/Weights(end);
                total_length=total_length+edg_length;
            end
            Mults(1)=Mults(1)+1;

            % connect
            Knots=Knots/total_length;
            Knots=Knots/max(Knots);

            edg=EdgeNURBS(name,Poles,Degree,Mults,Knots,Weights);
        end

        function [edg_1,edg_2,Degree,Mults,Knots]=matchEdge(edg_1,edg_2)
            Degree=max(edg_1.Degree,edg_2.Degree);
            edg_1.addDegree(Degree);
            edg_2.addDegree(Degree);
            U1=edg_1.u_list;
            U2=edg_2.u_list;

            % merge the knot vectors of u
            UC=unique([U1,U2]);
            U1_ins=[];
            U2_ins=[];
            for i=1:length(UC)
                i1=sum(U1 == UC(i));
                i2=sum(U2 == UC(i));
                m=max(i1,i2);
                U1_ins=[U1_ins,UC(i)*ones(1,m-i1)];
                U2_ins=[U2_ins,UC(i)*ones(1,m-i2)];
            end

            if ~isempty(U1_ins),edg_1.insertKnot(U1_ins);end
            if ~isempty(U2_ins),edg_2.insertKnot(U2_ins);end
            Mults=edg_1.Mults;
            Knots=edg_1.Knots;
        end

        function Point=MapGrid(line_u0,line_u1,line_0v,line_1v,u_list,v_list)
            % generate discrete area by mapping
            % using Lagrange polynomial mapping method
            % local parameter u and v is equispaced
            %
            % input:
            % line_u0 (matrix): u_num x dimension
            % line_u1 (matrix): u_num x dimension
            % line_v0 (matrix): v_num x dimension
            % line_v1 (matrix): v_num x dimension
            %
            % output:
            % Point (matrix): v_num x u_num x dimension
            %
            if nargin < 6,v_list=[];if nargin < 5,u_list=[];end;end

            geom_torl=100*eps;
            line_list={line_u0,line_1v,flipud(line_u1),flipud(line_0v)};
            line_list=GeomApp.correctLine(line_list,geom_torl);
            if ~any(norm(line_list{4}(end,:)-line_list{1}(1,:)) < geom_torl)
                error('GeomApp.MapGrid: line not connect');
            end

            line_u0=line_list{1};
            line_1v=line_list{2};
            line_u1=flipud(line_list{3});
            line_0v=flipud(line_list{4});

            u_num=size(line_u0,1);v_num=size(line_0v,1);
            if u_num ~= size(line_u1,1) || v_num ~= size(line_1v,1)
                error('GeomApp.MapGrid: opposite line of boundary discrete number no equal')
            end
            dimension=size(line_u0,2);
            if isempty(u_list),u_list=linspace(0,1,u_num);end,u_list=u_list(:)';
            if isempty(v_list),v_list=linspace(0,1,v_num);end,v_list=v_list(:)';

            Point=zeros(v_num,u_num,dimension);

            % preproces boundary
            Point(1,:,:)=line_u0(:,:);
            Point(end,:,:)=line_u1(:,:);
            Point(:,1,:)=line_0v(:,:);
            Point(:,end,:)=line_1v(:,:);

            if u_num > 2 && v_num > 2
                H=diff(u_list);G=diff(v_list);

                % solve prepare
                H_ipj=H(1:end-1)+H(2:end);
                H_ij=H(1:end-1).*H(2:end).*H_ipj/2;
                G_ipj=G(1:end-1)+G(2:end);
                G_ij=G(1:end-1).*G(2:end).*G_ipj/2;

                % D_xx
                d_xx=spdiags([H(1:end-1)./H_ij;-H_ipj./H_ij;H(2:end)./H_ij]',-1:1,u_num-2,u_num-2)';
                I_M=speye(v_num-2,v_num-2);
                D_xx=kron(d_xx,I_M);

                % D_yy
                d_yy=spdiags([G(1:end-1)./G_ij;-G_ipj./G_ij;G(2:end)./G_ij]',-1:1,v_num-2,v_num-2)';
                I_K=speye(u_num-2,u_num-2);
                D_yy=kron(I_K,d_yy);

                % Laplace matrix
                Lap=(D_xx+D_yy);

                Point_B=zeros(v_num-2,u_num-2,dimension);
                Point_B(1,:,:)=Point_B(1,:,:)-Point(1,2:u_num-1,:)*(G(2)/G_ij(1));
                Point_B(end,:,:)=Point_B(end,:,:)-Point(end,2:u_num-1,:)*(G(end-1)/G_ij(end));
                Point_B(:,1,:)=Point_B(:,1,:)-Point(2:v_num-1,1,:)*(H(2)/H_ij(1));
                Point_B(:,end,:)=Point_B(:,end,:)-Point(2:v_num-1,end,:)*(H(end-1)/H_ij(end));
                for dim_idx=1:dimension
                    Point_B_page=Point_B(:,:,dim_idx);
                    Point(2:v_num-1,2:u_num-1,dim_idx)=reshape(Lap\Point_B_page(:),v_num-2,u_num-2);
                end
            end
        end

        function [line_list,map_list,order_list]=correctLine(line_list,geom_torl)
            % correct line point order
            % order should be anti clockwise and start from first line
            %
            line_num=length(line_list);
            map_list=1:line_num;
            order_list=false(line_num,1);

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
                    overlap_idx=find(dist < geom_torl,1);
                    if ~any(overlap_idx)
                        error('GeomApp.correctLine: line not connect');
                    end
                    exchange_idx=ceil(overlap_idx/2)+line_idx;

                    % exchange line
                    line_temp=line_list{exchange_idx};
                    line_list{exchange_idx}=line_next;
                    line_list{line_idx+1}=line_temp;
                    map_list(line_idx+1)=exchange_idx;

                    % reorder point in line order
                    if mod(overlap_idx,2) == 0
                        line_list{line_idx+1}=flipud(line_list{line_idx+1});
                        order_list(line_idx+1)=true;
                    end
                end
            end
        end
    end

    methods(Static)
        function [U,Fval,node_list]=meshAdapt1D(fcn,low_bou,up_bou,torl,min_level,max_level,fval_num)
            % Binary-tree
            % adapt capture function value
            % ensure error of linear interplation will less than torl
            %
            if nargin < 7,fval_num=1;end

            % node_list which is a matrix store all node
            % a node is a array,contain level,index_1,index_2,index_c,node_index_1,node_index_2
            % place:
            % 1-c-2
            % cell:
            % 1-2
            % if children node is empty,left_index or right_index will be zero
            list_add_num=50; % node_list will be extend only when node_list is not enough
            node_list=zeros(list_add_num,6,'int64');

            % data_list use to sort all float data include coodinate,function value
            data_list=zeros(list_add_num,fval_num+1);

            % add vertex of cell into data_list first
            data_list(1,:)=[low_bou,fcn(low_bou)];
            data_list(2,:)=[up_bou,fcn(up_bou)];

            % create base root
            node_list(1,:)=[0,1,2,0,0,0];

            node_num=createNodeTree(1,2); % create node tree from root
            node_list=node_list(1:node_num,:);
            [U,Fval]=traversalInorder(1); % from small to large get list

            % add boundary info
            U=[data_list(1,1);U;data_list(2,1)];
            Fval=[data_list(1,2:end);Fval;data_list(2,2:end)];

            function node_num=createNodeTree(root_idx,data_num)
                % create node tree from root
                %
                stack=root_idx;
                node_num=root_idx;

                while ~isempty(stack)
                    % current node information
                    node_idx=stack(end);
                    node=node_list(node_idx,:);
                    stack=stack(1:end-1);

                    if node(1) < max_level
                        % judge if linear predict if accptable
                        % if not,create children cell
                        %
                        coord_c=(data_list(node(2),1)+data_list(node(3),1))/2;
                        fval_c=fcn(coord_c);
                        fval_c_pred=(data_list(node(2),2:end)+data_list(node(3),2:end))/2;

                        % add 1 new data into data_list
                        data_new_idx=data_num+1;
                        if data_num+1 > size(data_list,1)
                            data_list=[data_list;zeros(list_add_num,fval_num+1)];
                        end
                        data_list(data_new_idx,:)=[coord_c,fval_c];
                        node(4)=data_new_idx;
                        data_num=data_num+1;

                        % add 2 new node to node_list
                        node_new_idx=node_num+[1,2];
                        if node_num+2 > size(node_list,1)
                            node_list=[node_list;zeros(list_add_num,6)];
                        end
                        node_list(node_new_idx,:)=[
                            node(1)+1,node(2),node(4),0,0,0;
                            node(1)+1,node(4),node(3),0,0,0];
                        node([5,6])=node_new_idx;
                        node_num=node_num+2;

                        if any(abs(fval_c-fval_c_pred) > torl) || node(1) <= min_level
                            % add to stack
                            stack=[stack,node_new_idx];
                            node_list(node_idx,:)=node;
                        end
                    end
                end
            end

            function [x_list,fval_list]=traversalInorder(root_idx)
                % inorder traversal node tree to obtain data
                % inorder will make sure x_list is from small to large
                %
                stack=[];
                x_list=[];
                fval_list=[];
                node_idx=root_idx;

                while ~isempty(stack) || node_idx ~= 0
                    while node_idx ~= 0
                        stack=[stack,node_idx];

                        % node=node.left;
                        node_idx=node_list(node_idx,5);
                    end

                    node_idx=stack(end);
                    stack=stack(1:end-1);
                    data_idx=node_list(node_idx,4);
                    if data_idx ~= 0
                        x_list=[x_list;data_list(data_idx,1)];
                        fval_list=[fval_list;data_list(data_idx,2:end)];
                    end

                    % node=node.right;
                    node_idx=node_list(node_idx,6);
                end
            end
        end

        function [U,V,Fval,node_list]=meshAdapt2DUV(fcn,low_bou,up_bou,torl,min_level,max_level,fval_num)
            % 2D Omni-Tree
            % adapt capture 2 dimemsion function value
            % ensure error of linear interplation will less than torl
            %
            if nargin < 7,fval_num=1;end

            % node_list which is a matrix store all node
            % a node is a array,contain:
            % level,idx_1-8(index of data_list),idx_c,node_index_1-4
            % place:
            % 3-8-4 or 3-4 or 3-8-4
            % 1-5-2    6-7    6-c-7
            %          1-2    1-5-2
            % node:
            % 1-2 or 3 or 3-4
            %        1    1-2
            % if children node is empty,left_index or right_index will be zero
            list_add_num=50; % list will be extend only when list is not enough
            node_list=zeros(list_add_num,14,'int64');
            % data_list use to sort all float data include coodinate,function value
            data_list=zeros(list_add_num,fval_num+2);

            % add vertex of cell into data_list first
            bou_1=[low_bou(1),low_bou(2)];
            data_list(1,:)=[bou_1,fcn(bou_1)];
            bou_2=[up_bou(1),low_bou(2)];
            data_list(2,:)=[bou_2,fcn(bou_2)];
            bou_3=[low_bou(1),up_bou(2)];
            data_list(3,:)=[bou_3,fcn(bou_3)];
            bou_4=[up_bou(1),up_bou(2)];
            data_list(4,:)=[bou_4,fcn(bou_4)];

            % create base root
            node_list(1,:)=[0,1,2,3,4,0,0,0,0,0,0,0,0,0];

            % create node tree from root
            [node_num,data_num]=createNodeTree(1,4);
            node_list=node_list(1:node_num,:);
            data_list=data_list(1:data_num,:);

            % generate U and V
            [u_list,~,u_idx]=unique(data_list(:,1));
            [v_list,~,v_idx]=unique(data_list(:,2));
            idx_list=[u_idx,v_idx];

            [U,V]=meshgrid(u_list,v_list);

            % local data to Fval
            Fval=nan(length(v_list),length(u_list),fval_num);

            for data_index=1:data_num
                if all([idx_list(data_index,2),idx_list(data_index,1)] ~= 0)
                    Fval(idx_list(data_index,2),idx_list(data_index,1),:)=data_list(data_index,3:end);
                end
            end

            % fit nan data
            idx=find(isnan(Fval(:,:,1)));
            Fval_sub=fcn([U(idx),V(idx)]);
            for favl_idx=1:fval_num
                Fval(idx+(favl_idx-1)*(length(v_list)*length(u_list)))=Fval_sub(:,favl_idx);
            end

            function [node_num,data_num]=createNodeTree(root_idx,data_num)
                % create quad tree
                %
                stack=root_idx;
                node_num=root_idx;

                while ~isempty(stack)
                    % current node information
                    node_idx=stack(end);
                    node=node_list(node_idx,:);
                    stack=stack(1:end-1);

                    if node(1) < max_level
                        % judge if linear predict if accptable
                        % if not,create children cell
                        %
                        [coord_c,coord_5,coord_6,coord_7,coord_8,...
                            fval_c,fval_5,fval_6,fval_7,fval_8]=calCell(node(2),node(3),node(4),node(5));
                        [fval_pred_c,fval_pred_5,fval_pred_6,...
                            fval_pred_7,fval_pred_8]=calCellPred(node(2),node(3),node(4),node(5));

                        % check u direction
                        if any(abs(fval_5-fval_pred_5) > torl) || any(abs(fval_8-fval_pred_8) > torl)
                            add_u_flag=true(1);
                        else
                            add_u_flag=false(1);
                        end

                        % check v direction
                        if any(abs(fval_6-fval_pred_6) > torl) || any(abs(fval_7-fval_pred_7) > torl)
                            add_v_flag=true(1);
                        else
                            add_v_flag=false(1);
                        end

                        % check center place
                        if any(abs(fval_c-fval_pred_c) > torl)
                            add_u_flag=true(1);
                            add_v_flag=true(1);
                        end

                        if (add_u_flag && add_v_flag) || node(1) < min_level
                            % add 5 data into data_list
                            data_new_idx=data_num+(1:5);
                            if data_num+5 > size(data_list,1)
                                data_list=[data_list;zeros(list_add_num,fval_num+2)];
                            end
                            data_list(data_new_idx,:)=[
                                coord_5,fval_5;
                                coord_6,fval_6;
                                coord_7,fval_7;
                                coord_8,fval_8;
                                coord_c,fval_c;];
                            node(6:10)=data_new_idx;
                            data_num=data_num+5;

                            % add 4 new node to node_list
                            node_new_idx=node_num+(1:4);
                            if node_num+4 > size(node_list,1)
                                node_list=[node_list;zeros(list_add_num,14)];
                            end
                            node_list(node_new_idx,:)=[...
                                node(1)+1,node(2),node(6),node(7),node(10),0,0,0,0,0,0,0,0,0;
                                node(1)+1,node(6),node(3),node(10),node(8),0,0,0,0,0,0,0,0,0;
                                node(1)+1,node(7),node(10),node(4),node(9),0,0,0,0,0,0,0,0,0;
                                node(1)+1,node(10),node(8),node(9),node(5),0,0,0,0,0,0,0,0,0;];
                            node(11:14)=node_new_idx;
                            node_num=node_num+4;

                        elseif add_u_flag
                            % add 2 data into data_list
                            data_new_idx=data_num+(1:2);
                            if data_num+2 > size(data_list,1)
                                data_list=[data_list;zeros(list_add_num,fval_num+2)];
                            end
                            data_list(data_new_idx,:)=[
                                coord_5,fval_5;
                                coord_8,fval_8;];
                            node([6,9])=data_new_idx;
                            data_num=data_num+2;

                            % add 2 new node to node_list
                            node_new_idx=node_num+(1:2);
                            if node_num+2 > size(node_list,1)
                                node_list=[node_list;zeros(list_add_num,14)];
                            end
                            node_list(node_new_idx,:)=[...
                                node(1)+1,node(2),node(6),node(4),node(9),0,0,0,0,0,0,0,0,0;
                                node(1)+1,node(6),node(3),node(9),node(5),0,0,0,0,0,0,0,0,0;];
                            node([11,12])=node_new_idx;
                            node_num=node_num+2;

                        elseif add_v_flag
                            % add 2 data into data_list
                            data_new_idx=data_num+(1:2);
                            if data_num+2 > size(data_list,1)
                                data_list=[data_list;zeros(list_add_num,fval_num+2)];
                            end
                            data_list(data_new_idx,:)=[
                                coord_6,fval_6;
                                coord_7,fval_7;];
                            node([7,8])=data_new_idx;
                            data_num=data_num+2;

                            % add 2 new node to node_list
                            node_new_idx=node_num+(1:2);
                            if node_num+2 > size(node_list,1)
                                node_list=[node_list;zeros(list_add_num,14)];
                            end
                            node_list(node_num+(1:2),:)=[...
                                node(1)+1,node(2),node(3),node(7),node(8),0,0,0,0,0,0,0,0,0;
                                node(1)+1,node(7),node(8),node(4),node(5),0,0,0,0,0,0,0,0,0;];
                            node([11,13])=node_new_idx;
                            node_num=node_num+2;
                        else
                            node_new_idx=[];
                        end

                        % add to stack
                        stack=[stack,node_new_idx];
                        node_list(node_idx,:)=node;
                    end
                end
            end

            function [coord_c,coord_5,coord_6,coord_7,coord_8,...
                    fval_c,fval_5,fval_6,fval_7,fval_8]=calCell(vidx_1,vidx_2,vidx_3,vidx_4)
                % calculate index of c,5,6,7,8 place function value
                % abbreviation:
                % vidx: vertex index
                %
                coord_c=(data_list(vidx_1,[1,2])+data_list(vidx_4,[1,2]))/2;
                fval_c=fcn(coord_c);

                coord_5=(data_list(vidx_1,[1,2])+data_list(vidx_2,[1,2]))/2;
                fval_5=fcn(coord_5);
                coord_6=(data_list(vidx_1,[1,2])+data_list(vidx_3,[1,2]))/2;
                fval_6=fcn(coord_6);
                coord_7=(data_list(vidx_2,[1,2])+data_list(vidx_4,[1,2]))/2;
                fval_7=fcn(coord_7);
                coord_8=(data_list(vidx_3,[1,2])+data_list(vidx_4,[1,2]))/2;
                fval_8=fcn(coord_8);
            end

            function [fval_pred_c,fval_pred_5,fval_pred_6,...
                    fval_pred_7,fval_pred_8]=calCellPred(vidx_1,vidx_2,vidx_3,vidx_4)
                % calculate index of c,5,6,7,8 place linear predict function value
                %
                fval_pred_c=(data_list(vidx_1,3:end)+data_list(vidx_2,3:end)+data_list(vidx_3,3:end)+data_list(vidx_4,3:end))/4;
                fval_pred_5=(data_list(vidx_1,3:end)+data_list(vidx_2,3:end))/2;
                fval_pred_6=(data_list(vidx_1,3:end)+data_list(vidx_3,3:end))/2;
                fval_pred_7=(data_list(vidx_2,3:end)+data_list(vidx_4,3:end))/2;
                fval_pred_8=(data_list(vidx_3,3:end)+data_list(vidx_4,3:end))/2;
            end

        end

    end
end