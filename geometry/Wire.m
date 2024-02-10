classdef Wire < handle & matlab.mixin.Copyable
    % Topology Entity Wire
    %
    properties
        name=''; % name of wire
        edge_list=Edge.empty(); % edge list
        vertex_list=[]; % vertex list
        looped=false; % if wire is looped
    end

    properties
        topology=struct(); % association with edge

        bound_box=[]; % geometry bounding box
        dimension=[]; % dimension of curve
    end

    % define wire
    methods
        function self=Wire(name,edg_list,geom_torl)
            % generate wire entity
            %
            if nargin < 3, geom_torl=[];end
            if isempty(geom_torl), geom_torl=1e-6;end
            self.name=name;
            edg_list=edg_list(:);

            if isempty(edg_list)
                return;
            end

            % check dimension of each edge
            edg_num=length(edg_list);
            dimension=edg_list(1).dimension;
            for edg_idx=2:edg_num
                if edg_list(edg_idx).dimension ~= dimension
                    error('Wire: dimension of each edge is not equal');
                end
            end
            self.edge_list=edg_list;

            % generate bounding box
            bound_box=edg_list(1).bound_box;
            for edg_idx=2:edg_num
                bound_box(1,:)=min([bound_box(1,:);edg_list(edg_idx).bound_box(1,:)],[],1);
                bound_box(2,:)=max([bound_box(2,:);edg_list(edg_idx).bound_box(2,:)],[],1);
            end
            self.bound_box=bound_box;
            self.dimension=dimension;

            % calculate topology of wire
            self.calTopology(geom_torl);
        end

        function self=calTopology(self,geom_torl)
            % calculate topology property of wire
            %
            edg_list=self.edge_list;

            edg_num=length(edg_list);
            if edg_num > 1
                % check whether edge is connect
                % load all boundary vertex of edge
                % and create line list for correctLine
                line_list=cell(1,edg_num);
                for edg_idx=1:edg_num

                    line_list{edg_idx}=edg_list(edg_idx).vertex_list;
                end
                [line_list,map_list,reverse_list]=GeomApp.correctLine(line_list,geom_torl);
                edg_list=edg_list(map_list);
                for edg_idx=1:edg_num
                    if reverse_list(edg_idx)
                        edg_list(edg_idx).reverse();
                    end
                end

                % load all vertex
                vtx_list=zeros(2*edg_num,self.dimension);
                for edg_idx=1:edg_num
                    vtx_list((2*edg_idx-1):(2*edg_idx),:)=line_list{edg_idx};
                end
            else
                vtx_list=edg_list(1).vertex_list;
            end

            % check if wire is loop
            if norm(vtx_list(1,:)-vtx_list(end,:),2) < geom_torl
                lopd=true;
                vtx_list=vtx_list(1:2:(2*edg_num-1),:);
            else
                lopd=false;
                vtx_list=vtx_list([1:2:(2*edg_num-1),2*edg_num],:);
            end

            % generate association between edge and vertex
            vtx_num=size(vtx_list,1);
            vtx_topo=repmat(struct('edge',[]),vtx_num,1);
            if lopd
                vtx_topo(1).edge=[edg_num,1];
                vtx_topo(end).edge=[edg_num-1,edg_num];
            else
                vtx_topo(1).edge=1;
                vtx_topo(end).edge=edg_num;
            end
            for vtx_idx=2:vtx_num-1
                vtx_topo(vtx_idx).edge=[vtx_idx-1,vtx_idx];
            end

            self.edge_list=edg_list;
            self.vertex_list=vtx_list;
            self.looped=lopd;
            
            topo.vertex=vtx_topo;
            self.topology=topo;
        end

        function self=addEdge(self,edg)
            % add edge into wire
            %
            if self.looped
                error('Wire.addEdge: can not add new edge to a looped wire')
            end

            % check whether edge is connect
            line_list={self.vertex_list,edg.vertex_list};
            [line_list,map_list,reverse_list]=GeomApp.correctLine(line_list,geom_torl);
            
            if reverse_list(1), self.reverse();end
            if reverse_list(2), edg.reverse();end

            % updata properties
            if map_list(1) == 1 && map_list(2) == 2
                self.edge_list=[self.edge_list;edg];
            elseif map_list(1) == 2 && map_list(2) == 1
                self.edge_list=[edg;self.edge_list];
            end
            self.vertex_list=[line_list{1};line_list{2}(2,:)];

            % check if wire is loop
            if norm(self.vertex_list(1,:)-self.vertex_list(end,:),2) < geom_torl
                self.looped=true;
                self.vertex_list(end,:)=[];
            end

            % updata bounding box
            self.bound_box(1,:)=min([self.bound_box(1,:);edg.bound_box(1,:)],[],1);
            self.bound_box(2,:)=max([self.bound_box(2,:);edg.bound_box(2,:)],[],1);
        end

        function [edg,edg_idx]=getEdge(self,edg_name)
            % load edge from edge_list base on input edg_name
            %
            edg=[];
            edg_idx=[];
            for edg_i=1:length(self.edge_list)
                edg_u=self.edge_list(edg_i);
                if strcmp(edg_u.name,edg_name)
                    edg=[edg,edg_u];
                    edg_idx=[edg_idx,edg_i];
                end
            end
        end

        function self=splitWire(self,edg_idx,u_b)
            % split curve at ub
            %

            % split curve
            vtx_b=self.edge_list(edg_idx).calPoint(u_b);
            [edg_1,edg_2]=self.edge_list(edg_idx).splitEdge(u_b);

            % modify wire
            self.vertex_list=[self.vertex_list(1:edg_idx,:);vtx_b;self.vertex_list(edg_idx+1:end,:)];
            self.edge_list=[self.edge_list(1:edg_idx-1,:);edg_1;edg_2;self.edge_list(edg_idx+1:end,:)];
        end

        function [ovlp,topo_wir]=checkOverlap(self,wir,geom_torl)
            % check if two edge is overlap
            %
            % output:
            % ovlp: whether have overlap edge
            % topo_edg: cell matrix, overlap check of each pair edge
            %
            if nargin < 3, geom_torl=[];end
            if isempty(geom_torl), geom_torl=1e-6;end

            % check if bounding box is overlap
            if ~GeomApp.checkBoxOverlap(self,wir,geom_torl)
                % if no, pass
                ovlp=false;topo_wir={};
                return;
            end

            % check each edge if overlap
            ovlp=false;
            ovlp_mat=false(length(self.edge_list),length(wir.edge_list));
            U_bas=cell(length(self.edge_list),length(wir.edge_list));
            U_ref=cell(length(self.edge_list),length(wir.edge_list));
            for edg_bas_idx=1:length(self.edge_list)
                edg_bas=self.edge_list(edg_bas_idx);
                for edg_ref_idx=1:length(wir.edge_list)
                    edg_ref=wir.edge_list(edg_ref_idx);
                    
                    % check edge overlap
                    [ovlp_edg,self_u_proj,edg_u_proj]=edg_bas.checkOverlap(edg_ref,geom_torl);
                    ovlp=ovlp | ovlp_edg;
                    ovlp_mat(edg_bas_idx,edg_ref_idx)=ovlp_edg;
                    U_bas(edg_bas_idx,edg_ref_idx)={self_u_proj};
                    U_ref(edg_bas_idx,edg_ref_idx)={edg_u_proj};
                end
            end
            topo_wir.overlap=ovlp_mat;
            topo_wir.U_base=U_bas;
            topo_wir.U_ref=U_ref;
        end

        function self=reverse(self)
            % reverse wire direction
            %
            for edg_idx=1:length(self.edge_list)
                self.edge_list(edg_idx).reverse();
            end
            self.edge_list=self.edge_list(end:-1:1);
            self.vertex_list=self.vertex_list(end:-1:1,:);
        end

        function [Points_list,U_list]=calGeom(self,u_param)
            % calculate point on all curve
            %
            if nargin < 2
                u_param=[];
            end

            % using bound box to auto calculate capture precision
            if isempty(u_param), u_param=2^-6*mean(self.bound_box(2,:)-self.bound_box(1,:));end

            edg_num=length(self.edge_list);
            Points_list=cell(edg_num,1);
            U_list=cell(edg_num,1);
            for edg_idx=1:edg_num
                [Pnts,U]=self.edge_list(edg_idx).calGeom(u_param);
                Points_list{edg_idx}=Pnts;
                U_list{edg_idx}=U;
            end
        end

        function [ln_hdl_list,sctr_hdl]=displayModel(self,axe_hdl,u_param)
            % display wire on axes
            %
            if nargin < 3
                u_param=[];
                if nargin < 2
                    axe_hdl=[];
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end

            % calculate point on edge
            [Points_list]=self.calGeom(u_param);

            % plot edge
            edg_num=length(self.edge_list);
            for edg_idx=1:edg_num
                ln_hdl_list(edg_idx)=linePoints(axe_hdl,Points_list{edg_idx},self.dimension);
            end
            sctr_hdl=scatterPoints(axe_hdl,self.vertex_list,self.dimension);

            if self.dimension == 2
            else
                zlabel('z');
                view(3);
            end
            xlabel('x');
            ylabel('y');

            function ln_hdl=linePoints(axe_hdl,Points,dimension)
                % line point on different dimension
                %
                if dimension == 2
                    ln_hdl=line(axe_hdl,Points(:,1),Points(:,2));
                else
                    ln_hdl=line(axe_hdl,Points(:,1),Points(:,2),Points(:,3));
                end
            end

            function sctr_hdl=scatterPoints(axe_hdl,Points,dimension)
                % scatter point on different dimension
                %
                hold on;
                if dimension == 2
                    sctr_hdl=scatter(axe_hdl,Points(:,1),Points(:,2));
                else
                    sctr_hdl=scatter3(axe_hdl,Points(:,1),Points(:,2),Points(:,3));
                end
                hold off;
            end
        end
    end

    % parameterized function
    methods

    end
end