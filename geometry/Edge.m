classdef Edge < handle & matlab.mixin.Copyable
    % Topology Entity Edge with Bound
    %
    properties
        name=''; % name of edge
        curve; % Curve handle
        vertex_list=[]; % vertex list of boundary
        bound=[]; % boundary of outer
    end

    properties
        topology=struct(); % association with outer_bound and vertex_list

        bound_box=[]; % geometry bounding box
        dimension=[]; % dimension of curve
    end

    % define edge
    methods
        function self=Edge(name,crv,vtx_list)
            % generate edge entity
            %
            if nargin < 3, vtx_list=[];end
            self.name=name;
            self.curve=crv;

            if isempty(crv)
                return;
            end

            if isa(crv,'Curve'), region=[min(crv.Knots),max(crv.Knots)];
            elseif isa(crv,'CurveCST'), region=[0,1];end
            if isempty(vtx_list)
                bou=region;
                vtx_list=crv.calPoint(bou);
            else
                % project vertex list to get outer_bound
                if size(vtx_list,1) ~= 2
                    error('Edge: only can have a pair of vertex');
                end
                bou=crv.projectPoint(vtx_list);
                if bou(2) < bou(1)
                    bou=bou(end:-1:1);
                    vtx_list=vtx_list(end:-1:1,:);
                end
                bou(1)=max(bou(1),region(1));
                bou(2)=min(bou(2),region(2));
            end
            self.vertex_list=vtx_list;
            self.bound=bou;

            if self.bound(1) < min(crv.Knots) || self.bound(2) > max(crv.Knots)
                disp('?');
            end

            topo.vertex=[1,2];
            self.topology=topo;

            U_x=linspace(bou(1),bou(2),101);
            Pnts=crv.calPoint(U_x);
            self.bound_box=[min(Pnts,[],1);max(Pnts,[],1)];
            self.dimension=size(Pnts,2);
        end

        function self=reverse(self)
            % revese direction of edge
            %
            self.curve.reverse();
            self.vertex_list=flipud(self.vertex_list);

            if isa(self.curve,'Curve'), bou_max=max(self.curve.Knots);
            elseif isa(self.curve,'CurveCST'), bou_max=1;end
            self.bound=bou_max-self.bound(end:-1:1);
        end

        function [edg_1,edg_2]=splitEdge(self,u_b)
            % split curve at ub
            %
            if u_b <= self.bound(1) || u_b >= self.bound(2)
                error('Edge.splitEdge: ub is out of boundary of curve');
            end

            % split curve
            vtx_b=self.curve.calPoint(u_b);
            [crv_1,crv_2]=self.curve.splitCurve(u_b);

            % generate new edge
            vtx_list_1=[self.vertex_list(1,:);vtx_b];
            edg_1=Edge(self.name,crv_1,vtx_list_1);

            vtx_list_2=[vtx_b;self.vertex_list(2,:)];
            edg_2=Edge(self.name,crv_2,vtx_list_2);
        end

        function [ovlp,self_u_proj,edg_u_proj]=checkOverlap(self,edg,geom_torl)
            % check if two edge is overlap
            %
            % output:
            % ovlp: whether edge is overlap
            % self_u_proj: u of self vertex project to edg
            % edg_u_proj: u of edg vertex project to self
            %
            if nargin < 3, geom_torl=[];end
            if isempty(geom_torl), geom_torl=1e-6;end

            ovlp=false;self_u_proj=[];edg_u_proj=[];
            % check if bounding box is overlap
            if ~GeomApp.checkBoxOverlap(self,edg,geom_torl)
                % if no, pass
                return;
            end

            if norm(self.vertex_list-edg.vertex_list) < geom_torl
                Pnts_self=self.calGeom(21);Pnts_ref=edg.calGeom(21);
                if norm(Pnts_self-Pnts_ref) < geom_torl
                    ovlp=true;self_u_proj=[0,1];edg_u_proj=[0,1];
                    return;
                end
            elseif norm(self.vertex_list-flipud(edg.vertex_list)) < geom_torl
                Pnts_self=self.calGeom(21);Pnts_ref=edg.calGeom(21);
                if norm(Pnts_self-flipud(Pnts_ref)) < geom_torl
                    ovlp=true;self_u_proj=[1,0];edg_u_proj=[1,0];
                    return;
                end
            end
            
            % project vertex to curve, if distanc too far, mean no overlap

        end

        function point=calPoint(self,u_x)
            % warp of Curve.calPoint
            %
            point=self.curve.calPoint(u_x);
        end

        function [Points,U]=calGeom(self,u_param)
            % generate point matrix on edge by u_param
            % if u_param is float, adapt calculate point on curve to mean u_param
            %
            if nargin < 2
                u_param=[];
            end

            % using bound box to auto calculate capture precision
            if isempty(u_param), u_param=2^-6*mean(self.bound_box(2,:)-self.bound_box(1,:));end

            if isa(self.curve,'Curve'), region=[min(self.curve.Knots),max(self.curve.Knots)];
            elseif isa(self.curve,'CurveCST'), region=[0,1];end

            low_bou=max(self.bound(1),region(1));
            up_bou=min(self.bound(2),region(2));
            if length(u_param) == 1 && u_param ~= fix(u_param)
                % auto capture Points
                value_torl=u_param;min_level=2;max_level=16;

                [U,Points,~]=GeomApp.meshAdapt1D(@(x) self.curve.calPoint(x),low_bou,up_bou,value_torl,min_level,max_level);
                [U,map]=sort(U);
                Points=Points(map,:);
            else
                % manual capture Points
                if length(u_param) == 1, u_param=linspace(low_bou,up_bou,u_param);end
                U=u_param;
                Points=self.curve.calPoint(u_param);
            end
        end

        function [ln_hdl,sctr_hdl]=displayModel(self,axe_hdl,u_param)
            % display edge on axes
            %
            if nargin < 3
                u_param=[];
                if nargin < 2
                    axe_hdl=[];
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end

            % calculate point on curve
            Points_edg=self.calGeom(u_param);

            % plot boundary
            sctr_hdl=scatterPoints(axe_hdl,self.vertex_list,self.dimension);

            % plot edge
            if self.dimension == 2
                ln_hdl=line(axe_hdl,Points_edg(:,1),Points_edg(:,2));
            else
                ln_hdl=line(axe_hdl,Points_edg(:,1),Points_edg(:,2),Points_edg(:,3));
                zlabel('z');
                view(3);
            end
            xlabel('x');
            ylabel('y');

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