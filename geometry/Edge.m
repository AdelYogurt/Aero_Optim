classdef Edge < handle
    % Topology Entity Edge with Bound
    %
    properties
        name=''; % name of edge
        curve; % Curve handle
        vertex_list=[]; % vertex list of boundary
        outer_bound=[]; % boundary of outer
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

            if isempty(vtx_list)
                if isa(crv,'Curve'), otr_bou=[min(crv.Knots),max(crv.Knots)];
                elseif isa(crv,'CurveCST'), otr_bou=[0,1];end
                vtx_list=crv.calPoint(otr_bou);
            else
                % project vertex list to get outer_bound
                if size(vtx_list,2) ~= 2
                    error('Edge: only can have a pair of vertex');
                end
                [~,otr_bou]=crv.projectPoint(vtx_list);
                if otr_bou(2) < otr_bou(1)
                    otr_bou=otr_bou(end:-1:1);
                    vtx_list=vtx_list(end:-1:1,:);
                end
            end
            self.vertex_list=vtx_list;
            self.outer_bound=otr_bou;

            topo.vertex=[1,2];
            self.topology=topo;

            U_x=linspace(otr_bou(1),otr_bou(2),101);
            Pnts=crv.calPoint(U_x);
            self.bound_box=[min(Pnts,[],1);max(Pnts,[],1)];
            self.dimension=size(Pnts,2);
        end

        function self=reverse(self)
            % revese direction of edge
            %
            self.curve.reverse;
            self.vertex_list=flipud(self.vertex_list);
            
            if isa(self.curve,'Curve'), bou_max=max(self.curve.Knots);
            elseif isa(self.curve,'CurveCST'), bou_max=1;end
            self.outer_bound=bou_max-self.outer_bound(end:-1:1);
        end

        function point=calPoint(self,u_x)
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
            if isempty(u_param), u_param=2^-5*mean(self.bound_box(2,:)-self.bound_box(1,:));end

            low_bou=self.outer_bound(1);
            up_bou=self.outer_bound(2);
            if length(u_param) == 1 && u_param ~= fix(u_param)
                % auto capture Points
                value_torl=u_param;min_level=1;max_level=16;
                dim=self.dimension;

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