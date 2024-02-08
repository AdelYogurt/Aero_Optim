classdef Face < handle
    % Topology Entity Face with Bound
    %
    properties
        name=''; % name of face
        surface; % Surface handle
        outer_bound; % boundary of outer
        bound; % boundary of inner
        wire_list=Wire.empty(); % vertex list of boundary
    end

    properties
        topology=struct(); % association with outer_bound, bound and wire_list

        bound_box=[]; % geometry bounding box
        dimension=[]; % dimension of curve
    end

    % define face
    methods
        function self=Face(name,srf,otr_bou,bou)
            if nargin < 4, bou=[];if nargin < 3, otr_bou=[];end;end
            if isempty(otr_bou)
                if isa(srf,'Surface'), otr_bou=[min(srf.UKnots),min(srf.VKnots);max(srf.UKnots),max(srf.VKnots)];
                elseif isa(srf,'SurfaceCST'), otr_bou=[0,0;1,1];end
            end
            self.name=name;
            self.surface=srf;
            self.outer_bound=otr_bou;
            self.bound=bou;

            if isnumeric(otr_bou)
                % generate boundary edge
                edg_1=Edge('',srf.getBoundCurve('u0'));
                edg_2=Edge('',srf.getBoundCurve('1v'));
                edg_3=Edge('',srf.getBoundCurve('u1').reverse);
                edg_4=Edge('',srf.getBoundCurve('0v').reverse);
                self.wire_list=Wire('',[edg_1;edg_2;edg_3;edg_4]);
                topo.edge={[1,2,3,4]};
            else
                % fit point on boundary to get edge
                topo.edge={};
            end
            if ~isempty(bou)
                % add inner edge
                topo.edge=[topo.edge,{}];
            end
            
            self.topology=topo;

            U_x=linspace(otr_bou(1,1),otr_bou(2,1),21);
            V_x=linspace(otr_bou(1,2),otr_bou(2,2),21);
            [U_x,V_x]=meshgrid(U_x,V_x);
            Pnts=reshape(srf.calPoint(U_x(:),V_x(:)),numel(U_x),[]);
            self.bound_box=[min(Pnts,[],1);max(Pnts,[],1)];
            self.dimension=size(Pnts,2);
        end

        function point=calPoint(self,u_x,v_x)
            point=self.surface.calPoint(u_x,v_x);
        end

        function [Points,U,V]=calGeom(self,u_param,v_param)
            % generate point matrix on face by u_param, v_param
            % if u_param is float, adapt calculate point on curve to mean u_param
            %
            if nargin < 3
                v_param=[];
                if nargin < 2
                    u_param=[];
                end
            end

            % using bound box to auto calculate capture precision
            if isempty(u_param), u_param=2^-5*mean(self.bound_box(2,:)-self.bound_box(1,:));end

            low_bou=self.outer_bound(1,:);
            up_bou=self.outer_bound(2,:);
            if length(u_param) == 1 && u_param ~= fix(u_param)
                % auto capture Points
                value_torl=u_param;min_level=4;max_level=16;
                dim=self.dimension;

                [U,V,Points,~]=GeomApp.meshAdapt2DUV(@(x) reshape(self.surface.calPoint(x(:,1),x(:,2)),[],dim),low_bou,up_bou,value_torl,min_level,max_level,dim);
            else
                % manual capture Points
                if isempty(u_param), u_param=21;end
                if isempty(v_param), v_param=21;end

                if length(u_param) == 1, u_param=linspace(low_bou,up_bou,u_param);end
                if length(v_param) == 1, v_param=linspace(low_bou,up_bou,v_param);end

                % calculate local coordinate matrix
                [u_param,v_param]=meshgrid(u_param(:),v_param(:));
                U=u_param;V=v_param;
                [Points]=self.calPoint(u_param,v_param);
            end
        end

        function [srf_hdl,ln_hdl_list,sctr_hdl_list]=displayModel(self,axe_hdl,u_param,v_param)
            % display surface on axes
            %
            if nargin < 4
                v_param=[];
                if nargin < 3
                    u_param=[];
                    if nargin < 2
                        axe_hdl=[];
                    end
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end
            srf_option=struct('LineStyle','none');

            % calculate point on surface
            [Points_fac,U,V]=self.calGeom(u_param,v_param);

            % plot outer boundary
            wir=self.wire_list(1);
            if length(wir.edge_list) == 4
                Points=wir.edge_list(1).calPoint(U(1,:));
                ln_hdl_list{1}(1)=linePoints(axe_hdl,Points,self.dimension);
                Points=wir.edge_list(2).calPoint(V(:,1));
                ln_hdl_list{1}(2)=linePoints(axe_hdl,Points,self.dimension);
                Points=wir.edge_list(3).calPoint(U(1,end:-1:1));
                ln_hdl_list{1}(3)=linePoints(axe_hdl,Points,self.dimension);
                Points=wir.edge_list(4).calPoint(V(end:-1:1,1));
                ln_hdl_list{1}(4)=linePoints(axe_hdl,Points,self.dimension);
            else
                for edg_idx=1:length(wir.edge_list)
                    Points=self.wire_list(edg_idx).calGeom(u_param);
                    ln_hdl_list{1}(edg_idx)=linePoints(axe_hdl,Points,self.dimension);
                end
            end
            sctr_hdl_list{1}=scatterPoints(axe_hdl,wir.vertex_list,self.dimension);

            % plot inner boundary
            bou_num=length(self.wire_list);
            for bou_idx=2:bou_num
                wir=self.wire_list(bou_idx);
                for edg_idx=1:length(wir.edge_list)
                    Points=self.wire_list(edg_idx).calGeom(u_param);
                    ln_hdl_list{bou_idx}(edg_idx)=linePoints(axe_hdl,Points,self.dimension);
                end
                sctr_hdl_list{bou_idx}=scatterPoints(axe_hdl,wir.vertex_list,self.dimension);
            end

            % plot face
            hold on;
            if self.dimension == 2
                srf_hdl=surf(axe_hdl,Points_fac(:,:,1),Points_fac(:,:,2),srf_option);
            else
                srf_hdl=surf(axe_hdl,Points_fac(:,:,1),Points_fac(:,:,2),Points_fac(:,:,3),srf_option);
                zlabel('z');
                view(3);
            end
            hold off;
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