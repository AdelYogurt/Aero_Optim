classdef Face < handle & matlab.mixin.Copyable
    % Topology Entity Face with Bound
    %
    properties
        name=''; % name of face
        surface; % Surface handle
        wire_list=Wire.empty(); % wire list of boundary
        bound=Curve.empty(); % boundary of outer and inner
    end

    properties
        topology=struct(); % association with outer_bound, bound and wire_list

        bound_box=[]; % geometry bounding box
        dimension=[]; % dimension of curve
    end

    % define face
    methods
        function self=Face(name,srf,wire_list)
            % generate face entity
            %
            if nargin < 3, wire_list=[];end
            self.name=name;
            self.surface=srf;

            if isempty(srf)
                return;
            end
            
            if isa(srf,'Surface'), region=[min(srf.UKnots),min(srf.VKnots);max(srf.UKnots),max(srf.VKnots)];
            elseif isa(srf,'SurfaceCST'), region=[0,0;1,1];end
            if isempty(wire_list)
                % generate boundary edge
                edg_1=Edge('',srf.getBoundCurve('u0'));
                edg_2=Edge('',srf.getBoundCurve('1v'));
                edg_3=Edge('',srf.getBoundCurve('u1').reverse);
                edg_4=Edge('',srf.getBoundCurve('0v').reverse);
                wire_list=Wire('',[edg_1;edg_2;edg_3;edg_4]);

                [U_point,V_point]=meshgrid([region(1,1),region(2,1)],[region(1,2),region(2,2)]);
                U_point=U_point';V_point=V_point';
                point=[U_point(:),V_point(:)];
                bou=[
                    Curve([point(1,:);point(2,:)]);
                    Curve([point(2,:);point(4,:)]);
                    Curve([point(4,:);point(3,:)]);
                    Curve([point(3,:);point(1,:)])];
                topo.edge={[1,2,3,4]};
            else
                geom_torl=1e-6;
                
                % fit point on boundary to get outer or inner bound
                wir_num=length(wire_list);
                topo.edge={};bou=[];

                region=[];
                for wir_idx=1:wir_num
                    % project edge to surface
                    wir=wire_list(wir_idx);
                    edg_num=length(wir.edge_list);

                    if wir_idx > 1
                        disp('');
                    end

                    % add edge topology
                    topo.edge{wir_idx}=(1:edg_num)+length(bou);
                    for edg_idx=1:edg_num
                        % project point to surface
                        edg=wir.edge_list(edg_idx);
                        Pnts=edg.calGeom();
                        [U_proj,V_proj]=srf.projectPoint(reshape(Pnts,[],1,size(Pnts,2)));
                        idx=find(vecnorm(squeeze(srf.calPoint(U_proj,V_proj))-Pnts,2,2) < geom_torl);
                        UV_proj=[U_proj,V_proj];

                        % interplote to get boundary curve 
                        bou=[bou;GeomApp.interpPointToCurve(UV_proj(idx,:),edg.curve.Degree)];

                        if isempty(region)
                            region=[min(UV_proj,[],1);max(UV_proj,[],1)];
                        else
                            region=[min([region;UV_proj],[],1);max([region;UV_proj],[],1)];
                        end
                    end
                end
            end
            self.wire_list=wire_list;
            self.bound=bou;
            
            self.topology=topo;

            U_x=linspace(region(1,1),region(2,1),21);
            V_x=linspace(region(1,2),region(2,2),21);
            [U_x,V_x]=meshgrid(U_x,V_x);
            Pnts=reshape(srf.calPoint(U_x(:),V_x(:)),numel(U_x),[]);
            self.bound_box=[min(Pnts,[],1);max(Pnts,[],1)];
            self.dimension=size(Pnts,2);
        end

        function point=calPoint(self,u_x,v_x)
            % warp of Surface.calPoint
            %
            point=self.surface.calPoint(u_x,v_x);
        end

        function [ctnt,topo_fca]=checkFaceContinuity(self,fce,geom_torl)
            % find out if connect between face
            %
            if nargin < 3, geom_torl=[];end
            if isempty(geom_torl), geom_torl=1e-6;end

            % check if bounding box overlap
            if ~GeomApp.checkBoxOverlap(self,fce,geom_torl)
                % if no, pass
                ctnt=false;topo_fca={};
                return;
            end

            % check each wire if overlap
            ctnt=false;
            ovlp_mat=false(length(self.wire_list),length(fce.wire_list));
            ovlp_topo=cell(length(self.wire_list),length(fce.wire_list));
            for wir_ba_idx=1:length(self.wire_list)
                wir_ba=self.wire_list(wir_ba_idx);
                for wir_idx=1:length(fce.wire_list)
                    wir=fce.wire_list(wir_idx);
                    % check wire overlap
                    [ovlp,topo_wir]=wir_ba.checkOverlap(wir,geom_torl);
                    ctnt=ctnt | ovlp;
                    ovlp_mat(wir_ba_idx,wir_idx)=ovlp;
                    ovlp_topo(wir_ba_idx,wir_idx)={topo_wir};
                end
            end
            topo_fca.overlap=ovlp_mat;
            topo_fca.wire=ovlp_topo;
        end

        function [Pnts,UV,Elems,Pnts_bous,U_bous]=calGeom(self,u_param)
            % generate point matrix on face by u_param
            % if u_param is float, adapt calculate point on curve to mean u_param
            %
            if nargin < 2
                u_param=[];
            end

            % using bound box to auto calculate capture precision
            if isempty(u_param), u_param=2^-6*mean(self.bound_box(2,:)-self.bound_box(1,:));end

            srf=self.surface;
            if isa(srf,'Surface'), region=[min(srf.UKnots),min(srf.VKnots);max(srf.UKnots),max(srf.VKnots)];
            elseif isa(srf,'SurfaceCST'), region=[0,0;1,1];end

            % load outer bound curve pole to get region
            Poles_otr_bou=cell(1,length(self.bound));
            poly_otr_bou=[];
            for crv_idx=1:length(self.bound)
                Poles_otr_bou{crv_idx}=self.bound(crv_idx).Poles;
                poly_otr_bou=[poly_otr_bou;self.bound(crv_idx).Poles];
            end
            low_bou=min(poly_otr_bou,[],1);low_bou=max(low_bou,region(1,:));
            up_bou=max(poly_otr_bou,[],1);up_bou=min(up_bou,region(2,:));

            if length(u_param) == 1 && u_param ~= fix(u_param)
                % auto capture Points
                value_torl=u_param;min_level=1;max_level=16;
                dim=self.dimension;

                [UV_srf,Points_srf]=GeomApp.meshAdapt2D(@(x) reshape(self.surface.calPoint(x(:,1),x(:,2)),[],dim),low_bou,up_bou,value_torl,min_level,max_level);
            else
                % manual capture Points
                if isempty(u_param), u_param=21;end

                if length(u_param) == 1, v_param=linspace(low_bou(2),up_bou(2),u_param);end
                if length(u_param) == 1, u_param=linspace(low_bou(1),up_bou(1),u_param);end

                % calculate local coordinate matrix
                [u_param,v_param]=meshgrid(u_param(:),v_param(:));
                UV_srf=[u_param(:),v_param(:)];
                Points_srf=self.calPoint(u_param,v_param);
                Points_srf=reshape(Points_srf,[],self.dimension);
            end

            Pnts_bous=[];
            U_bous=[];

            % calculate outer bound
            wir_num=length(self.wire_list);
            [Pnts_otr_bou,U_otr_bou]=self.wire_list(1).calGeom(u_param);
            Pnts_bous=[Pnts_bous;Pnts_otr_bou];
            U_bous=[U_bous;U_otr_bou];
            
            Pnts_otr_bou=cell2mat(Pnts_otr_bou);poly_otr_bou=[];
            for crv_idx=1:length(U_otr_bou)
                poly_otr_bou=[poly_otr_bou;self.bound(crv_idx).calPoint(U_otr_bou{crv_idx})];
            end

            if wir_num > 1
                % calculate inner bound
                topo_edg=self.topology.edge;
                Pnts_bou_list=cell(wir_num-1,1);
                U_bou_list=cell(wir_num-1,1);
                poly_bou_list=cell(wir_num-1,1);

                % calculate each loop
                for wir_idx=2:wir_num
                    [Pnts_bou,U_bou]=self.wire_list(wir_idx).calGeom(u_param);
                    Pnts_bous=[Pnts_bous;Pnts_bou];
                    U_bous=[U_bous;U_bou];

                    crv_inner_idx=topo_edg{wir_idx};
                    poly_bou=[];
                    for crv_idx=1:length(crv_inner_idx)
                        poly_bou=[poly_bou;self.bound(crv_inner_idx(crv_idx)).calPoint(U_bou{crv_idx})];
                    end

                    Pnts_bou=cell2mat(Pnts_bou);
                    Pnts_bou_list{wir_idx-1}=Pnts_bou;
                    U_bou_list{wir_idx-1}=U_bou;
                    poly_bou_list{wir_idx-1}=poly_bou;
                end
            else
                Pnts_bou_list=[];
                U_bou_list=[];
                poly_bou_list=[];
            end

            % cutting outer area
            [in_region,on_region]=inpolygon(UV_srf(:,1),UV_srf(:,2),poly_otr_bou(:,1),poly_otr_bou(:,2));
            Bool_region=in_region & ~on_region;

            % cutting inner area
            Pnts_bou_all=[];poly_bou_all=[];
            for wir_idx=2:wir_num
                Pnts_bou=Pnts_bou_list{wir_idx-1};
                poly_bou=poly_bou_list{wir_idx-1};

                [in_inner,on_inner]=inpolygon(UV_srf(:,1),UV_srf(:,2),poly_bou(:,1),poly_bou(:,2));
                Bool_region=Bool_region & ~in_inner & ~on_inner;

                Pnts_bou_all=[Pnts_bou_all;Pnts_bou];
                poly_bou_all=[poly_bou_all;poly_bou];
            end

            UV=[poly_bou_all;poly_otr_bou;UV_srf(Bool_region,:)];
            Pnts=[Pnts_bou_all;Pnts_otr_bou;Points_srf(Bool_region,:)];

            % delete duplicate data
            [~,idx,~]=unique(UV,'first','rows');idx=sort(idx);
            UV=UV(idx,:);
            Pnts=Pnts(idx,:);

            % delaunay triangle
            if isempty(poly_bou_all)
                Elems=delaunay(UV);
            else
                Con_DT=[(1:(size(poly_bou_all,1)-1))',(2:size(poly_bou_all,1))';size(poly_bou_all,1),1];
                DT=delaunayTriangulation(UV,Con_DT);
                Elems=DT.ConnectivityList(~isInterior(DT),:);
            end
        end

        function [pth_hdl,ln_hdl_list,sctr_hdl]=displayModel(self,axe_hdl,u_param)
            % display surface on axes
            %
            if nargin < 3
                u_param=[];
                if nargin < 2
                    axe_hdl=[];
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end
            srf_option=struct('LineStyle','none');

            % calculate point on surface
            [Pnts,~,Elems,Pnts_bous,~]=self.calGeom(u_param);

            % plot outer boundary
            bou_num=length(Pnts_bous);
            vtx_list=zeros(bou_num,self.dimension);
            for bou_idx=1:bou_num
                Pnts_bou=Pnts_bous{bou_idx};
                vtx_list(bou_idx,:)=Pnts_bou(1,:);
                ln_hdl_list(bou_idx)=linePoints(axe_hdl,Pnts_bou,self.dimension);
            end
            sctr_hdl=scatterPoints(axe_hdl,vtx_list,self.dimension);

            % plot face
            hold on;
            if self.dimension == 2
                pth_hdl=patch(axe_hdl,'faces',Elems,'vertices',Pnts,'facevertexcdata',Pnts(:,3),...
                    'facecolor',get(axe_hdl,'DefaultSurfaceFaceColor'), ...
                    'edgecolor',get(axe_hdl,'DefaultSurfaceEdgeColor'),srf_option);
            else
                pth_hdl=patch(axe_hdl,'faces',Elems,'vertices',Pnts,'facevertexcdata',Pnts(:,3),...
                    'facecolor',get(axe_hdl,'DefaultSurfaceFaceColor'), ...
                    'edgecolor',get(axe_hdl,'DefaultSurfaceEdgeColor'),srf_option);
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