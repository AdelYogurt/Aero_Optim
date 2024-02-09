classdef Face < handle
    % Topology Entity Face with Bound
    %
    properties
        name=''; % name of face
        surface; % Surface handle
        wire_list=Wire.empty(); % wire list of boundary
        outer_bound=Curve.empty(); % boundary of outer and inner
        bound={}; % boundary of inner
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

            if isa(srf,'Surface'), region=[min(srf.UKnots),min(srf.VKnots);max(srf.UKnots),max(srf.VKnots)];
            elseif isa(srf,'SurfaceCST'), region=[0,0;1,1];end
            if isempty(wire_list)
                % generate boundary edge
                edg_1=Edge('',srf.getBoundCurve('u0'));
                edg_2=Edge('',srf.getBoundCurve('1v'));
                edg_3=Edge('',srf.getBoundCurve('u1').reverse);
                edg_4=Edge('',srf.getBoundCurve('0v').reverse);
                self.wire_list=Wire('',[edg_1;edg_2;edg_3;edg_4]);

                [U_point,V_point]=meshgrid([region(1,1),region(2,1)],[region(1,2),region(2,2)]);
                U_point=U_point';V_point=V_point';
                point=[U_point(:),V_point(:)];
                otr_bou=[
                    Curve([point(1,:);point(2,:)]);
                    Curve([point(2,:);point(4,:)]);
                    Curve([point(4,:);point(3,:)]);
                    Curve([point(3,:);point(1,:)])];
                bou={};
                topo.edge={[1,2,3,4]};
            else
                % fit point on boundary to get outer bound
                wir_num=length(wire_list);
                topo.edge={};

                % if exist inner boundary
                if wir_num > 1
                end
            end
            self.outer_bound=otr_bou;
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
            point=self.surface.calPoint(u_x,v_x);
        end

        function [Points,UV,Elements,Points_bound_list,U_bound_list]=calGeom(self,u_param)
            % generate point matrix on face by u_param
            % if u_param is float, adapt calculate point on curve to mean u_param
            %
            if nargin < 2
                u_param=[];
            end
            geom_torl=1e-6;

            % using bound box to auto calculate capture precision
            if isempty(u_param), u_param=2^-5*mean(self.bound_box(2,:)-self.bound_box(1,:));end

            % load outer bound curve pole to get region
            Poles_otr_bou=cell(1,length(self.outer_bound));
            poly_otr_bou=[];
            for crv_idx=1:length(self.outer_bound)
                Poles_otr_bou{crv_idx}=self.outer_bound(crv_idx).Poles;
                poly_otr_bou=[poly_otr_bou;self.outer_bound(crv_idx).Poles];
            end
            low_bou=min(poly_otr_bou,[],1);
            up_bou=max(poly_otr_bou,[],1);

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
                UV_srf=cat(3,u_param,v_param);
                Points_srf=self.calPoint(u_param,v_param);

                if length(self.wire_list) == 1 && length(self.outer_bound) == 4
                    % check if just cutting plane
                    Poles=Poles_otr_bou{1};
                    if all(abs(Poles(:,1)-low_bou(1)) < geom_torl)
                        just_cut=true;check_order=[2,3,4,1];
                    elseif all(abs(Poles(:,1)-up_bou(1)) < geom_torl)
                        just_cut=true;check_order=[4,1,2,3];
                    elseif all(abs(Poles(:,2)-low_bou(2)) < geom_torl)
                        just_cut=true;check_order=[1,2,3,4];
                    elseif all(abs(Poles(:,2)-up_bou(2)) < geom_torl)
                        just_cut=true;check_order=[3,4,1,2];
                    else
                        just_cut=false;
                    end

                    if just_cut
                        Poles_otr_bou=Poles_otr_bou(check_order);
                        % check each poles
                        if all(abs(Poles_otr_bou{1}(:,2)-low_bou(2)) < geom_torl) && ...
                            all(abs(Poles_otr_bou{2}(:,1)-up_bou(1)) < geom_torl) && ...
                            all(abs(Poles_otr_bou{3}(:,2)-up_bou(2)) < geom_torl) && ...
                            all(abs(Poles_otr_bou{4}(:,1)-low_bou(1)) < geom_torl)
                            Elements=[];
                            Points_bound_list={
                                reshape(Points_srf(1,:,:),[],self.dimension);
                                reshape(Points_srf(:,1,:),[],self.dimension);
                                reshape(Points_srf(end,end:-1:1,:),[],self.dimension);
                                reshape(Points_srf(end:-1:1,1,:),[],self.dimension)};
                            U_bound_list={
                                reshape(UV_srf(1,:,:),[],self.dimension);
                                reshape(UV_srf(:,1,:),[],self.dimension);
                                reshape(UV_srf(end,end:-1:1,:),[],self.dimension);
                                reshape(UV_srf(end:-1:1,1,:),[],self.dimension)};
                            return;
                        end
                    else
                        % convert into list for generta triangle element
                        UV_srf=reshape(UV_srf,[],2);
                        Points_srf=reshape(Points_srf,[],self.dimension);
                    end
                end
            end

            Points_bound_list=[];
            U_bound_list=[];

            % calculate outer bound
            wir_num=length(self.wire_list);
            [Pnts_otr_bou,U_otr_bou]=self.wire_list(1).calGeom(u_param);
            Points_bound_list=[Points_bound_list;Pnts_otr_bou];
            U_bound_list=[U_bound_list;U_otr_bou];
            
            Pnts_otr_bou=connectCellToMat(Pnts_otr_bou);poly_otr_bou=[];
            for crv_idx=1:length(self.outer_bound)
                poly_otr_bou=[poly_otr_bou;self.outer_bound(crv_idx).calPoint(U_otr_bou{crv_idx})];
            end

            if wir_num > 1
                % calculate inner bound
                topo_edg=self.topology.edge;
                Pnts_bou_list=cell(wir_num-1,1);
                U_bou_list=cell(wir_num-1,1);
                poly_bou_list=cell(wir_num-1,1);

                % calculate each loop
                for wir_idx=2:wir_num
                    [Pnts_bou,U_bou]=self.wire_list(wir_idx-1).calGeom(u_param);
                    Points_bound_list=[Points_bound_list;Pnts_bou];
                    U_bound_list=[U_bound_list;U_bou];

                    crv_inner_idx=topo_edg{wir_idx};
                    poly_bou=[];
                    for crv_idx=1:length(crv_inner_idx)
                        poly_bou=[poly_bou;self.bound(crv_inner_idx(crv_idx)).calPoint(U_otr_bou{crv_idx})];
                    end

                    Pnts_bou=connectCellToMat(Pnts_bou);
                    Pnts_bou_list{wir_idx-1}=Pnts_bou;
                    Pnts_bou_list{wir_idx-1}=U_bou;
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
            Points=[Pnts_bou_all;Pnts_otr_bou;Points_srf(Bool_region,:)];

            % delete duplicate data
            [~,I,~]=unique(UV,'first','rows');I=sort(I);
            UV=UV(I,:);
            Points=Points(I,:);

            % delaunay triangle
            if isempty(poly_bou_all)
                Elements=delaunay(UV);
            else
                Con_DT=[(1:(size(poly_bou_all,1)-1))',(2:size(poly_bou_all,1))';size(poly_bou_all,1),1];
                DT=delaunayTriangulation(UV,Con_DT);
                Elements=DT.ConnectivityList(~isInterior(DT),:);
            end

            function M=connectCellToMat(C)
                % connect matrix in cell to matrix
                M=[];
                for C_idx=1:length(C)
                    M=[M;C{C_idx}];
                end
            end
        end

        function [srf_hdl,ln_hdl_list,sctr_hdl]=displayModel(self,axe_hdl,u_param)
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
            [Pnts,~,Elems,Pnts_bou_list,~]=self.calGeom(u_param);

            % plot outer boundary
            bou_num=length(Pnts_bou_list);
            vtx_list=zeros(bou_num,self.dimension);
            for bou_idx=1:bou_num
                Pnts_bou=Pnts_bou_list{bou_idx};
                vtx_list(bou_idx,:)=Pnts_bou(1,:);
                ln_hdl_list(bou_idx)=linePoints(axe_hdl,Pnts_bou,self.dimension);
            end
            sctr_hdl=scatterPoints(axe_hdl,vtx_list,self.dimension);

            % plot face
            hold on;
            if isempty(Elems)
                % output is matrix
                if self.dimension == 2
                    srf_hdl=surf(axe_hdl,Pnts(:,:,1),Pnts(:,:,2),srf_option);
                else
                    srf_hdl=surf(axe_hdl,Pnts(:,:,1),Pnts(:,:,2),Pnts(:,:,3),srf_option);
                    zlabel('z');
                    view(3);
                end
            else
                % output is triangle element
                if self.dimension == 2
                    srf_hdl=trisurf(axe_hdl,Elems,Pnts(:,1),Pnts(:,2),srf_option);
                else
                    srf_hdl=patch(axe_hdl,'faces',Elems,'vertices',Pnts,'facevertexcdata',Pnts(:,3),...
                        'facecolor',get(axe_hdl,'DefaultSurfaceFaceColor'), ...
                        'edgecolor',get(axe_hdl,'DefaultSurfaceEdgeColor'),srf_option);
%                     srf_hdl=trisurf(axe_hdl,Elems,Pnts(:,1),Pnts(:,2),Pnts(:,3),srf_option);
                    zlabel('z');
                    view(3);
                end
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