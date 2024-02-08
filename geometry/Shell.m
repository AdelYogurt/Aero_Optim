classdef Shell < handle
    % Topology Entity Shell
    %
    properties
        name=''; % name of shell
        face_list=Face.empty(); % face list
        edge_list=Edge.empty(); % edge list
        vertex_list=[]; % vertex list
        closed=false; % if shell is closed
    end

    properties
        topology=struct(); % association with face

        bound_box=[]; % geometry bounding box
        dimension=[]; % dimension of curve
    end

    % define shell
    methods
        function self=Shell(name,fce_list,geom_torl)
            % generate shell entity
            %
            if nargin < 3, geom_torl=[];end
            if isempty(geom_torl), geom_torl=1e-6;end
            self.name=name;

            % check dimension of each face
            fce_num=length(fce_list);
            dimension=fce_list(1).dimension;
            for fce_idx=2:fce_num
                if fce_list(fce_idx).dimension ~= dimension
                    error('Shell: dimension of each face is not equal');
                end
            end
            self.face_list=fce_list;

            % generate bounding box
            bound_box=fce_list(1).bound_box;
            for fce_idx=2:fce_num
                bound_box(1,:)=min([bound_box(1,:);fce_list(fce_idx).bound_box(1,:)],[],1);
                bound_box(2,:)=max([bound_box(2,:);fce_list(fce_idx).bound_box(2,:)],[],1);
            end
            self.bound_box=bound_box;
            self.dimension=dimension;

            % calculate topology of shell
            self.calTopology(geom_torl);
        end

        function self=calTopology(self,geom_torl)
            % calculate topology property of wire
            %
            fce_list=self.face_list;

            fce_num=length(fce_list);
            if fce_num > 1
                % check whether face is connect
                edg_list=[];vtx_list=[];
                for fce_idx=1:fce_num
                    edg_list=[edg_list;fce_list(fce_idx).wire_list(1).edge_list];
                    vtx_list=[vtx_list;fce_list(fce_idx).wire_list(1).vertex_list];
                end

                % generate association between edge and vertex
                vtx_num=size(vtx_list,1);
                vtx_topo=repmat(struct('edge',[],'face',[]),vtx_num,1);

                edg_num=size(edg_list,1);
                edg_topo=repmat(struct('face',[],'map',[]),edg_num,1);
            else
                edg_list=[];
                vtx_list=[];

                % load all edge and vertex
                wir_list=fce_list(1).wire_list;
                for bou_idx=1:length(wir_list)
                    edg_list=[edg_list;wir_list(bou_idx).edge_list];
                    vtx_list=[vtx_list;wir_list(bou_idx).vertex_list];
                end

                % generate association between edge and vertex
                vtx_num=size(vtx_list,1);
                vtx_topo=repmat(struct('edge',1,'face',1),vtx_num,1);

                edg_num=size(edg_list,1);
                edg_topo=repmat(struct('face',1,'map',[0,1]),edg_num,1);
            end

            % check if shell is closed
            if length([edg_topo.face])/edg_num == 2
                clsd=true;
            else
                clsd=false;
            end

            self.face_list=fce_list;
            self.edge_list=edg_list;
            self.vertex_list=vtx_list;
            self.closed=clsd;

            topo.vertex=vtx_topo;
            topo.edge=edg_topo;
            self.topology=topo;
        end

        function self=addFace(self,fce)
            % add face into shell
            %
            self.face_list=[self.face_list;fce];
        end

        function [fce,fce_idx]=getFace(self,fce_name)
            % load face from face_list base on input fce_name
            %
            fce=[];
            fce_idx=[];
            for fce_i=1:length(self.face_list)
                fce_u=self.face_list(fce_i);
                if strcmp(fce_u.name,fce_name)
                    fce=[fce;fce_u];
                    fce_idx=[fce_idx;fce_i];
                end
            end
        end

        function fce_list=calGeom(self,u_param,v_param)
            % calculate point on all face
            %
            if nargin < 3
                v_param=[];
                if nargin < 2
                    u_param=[];
                end
            end

            if length(u_param) == 1 && u_param ~= fix(u_param)
                % auto capture Points

            else
                % manual capture Points
                fce_list=[];fce_num=length(self.face_list);
                for fce_idx=1:fce_num
                    fce=self.face_list(fce_idx);
                    fce_list=[fce_list;{fce.calGeom(u_param,v_param)}];
                end
            end
        end

        function [srf_hdl_list,ln_hdl_list,sctr_hdl]=displayModel(self,axe_hdl,u_param,v_param)
            % display shell on axes
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
            if isempty(axe_hdl),axe_hdl=axes(figure());end

            fce_num=length(self.face_list);
            for fce_idx=1:fce_num
                [Points,U,V]=self.face_list(fce_idx).calGeom(u_param,v_param);
                srf_hdl_list(fce_idx)=surfacePoints(axe_hdl,Points,self.dimension);

                % plot outer boundary
                wir=self.face_list(fce_idx).wire_list(1);
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
                        Points=self.face_list(fce_idx).wire_list(edg_idx).calGeom(u_param);
                        ln_hdl_list{1}(edg_idx)=linePoints(axe_hdl,Points,self.dimension);
                    end
                end

                % plot inner boundary
                bou_num=length(self.face_list(fce_idx).wire_list);
                for bou_idx=2:bou_num
                    wir=self.face_list(fce_idx).wire_list(bou_idx);
                    for edg_idx=1:length(wir.edge_list)
                        Points=self.face_list(fce_idx).wire_list(edg_idx).calGeom(u_param);
                        ln_hdl_list{bou_idx}(edg_idx)=linePoints(axe_hdl,Points,self.dimension);
                    end
                end
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

            function srf_hdl=surfacePoints(axe_hdl,Points,dimension)
                % surface point on different dimension
                %
                srf_option=struct('LineStyle','none');
                hold on;
                if dimension == 2
                    srf_hdl=surf(axe_hdl,Points(:,:,1),Points(:,:,2),srf_option);
                else
                    srf_hdl=surf(axe_hdl,Points(:,:,1),Points(:,:,2),Points(:,:,3),srf_option);
                    zlabel('z');
                    view(3);
                end
                hold off;
            end
        end

    end

    % parameterized function
    methods

    end
end