classdef Shell < handle & matlab.mixin.Copyable
    % Topology Entity Shell
    %
    properties
        name=''; % name of shell
        face_list=Face.empty(); % face list
        wire_list=Wire.empty(); % wire list of boundary
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
            self.face_list=fce_list;

            if isempty(fce_list)
                return;
            end

            % check dimension of each face
            fce_num=length(fce_list);
            dimension=fce_list(1).dimension;
            for fce_idx=2:fce_num
                if fce_list(fce_idx).dimension ~= dimension
                    error('Shell: dimension of each face is not equal');
                end
            end
            
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
                % find connect face
                ctnt_mat=false(fce_num,fce_num);
                ctnt_topo=cell(fce_num,fce_num);
                for fce_bas_idx=1:fce_num-1
                    fce_bas=fce_list(fce_bas_idx);
                    for fce_ref_idx=fce_bas_idx+1:length(fce_list)
                        fce_ref=fce_list(fce_ref_idx);

                        % check face connect
                        [ctnt,topo_fca]=fce_bas.checkFaceContinuity(fce_ref,geom_torl);
                        
                        ctnt_mat(fce_bas_idx,fce_ref_idx)=ctnt;
                        ctnt_topo(fce_bas_idx,fce_ref_idx)={topo_fca};
                    end
                end

                % load all wire
                fce_num=length(fce_list);
                wir_list=[];
                for fce_idx=1:fce_num
                    fce=fce_list(fce_idx);
                    wir_list=[wir_list;fce.wire_list];
                end

                topo_shl.continuity=ctnt_mat;
                topo_shl.face=ctnt_topo;
            else
                % load wire_list
                wir_list=fce_list(1).wire_list;
                topo_shl.continuity=[];
                topo_shl.face=[];
            end

            % check if shell is closed
%             if length([edg_topo.face])/edg_num == 2
%                 clsd=true;
%             else
                clsd=false;
%             end

            self.face_list=fce_list;
            self.wire_list=wir_list;
            self.closed=clsd;

            self.topology=topo_shl;
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

        function [Pnts_list,UV_list,Elems_list,Pnts_bous_list,U_bous_list]=calGeom(self,u_param)
            % calculate point on all face
            %
            if nargin < 2
                u_param=[];
            end

            % using bound box to auto calculate capture precision
            if isempty(u_param), u_param=2^-6*mean(self.bound_box(2,:)-self.bound_box(1,:));end

            fce_num=length(self.face_list);
            Pnts_list=cell(fce_num,1);
            UV_list=cell(fce_num,1);
            Elems_list=cell(fce_num,1);
            Pnts_bous_list=cell(fce_num,1);
            U_bous_list=cell(fce_num,1);
            for fce_idx=1:fce_num
                [Pnts,UV,Elems,Pnts_bous,U_bous]=self.face_list(fce_idx).calGeom(u_param);
                Pnts_list{fce_idx}=Pnts;
                UV_list{fce_idx}=UV;
                Elems_list{fce_idx}=Elems;
                Pnts_bous_list{fce_idx}=Pnts_bous;
                U_bous_list{fce_idx}=U_bous;
            end
        end

        function [srf_hdl_list,ln_hdl_list,sctr_hdl]=displayModel(self,axe_hdl,u_param)
            % display shell on axes
            %
            if nargin < 3
                u_param=[];
                if nargin < 2
                    axe_hdl=[];
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end

            % calculate point on face
            [Pnts_list,~,Elems_list,Pnts_bous_list,~]=calGeom(self,u_param);

            % plot edge
            fce_num=length(self.face_list);
            vtx_list=[];
            for fce_idx=1:fce_num
                Pnts=Pnts_list{fce_idx};
                Elems=Elems_list{fce_idx};
                Pnts_bous=Pnts_bous_list{fce_idx};

                pth_hdl(fce_idx)=patchPoints(axe_hdl,Elems,Pnts);
                Pnts_bou=[];bou_num=length(Pnts_bous);
                for bou_idx=1:bou_num
                    Pnts_bou=[Pnts_bou;Pnts_bous{bou_idx}];
                    vtx_list=[vtx_list;Pnts_bous{bou_idx}(1,:)];
                end
                ln_hdl_list(bou_idx)=linePoints(axe_hdl,Pnts_bou,self.dimension);
            end
            sctr_hdl=scatterPoints(axe_hdl,vtx_list,self.dimension);

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

            function pth_hdl=patchPoints(axe_hdl,Elems,Pnts)
                % surface point on different dimension
                %
                srf_option=struct('LineStyle','none');
                hold on;
                pth_hdl=patch(axe_hdl,'faces',Elems,'vertices',Pnts,'facevertexcdata',Pnts(:,3),...
                    'facecolor',get(axe_hdl,'DefaultSurfaceFaceColor'), ...
                    'edgecolor',get(axe_hdl,'DefaultSurfaceEdgeColor'),srf_option);
                hold off;
            end
        end

    end

    % parameterized function
    methods

    end
end