classdef Shell < handle
    % shell
    %
    properties
        name='';
        face_list={}; % cell, can be surface CST3D or surface BSpline
        topo_list=[]; % define the geometry topology
    end

    % define Shell
    methods
        function self=Shell(name,face_list)
            if nargin < 2,face_list={};end
            self.name=name;
            self.face_list=face_list;
        end

        function [surf,surf_idx]=getFace(self,surf_name)
            % load surface from face_list base on input surface name
            %
            for surf_idx=1:length(self.face_list)
                surf=self.face_list{surf_idx};
                if strcmp(surf.name,surf_name)
                    return;
                end
            end
            surf=[];
            surf_idx=0;
        end

        function addFace(self,fce_new)
            % add face into shell
            %
            if ~iscell(fce_new),fce_new={fce_new};end
            
            for fce_idx=1:length(fce_new)
                self.face_list=[self.face_list;fce_new(fce_idx)];
                self.topo_list=[self.topo_list;zeros(1,1)];
            end
        end

        function addTopo(self,fceA,fceB,idx)
            % add new topology into shell
            %
            idxA=find(self.face_list == fceA);
            idxB=find(self.face_list == fceB);
            data=[idxA,idxB,idx];
            self.topo_list=[self.topo_list;data];
        end
    end

    % control Face
    methods
        function self=translate(self,tran_vctr)
            for fce_idx=1:length(self.face_list)
                self.face_list{fce_idx}.translate(tran_vctr);
            end
        end

        function self=rotate(self,rotate_matrix)
            for fce_idx=1:length(self.face_list)
                self.face_list{fce_idx}.rotate(rotate_matrix);
            end
        end
    end

    % calculate point
    methods
        function [srf_list,U,V]=calGeom(self,u_param,v_param)
            % calculate all face point
            %
            if nargin < 3
                v_param=[];
                if nargin < 2
                    u_param=[];
                end
            end

            fce_num=length(self.face_list);
            srf_list=cell(1,fce_num);
            U_list=cell(1,fce_num);
            V_list=cell(1,fce_num);
            for fce_idx=1:fce_num
                fce=self.face_list{fce_idx};
                [Pnt,U,V]=fce.calGeom(u_param,v_param);
                srf_list{fce_idx}=Pnt;
                U_list{fce_idx}=U;
                V_list{fce_idx}=V;
            end

            % adjust U and V to satisfy topology tight
            U=[];V=[];   
            topo_num=size(self.topo_list,1);
            for topo_idx=topo_num
                data=self.topo_list(topo_idx,:);
            end

        end

        function mesh_data=getMeshWGS(self,u_param,v_param)
            % generate LaWGS mesh data
            %
            if nargin < 3
                v_param=[];
                if nargin < 2
                    u_param=[];
                end
            end

            fce_list=self.calGeom(u_param,v_param);
            mesh_data=struct();
            unknown_idx=1;fce_num=length(self.face_list);
            for fce_idx=1:fce_num
                fce=self.face_list{fce_idx};
                if ~isempty(fce)
                    Pnt=fce_list{fce_idx};
                    srf.X=Pnt(:,:,1);
                    srf.Y=Pnt(:,:,2);
                    srf.Z=Pnt(:,:,3);
                    if fce.reverse
                        srf.X=flipud(srf.X);
                        srf.Y=flipud(srf.Y);
                        srf.Z=flipud(srf.Z);
                    end
                    srf.type='wgs';
                    srf_name=fce.name;
                    if isempty(srf_name)
                        srf_name=['srf_',num2str(unknown_idx)];
                        unknown_idx=unknown_idx+1;
                    end
                    mesh_data.(srf_name)=srf;
                end
            end
        end
    end

    % parameterized function
    methods
        function coord=calCoord(self,point_list,fce_index_list)
            % calculate all input point local coordinate
            %
            coord=struct();

            fce_name_list=fieldnames(fce_index_list);
            for fce_idx=1:length(fce_name_list)
                % get surf
                fce_name=fce_name_list{fce_idx};
                fce=self.getFace(fce_name);
                if isempty(fce)
                    continue;
                end
                point_idx=fce_index_list.(fce_name);

                % calculate coordinate
                point=point_list(point_idx,1:3);
                [U,V,~]=fce.calCoord(reshape(point,[],1,3));
                coord.(fce_name).index=point_idx;
                coord.(fce_name).U=U;
                coord.(fce_name).V=V;
            end

        end

        function mesh_data=calMeshPoint(self,coord)
            % calculate all mesh point by input W coord
            %
            mesh_data=struct();

            surf_name_list=fieldnames(coord);
            for surf_idx=1:length(surf_name_list)
                % get surface
                surf_name=surf_name_list{surf_idx};
                surf=self.getFace(surf_name);
                if isempty(surf)
                    continue;
                end
                point_idx=coord.(surf_name).index;
                U=coord.(surf_name).U;
                V=coord.(surf_name).V;

                % calculate coordinate
                Pnts=surf.calPoint(U,V);
                mesh_data.(surf_name).type='scatter';
                mesh_data.(surf_name).index=point_idx;
                mesh_data.(surf_name).X=Pnts(:,:,1);
                mesh_data.(surf_name).Y=Pnts(:,:,2);
                mesh_data.(surf_name).Z=Pnts(:,:,3);
            end

        end
    end

    % visualizate function
    methods
        function plotGeom(self,axe_hdl,u_param,v_param,crv_option,ctrl_option)
            % draw all surface of shell
            % wrapper of plotGeom
            %
            if nargin < 6
                ctrl_option=[];
                if nargin < 5
                    crv_option=[];
                    if nargin < 4
                        v_param=[];
                        if nargin < 3
                            u_param=[];
                            if nargin < 2
                                axe_hdl=[];
                            end
                        end
                    end
                end
            end

            if isempty(axe_hdl),axe_hdl=axes(figure());end

            fce_num=length(self.face_list);
            for fce_idx=1:fce_num
                fce=self.face_list{fce_idx};
                if ~isempty(fce)
                    fce.plotGeom(axe_hdl,u_param,v_param,crv_option,ctrl_option);
                end
            end
            
            xlabel('x');
            ylabel('y');
            zlabel('z');
            view(3);

            % axis equal;
            % x_range=xlim();
            % y_range=ylim();
            % z_range=zlim();
            % center=[mean(x_range),mean(y_range),mean(z_range)];
            % range=max([x_range(2)-x_range(1),y_range(2)-y_range(1),z_range(2)-z_range(1)])/2;
            % xlim([center(1)-range,center(1)+range]);
            % ylim([center(2)-range,center(2)+range]);
            % zlim([center(3)-range,center(3)+range]);
        end
    end
end