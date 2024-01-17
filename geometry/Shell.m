classdef Shell < handle
    % shell
    %
    properties
        name;
        face_list; % cell, can be surface CST3D or surface BSpline
    end

    % define Shell
    methods
        function self=Shell(name,face_list)
            if nargin < 2,face_list={};end
            self.name=name;
            self.face_list=face_list;
        end

        function addFace(self,fce)
            % add face into shell
            %
            self.face_list=[self.face_list;{fce}];
        end

        function fce_list=calShell(self,u_param,v_param)
            % calculate all face point
            %
            if nargin < 3
                v_param=[];
                if nargin < 2
                    u_param=[];
                end
            end

            fce_list=[];fce_num=length(self.face_list);
            for fce_idx=1:fce_num
                fce=self.face_list{fce_idx};
                if ~isempty(fce)
                    fce_list=[fce_list,{fce.calFace(u_param,v_param)}];
                end
            end
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
    end

    % generate mesh
    methods
        function mesh_data=getMeshWGS(self,u_param,v_param)
            % generate LaWGS mesh data
            %
            if nargin < 3
                v_param=[];
                if nargin < 2
                    u_param=[];
                end
            end

            fce_list=self.calShell(u_param,v_param);
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

    % visualizate function
    methods
        function gplot(self,axe_hdl,u_param,v_param,crv_option,ctrl_option)
            % draw all surface of shell
            % wrapper of gplot
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

            surf_num=length(self.face_list);
            for surf_idx=1:surf_num
                surf=self.face_list{surf_idx};
                if ~isempty(surf)
                    surf.gplot(axe_hdl,u_param,v_param,crv_option,ctrl_option);
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