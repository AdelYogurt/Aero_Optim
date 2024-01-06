classdef AirfoilModel < handle
    % basical problem parameter
    properties
        % run paramter
        partitions=[];

        % temp file parameter
        dir_temp='';
        REMOVE_TEMP=false(1);

        % info parameter
        run_desc=[];
        out_logger=[];
    end

    % initial input parameter, include software parameter
    properties
        geom_param;
        mesh_param; 
        CFD_param
    end

    % main function
    methods
        function self=AirfoilModel(partitions,geom_param,mesh_param,CFD_param,...
                dir_temp,REMOVE_TEMP,run_desc,out_logger)
            % airfoil model setup function
            %
            % notice:
            % airfoil aerodynamic model need three part:
            % geometry, mesh generater, CFD solver
            %
            % input:
            % partitions, geometry_varargin, mesh_varargin, CFD_varargin
            % dir_temp,REMOVE_TEMP,run_desc,out_logger
            %
            if nargin < 8
                out_logger=[];
                if nargin < 7
                    run_desc=[];
                    if nargin < 6
                        REMOVE_TEMP=[];
                        if nargin < 5
                            dir_temp=[];
                        end
                    end
                end
            end

            % run paramter
            self.partitions=partitions;
            if ~isempty(dir_temp), self.out_logger=dir_temp;end
            if ~isempty(REMOVE_TEMP), self.REMOVE_TEMP=REMOVE_TEMP;end
            if ~isempty(run_desc), self.run_desc=run_desc;end
            if ~isempty(out_logger), self.out_logger=out_logger;end

            % geomertry module setup
            if ~isempty(geom_param)
                switch geom_param.solver
                    case 'CST'
                        if ~isfield(geom_param,'mesh_coord'), error('AirfoilModel: lack mesh_coord');end
                    case 'FFD'

                    case 'HICKS_HENNE'
                        if ~isfield(geom_param,'mesh_point'), error('AirfoilModel: lack mesh_point');end
                        if ~isfield(geom_param,'hicks_hemme_define'), error('AirfoilModel: lack hicks_hemme_define');end
                        if ~isfield(geom_param.hicks_hemme_define,'SCALE'),...
                                geom_param.hicks_hemme_define.SCALE=ones(size(geom_param.hicks_hemme_define.MARKER));end
                        if ~isfield(geom_param.hicks_hemme_define,'MARKER'), error('AirfoilModel: lack MARKER');end
                        if ~isfield(geom_param.hicks_hemme_define,'PARAM'), error('AirfoilModel: lack PARAM');end
                    otherwise
                        error('AirfoilModel: unsupported geometry solver')
                end
            end
            self.geom_param=geom_param;

            % mesh generate module setup
            if ~isempty(mesh_param)
                switch mesh_param.solver
                    case 'SU2_DEF'
                        if ~isfield(mesh_param,'initial_mesh_filestr'), error('AirfoilModel: lack initial_mesh_filestr');end
                        if ~isfield(mesh_param,'SU2_DEF_param'), error('AirfoilModel: lack SU2_DEF_param');end
                        if ~isfield(mesh_param,'dat_filestr'), mesh_param.dat_filestr='model/airfoil.dat';end
                        if ~isfield(mesh_param,'mesh_filestr'), mesh_param.mesh_filestr='model/airfoil.su2';end
                    otherwise
                        error('AirfoilModel: unsupported mesh solver')
                end
            end
            self.mesh_param=mesh_param;

            % CFD solver module setup
            if ~isempty(CFD_param)
                switch CFD_param.solver
                    case 'SU2_CFD'
                        if ~isfield(CFD_param,'SU2_CFD_param'), error('AirfoilModel: lack SU2_CFD_param');end
                        if ~isfield(CFD_param,'restart_filestr'), CFD_param.restart_filestr='';end
                    case 'Fluent'
                        if ~isfield(CFD_param,'fluent_jou_filestr'), error('AirfoilModel: lack fluent_jou_filestr');end
                        if ~isfield(CFD_param,'fluent_dir'), error('AirfoilModel: lack fluent_dir');end
                        if ~isfield(CFD_param,'solver_dimension'), error('AirfoilModel: lack solver_dimension');end
                        if ~isfield(CFD_param,'out_filename'), error('AirfoilModel: lack out_filename');end
                    otherwise
                        error('AirfoilModel: unsupported CFD solver')
                end
            end
            self.CFD_param=CFD_param;
        end

        % solve airfoil aerodynamic model function
        function [geo_out,mesh_out,CFD_out]=solveAirfoil(self,geo_in)
            % call each module to solve aerodynamic model
            %
            % solve flow:
            % (geo_in)>geo_module>(surface_mesh)>mesh_module>(mesh)>CFD_module
            %

            % step 1: geometry module
            geo_out=self.geoModule(geo_in);

            % step 2: mesh module
            % run SU2_DEF deform mesh
            writePoint(geo_out.mesh_point,self.mesh_param.dat_filestr,2)
            if isa(self.mesh_param.SU2_DEF_param,'SU2Config')
                config=self.mesh_param.SU2_DEF_param;
            else
                config=SU2Config(self.mesh_param.SU2_DEF_param);
            end
            config.setParameter('MESH_FILENAME',self.mesh_param.initial_mesh_filestr);
            config.setParameter('DV_FILENAME',self.mesh_param.dat_filestr);
            config.setParameter('MESH_OUT_FILENAME',self.mesh_param.mesh_filestr);
            [mesh_out.mesh_filestr,mesh_out.SU2_DEF_info]=runSU2DEF...
                (config,self.partitions,self.dir_temp,self.run_desc,self.REMOVE_TEMP,self.out_logger);

            % step 3: CFD module
            if strcmp(self.CFD_param.solver,'SU2_CFD')
                % run SU2_CFD to get result
                if isa(self.CFD_param.SU2_CFD_param,'SU2Config')
                    config=self.CFD_param.SU2_CFD_param;
                else
                    config=SU2Config(self.CFD_param.SU2_CFD_param);
                end
                config.setParameter('MESH_FILENAME',mesh_out.mesh_filestr);
                config.setParameter('RESTART_FILENAME',self.CFD_param.restart_filestr);
                [CFD_out.SU2_data,CFD_out.SU2_history,CFD_out.SU2_surface,CFD_out.SU2_CFD_info]=runSU2CFD...
                    (config,self.partitions,self.dir_temp,self.run_desc,self.REMOVE_TEMP,self.out_logger);
            elseif strcmp(self.CFD_param.solver,'Fluent')
                % convert SU2 to cgns
                mesh_data=readMeshSU2(mesh_out.mesh_filestr);
                [airfoil_mesh_dir,airfoil_mesh_filename,~]=fileparts(mesh_out.mesh_filestr);
                if isfield(mesh_data,airfoil_mesh_filename)
                    mesh_data.FLUID=mesh_data.(airfoil_mesh_filename);
                    mesh_data=rmfield(mesh_data,airfoil_mesh_filename);
                end
                mesh_filestr=fullfile(airfoil_mesh_dir,[airfoil_mesh_filename,'.cgns']);
                writeMeshCGNS(mesh_data,mesh_filestr)

                % run fluent to get result
                [CFD_out.fluent_out,CFD_out.fluent_info]=runFluentCFD...
                    (mesh_filestr,self.CFD_param.fluent_jou_filestr,self.partitions,self.CFD_param.fluent_dir,self.CFD_param.solver_dimension,self.CFD_param.out_filename,...
                    self.dir_temp,self.run_desc);
            end
        end

        function drawGeo(self,geo_in,fig_hdl,draw_option)
            % draw specific geometry
            %
            if nargin < 4
                draw_option=struct();
                if nargin < 3
                    fig_hdl=[];
                end
            end

            geo_out=self.geoModule(geo_in);
            displayMesh(geo_out.mesh_point,[],fig_hdl,draw_option)
        end

        function geo_out=geoModule(self,geo_in)
            % calculate geo out
            %
            if strcmp(self.geom_param.solver,'CST')
                % rebulid surface mesh point
                airfoil=AirfoilGeom(geo_in.C_par_low,geo_in.Poles_low,geo_in.C_par_up,geo_in.Poles_up);
                geo_out.mesh_point=airfoil.calMeshPoint(self.geom_param.mesh_coord);
            elseif strcmp(self.geom_param.solver,'FFD')

            elseif strcmp(self.geom_param.solver,'HICKS_HENNE')
                % deform surface mesh point
                self.geom_param.hicks_hemme_define.SIZE=geo_in.SIZE;
                geo_out.mesh_point=AirfoilModel.geoHicksHemme...
                    (self.geom_param.mesh_point,self.geom_param.hicks_hemme_define);
            end
        end
    end

    % parameterization function
    methods(Static)
        function mesh_point=geoHicksHemme(mesh_point,hicks_hemme_define)
            % mesh_point_surface.(marker)
            % marker.index, marker.X, marker.Y
            %
            % hicks_hemme_define.SCALE, hicks_hemme_define.MARKER,...
            % hicks_hemme_define.PARAM, hicks_hemme_define.SIZE
            %
            bump_scale=hicks_hemme_define.SCALE;
            bump_marker=hicks_hemme_define.MARKER;
            bump_direction=hicks_hemme_define.PARAM(:,1);
            bump_location=hicks_hemme_define.PARAM(:,2);
            bump_size=hicks_hemme_define.SIZE;

            bump_E=log(0.5)./log(bump_location);

            % calculate the new coordinates
            for idx=1:length(bump_size)
                marker_name=bump_marker{idx};

                direction=bump_direction(idx);
                if direction == 0, direction=-1;end

                mesh_point.(marker_name).Y=mesh_point.(marker_name).Y+...
                    direction*bump_size(idx)*bump_scale(idx)*sin(pi*mesh_point.(marker_name).X.^bump_E(idx)).^3;
            end
        end
    end
end
