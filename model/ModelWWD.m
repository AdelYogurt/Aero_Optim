classdef ModelWWD < handle
    % basical problem parameter
    properties
        % run paramter
        partitions=[];

        % temp file parameter
        dir_temp='';
        REMOVE_TEMP=false(1);

        % info parameter
        run_description=[];
        out_logger=[];
    end

    % initial input parameter, include software parameter
    properties
        geo_param;
        mesh_param; 
        CFD_param
    end

    % initial problem parameter
    properties
        % CFD
        CFD_config_filestr;
        flight_condition=[];
        restart_filestr=[];
    end

    methods
        % problem setup function
        function self=ModelWWD(partitions,geo_param,mesh_param,CFD_param,...
                dir_temp,REMOVE_TEMP,run_description,out_logger)
            % WWD model setup function
            %
            % notice:
            % airfoil aerodynamic model need three part:
            % geometry, mesh generater, CFD solver
            %
            % input:
            % partitions, geometry_varargin, mesh_varargin, CFD_varargin
            % dir_temp,REMOVE_TEMP,run_description,out_logger
            %
            if nargin < 8
                out_logger=[];
                if nargin < 7
                    run_description=[];
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
            if ~isempty(run_description), self.run_description=run_description;end
            if ~isempty(out_logger), self.out_logger=out_logger;end

            % geomertry module setup
            if ~isempty(geo_param)
                switch geo_param.solver
                    case 'CST'
                        if ~isfield(geo_param,'mesh_coord'), error('ModelAirfoil: lack mesh_coord');end
                    case 'FFD'

                    case 'BSpline'

                    otherwise
                        error('ModelAirfoil: unsupported geometry solver')
                end
            end
            self.geo_param=geo_param;
            
            % mesh generate module setup
            if ~isempty(mesh_param)
                switch mesh_param.solver
                    case 'SU2_DEF'
                        if ~isfield(mesh_param,'initial_mesh_filestr'), error('ModelAirfoil: lack initial_mesh_filestr');end
                        if ~isfield(mesh_param,'SU2_DEF_param'), error('ModelAirfoil: lack SU2_DEF_param');end
                        if ~isfield(mesh_param,'dat_filestr'), mesh_param.dat_filestr='model/WWD.dat';end
                        if ~isfield(mesh_param,'mesh_filestr'), mesh_param.mesh_filestr='model/WWD.su2';end
                    otherwise
                        error('ModelAirfoil: unsupported mesh solver')
                end
            end
            self.mesh_param=mesh_param;
            
            % CFD solver module setup
            if ~isempty(CFD_param)
                switch CFD_param.solver
                    case 'SU2_CFD'
                        if ~isfield(CFD_param,'SU2_CFD_param'), error('ModelAirfoil: lack SU2_CFD_param');end
                        if ~isfield(CFD_param,'restart_filestr'), CFD_param.restart_filestr='';end
                    otherwise
                        error('ModelAirfoil: unsupported CFD solver')
                end
            end
            self.CFD_param=CFD_param;
        end

        % solve WWD aerodynamic model function
        function [geo_out,mesh_out,CFD_out]=solveWWD(self,geo_in)
            % call each module to solve aerodynamic model
            %
            % solve flow:
            % (geo_in)>geo_module>(surface_mesh)>mesh_module>(mesh)>CFD_module
            %

            % step 1: geometry module
            geo_out=self.geoModule(geo_in);

            % step 2: mesh module
            % run SU2_DEF deform mesh
            writePoint(geo_out.mesh_point,self.mesh_param.dat_filestr,3)
            [mesh_out.mesh_filestr,mesh_out.SU2_DEF_info]=runSU2DEF...
                (self.mesh_param.initial_mesh_filestr,self.mesh_param.SU2_DEF_param,self.mesh_param.dat_filestr,self.partitions,...
                self.mesh_param.mesh_filestr,self.dir_temp,self.run_description,self.REMOVE_TEMP,self.out_logger);

            % step 3: CAE part
            % run SU2_CFD to get result
            [CFD_out.SU2_history,CFD_out.SU2_surface,CFD_out.SU2_CFD_info]=runSU2CFD...
                (mesh_out.mesh_filestr,self.CFD_param.SU2_CFD_param,self.partitions,...
                self.CFD_param.restart_filestr,self.dir_temp,self.run_description,self.REMOVE_TEMP,self.out_logger);
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

            if ispc()
                geo_out=self.geoModule(geo_in);
                displayMesh(geo_out.mesh_point,[],fig_hdl,draw_option);
            end
        end

        function geo_out=geoModule(self,geo_in)
            % calculate geo out
            %
            WWD=WaveriderWingDia(geo_in.total_length,geo_in.param);
            geo_out.mesh_point=WWD.calMeshPoint(self.geo_param.mesh_coord);

            mesh_data=WWD.getWGSMesh();
            mesh_data=convertWGSToSTL(mesh_data);
            [area,volume]=calSTLGeometry(mesh_data);
            area=area*2;volume=volume*2;

            geo_out.area=area;
            geo_out.volume=volume;
            geo_out.V_eff=6*sqrt(pi)*volume/area^(3/2);
        end
    end
    
end

