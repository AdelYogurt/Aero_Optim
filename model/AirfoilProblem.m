classdef AirfoilProblem < handle
    % initial problem parameter
    properties
        % run paramter
        partitions;

        % geomertry
        class_par_Y=[0.5,1];
        dat_filestr='model/airfoil_deform.dat';

        % DEF parameter
        airfoil_coord;
        mesh_filestr;
        DEF_config_filestr;
        mesh_out_filestr='model/airfoil_deform.su2';

        % CFD
        CFD_config_filestr;
        flight_condition=[];
        restart_filestr=[];

        % info paramete
        out_logger=[];
        run_description=[];
    end

    methods
        % problem setup function
        function self=AirfoilProblem(partitions,...
                airfoil_coord,mesh_filestr,DEF_config_filestr,CFD_config_filestr,...
                dat_filestr,mesh_out_filestr,flight_condition,restart_filestr,...
                out_logger,run_description)
            if nargin < 11
                run_description=[];
                if nargin < 10
                    out_logger=[];
                    if nargin < 9
                        restart_filestr=[];
                        if nargin < 8
                            flight_condition=[];
                            if nargin < 7
                                mesh_out_filestr=[];
                                if nargin < 6
                                    dat_filestr=[];
                                end
                            end
                        end
                    end
                end
            end
                        
            if ~isempty(dat_filestr),self.dat_filestr=dat_filestr;end
            if ~isempty(mesh_out_filestr),self.mesh_out_filestr=mesh_out_filestr;end
            if isempty(flight_condition)
                flight_condition.MACH_NUMBER=0.63;
                flight_condition.SIDESLIP_ANGLE=0.0;
                flight_condition.AOA=2.0;
                flight_condition.FREESTREAM_TEMPERATURE=273.0;
                flight_condition.FREESTREAM_PRESSURE=101325;
            end

            % run paramter
            self.partitions=partitions;

            % geomertry

            % DEF parameter
            self.airfoil_coord=airfoil_coord;
            self.mesh_filestr=mesh_filestr;
            self.DEF_config_filestr=DEF_config_filestr;

            % CFD
            self.CFD_config_filestr=CFD_config_filestr;
            self.flight_condition=flight_condition;
            self.restart_filestr=restart_filestr;

            % info paramete
            self.out_logger=out_logger;
            self.run_description=run_description;
        end

        % CAD CAE function
        function [mesh_point,mesh_out_filestr,SU2_DEF_info,SU2_history,SU2_surface,SU2_CFD_info]=solveAirfoilSU2...
                (self,control_point,airfoil_coord,...
                mesh_filestr,DEF_config_filestr,...
                CFD_config_filestr,flight_condition,...
                dat_filestr,mesh_out_filestr,restart_filestr,run_description)
            % call SU2 program to obtain airfoil aerodynamic parameter
            %
            if nargin < 11
                run_description=[];
                if nargin < 10
                    restart_filestr=[];
                    if nargin < 9
                        mesh_out_filestr=[];
                        if nargin < 8
                            dat_filestr=[];
                        end
                    end
                end
            end

            if isempty(dat_filestr),dat_filestr=self.dat_filestr;end
            if isempty(mesh_out_filestr),mesh_out_filestr=self.mesh_out_filestr;end
            if isempty(flight_condition),flight_condition=self.flight_condition;end
            if isempty(restart_filestr),restart_filestr=self.restart_filestr;end

            % step 1: CAD part
            % rebulid surface mesh point
            airfoil=Airfoil2DCST(control_point.up,control_point.low);
            mesh_point=airfoil.calMeshPoint(airfoil_coord);
            writePoint(mesh_point,dat_filestr,2)

            % step 2: CAD/CAE part
            % deform mesh
            SU2_DEF_info=runSU2DEF...
                (mesh_filestr,DEF_config_filestr,dat_filestr,mesh_out_filestr,self.partitions,[],...
                [],run_description,[],self.out_logger);

            % step 3: CAE part
            % run SU2_CFD to get result
            mesh_filestr=mesh_out_filestr;
            [SU2_history,SU2_surface,SU2_CFD_info]=runSU2CFD...
                (mesh_filestr,CFD_config_filestr,self.partitions,flight_condition,restart_filestr,...
                [],run_description,[],self.out_logger);
        end
    end
end

