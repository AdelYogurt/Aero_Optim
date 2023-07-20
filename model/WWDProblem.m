classdef WWDProblem < handle
    % initial problem parameter
    properties
        % run paramter
        partitions;

        % geomertry
        total_length=4;
        dat_filestr='model/WWD_deform.dat';

        % DEF parameter
        WWD_coord;
        mesh_filestr;
        DEF_config_filestr;
        mesh_out_filestr='model/WWD_deform.su2';

        % CFD
        CFD_config_filestr;
        flight_condition=[];
        restart_filestr=[];

        % info paramete
        out_logger=[];
        run_description=[];
    end

    % optimal problem define
    properties
        % basical parameter
        vari_num=16;
        con_num=4;
        coneq_num=0;
        low_bou=[0.7,0.7, ...
            2.4,0.4,0.5,0.5, ...
            0.1,0.1,0.005, ...
            0.4,0.4,0.4,0.6,0.5,0.01,0.01];
        up_bou=[1.0,1.0, ...
            3.2,0.8,2.0,2.0, ...
            0.5,0.5,0.02, ...
            0.6,0.6,0.6,0.8,0.7,0.05,0.05];

        % problem constraint define
        CL_min=0.15;
        CD_max=0.03;
        V_eff_min=0.23;
        HF_max=6e6;
    end

    methods
        % problem setup function
        function self=WWDProblem(partitions,...
                WWD_coord,mesh_filestr,DEF_config_filestr,CFD_config_filestr,...
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
                flight_condition.MACH_NUMBER=13.8;
                flight_condition.SIDESLIP_ANGLE=0.0;
                flight_condition.AOA=10.0;
                flight_condition.FREESTREAM_TEMPERATURE=264.9206;
                flight_condition.FREESTREAM_PRESSURE=143.9885;
            end

            % run paramter
            self.partitions=partitions;

            % geomertry
            
            % DEF parameter
            self.WWD_coord=WWD_coord;
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
        function [mesh_point,mesh_out_filestr,SU2_DEF_info,SU2_history,SU2_surface,SU2_CFD_info]=solveWWDSU2...
                (self,geo_par,WWD_coord,...
                mesh_filestr,DEF_config_filestr,...
                CFD_config_filestr,flight_condition,...
                dat_filestr,mesh_out_filestr,restart_filestr,run_description)
            % call SU2 program to obtain WWD aerodynamic parameter
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
            WWD=WaveriderWingDia(self.total_length,geo_par);
            mesh_point=WWD.calMeshPoint(WWD_coord);
            writePoint(mesh_point,dat_filestr,3);

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

        % optimization funciton
        function [obj,con,coneq]=objcon_fcn(self,x)
            % objcon function for optimize
            %
            % input:
            % x
            %
            % output:
            % obj, con, coneq
            %

            % analysis CFD parameter
            geo_par=x;
            [mesh_point,mesh_out_filestr,SU2_DEF_info,SU2_history,SU2_surface,SU2_CFD_info]=self.solveWWDSU2...
                (geo_par,self.WWD_coord,...
                self.mesh_filestr,self.DEF_config_filestr,...
                self.CFD_config_filestr,self.flight_condition,...
                self.dat_filestr,self.mesh_out_filestr,self.restart_filestr,self.run_description);

            % objective
            obj=-SU2_history.('CEff')(end);

            % constaints
            g1=self.CL_min-SU2_history.('CL')(end);
            g2=SU2_history.('CD')(end)-self.CD_max;
            g3=self.V_eff_min-V_eff;
            g4=max(SU2_surface.('HF'))-self.HF_max;
            con=[g1,g2,g3,g4];
            coneq=[];

        end
    end
    
end

