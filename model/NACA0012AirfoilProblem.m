classdef NACA0012AirfoilProblem < AirfoilProblem
    % optimal problem define
    properties
        % basical parameter
        vari_num=10;
        con_num=2;
        coneq_num=0;
        low_bou=[0.5040,2.1551,-0.4582,1.1015,0.1962,...
            0.5040,2.1551,-0.4582,1.1015,0.1962]*1e-2;
        up_bou=[10.5040,12.1551,9.5418,11.1015,10.1962,...
            10.5040,12.1551,9.5418,11.1015,10.1962]*1e-2;

        % problem constraint define
        t_mim=0.096;
        CL_min=0.30;
    end

    methods
        % problem setup function
        function self=NACA0012AirfoilProblem(partitions,...
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

            % setup initial problem parameter
            self@AirfoilProblem(partitions,...
                airfoil_coord,mesh_filestr,DEF_config_filestr,CFD_config_filestr,...
                dat_filestr,mesh_out_filestr,flight_condition,restart_filestr,...
                out_logger,run_description);
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
            x_control_number=length(x)/2;

            % generate Bezier controll point
            control_point_up=zeros(x_control_number+2,2);
            control_point_up(1,:)=[0.000000000000,0.067924153974];
            control_point_up(end,:)=[1.000000000000,0.055565107922];
            control_point_up(:,1)=linspace(0,1,x_control_number+2);
            control_point_low=zeros(x_control_number+2,2);
            control_point_low(1,:)=[0.000000000000,0.067924153974];
            control_point_low(end,:)=[1.000000000000,0.055565107922];
            control_point_low(:,1)=linspace(0,1,x_control_number+2);

            % add x
            control_point_up(2:end-1,2)=x(1:x_control_number);
            control_point_low(2:end-1,2)=x(x_control_number+1:end);

            % analysis CFD parameter
            control_point.up=control_point_up;
            control_point.low=control_point_low;
            [mesh_point,mesh_out_filestr,SU2_DEF_info,SU2_history,SU2_surface,SU2_CFD_info]=self.solveAirfoilSU2...
                (control_point,self.airfoil_coord,...
                self.mesh_filestr,self.DEF_config_filestr,...
                self.CFD_config_filestr,self.flight_condition,...
                self.dat_filestr,self.mesh_out_filestr,self.restart_filestr,self.run_description);

            % evaluate airfoil max thickness
            point_up=mesh_point.AIRFOIL_UP;
            point_low=mesh_point.AIRFOIL_LOW;
            t=max(point_up.Y-point_low.Y);

            % objective
            obj=-SU2_history.('CEff')(end);

            % constaints
            g1=self.CL_min-SU2_history.('CL')(end);
            g2=self.t_mim-t;
            con=[g1,g2];
            coneq=[];

        end
    end
end
