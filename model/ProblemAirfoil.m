classdef ProblemAirfoil < ModelAirfoil
    % optimal problem define
    properties
        % basical parameter
        vari_num=10;
        con_num=2;
        coneq_num=0;
        low_bou=0.01*ones(1,10);
        up_bou=0.1*ones(1,10);

        % problem constraint define
        t_mim=0.12;
        CL_min=0.28;
    end

    % main function
    methods
        % problem setup function
        function self=ProblemAirfoil(partitions,geo_param,mesh_param,CFD_param,...
                dir_temp,REMOVE_TEMP,run_description,out_logger)
            % setup airfoil problem function
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
            self@ModelAirfoil(partitions,geo_param,mesh_param,CFD_param,...
                dir_temp,REMOVE_TEMP,run_description,out_logger);

            if strcmp(geo_param.solver,'CST')
                self.low_bou=0.01*ones(1,10);
                self.up_bou=0.1*ones(1,10);
            elseif strcmp(geo_param.solver,'FFD')
                self.low_bou=-0.1*ones(1,10);
                self.up_bou=0.1*ones(1,10);
            elseif strcmp(geo_param.solver,'HICKS_HENNE')
                self.low_bou=-0.01*ones(1,10);
                self.up_bou=0.01*ones(1,10);
            end
        end
    end

    % optimization interface funciton
    methods
        function [obj,con,coneq]=objcon_fcn(self,x)
            % objcon function for optimize
            %
            % input:
            % x
            %
            % output:
            % obj, con, coneq
            %
            geo_in=self.decode(x);

            [geo_out,mesh_out,CFD_out]=self.solveAirfoil(geo_in);
            mesh_point=geo_out.mesh_point;
            if strcmp(self.CFD_param.solver,'SU2_CFD')
                CEff=CFD_out.SU2_history.('CEff')(end);
                CL=CFD_out.SU2_history.('CL')(end);
            elseif strcmp(self.CFD_param.solver,'Fluent')
                CD=CFD_out.fluent_out.data(end,2);
                CL=CFD_out.fluent_out.data(end,3);
                CEff=CL/CD;
            end

            % evaluate airfoil max thickness
            point_up=mesh_point.AIRFOIL_UP;
            point_low=mesh_point.AIRFOIL_LOW;
            t=max(point_up.Y-point_low.Y);

            % objective
            obj=-CEff;

            % constaints
            g1=self.CL_min-CL;
            g2=self.t_mim-t;
            con=[g1,g2];
            coneq=[];
        end
    end

    % parameterization function
    methods
        function drawX(self,x,fig_hdl,draw_option)
            % draw x geometry
            %
            if nargin < 4
                draw_option=struct();
                if nargin < 3
                    fig_hdl=[];
                end
            end

            geo_in=self.decode(x);
            self.drawGeo(geo_in);
        end

        function geo_in=decode(self,x)
            % decode x into geo_in
            %
            if strcmp(self.geo_param.solver,'CST')
                x_number=int64(length(x));
                geo_in.low=[linspace(0,1,x_number/2);x(1:x_number/2)]';
                geo_in.up=[linspace(0,1,x_number/2);x(x_number/2+1:end)]';
            elseif strcmp(self.geo_param.solver,'FFD')

            elseif strcmp(self.geo_param.solver,'HICKS_HENNE')
                geo_in.SIZE=x;
            end
        end
    end
end
