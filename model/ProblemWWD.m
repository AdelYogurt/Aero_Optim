classdef ProblemWWD < ModelWWD
    % optimal problem define
    properties
        % basical parameter
        vari_num=16;
        con_num=4;
        coneq_num=0;
        low_bou=[0.7,0.7, ...
            2.4,0.4,0.5,0.5, ...
            0.1,0.1,0.01, ...
            0.4,0.4,0.4,0.6,0.5,0.01,0.01];
        up_bou=[1.0,1.0, ...
            3.2,0.8,2.0,2.0, ...
            0.5,0.5,0.04, ...
            0.6,0.6,0.6,0.8,0.7,0.05,0.05];

        % problem constraint define
        CD_max=0.050;
        CL_min=0.125;
        HF_max=8e6;
        V_eff_min=0.26;
    end

    methods
        % problem setup function
        function self=ProblemWWD(partitions,geo_param,mesh_param,CFD_param,...
                dir_temp,REMOVE_TEMP,run_description,out_logger)
            % setup WWD problem function
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
            self@ModelWWD(partitions,geo_param,mesh_param,CFD_param,...
                dir_temp,REMOVE_TEMP,run_description,out_logger);
        end

        % optimization interface funciton
        function [obj,con,coneq]=objcon_fcn(self,x)
            % interface of objcon function for optimize
            %
            % input:
            % x
            %
            % output:
            % obj, con, coneq
            %

            geo_in=self.decode(x);

            % analysis CFD result
            try
                [geo_out,mesh_out,CFD_out]=solveWWD(self,geo_in);
            catch
                error_massage=['fatal solve WWD SU2 at ',num2str(x(:)','%f '),'\n','error message: '];
                if ~isempty(self.out_logger)
                    self.out_logger.error(error_massage);
                end
                disp(error_massage);
                geo_out=geoModule(self,geo_in);
                CFD_out.SU2_history.('CEff')=2.4;
                CFD_out.SU2_history.('CD')=0.055;
                CFD_out.SU2_history.('CL')=0.130;
                CFD_out.SU2_surface.('Heat_Flux')=9.5e+06;
            end

            % generate value
            value.('CEff')=CFD_out.SU2_history.('CEff')(end);
            value.('CD')=CFD_out.SU2_history.('CD')(end);
            value.('CL')=CFD_out.SU2_history.('CL')(end);
            value.('maxHF')=max(CFD_out.SU2_surface.('Heat_Flux'));
            value.('V_eff')=geo_out.V_eff;

            [obj,con,coneq]=self.objconValue(value);
        end

        function [Obj,Con,Coneq]=objconValue(self,Value)
            % base on input value to calculate obj, con and coneq
            %

            % objective
            Obj=-[Value.('CEff')];
            Obj=Obj(:);

            % constaints
            G1=[Value.('CD')]-self.CD_max;
            G2=self.CL_min-[Value.('CL')];
            G3=([Value.('maxHF')]-self.HF_max)/1e6;
            G4=self.V_eff_min-[Value.('V_eff')];
            Con=[G1(:),G2(:),G3(:),G4(:)];
            Coneq=[];
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

            if ispc()
                geo_in=self.decode(x);
                self.drawGeo(geo_in,fig_hdl,draw_option);
            end
        end

        function geo_in=decode(~,x)
            % decode x into geo_in
            %
            geo_in.total_length=4;
            geo_in.param=x;
        end
    end
end

