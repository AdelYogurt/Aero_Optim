classdef AirfoilProblem < AirfoilModel
    % optimal problem define
    properties
        % basical parameter
        vari_num=10;
        con_num=2;
        coneq_num=0;
        low_bou=0.01*ones(1,10);
        up_bou=0.1*ones(1,10);

        % constraint define
        CL_min=0.28;
        t_min=0.12;
    end

    % main function
    methods
        % problem setup function
        function self=AirfoilProblem(partitions,geom_param,mesh_param,CFD_param,...
                dir_temp,REMOVE_TEMP,run_desc,out_logger)
            % setup airfoil problem function
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
            self@AirfoilModel(partitions,geom_param,mesh_param,CFD_param,...
                dir_temp,REMOVE_TEMP,run_desc,out_logger);

            if strcmp(geom_param.solver,'CST')
                geom_low_bou=-0.03*ones(1,self.vari_num);
                geom_up_bou=0.1*ones(1,self.vari_num);
            elseif strcmp(geom_param.solver,'FFD')
                geom_low_bou=-0.1*ones(1,self.vari_num);
                geom_up_bou=0.1*ones(1,self.vari_num);
            elseif strcmp(geom_param.solver,'HICKS_HENNE')
                geom_low_bou=-0.01*ones(1,self.vari_num);
                geom_up_bou=0.01*ones(1,self.vari_num);
            end

            self.low_bou=geom_low_bou;
            self.up_bou=geom_up_bou;
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

            [geo_out,mesh_out,CFD_out]=self.solveModel(geo_in);
            mesh_point=geo_out.mesh_point;
            if strcmp(self.CFD_param.solver,'SU2_CFD')
                CEff=CFD_out.SU2_data.('CEff')(end);
                CL=CFD_out.SU2_data.('CL')(end);
            elseif strcmp(self.CFD_param.solver,'Fluent')
                CD=CFD_out.fluent_out.data(end,2);
                CL=CFD_out.fluent_out.data(end,3);
                CEff=CL/CD;
            end

            % evaluate airfoil max thickness
            point_up=mesh_point.AIRFOIL_UP;
            point_low=mesh_point.AIRFOIL_LOW;
            t=max(point_up.Y-interpLinear(point_low.X,point_low.Y,point_up.X));

            % objective
            obj=-CEff;

            % constaints
            g1=self.CL_min-CL;
            g2=self.t_min-t;
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
            self.drawGeo(geo_in,fig_hdl,draw_option);
        end

        function geo_in=decode(self,x)
            % decode x into geo_in
            %
            if strcmp(self.geom_param.solver,'CST')
                x_num=length(x);
                geo_in.Poles_low=[linspace(0,1,x_num/2);-x(1:x_num/2)]';
                geo_in.Poles_up=[linspace(0,1,x_num/2);x(x_num/2+1:end)]';

                if ~isfield(self.geom_param,'C_par_low'), geo_in.C_par_low=[0.5,1];
                else,geo_in.C_par_low=self.geom_param.C_par_low;end
                if ~isfield(self.geom_param,'C_par_up'), geo_in.C_par_up=[0.5,1];
                else,geo_in.C_par_up=self.geom_param.C_par_up;end

            elseif strcmp(self.geom_param.solver,'FFD')

            elseif strcmp(self.geom_param.solver,'HICKS_HENNE')
                geo_in.SIZE=x;
            end
            
        end
    end
end

function [Y_pred,X_pred_idx]=interpLinear(X,Y,X_pred)
% a simple linear interpolation to calculate data
% Y can be m x n matrix, n is variable
%
[X,idx]=sort(X);
Y=Y(idx,:);
Y_pred=zeros(length(X_pred),size(Y,2));
X_pred_idx=zeros(length(X_pred),1);
for x_idx=1:length(X_pred)
    x_pred=X_pred(x_idx);
    num=length(X);
    idx=num; % search start from last one, find out X samll than x
    while ((idx > 1) && (X(idx) > x_pred))
        idx=idx-1;
    end

    if (idx == num)
        % out interp
        Y_pred(x_idx,:)=(Y(num,:)-Y(num-1,:))/(X(num)-X(num-1))*(x_pred-X(num))+Y(num,:);
        X_pred_idx(x_idx)=idx;
    elseif (idx == 0)
        Y_pred(x_idx,:)=(Y(2,:)-Y(1,:))/(X(2)-X(1))*(x_pred-X(1))+Y(1,:);
        X_pred_idx(x_idx)=idx;
    else
        % linear interpolation
        Y_pred(x_idx,:)=Y(idx,:)+...
            (Y(idx+1,:)-Y(idx,:))*...
            (x_pred-X(idx,:))/...
            (X(idx+1)-X(idx));
        X_pred_idx(x_idx)=idx;
    end
end

end

