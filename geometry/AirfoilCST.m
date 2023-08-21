classdef AirfoilCST < handle
    % parameterized airfoil entity
    %
    properties
        class_par_Y=[0.5,1];

        curve_list; % cell
    end

    % define function
    methods
        function self=AirfoilCST(control_point_low,control_point_up,curve_type)
            % generate airfoil object
            % default use control point to generate Bezier curve
            %
            if nargin < 3
                curve_type='Bezier';
            end

            if strcmp(curve_type,'Bezier')
                bezier_up=CurveBezier('AIRFOIL_UP',control_point_up);
                bezier_low=CurveBezier('AIRFOIL_LOW',control_point_low);
            end

            shape_fun_up=@(x) shapeFunctionCurve(bezier_up,x);
            shape_fun_low=@(x) shapeFunctionCurve(bezier_low,x);

            airfoil_up=CurveCST2D('AIRFOIL_UP',1,1,shape_fun_up,self.class_par_Y);
            airfoil_low=CurveCST2D('AIRFOIL_LOW',1,-1,shape_fun_low,self.class_par_Y);

            self.curve_list={airfoil_up,airfoil_low};
        end

        function [curve,curve_idx]=getCurve(self,curve_name)
            % load curve from curve_list base on input curve name
            %
            for curve_idx=1:length(self.curve_list)
                curve=self.curve_list{curve_idx};
                if strcmp(curve.name,curve_name)
                    return;
                end
            end
            curve=[];
            curve_idx=0;
        end
    end

    % parameterized function
    methods
        function airfoil_coord=calCoord(self,point_list,curve_index_list)
            % calculate all input point local coordinate
            %
            airfoil_coord=struct();

            curve_name_list=fieldnames(curve_index_list);
            for curve_idx=1:length(curve_name_list)
                % get curve
                curve_name=curve_name_list{curve_idx};
                curve=self.getCurve(curve_name);
                point_idx=curve_index_list.(curve_name);

                % calculate coordinate
                point=point_list(point_idx,1:2);
                U=curve.calCoordinate(point(:,1));
                airfoil_coord.(curve_name).index=point_idx;
                airfoil_coord.(curve_name).U=U;
            end

        end

        function mesh_point=calMeshPoint(self,airfoil_coord)
            % calculate all mesh point by input curve coord
            %
            mesh_point=struct();

            curve_name_list=fieldnames(airfoil_coord);
            for curve_idx=1:length(curve_name_list)
                % get curve
                curve_name=curve_name_list{curve_idx};
                curve=self.getCurve(curve_name);
                point_idx=airfoil_coord.(curve_name).index;
                U=airfoil_coord.(curve_name).U;

                % calculate coordinate
                [X,Y]=curve.calPoint(U);
                mesh_point.(curve_name).index=point_idx;
                mesh_point.(curve_name).X=X;
                mesh_point.(curve_name).Y=Y;
                mesh_point.(curve_name).type='scatter';
            end 

        end
    end

    % visualization function
    methods
        function curve_total=calCurveMatrix(self,U_par)
            % obtain curve line
            %
            if nargin < 1
                U_par=[];
            end

            curve_total=cell(length(self.curve_list),1);
            for curve_idx=1:length(self.curve_list)
                [curve.X,curve.Y]=self.curve_list{curve_idx}.calCurve(U_par);
                curve.name=self.curve_list{curve_idx}.name;
                curve_total{curve_idx}=curve;
            end

        end

        function drawCurve(self,U_par,fig_hdl,draw_option)
            % draw curve on figure handle
            %
            if nargin < 4
                draw_option=struct();
                if nargin < 3
                    fig_hdl=[];
                    if nargin < 2
                        U_par=[];
                    end
                end
            end
            curve_total=calCurveMatrix(self,U_par);

            if isempty(fig_hdl)
                fig_hdl=figure(1);
            end
            axes_hdl=fig_hdl.CurrentAxes;
            if isempty(axes_hdl)
                axes_hdl=axes(fig_hdl);
            end
            
            for curve_idx=1:length(curve_total)
                curve=curve_total{curve_idx};
                line(axes_hdl,curve.X,curve.Y,draw_option);
            end

        end
    end
end

function s=shapeFunctionCurve(curve,u)
s=curve.calPoint(u);
s=s(:,2);
end
