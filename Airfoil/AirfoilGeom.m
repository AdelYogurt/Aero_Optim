classdef AirfoilGeom < handle
    % parameterized airfoil entity
    %
    properties
        class_par_Y=[0.5,1];

        curve_list; % cell
    end

    % define function
    methods
        function self=AirfoilGeom(control_point_low,control_point_up)
            % generate airfoil object
            % default use control point to generate Bezier curve
            % while degree equal to control point number -1,...
            % BSpline will be the same as Bezier curve
            %
            airfoil_low=CurveCST2D('AIRFOIL_LOW',1,-1,self.class_par_Y);
            airfoil_up=CurveCST2D('AIRFOIL_UP',1,1,self.class_par_Y);

            airfoil_low.addShapeBSpline(control_point_low(:,1),control_point_low(:,2));
            airfoil_up.addShapeBSpline(control_point_up(:,1),control_point_up(:,2));

            self.curve_list={airfoil_low,airfoil_up};
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
                U=curve.calCoord(point(:,1));
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

        function drawCurve(self,axe_hdl,U_par,draw_option)
            % draw curve on figure handle
            %
            if nargin < 4
                draw_option=struct();
                if nargin < 3
                    U_par=[];
                    if nargin < 2
                        axe_hdl=[];
                    end
                end
            end
            curve_total=calCurveMatrix(self,U_par);

            if isempty(axe_hdl),axe_hdl=axes(figure());end
            
            for curve_idx=1:length(curve_total)
                curve=curve_total{curve_idx};
                line(axe_hdl,curve.X,curve.Y,draw_option);
            end

        end
    end
end
