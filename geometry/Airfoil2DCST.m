classdef Airfoil2DCST < handle
    % parameterized airfoil entity
    %
    properties
        class_par_Y=[0.5,1];

        curve_list; % cell
    end

    % define function
    methods
        function self=Airfoil2DCST(control_point_up,control_point_low,curve_type)
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
            end

        end
    end

    % decode x function
    methods(Static)
        function [control_point_up,control_point_low]=decode(control_point_up,control_point_low,x)
            control_number=length(x)/2;

            control_point_up(2:end-1,2)=x(1:control_number)';
            control_point_low(2:end-1,2)=x(control_number+1:end)';
        end
    end

end

function s=shapeFunctionCurve(curve,u)
s=curve.calPoint(u);
s=s(:,2);
end
