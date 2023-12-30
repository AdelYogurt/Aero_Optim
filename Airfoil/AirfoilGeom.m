classdef AirfoilGeom < handle
    % parameterized airfoil entity
    %
    properties
        edge_list; % cell
    end

    % define function
    methods
        function self=AirfoilGeom(C_par_low,Poles_low,C_par_up,Poles_up)
            % generate airfoil object
            % default use control point to generate Bezier curve
            % while Degree equal to control point number -1,...
            % BSpline will be the same as Bezier curve
            %
            airfoil_low=EdgeCST2D('AIRFOIL_LOW',C_par_low,false,1,-1);
            airfoil_up=EdgeCST2D('AIRFOIL_UP',C_par_up,false,1,1);

            airfoil_low.addNURBS(Poles_low);
            airfoil_up.addNURBS(Poles_up);

            self.edge_list={airfoil_low,airfoil_up};
        end

        function [curve,curve_idx]=getEdge(self,curve_name)
            % load curve from curve_list base on input curve name
            %
            for curve_idx=1:length(self.edge_list)
                curve=self.edge_list{curve_idx};
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
                curve=self.getEdge(curve_name);
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
                curve=self.getEdge(curve_name);
                point_idx=airfoil_coord.(curve_name).index;
                U=airfoil_coord.(curve_name).U;

                % calculate coordinate
                Point=curve.calPoint(U);
                mesh_point.(curve_name).index=point_idx;
                mesh_point.(curve_name).X=Point(:,1);
                mesh_point.(curve_name).Y=Point(:,2);
                mesh_point.(curve_name).type='scatter';
            end 

        end
    end

    % visualization function
    methods
        function curve_total=calWire(self,U_par)
            % obtain curve line
            %
            if nargin < 1
                U_par=[];
            end

            curve_total=cell(length(self.edge_list),1);
            for curve_idx=1:length(self.edge_list)
                [curve.X,curve.Y]=self.edge_list{curve_idx}.calEdge(U_par);
                curve.name=self.edge_list{curve_idx}.name;
                curve_total{curve_idx}=curve;
            end

        end

        function drawEdge(self,axe_hdl,U_param,draw_option)
            % draw curve on figure handle
            %
            if nargin < 4
                draw_option=struct();
                if nargin < 3
                    U_param=[];
                    if nargin < 2
                        axe_hdl=[];
                    end
                end
            end
            curve_total=calWire(self,U_param);

            if isempty(axe_hdl),axe_hdl=axes(figure());end
            
            for curve_idx=1:length(curve_total)
                curve=curve_total{curve_idx};
                line(axe_hdl,curve.X,curve.Y,draw_option);
            end

        end
    end
end
