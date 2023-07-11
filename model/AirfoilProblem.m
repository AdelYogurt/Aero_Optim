classdef AirfoilProblem < handle
    properties
        partitions;
        out_logger;
        class_par_Y=[0.5,1];
    end
    methods
        function self=AirfoilProblem(SU2_partitions,out_logger)
            if nargin < 2
                out_logger=[];
            end
            self.partitions=SU2_partitions;
            self.out_logger=out_logger;
        end

        function total_point_list=prePoint(self,control_point_up,control_point_low,coord_data_dir,dat_filestr)
            if nargin < 5
                dat_filestr=[];
            end

            bezier_up=CurveBezier(control_point_up);
            bezier_low=CurveBezier(control_point_low);
            shape_fun_up=@(x) shapeFunctionCurve(bezier_up,x);
            shape_fun_low=@(x) shapeFunctionCurve(bezier_low,x);
            airfoil_up=CST2DCurve(1,1,shape_fun_up,self.class_par_Y);
            airfoil_low=CST2DCurve(1,-1,shape_fun_low,self.class_par_Y);

            airfoil_coord_up=importdata(fullfile(coord_data_dir,'airfoil_up_local_coord.txt'));
            [X_up,Y_up]=airfoil_up.calPoint(airfoil_coord_up(:,2)');

            airfoil_coord_low=importdata(fullfile(coord_data_dir,'airfoil_low_local_coord.txt'));
            [X_low,Y_low]=airfoil_low.calPoint(airfoil_coord_low(:,2)');

            total_point_list={[airfoil_coord_up(:,1),X_up',Y_up'],[airfoil_coord_low(:,1),X_low',Y_low']};

            if ~isempty(dat_filestr)
                file_out=fopen(dat_filestr,'w');
                for surface_index=1:length(total_point_list)
                    point_list=total_point_list{surface_index};

                    index=point_list(:,1);
                    X=point_list(:,2);
                    Y=point_list(:,3);

                    for point_index=1:size(index,1)
                        fprintf(file_out,'%d %f %f\n',index(point_index)-1,X(point_index),Y(point_index));
                    end
                end
                fclose(file_out);
            end
        end

        function [total_point_list,SU2_history,SU2_surface]=solveAirfoilSU2...
                (self,control_point_up,control_point_low,coord_data_dir,...
                AOA,SIDESLIP_ANGLE,Ma,T,P,T_w,...
                description,initial_data_dir,dir_out,restart_filestr)
            % call SU2 function to obtain airfoil aerodynamic
            %
            if nargin < 14
                restart_filestr=[];
                if nargin < 13
                    dir_out=[];
                    if nargin < 12
                        initial_data_dir=[];
                        if nargin < 11
                            description=[];
                            if nargin < 10
                                T_w=[];
                                if nargin < 9
                                    P=[];
                                    if nargin < 8
                                        T=[];
                                        if nargin < 7
                                            Ma=[];
                                            if nargin < 6
                                                SIDESLIP_ANGLE=[];
                                                if nargin < 5
                                                    AOA=[];
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end

            if isempty(initial_data_dir),initial_data_dir='initial_data';end
            if isempty(dir_out),dir_out='model';end
            if isempty(T_w),T_w=273.0;end
            if isempty(P),P=101325;end
            if isempty(T),T=273.0;end
            if isempty(Ma),Ma=0.63;end
            if isempty(SIDESLIP_ANGLE),SIDESLIP_ANGLE=0.0;end
            if isempty(AOA),AOA=2.0;end

            dat_filestr=fullfile(dir_out,'airfoil_deform.dat');
            mesh_out_filestr=fullfile(dir_out,'airfoil_deform.su2');

            % rebulid surface point
            mesh_filestr=fullfile(initial_data_dir,'NACA0012.cgns');
            cfg_filestr=fullfile(initial_data_dir,'NACA0012_deform.cfg');
            total_point_list=self.prePoint(control_point_up,control_point_low,coord_data_dir,dat_filestr);
            runSU2DEF(mesh_filestr,cfg_filestr,dat_filestr,mesh_out_filestr,self.partitions,[],...
                [],description,[],self.out_logger);

            % run SU2_CFD to get result
            cfg_filestr=fullfile(initial_data_dir,'airfoil.cfg');
            mesh_filestr=mesh_out_filestr;

            CFD_config.AOA=AOA;
            CFD_config.SIDESLIP_ANGLE=SIDESLIP_ANGLE;
            CFD_config.MACH_NUMBER=Ma;
            CFD_config.FREESTREAM_TEMPERATURE=T;
            CFD_config.FREESTREAM_PRESSURE=P;

            [SU2_history,SU2_surface]=runSU2CFD(mesh_filestr,cfg_filestr,self.partitions,CFD_config,...
                restart_filestr,[],description,[],self.out_logger);

        end
    end
end

function s=shapeFunctionCurve(curve,u)
s=curve.interpPoint(u);
s=s(:,2)';
end
