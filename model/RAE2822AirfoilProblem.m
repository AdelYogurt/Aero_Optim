classdef RAE2822AirfoilProblem < AirfoilProblem
    properties
        AOA;
        SIDESLIP_ANGLE;
        Ma;
        T;
        P;
        T_w;

        % problem basic problem define
        low_bou=[0.2467,1.0626,3.0626,-1.8196,-2.0042,...
            0.0370,1.7171,0.8304,3.7183,1.6115]*1e-2;
        up_bou=[10.2467,11.0626,13.0626, 8.1804, 7.9958,...
            10.0370,11.7171,10.8304,13.7183,11.6115]*1e-2;
        vari_num=10;
        con_num=2;
        coneq_num=0;

        % problem constraint define
        t_mim=0.096;
        Cl_min=0.268863;

        description;
        initial_data_dir;
        coord_data_dir;
        dir_out;
        restart_filestr;
    end
    methods
        function self=RAE2822AirfoilProblem...
                (SU2_partitions,AOA,SIDESLIP_ANGLE,Ma,T,P,T_w,...
                description,initial_data_dir,dir_out,restart_filestr,my_logger)
            if nargin < 12
                my_logger=[];
                if nargin < 11
                    restart_filestr=[];
                    if nargin < 10
                        dir_out=[];
                        if nargin < 9
                            initial_data_dir=[];
                            if nargin < 8
                                description=[];
                                if nargin < 7
                                    T_w=[];
                                    if nargin < 6
                                        P=[];
                                        if nargin < 5
                                            T=[];
                                            if nargin < 4
                                                Ma=[];
                                                if nargin < 3
                                                    SIDESLIP_ANGLE=[];
                                                    if nargin < 2
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
            end

            if isempty(initial_data_dir),initial_data_dir='initial_data';end
            if isempty(dir_out),dir_out='model';end
            if isempty(T_w),T_w=273.0;end
            if isempty(P),P=101325;end
            if isempty(T),T=273.0;end
            if isempty(Ma),Ma=0.63;end
            if isempty(SIDESLIP_ANGLE),SIDESLIP_ANGLE=0.0;end
            if isempty(AOA),AOA=2.0;end

            self@AirfoilProblem(SU2_partitions,my_logger);

            self.partitions=SU2_partitions;
            self.AOA=AOA;
            self.SIDESLIP_ANGLE=SIDESLIP_ANGLE;
            self.Ma=Ma;
            self.T=T;
            self.P=P;
            self.T_w=T_w;

            self.description=description;
            self.initial_data_dir=initial_data_dir;
            self.coord_data_dir='coord_data';
            self.dir_out=dir_out;
            self.restart_filestr=restart_filestr;

            % temp file
            safeMakeDirs('model/SU2_temp',my_logger);
            safeMakeDirs(dir_out,my_logger);
        end

        function [obj,con,coneq]=objcon_fcn(self,x)
            x_control_number=length(x)/2;

            % generate Bezier controll point
            control_point_up=zeros(x_control_number+2,2);
            control_point_up(1,:)=[0.000000000000,0.049430665250];
            control_point_up(end,:)=[1.000000000000,0.081613682039];
            control_point_up(:,1)=linspace(0,1,x_control_number+2);
            control_point_low=zeros(x_control_number+2,2);
            control_point_low(1,:)=[0.000000000000,0.049716816054];
            control_point_low(end,:)=[1.000000000000,-0.023031147026];
            control_point_low(:,1)=linspace(0,1,x_control_number+2);

            % add x
            control_point_up(2:end-1,2)=x(1:x_control_number);
            control_point_low(2:end-1,2)=x(x_control_number+1:end);

            [total_point_list,SU2_history,SU2_surface]=self.solveAirfoilSU2...
                (control_point_up,control_point_low,self.coord_data_dir,...
                self.AOA,self.SIDESLIP_ANGLE,self.Ma,self.T,self.P,self.T_w,...
                self.description,self.initial_data_dir,self.dir_out,self.restart_filestr);

            % evaluate geometry
            up=total_point_list{1};
            low=total_point_list{2};
            t=max(up(:,3)-low(:,3));

            Cl=SU2_history.('CL')(end);
            CEff=SU2_history.('CEff')(end);

            obj=-CEff;
            con=[self.Cl_min-Cl,self.t_mim-t];
            coneq=[];
        end
    end
end
