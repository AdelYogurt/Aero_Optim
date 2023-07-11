classdef NACA0012AirfoilProblem < AirfoilProblem
    properties
        AOA;
        SIDESLIP_ANGLE;
        Ma;
        T;
        P;
        T_w;

        % problem basic problem define
        low_bou=[0.5040,2.1551,-0.4582,1.1015,0.1962,...
            0.5040,2.1551,-0.4582,1.1015,0.1962]*1e-2;
        up_bou=[10.5040,12.1551,9.5418,11.1015,10.1962,...
            10.5040,12.1551,9.5418,11.1015,10.1962]*1e-2;
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
        function self=NACA0012AirfoilProblem...
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

            [total_point_list,SU2_out,SU2_history,SU2_surface]=self.solveAirfoilSU2...
                (control_point_up,control_point_low,self.coord_data_dir,...
                self.AOA,self.SIDESLIP_ANGLE,self.Ma,self.T,self.P,self.T_w,...
                self.description,self.initial_data_dir,self.dir_out,self.restart_filestr);

            % evaluate geometry
            up=total_point_list{1};
            low=total_point_list{2};
            t=max(up(:,3)-low(:,3));

            Cl=getValue(SU2_out,'CL');
            CD=getValue(SU2_out,'CD');
            CEff=Cl/CD;

            obj=-CEff;
            con=[self.Cl_min-Cl,self.t_mim-t];
            coneq=[];
        end
    end
end

function [value,idx]=getValue(SU2_out,parameter)
% get value from SU2_out
%
value=[];
for idx=1:length(SU2_out.parameter)
    if strcmp(SU2_out.parameter{idx},parameter)
        value=SU2_out.value{idx};
        return;
    end
end
if isempty(value)
    error('getValue: parameter do not exist');
end

end