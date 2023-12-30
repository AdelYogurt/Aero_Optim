classdef EdgeNURBS < handle
    % Non-Uniform Rational B-Splines Edge
    % define reference to step standard
    %
    properties % Explicit Attributes
        name='';
        Degree; % degree
        Poles; % control_points_list
        curve_form='.UNSPECIFIED.'; % curve_form
        Periodic='.F.'; % closed_curve (boolean)
        self_intersect='.F.'; % self_intersect (boolean)
        Weights; % weights_data

        Mults; % knot_multiplicities
        Knots; % Knots
        knot_spec='.UNSPECIFIED.';
    end

    properties % Derived Attributes
        pole_num;
        dimension;
        u_list;
    end

    % define Edge
    methods
        function self=EdgeNURBS(name,Poles,Degree,Mults,Knots,Weights)
            % generate Non-Uniform Rational B-Splines Edge
            %
            % input:
            % name (str):
            % Poles (matrix): control point, pole_num x dimension matrix
            % Degree (matrix): optional input
            % Mults (matrix): optional input
            % Knots (matrix): optional input
            % Weights (matrix): optional input
            %
            % output:
            % EdgeNURBS
            %
            % notice:
            % if input Degree is empty,...
            % Degree default is pole_num-1,...
            % which is Bezier curve.
            %
            if nargin < 6
                Weights=[];
                if nargin < 5
                    Knots = [];
                    if nargin < 4
                        Mults=[];
                        if nargin < 3
                            Degree=[];
                        end
                    end
                end
            end
            self.name=name;

            if nargin > 1
                [pole_num,dimension]=size(Poles);

                % default value
                if isempty(Degree),Degree=pole_num-1;end
                if isempty(Mults),Mults=[Degree+1,ones(1,pole_num-Degree-1),Degree+1];end
                if isempty(Knots),Knots=linspace(0,1,pole_num-Degree+1);end
                if isempty(Weights), Weights=ones(1,pole_num);end;Weights=Weights(:)';
                u_list=baseKnotVec(Mults,Knots);

                if pole_num < (Degree+1)
                    error('EdgeNURBS: pole_num less than Degree+1');
                end

                if length(u_list) ~= pole_num+Degree+1
                    error('EdgeNURBS: knot_num is not equal to pole_num+Degree+1');
                end

                % Explicit Attributes
                self.Degree=Degree;
                self.Poles=Poles;
                self.Mults=Mults;
                self.Knots=Knots;
                self.Weights=Weights;

                % Derived Attributes
                self.pole_num=pole_num;
                self.dimension=dimension;
                self.u_list=u_list;
            end
        end

        function Point=calNURBS(self,U)
            % calculate point on Non-Uniform Rational B-Splines Curve
            %
            U=U(:);point_num=length(U);

            Point=zeros(point_num,self.dimension);

            N_list=zeros(1,self.Degree+1);
            for point_idx=1:point_num
                % local u_x in u_list index
                u=U(point_idx);
                idx_end=self.pole_num; % is equal to the section index
                while idx_end > self.Degree+1 && u < self.u_list(idx_end)
                    idx_end=idx_end-1;
                end
                idx_start=idx_end-self.Degree;
                for N_idx=1:self.Degree+1
                    N_list(N_idx)=baseFcnN(self.u_list,u,N_idx+idx_start-1,self.Degree);
                end
                NW_list=N_list.*self.Weights(idx_start:idx_end);

                Point(point_idx,:)=NW_list*self.Poles(idx_start:idx_end,:)/sum(NW_list);
            end
        end
    end

    % control curve
    methods
        function reverse(self)
            % revese direction of curve
            %
            self.Poles=flipud(self.Poles);
            self.Mults=fliplr(self.Mults);
            self.Knots=1-fliplr(self.Knots);
            self.Weights=fliplr(self.Weights);
            self.u_list=1-fliplr(self.u_list);
        end

        function addDegree(self,degree_target)
            % increase curve degree
            %
            if degree_target > self.Degree
                Poles_old=self.Poles;
                Degree_old=self.Degree;
                Mults_old=self.Mults;
                Weights_old=self.Weights;
                pole_num_old=self.pole_num;
                u_list_old=self.u_list;
                
                for degree_temp=self.Degree+1:degree_target
                    % add repeat node
                    Degree_new=Degree_old+1;
                    pole_num_new=pole_num_old+1;
                    Mults_new=Mults_old+1;
                    u_list_new=baseKnotVec(Mults_new,self.Knots);

                    % calculate new ctrl
                    matrix=zeros(pole_num_new,pole_num_old);
                    for j=1:pole_num_new
                        for i=1:pole_num_old
                            matrix(j,i)=baseFcnL(u_list_old,u_list_new,i,j,Degree_old);
                        end
                    end
                    matrix=1/Degree_new*matrix;
                    Weights_new=(matrix*Weights_old')'; % different to BSpline
                    Poles_new=(matrix.*Weights_old*Poles_old)./Weights_new'; % different to BSpline

                    % sort old curve data
                    Poles_old=Poles_new;
                    Degree_old=Degree_new;
                    Mults_old=Mults_new;
                    Weights_old=Weights_new;
                    pole_num_old=pole_num_new;
                    u_list_old=u_list_new;
                end

                self.Poles=Poles_new;
                self.Degree=Degree_new;
                self.Mults=Mults_new;
                self.Weights=Weights_new;
                self.pole_num=pole_num_new;
                self.u_list=u_list_new;
            end
        end
    end

    % calculate point
    methods
        function [Point,U]=calEdge(self,u_param)
            % generate curve matrix by u_x_list or point_number
            %
            if nargin < 2 || isempty(u_param)
                u_param=1e-3;
            end

            low_bou=0;up_bou=1;
            if length(u_param) == 1 && u_param ~= fix(u_param)
                value_torl=u_param;min_level=2;max_level=50;
                [U,Point,~]=GeomApp.meshAdapt1D(@(x) self.calPoint(x),low_bou,up_bou,value_torl,min_level,max_level,self.dimension);
            else
                if length(u_param) == 1
                    U=[];
                    u_grid_num=u_param;
                else
                    % mean input U matrix
                    U=u_param;
                    u_grid_num=length(U)-1;
                end
                if isempty(U),U=linspace(low_bou,up_bou,u_grid_num+1);end
                Point=self.calPoint(U);
            end
        end

        function Point=calPoint(self,U)
            % according U to calculate point
            %
            Point=calNURBS(self,U);
        end
    end

    % calculate coord
    methods
        function [U,Point]=calCoord(self,Point)

        end

        function [Point,U]=calPorject(self,Point0,U,geo_torl)
            % adjust u by Jacobian transformation
            % also can project point to surface
            %
            if nargin < 6,geo_torl=sqrt(eps);end
            
        end

        function dPoint_dU=calGradient(self,U,step)
            % use differ to calculate gradient
            %
            if nargin < 3 || isempty(step),step=100*eps;end

            U=U(:);
            Point=self.calPoint(U);
            
            Point_UF=self.calPoint(U+step);
            Point_UB=self.calPoint(U-step);
            Bool_U=(U+step) > 1;
            Bool_D=(U-step) < 0;
            Bool(:,1)=Bool_U;Bool(:,2)=Bool_D;
            Bool_C=~any(Bool,2);
            dPoint_dU=zeros(size(U,1),self.dimension); % allocate memory
            dPoint_dU(Bool_C,:)=(Point_UF(Bool_C,:)-Point_UB(Bool_C,:))/2/step;
            dPoint_dU(Bool_U,:)=(Point(Bool_U,:)-Point_UB(Bool_U,:))/step;
            dPoint_dU(Bool_D,:)=(Point_UF(Bool_D,:)-Point(Bool_D,:))/step;
            dPoint_dU=real(dPoint_dU);
        end
    end

    % visualizate curve
    methods
        function drawEdge(self,axe_hdl,u_param,crv_option,ctrl_option)
            % draw curve on figure handle
            %
            if nargin < 5
                ctrl_option=[];
                if nargin < 4
                    crv_option=[];
                    if nargin < 3
                        u_param=[];
                        if nargin < 2
                            axe_hdl=[];
                        end
                    end
                end
            end

            if isempty(axe_hdl),axe_hdl=axes(figure());end

            % default draw option
            if isempty(crv_option)
                crv_option=struct();
            end
            if isempty(ctrl_option)
                ctrl_option=struct('Marker','s','LineStyle','--','Color','r');
            end

            % calculate point on curve
            Point=self.calEdge(u_param);

            % plot curve
            if self.dimension == 2
                line(axe_hdl,Point(:,1),Point(:,2),crv_option);
                if ~isempty(self.Poles)
                    line(axe_hdl,self.Poles(:,1),self.Poles(:,2),ctrl_option);
                end
            else
                line(axe_hdl,Point(:,1),Point(:,2),Point(:,3),crv_option);
                if ~isempty(self.Poles)
                    line(axe_hdl,self.Poles(:,1),self.Poles(:,2),self.Poles(:,3),ctrl_option);
                end
                zlabel('z');
                view(3);
            end
            xlabel('x');
            ylabel('y');

            % axis equal
            % x_range=xlim();
            % y_range=ylim();
            % z_range=zlim();
            % center=[mean(x_range),mean(y_range),mean(z_range)];
            % range=max([x_range(2)-x_range(1),y_range(2)-y_range(1),z_range(2)-z_range(1)])/2;
            % xlim([center(1)-range,center(1)+range]);
            % ylim([center(2)-range,center(2)+range]);
            % zlim([center(3)-range,center(3)+range]);
        end

        function [step_str,obj_idx]=getStep(self,obj_idx)
            % write BSpline into step file
            %
            if nargin < 2,obj_idx=1;end
            out_name=self.name;
            if isempty(out_name),out_name='NONE';end

            % generate CARTESIAN_POINT
            CARTESIAN_POINT=obj_idx;

            str_ctrl=[];Ctrl=[self.Poles,zeros(self.pole_num,3-self.dimension)];
            for ctrl_idx=1:self.pole_num
                str=[num2str(obj_idx,'#%d ='),' CARTESIAN_POINT',...
                    ' ( ',...
                    '''NONE''',', ',...
                    '( ',num2str(Ctrl(ctrl_idx,1),'%.16f'),', ',...
                    num2str(Ctrl(ctrl_idx,2),'%.16f'),', ',...
                    num2str(Ctrl(ctrl_idx,2),'%.16f'),' )',...
                    ' ) ;\n'];
                str_ctrl=[str_ctrl,str];
                obj_idx=obj_idx+1;
            end

            % generate curve
            point_idx=((1:self.pole_num)-1)+CARTESIAN_POINT;
            str_curve=[num2str(obj_idx,'#%d ='),' B_SPLINE_CURVE_WITH_KNOTS',...
                ' ( ',...
                '''',out_name,'''',', ',...
                num2str(self.Degree,'%d'),', ',...
                '( ',num2str(point_idx(1),'#%d'),num2str(point_idx(2:end),', #%d'),' ),\n',...
                self.curve_form,', ',self.Periodic,', ',self.self_intersect,',\n',...
                '( ',num2str(self.Mults(1),'%d'),num2str(self.Mults(2:end),', %d'),' ),\n',...
                '( ',num2str(self.Knots(1),'%.16f'),num2str(self.Knots(2:end),', %.16f'),' ),\n',...
                self.knot_spec,...
                ' ) ;\n'];

            step_str=[str_ctrl,'\n',str_curve,'\n'];

        end

    end
end
