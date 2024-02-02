classdef EdgeNURBS < handle
    % Non-Uniform Rational B-Splines Edge
    % define reference to step standard
    %
    properties % Explicit Attributes
        name='';
        Degree=[]; % degree
        Poles=[]; % control_points_list
        edge_form='.UNSPECIFIED.'; % edge_form
        Periodic='.F.'; % closed_curve (boolean)
        self_intersect='.F.'; % self_intersect (boolean)
        Weights=[]; % weights_data

        Mults=[]; % knot_multiplicities
        Knots=[]; % Knots
        knot_spec='.UNSPECIFIED.';

        u_list=[]; % knot_vector
    end

    properties(Access=private)
        Ctrls=[]; % control points with weight
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
                    Knots=[];
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
                [pole_num,~]=size(Poles);

                % default value
                if isempty(Degree),Degree=pole_num-1;end
                if isempty(Mults),Mults=[Degree+1,ones(1,pole_num-Degree-1),Degree+1];end
                if isempty(Knots),Knots=linspace(0,1,pole_num-Degree+1);end
                if isempty(Weights), Weights=ones(1,pole_num);end;Weights=Weights(:);
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
                self.u_list=u_list;
            end
        end

        function Pnts=calNURBS(self,U)
            % calculate point on Non-Uniform Rational B-Splines Curve
            %
            Pnts=GeomBSpline.bspeval(self.Ctrls,self.Degree,self.u_list,U);
            Pnts=Pnts(:,1:end-1)./Pnts(:,end);
            % U=U(:);point_num=length(U);
            % 
            % Point=zeros(point_num,size(self.Poles,2));
            % 
            % N_list=zeros(1,self.Degree+1);
            % for point_idx=1:point_num
            %     % local u_x in u_list index
            %     u=U(point_idx);
            %     idx_end=size(self.Poles,1); % is equal to the section index
            %     while idx_end > self.Degree+1 && u < self.u_list(idx_end)
            %         idx_end=idx_end-1;
            %     end
            %     idx_start=idx_end-self.Degree;
            %     for N_idx=1:self.Degree+1
            %         N_list(N_idx)=baseFcnN(self.u_list,u,N_idx+idx_start-1,self.Degree);
            %     end
            %     NW_list=N_list.*self.Weights(idx_start:idx_end)';
            % 
            %     Point(point_idx,:)=NW_list*self.Poles(idx_start:idx_end,:)/sum(NW_list);
            % end
        end
    end

    % control Edge
    methods
        function self=reverse(self)
            % revese direction of curve
            %
            self.Poles=flipud(self.Poles);
            self.Mults=fliplr(self.Mults);
            self.Knots=1-fliplr(self.Knots);
            self.Weights=flipud(self.Weights);
            self.u_list=1-fliplr(self.u_list);
        end

        function self=addDegree(self,deg_tar)
            % increase edge degree
            %
            if deg_tar <= self.Degree
                return;
            end

            Ctrl=[self.Poles,self.Weights];
            [Ctrl_new,Mults_new]=GeomApp.addDegree(self.Degree,Ctrl,self.Mults,self.Knots,deg_tar);
            % Ctrl_new=Ctrl_new./Ctrl_new(:,end);

            self.Poles=Ctrl_new(:,1:end-1);
            self.Degree=deg_tar;
            self.Mults=Mults_new;
            self.Weights=Ctrl_new(:,end);
            self.u_list=baseKnotVec(self.Mults,self.Knots);
        end

        function self=insertKnot(self,U_ins)
            % insert knot to edge
            %
            Ctrl=[self.Poles.*self.Weights,self.Weights];
            U=self.u_list;
            deg=self.Degree;

            [Ctrl,U]=GeomApp.insertKnot(deg,Ctrl,U,U_ins);

            self.Poles=Ctrl(:,1:end-1)./Ctrl(:,end);
            self.Weights=Ctrl(:,end);
            self.Knots=unique(U);
            self.Mults=ones(size(self.Knots));
            for k_idx=1:length(self.Knots)
                self.Mults(k_idx)=sum(U == self.Knots(k_idx));
            end
            self.u_list=U;
        end
    
        function self=translate(self,tran_vctr)
            self.Poles=self.Poles+tran_vctr;
        end

        function self=rotate(self,rotate_matrix)
            self.Poles=self.Poles*rotate_matrix';
        end

    end

    % calculate point
    methods
        function [Point,U]=calGeom(self,u_param)
            % generate curve matrix by u_x_list or point_number
            %
            if nargin < 2
                u_param=[];
            end

            low_bou=0;up_bou=1;
            if isempty(u_param)
                pnt=self.calPoint([0;1]);
                bou_min=min(pnt,[],1);bou_max=max(pnt,[],1);
                u_param=2^-12*mean(bou_max-bou_min);
            end
            if length(u_param) == 1 && u_param ~= fix(u_param)
                value_torl=u_param;min_level=2;max_level=50;
                dim=numel(self.calPoint(0));
                [U,Point,~]=GeomApp.meshAdapt1D(@(x) self.calPoint(x),low_bou,up_bou,value_torl,min_level,max_level,dim);
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
            dPoint_dU=zeros(size(U,1),dim); % allocate memory
            dPoint_dU(Bool_C,:)=(Point_UF(Bool_C,:)-Point_UB(Bool_C,:))/2/step;
            dPoint_dU(Bool_U,:)=(Point(Bool_U,:)-Point_UB(Bool_U,:))/step;
            dPoint_dU(Bool_D,:)=(Point_UF(Bool_D,:)-Point(Bool_D,:))/step;
            dPoint_dU=real(dPoint_dU);
        end
    end

    % visualizate curve
    methods
        function li_hdl=plotGeom(self,axe_hdl,u_param,crv_option,ctrl_option)
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

            if isempty(axe_hdl),axe_hdl=gca();end

            % default draw option
            if isempty(crv_option)
                crv_option=struct();
            end
            if isempty(ctrl_option)
                ctrl_option=struct('Marker','s','LineStyle','--','Color','r');
            end

            % calculate point on curve
            Point=self.calGeom(u_param);

            % plot curve
            dim=numel(self.calPoint(0));
            if dim == 2
                li_hdl=line(axe_hdl,Point(:,1),Point(:,2),crv_option);
%                 if ~isempty(self.Poles)
%                     line(axe_hdl,self.Poles(:,1),self.Poles(:,2),ctrl_option);
%                 end
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

    end

    methods
        function set.Poles(self,Poles)
            self.Poles=Poles;
            if isempty(self.Weights) || size(self.Poles,1) ~= size(self.Weights,1)
                self.Ctrls=[self.Poles,ones(size(self.Poles,1),1)];
            else
                self.Ctrls=[self.Poles.*self.Weights,self.Weights];
            end
        end

        function set.Weights(self,Weights)
            self.Weights=Weights;
            self.Ctrls=[self.Poles.*self.Weights,self.Weights];
        end
    end
end
