classdef Curve < handle & matlab.mixin.Copyable
    % Non-Uniform Rational B-Splines Curve
    % define reference to step standard
    %
    properties % Explicit Attributes
        name='';
        Degree=[]; % degree
        Poles=[]; % control_points_list
        curve_form='.UNSPECIFIED.'; % curve_form
        Periodic='.F.'; % closed_curve (boolean)
        self_intersect='.F. '; % self_intersect (boolean)
        Weights=[]; % weights_data

        Mults=[]; % knot_multiplicities
        Knots=[]; % Knots
        knot_spec='.UNSPECIFIED.';

        u_list=[]; % knot_vector
    end

    properties % Derivate Attributes
        deriv_curve;
    end

    methods % define curve
        function self=Curve(Poles,Degree,Mults,Knots,Weights)
            % generate Non-Uniform Rational B-Splines Curve
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
            % Curve
            %
            % notice:
            % Degree default is pole_num-1 which will be Bezier curve
            %
            if nargin < 5
                Weights=[];
                if nargin < 4
                    Knots=[];
                    if nargin < 3
                        Mults=[];
                        if nargin < 2
                            Degree=[];
                        end
                    end
                end
            end

            [pole_num,~]=size(Poles);

            % default value
            if isempty(Degree),Degree=pole_num-1;end
            if isempty(Mults) && isempty(Knots)
                Mults=[Degree+1,ones(1,pole_num-Degree-1),Degree+1];
                Knots=linspace(0,1,pole_num-Degree+1);
            elseif ~isempty(Mults) && isempty(Knots)
                Knots=linspace(0,1,length(Mults));
            elseif isempty(Mults) && ~isempty(Knots)
                error('Curve: need Mults inpt');
            end
            if ~isempty(Weights), Weights=Weights(:);end
            u_list=baseKnotVec(Mults,Knots);

            if pole_num < (Degree+1)
                error('Curve: pole_num less than Degree+1');
            end

            if length(u_list) ~= pole_num+Degree+1
                error('Curve: knot_num is not equal to pole_num+Degree+1');
            end

            % Explicit Attributes
            self.Degree=Degree;
            self.Poles=Poles;
            self.Mults=Mults;
            self.Knots=Knots;
            self.Weights=Weights;
            self.u_list=u_list;
        end

        function [Points,Weights]=calPoint(self,U_x)
            % calculate point on Non-Uniform Rational B-Splines Curve
            %
            U_x=U_x(:);
            if isempty(self.Weights), Ctrls=self.Poles;
            else, Ctrls=[self.Poles.*self.Weights,self.Weights]; end

            % evaluate along the u direction
            Points=calBSpline(Ctrls,self.Degree,self.u_list,U_x);
            
            if isempty(self.Weights), Weights=[];
            else
                Weights=Points(:,end);
                Points=Points(:,1:end-1);
                if nargout < 2
                    Points=Points./Weights;
                end
            end

            function Points=calBSpline(Ctrls,Degree,u_list,U_x)
                [N_list,idx_srt,idx_end]=baseFcnN(U_x,Degree,u_list);
                Points=zeros(length(U_x),size(Ctrls,2));
                for deg_idx=1:Degree+1
                    Points=Points+N_list(:,deg_idx).*Ctrls(idx_srt+(deg_idx-1),:);
                end
            end

            % U=U(:);point_num=length(U);
            % Point=zeros(point_num,size(self.Poles,2));
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
        
        function [Points,dPoints_dU,dPoints_dU2]=calGradient(self,U_x)
            % calculate gradient of Curve
            %
            deriv_time=nargout-1;
            self.preGradient(deriv_time);

            % calculate initial curve point and weight
            [Points,Weits]=self.calPoint(U_x);
            if ~isempty(Weits), Points=Points./Weits;end

            % calculate first derivative
            [Points_deriv,Weits_deriv]=self.deriv_curve.calPoint(U_x);
            if isempty(Weits)
                dPoints_dU=Points_deriv;
            else
                Weits_deriv=Points_deriv(:,end);
                Points_deriv=Points_deriv(:,1:end-1);
                dPoints_dU=(Points_deriv-Weits_deriv.*Points)./Weits;
            end

            if nargout >= 3
                % calculate second derivative
                [Points_deriv2,~]=self.deriv_curve.deriv_curve.calPoint(U_x);
                Weits_deriv2=Points_deriv2(:,end);
                Points_deriv2=Points_deriv2(:,1:end-1);
                dPoints_dU2=(Points_deriv2-(2*Points_deriv.*Weits_deriv+Points.*Weits_deriv2)./Weits+2*Points.*Weits_deriv.^2./Weits.^2)./Weits;
            end
        end

        function ln_hdl=displayGeom(self,axe_hdl,crv_option,u_param)
            % display curve on axes
            %
            if nargin < 4
                u_param=[];
                if nargin < 3
                    crv_option=[];
                    if nargin < 2
                        axe_hdl=[];
                    end
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end

            if isempty(crv_option), crv_option=struct();end

            % calculate Points on curve
            if isempty(u_param), u_param=101;end
            if length(u_param) == 1, u_param=linspace(min(self.Knots),max(self.Knots),u_param);end
            Points=self.calPoint(u_param);

            % draw Points on axe_hdl
            if size(self.Poles,2) == 2
                ln_hdl=line(axe_hdl,Points(:,1),Points(:,2),crv_option);
            else
                ln_hdl=line(axe_hdl,Points(:,1),Points(:,2),Points(:,3),crv_option);
                zlabel('z');
                view(3);
            end
            xlabel('x');
            ylabel('y');
        end

        function ln_hdl=displayPoles(self,axe_hdl,pole_option)
            % draw curve on figure handle
            %
            if nargin < 3
                pole_option=[];
                if nargin < 2
                    axe_hdl=[];
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end

            % default draw option
            if isempty(pole_option)
                pole_option=struct('Marker','s','LineStyle','--','Color','r');
            end

            % draw Poles on axe_hdl
            if size(self.Poles,2) == 2
                ln_hdl=line(axe_hdl,self.Poles(:,1),self.Poles(:,2),pole_option);
            else
                ln_hdl=line(axe_hdl,self.Poles(:,1),self.Poles(:,2),self.Poles(:,3),pole_option);
                zlabel('z');
                view(3);
            end
            xlabel('x');
            ylabel('y');
        end
    end

    methods % control Edge
        function self=reverse(self)
            % revese direction of curve
            %
            self.Poles=self.Poles(end:-1:1,:);
            self.Mults=self.Mults(end:-1:1);
            min_k=min(self.Knots);max_k=max(self.Knots);dk=max_k-min_k;
            self.Knots=dk-(self.Knots(end:-1:1)-min_k)+min_k;
            self.Weights=self.Weights(end:-1:1,:);
            self.u_list=baseKnotVec(self.Mults,self.Knots);
        end

        function self=addDegree(self,Degree_target)
            % increase curve degree
            %
            if Degree_target <= self.Degree, return;end

            % modify Degree, Mults and Ctrls
            if isempty(self.Weights), Ctrls=[self.Poles];
            else, Ctrls=[self.Poles.*self.Weights,self.Weights];end
            [Ctrls,self.Degree,self.Mults]=GeomApp.addDegree(Ctrls,self.Degree,self.Mults,self.Knots,Degree_target);
            self.u_list=baseKnotVec(self.Mults,self.Knots);

            % modify Poles and Weights
            if isempty(self.Weights), self.Poles=Ctrls;
            else
                self.Poles=Ctrls(:,1:end-1)./Ctrls(:,end);
                self.Weights=Ctrls(:,end);
            end
        end

        function self=insertKnot(self,U_ins)
            % insert knot to curve
            %
            if isempty(U_ins), return;end

            % modify u_list, Mults, Knots and Ctrls
            if isempty(self.Weights), Ctrls=[self.Poles];
            else, Ctrls=[self.Poles.*self.Weights,self.Weights];end
            [Ctrls,~,self.u_list]=GeomApp.insertKnot(Ctrls,self.Degree,self.u_list,U_ins);
            self.Knots=unique(self.u_list);
            self.Mults=ones(size(self.Knots));
            for k_idx=1:length(self.Knots)
                self.Mults(k_idx)=sum(self.u_list == self.Knots(k_idx));
            end

            % modify Poles and Weights
            if isempty(self.Weights), self.Poles=Ctrls;
            else
                self.Poles=Ctrls(:,1:end-1)./Ctrls(:,end);
                self.Weights=Ctrls(:,end);
            end
        end
    
        function [crv_1,crv_2]=splitCurve(self,u_b)
            % split curve at ub
            %
            if u_b <= min(self.Knots) || u_b >= max(self.Knots)
                error('Curve.splitCurve: ub is out of boundary of curve');
            end

            % insert ub to modify poles
            rep_tim=sum(find(self.u_list == u_b));
            self.insertKnot(repmat(u_b,1,self.Degree-rep_tim));

            % locate u place
            u_num=find(self.u_list == u_b,1,'last');
            pole_num=u_num-self.Degree;
            k_num=find(self.Knots == u_b,1,'last');

            % generate new curve
            Poles_1=self.Poles(1:pole_num,:);
            Degree_1=self.Degree;
            Mults_1=self.Mults(1:k_num);Mults_1(end)=Mults_1(end)+1;
            Knots_1=self.Knots(1:k_num);
            Knots_1=(Knots_1-min(Knots_1))/(max(Knots_1)-min(Knots_1));
            if isempty(self.Weights), Weights_1=[];
            else, Weights_1=self.Weights(1:pole_num,:);end
            crv_1=Curve(Poles_1,Degree_1,Mults_1,Knots_1,Weights_1);

            Poles_2=self.Poles(pole_num:end,:);
            Degree_2=self.Degree;
            Mults_2=self.Mults(k_num:end);Mults_2(1)=Mults_2(1)+1;
            Knots_2=self.Knots(k_num:end);
            Knots_2=(Knots_2-min(Knots_2))/(max(Knots_2)-min(Knots_2));
            if isempty(self.Weights), Weights_2=[];
            else, Weights_2=self.Weights(pole_num:end,:);end
            crv_2=Curve(Poles_2,Degree_2,Mults_2,Knots_2,Weights_2);
        end

        function self=preGradient(self,deriv_time)
            % generate derivate BSpline for calculate gradient
            %
            if deriv_time > 0
                if isempty(self.deriv_curve)
                    % generate derivate curve
                    if isempty(self.Weights), Ctrls=[self.Poles];
                    else, Ctrls=[self.Poles.*self.Weights,self.Weights];end
                    [Ctrls_deriv,Deg_deriv,Mults_deriv]=GeomApp.calGradient(Ctrls,self.Degree,self.Mults,self.Knots);

                    % modify Poles and Weights
                    self.deriv_curve=Curve(Ctrls_deriv,Deg_deriv,Mults_deriv,self.Knots);
                end

                self.deriv_curve.preGradient(deriv_time-1);
            end
        end

        function self=translate(self,tran_vctr)
            self.Poles=self.Poles+tran_vctr;
        end

        function self=rotate(self,rotate_matrix)
            self.Poles=self.Poles*rotate_matrix';
        end

    end

    methods % calculate coord
        function U=calCoordinate(self,Points,geom_torl)
            % base on X, Y, Z calculate local coordinate in curve
            %
            if nargin < 3, geom_torl=[];end
            if isempty(geom_torl), geom_torl=sqrt(eps);end

            % find point to start
            U=findNearest(self,Points,20);

            % use project function to adjust parameter
            U=self.projectPoint(Points,geom_torl,U);
        end

        function U=projectPoint(self,Points_init,geom_torl,U)
            % adjust U by Jacobian transformation
            % also can project point to curve
            %
            if nargin < 4
                U=[];
                if nargin < 3
                    geom_torl=[];
                end
            end
            if isempty(geom_torl), geom_torl=100*eps;end

            % find point to start
            if isempty(U)
                U=findNearest(self,Points_init,20);
            end
            [num,~]=size(Points_init);U=U(:);

            iter=0;iter_max=50;
            done=false;idx=1:num;
            while ~done
                [Points,dPoints_dU]=self.calGradient(U(idx));
                dPoint=Points_init(idx,:)-Points;

                RU_RU=sum(dPoints_dU.^2,2);
                RU_D=sum(dPoints_dU.*dPoint,2);

                dU=RU_D./(RU_RU);
                dU(isnan(dU) | isinf(dU))=0;

                U(idx)=U(idx)+dU;
                U=max(U,min(self.Knots));U=min(U,max(self.Knots));

                idx=find(abs(RU_D) > geom_torl & abs(dU) > geom_torl);

                % Points_inv=self.calPoint(U);
                % scatter3(Points_inv(:,1),Points_inv(:,2),Points_inv(:,3));

                iter=iter+1;
                if isempty(idx) || iter >= iter_max
                    done=true;
                end
            end
        end

        function U=findNearest(self,Points,param)
            % find nearest U in grid
            %
            if nargin < 3, param=[];end
            if isempty(param), param=20;end
            num=size(Points,1);

            % generate rough mesh to initialize pre coord
            U_base=((0:(param-1))+(1:param))/2/param;
            U_base=U_base.*(max(self.Knots)-min(self.Knots))+min(self.Knots);
            Point_base=self.calPoint(U_base);
            U=zeros(num,1);
            for p_idx=1:num
                point=Points(p_idx,:);
                dis_abs=sum(abs(Point_base-point),2);
                [~,idx]=min(dis_abs,[],"all");
                U(p_idx)=U_base(idx(1));
            end
        end
    
    end
end
