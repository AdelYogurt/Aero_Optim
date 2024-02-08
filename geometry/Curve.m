classdef Curve < handle
    % Non-Uniform Rational B-Splines Curve
    % define reference to step standard
    %
    properties % Explicit Attributes
        name='';
        Degree=[]; % degree
        Poles=[]; % control_points_list
        edge_form='.UNSPECIFIED.'; % edge_form
        Periodic='.F.'; % closed_curve (boolean)
        self_intersect='.F. '; % self_intersect (boolean)
        Weights=[]; % weights_data

        Mults=[]; % knot_multiplicities
        Knots=[]; % Knots
        knot_spec='.UNSPECIFIED.';

        u_list=[]; % knot_vector
    end

    % define curve
    methods
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
            % if input Degree is empty,...
            % Degree default is pole_num-1,...
            % which is Bezier curve.
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

        function Points=calPoint(self,U_x)
            % calculate point on Non-Uniform Rational B-Splines Curve
            %
            U_x=U_x(:);
            Ctrls=[self.Poles.*self.Weights,self.Weights];

            % evaluate along the u direction
            Points=calBSpline(Ctrls,self.Degree,self.u_list,U_x);
            Points=Points(:,1:end-1)./Points(:,end);

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

    % control Edge
    methods
        function self=reverse(self)
            % revese direction of curve
            %
            self.Poles=self.Poles(end:-1:1,:);
            self.Mults=self.Mults(end:-1:1);
            self.Knots=max(self.Knots)-self.Knots(end:-1:1);
            self.Weights=self.Weights(end:-1:1,:);
            self.u_list=max(self.u_list)-self.u_list(end:-1:1);
        end

        function self=addDegree(self,deg_tar)
            % increase curve degree
            %
            if deg_tar <= self.Degree, return;end

            % modify Degree, Mults and Ctrls
            Ctrls=[self.Poles.*self.Weights,self.Weights];
            [Ctrls,self.Mults]=GeomApp.addDegree(self.Degree,Ctrls,self.Mults,self.Knots,deg_tar);
            self.Degree=deg_tar;
            self.u_list=baseKnotVec(self.Mults,self.Knots);

            % modify Poles and Weights
            self.Poles=Ctrls(:,1:end-1)./Ctrls(:,end);
            self.Weights=Ctrls(:,end);
        end

        function self=insertKnot(self,U_ins)
            % insert knot to curve
            %
            if isempty(U_ins), return;end

            % modify u_list, Mults, Knots and Ctrls
            Ctrls=[self.Poles.*self.Weights,self.Weights];
            [Ctrls,self.u_list]=GeomApp.insertKnot(self.Degree,Ctrls,self.u_list,U_ins);
            self.Knots=unique(self.u_list);
            self.Mults=ones(size(self.Knots));
            for k_idx=1:length(self.Knots)
                self.Mults(k_idx)=sum(self.u_list == self.Knots(k_idx));
            end

            % modify Poles and Weights
            self.Poles=Ctrls(:,1:end-1)./Ctrls(:,end);
            self.Weights=Ctrls(:,end);
        end
    
        function self=translate(self,tran_vctr)
            self.Poles=self.Poles+tran_vctr;
        end

        function self=rotate(self,rotate_matrix)
            self.Poles=self.Poles*rotate_matrix';
        end

    end
end
