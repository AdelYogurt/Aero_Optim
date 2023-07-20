classdef CurveBSpline < handle
    % B-spline curve 
    % define reference to step standard
    %
    properties
        name='';
        control_list;
        control_number;
        degree;
    end

    properties
        u_list;
        node_list;
        node_number;
        dimension;
    end

    % define function
    methods
        function self=CurveBSpline(name,control_list,node_list,degree,u_list)
            % generate B spline
            %
            % input:
            % [], node_list(Bspline will across), degree, u_list(optional)
            % control_list, [], degree, u_list(optional)
            %
            % notice:
            % u_list default is [zeros(1,degree),u_node_list,ones(1,degree)];
            %
            if nargin < 5
                u_list = [];
            end
            self.name=name;

            if isempty(control_list)
                % base on node point list inverse control point list
                [node_number,dimension]=size(node_list);
                if node_number == 2
                    degree=1;
                    self.degree=degree;
                    self.dimension=dimension;

                    u_list=[0,0,1,1];
                    
                    self.node_list=node_list;
                    self.node_number=node_number;
                    self.u_list=u_list;

                    control_list=node_list;
                    control_number=node_number;

                    self.control_list=control_list;
                    self.control_number=control_number;
                else
                    degree=3;
                    control_number=node_number+degree-1;

                    self.degree=degree;
                    self.dimension=dimension;

                    % base on curve curvature to calculate correct k of u_list
                    % modified sine length parameterization method
                    line_list=node_list(2:end,:)-node_list(1:(end-1),:);
                    length_list=sqrt(sum(line_list.^2,2));
                    line_vector_list=line_list./length_list;

                    % calculate angle min(pi-theta,pi/2);
                    theta_list=zeros(node_number-2,1);
                    for theta_index=1:node_number-2
                        cos_theta=line_vector_list(theta_index,:)*line_vector_list(theta_index+1,:)';
                        if cos_theta > 1
                            theta=0;
                        elseif cos_theta < -1
                            theta=pi;
                        else
                            theta=acos(cos_theta);
                        end
                        theta_list(theta_index)=min(pi-theta,pi/2);
                    end

                    % calcuate curvature
                    curvature_list=zeros(node_number-1,1);
                    curvature_list(1)=1+3*(theta_list(1)*length_list(1)/(length_list(1)+length_list(2)));
                    for curvature_index=2:node_number-2
                        curvature_list(curvature_index)=1+3/2*...
                            (theta_list(curvature_index-1)*length_list(curvature_index-1)/(length_list(curvature_index)+length_list(curvature_index-1))+...
                            theta_list(curvature_index)*length_list(curvature_index+1)/(length_list(curvature_index)+length_list(curvature_index+1)));
                    end
                    curvature_list(end)=1+3*(theta_list(end)*length_list(end)/(length_list(end)+length_list(end-1)));

                    %                 length_total=sum(length_list); % if just use
                    length_curvature_total=sum(length_list.*curvature_list);

                    % calculate u list
                    if isempty(u_list)
                        u_list=zeros(1,node_number+2*degree);
                        for i=node_number+degree+1:node_number+2*degree
                            u_list(i)=1;
                        end

                        for i=degree+2:node_number+degree
                            %                     u_list(i)=u_list(i-1)+length_list(i-degree-1)/length_total;
                            u_list(i)=u_list(i-1)+curvature_list(i-degree-1)*length_list(i-degree-1)/length_curvature_total;
                        end
                    end

                    if length(u_list) ~= control_number+degree+1
                        error('CurveBSpline: length of do not equal to control_number+degree+1');
                    end
                    u_delta_list=u_list(2:end)-u_list(1:end-1);

                    self.node_list=node_list;
                    self.node_number=node_number;
                    self.u_list=u_list;

                    % the start and end tangent vector is undefined
                    % use predict method to calculate
                    % concrete content see the reserve calculate of B-Spline

                    tangent_vector_list=zeros(node_number,dimension);
                    % % calculate tangent_vector beside start and end Fmill method
                    % for tangent_vector_index=2:node_number-1
                    %     vector_temp=node_list(tangent_vector_index+1,:)-node_list(tangent_vector_index-1,:);
                    %     tangent_vector_list(tangent_vector_index,:)=vector_temp/norm(vector_temp,2);
                    % end

                    % calculate tangent_vector beside start and end Bessell method
                    % see [1]王飞.计算机图形学[M].北京邮电大学出版社,2011,124.
                    for tangent_vector_index=2:node_number-1
                        u_length=u_delta_list(degree+tangent_vector_index-1)+u_delta_list(degree+tangent_vector_index);
                        tangent_vector_list(tangent_vector_index,:)=...
                            u_delta_list(degree+tangent_vector_index)/u_length*...
                            line_list(tangent_vector_index-1,:)/u_delta_list(degree-1+tangent_vector_index)+...
                            u_delta_list(degree-1+tangent_vector_index)/u_length*...
                            line_list(tangent_vector_index,:)/u_delta_list(degree+tangent_vector_index);
                    end

                    % calculate tangent_vector start and end
                    % use Parabolic end Boundary condition
                    % see [1]王飞.计算机图形学[M].北京邮电大学出版社,2011,125/128.
                    tangent_vector_list(1,:)=2*line_list(1,:)/u_delta_list(degree+1)-tangent_vector_list(2,:);
                    tangent_vector_list(end,:)=2*line_list(end,:)/u_delta_list(node_number+degree-1)-tangent_vector_list(end-1,:);

                    % tangent_vector_list(1,:)=(3*line_list(1,:)/u_delta_list(degree+1)-tangent_vector_list(2,:))/2;
                    % tangent_vector_list(end,:)=(3*line_list(end,:)/u_delta_list(node_number+degree-1)-tangent_vector_list(end-1,:))/2;

                    matrix=zeros(control_number);
                    matrix(1,1)=1;
                    matrix(end,end)=1;
                    for rank_index=3:control_number-2
                        u_x=u_list(rank_index+degree-1);
                        for k_index=1:degree
                            matrix(rank_index,k_index+rank_index-2)=baseFunction(u_list,u_x,k_index+rank_index-2,degree)/2;
                        end
                    end
                    matrix(2,1)=-1;matrix(2,2)=1;
                    matrix(end-1,end-1)=-1;matrix(end-1,end)=1;
                    colume=[node_list(1,:);
                        tangent_vector_list(1,:)*u_delta_list(degree+1)/3;
                        node_list(2:end-1,:);
                        tangent_vector_list(end,:)*u_delta_list(node_number+degree-1)/3;
                        node_list(end,:);];
                    control_list=matrix\colume;

                    self.control_list=control_list;
                    self.control_number=control_number;
                end
            else
                % generate B spline by control point
                % input control_list, degree, u_list(optional)
                %
                [control_number,dimension]=size(control_list);
                if control_number < (degree+1)
                    error('CurveBSpline: control_number less than degree+1');
                end
                self.control_list=control_list;
                self.control_number=control_number;
                self.degree=degree;
                self.dimension=dimension;

                if isempty(u_list)
                    u_list=[zeros(1,degree),linspace(0,1,control_number-degree+1),ones(1,degree)];
                end

                if length(u_list) ~= control_number+degree+1
                    error('CurveBSpline: length of do not equal to control_number+degree+1');
                end

                self.u_list=u_list;

                % calculate node point list
                node_number=control_number-degree+1;
                node_list=self.calPoint(u_list(degree+1:control_number+1));

                self.node_list=node_list;
                self.node_number=node_number;
            end
        end
    end

    % calculate point function
    methods
        function point_list=calCurve(self,u_x_list)
            % generate curve matrix by u_x_list or point_number
            %
            if nargin < 2 || isempty(u_x_list)
                u_x_list=100;
            end

            if length(u_x_list) == 1
                u_x_list=linspace(0,1,u_x_list);
            end

            point_list=calPoint(self,u_x_list);
        end

        function point_list=calPoint(self,u_x_list)
            % according u_x to calculate point
            % u_x_list is u_x_number x 1 matrix
            % point_list is point_number x dimension matrix
            %
            point_number=length(u_x_list);
            point_list=zeros(point_number,self.dimension);

            N_u_list=zeros(1,self.degree+1);
            for point_index=1:point_number

                % local u_x in u_list index
                u_x=u_x_list(point_index);
                index_end=self.control_number; % is equal to the section index
                while index_end > self.degree+1 && u_x < self.u_list(index_end)
                    index_end=index_end-1;
                end

                %                 for index_end=self.degree+1:self.control_number % is equal to the section index
                %                     if u_x >= self.u_list(index_end)
                %                         break;
                %                     end
                %                 end

                index_start=index_end-self.degree;

                for N_index=1:self.degree+1
                    N_u_list(N_index)=baseFunction(self.u_list,u_x,N_index+index_start-1,self.degree);
                end

                %                 if any(u_x == self.u_list(self.degree+1:self.control_number+1))
                %                     N_u_list=N_u_list/2;
                %                 end

                point_list(point_index,:)=N_u_list*self.control_list(index_start:index_end,:)/sum(N_u_list);
            end
        end

        function drawCurve(self,line_option,control_option,figure_handle,u_x_list,node_option)
            % draw curve on figure handle
            %
            if nargin < 6
                node_option=[];
                if nargin < 5
                    u_x_list=[];
                    if nargin < 4
                        figure_handle=[];
                        if nargin < 3
                            control_option=[];
                            if nargin < 2
                                line_option=[];
                            end
                        end
                    end
                end
            end

            if isempty(figure_handle)
                figure_handle=figure(101);
            end
            axes_handle=figure_handle.CurrentAxes;
            if isempty(axes_handle)
                axes_handle=axes(figure_handle);
            end

            % default draw option
            if isempty(line_option)
                line_option=struct();
            end
            if isempty(node_option)
                node_option=struct('Marker','o','LineStyle','none');
            end
            if isempty(control_option)
                control_option=struct('Marker','s','LineStyle','--','Color','r');
            end

            % calculate point on curve
            point_list=calCurve(self,u_x_list);

            % plot line
            if self.dimension == 2
                line(axes_handle,point_list(:,1),point_list(:,2),line_option);
                line(axes_handle,self.node_list(:,1),self.node_list(:,2),node_option);
                line(axes_handle,self.control_list(:,1),self.control_list(:,2),control_option);

            elseif self.dimension == 3
                line(axes_handle,point_list(:,1),point_list(:,2),point_list(:,3),line_option);
                line(axes_handle,self.node_list(:,1),self.node_list(:,2),self.node_list(:,3),node_option);
                line(axes_handle,self.control_list(:,1),self.control_list(:,2),self.control_list(:,3),control_option);
                view(3);

            end

        end
    end
end

% base function
function N=baseFunction(u_list,u_x,i,k)
if k==0
    if ((u_list(i) <= u_x) && (u_x <= u_list(i+1)))
        N=1;
    else
        N=0;
    end
else
    if u_list(i+k) == u_list(i)
        A=0;
    else
        A=(u_x-u_list(i))/(u_list(i+k)-u_list(i));
    end

    if u_list(i+k+1) == u_list(i+1)
        B=0;
    else
        B=(u_list(i+k+1)-u_x)/(u_list(i+k+1)-u_list(i+1));
    end

    N=A*baseFunction(u_list,u_x,i,k-1)+B*baseFunction(u_list,u_x,i+1,k-1);
end
end
