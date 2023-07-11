classdef CurveBezier < handle
    properties
        control_list;
        control_number;
        u_list;
        order;
        dimension;
    end
    methods
        function self=CurveBezier(control_list,node_list,order,u_list)
            % generate Bezier
            %
            % input:
            % [], node_list(Bezier will across), order, u_list(optional)
            % control_list
            %
            % notice:
            % u_list default is linspace(0,1,node_number)
            %
            if nargin < 4
                u_list = [];
            end

            if isempty(control_list)
                % fitting Bezier by node point
                % input node_list, order
                % fitting method is least square method
                %
                [node_number,dimension]=size(node_list);

                if order > node_number-1
                    error('CurveBezier: curve order large than node number-1')
                end

                self.dimension=dimension;
                self.order=order;

                if isempty(u_list)
                    u_list=linspace(0,1,node_number);
                end

                % generate fitting matrix
                matrix=zeros(node_number,order+1);
                for node_index=1:node_number
                    matrix(node_index,:)=self.baseFunction(u_list(node_index));
                end

                control_list=matrix\node_list;
                control_number=order+1;

                self.control_list=control_list;
                self.control_number=control_number;
            else
                % generate Bezier by control point
                % input control_list

                [control_number,dimension]=size(control_list);
                if control_number < 2
                    error('CurveBezier: control_number less than 2');
                end

                order=control_number-1;

                self.control_list=control_list;
                self.control_number=control_number;
                self.order=order;
                self.dimension=dimension;
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
            % according u_x calculate point
            % u_x_list is u_x_number x 1 matrix
            % point_list is point_number x dimension matrix
            %
            u_x_number=length(u_x_list);
            point_list=zeros(u_x_number,self.dimension);
            for u_x_index=1:u_x_number
                u_x=u_x_list(u_x_index);
                P_list=self.baseFunction(u_x);
                point_list(u_x_index,:)=P_list*self.control_list;
            end
        end

        function drawCurve(self,line_option,control_option,figure_handle,u_x_list)
            % draw curve on figure handle
            %
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
            if isempty(control_option)
                control_option=struct('Marker','s','LineStyle','--','Color','r');
            end

            % calculate point on curve
            point_list=calCurve(self,u_x_list);

            % plot line
            if self.dimension == 2
                line(axes_handle,point_list(:,1),point_list(:,2),line_option);
                line(axes_handle,self.control_list(:,1),self.control_list(:,2),control_option);

            elseif self.dimension == 3
                line(axes_handle,point_list(:,1),point_list(:,2),point_list(:,3),line_option);
                line(axes_handle,self.control_list(:,1),self.control_list(:,2),self.control_list(:,3),control_option);
                view(3);

            end

        end
    end

    % base function
    methods
        function P_list=baseFunction(self,u_x)
            order=self.order;
            b_list=zeros(1,order+1);
            x_list=zeros(1,order+1);
            middle_index=ceil((order+1)/2);

            % calculate binom
            for b_index=0:middle_index-1
                b_list(b_index+1)=nchoosek(order,b_index);
            end

            % symmetry
            b_list(middle_index+1:end)=b_list(floor((order+1)/2):-1:1);

            % calculate x
            for x_index=0:order
                x_list(x_index+1)=u_x^(x_index)*(1-u_x)^(order-x_index);
            end

            P_list=b_list.*x_list;
        end
    end
end
