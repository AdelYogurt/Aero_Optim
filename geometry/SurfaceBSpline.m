classdef SurfaceBSpline
    properties
        control_X;
        control_Y;
        control_Z;
        control_rank_number;
        control_colume_number;
        u_list;
        v_list;
        node_X;
        node_Y;
        node_Z;
        node_rank_number;
        node_colume_number;
        order;
        dimension=3;
    end

    methods
        function self=SurfaceBSpline(control_X,control_Y,control_Z,...
                node_X,node_Y,node_Z,order,u_list,v_list)
            % generate BSpline surface by input control point or node point
            %
            % notice:
            % colume of X, Y, Z is LaWGS format
            % colume of node_X, node_Y, node_Z will by reseve calculate
            %
            % input:
            %
            %
            if nargin < 9
                v_list = [];
                if nargin < 8
                    u_list = [];
                end
            end

            if isempty(control_X) && isempty(control_Y) && isempty(control_Z)
                order=3;
                [node_rank_number,node_colume_number]=size(node_X);
                if any(size(node_Y) ~= [node_rank_number,node_colume_number]) ||...
                        any(size(node_Z) ~= [node_rank_number,node_colume_number])
                    error('SurfaceBSpline: size of control_X, control_Y, control_Z do not equal');
                end
                self.node_X=node_X;
                self.node_Y=node_Y;
                self.node_Z=node_Z;
                self.node_rank_number=node_rank_number;
                self.node_colume_number=node_colume_number;
                self.order=order;

                control_rank_number=node_rank_number+order-1;
                control_colume_number=node_colume_number+order-1;
                control_X=zeros(control_rank_number,control_colume_number);
                control_Y=zeros(control_rank_number,control_colume_number);
                control_Z=zeros(control_rank_number,control_colume_number);

                % create control point by colume of node_X, node_Y, node_Z
                temp_control_X=zeros(control_rank_number,control_colume_number-order+1);
                temp_control_Y=zeros(control_rank_number,control_colume_number-order+1);
                temp_control_Z=zeros(control_rank_number,control_colume_number-order+1);
                if isempty(v_list)
                    v_list=[zeros(1,order),linspace(0,1,control_rank_number-order+1),ones(1,order)];
                end
                for colume_index=1:node_colume_number
                    curve=CurveBSpline([],[node_X(:,colume_index),node_Y(:,colume_index),node_Z(:,colume_index)],order,v_list);
                    temp_control_X(:,colume_index)=curve.control_list(:,1);
                    temp_control_Y(:,colume_index)=curve.control_list(:,2);
                    temp_control_Z(:,colume_index)=curve.control_list(:,3);
                end

                % create control point by colume of node_X, node_Y, node_Z
                if isempty(u_list)
                    u_list=[zeros(1,order),linspace(0,1,control_colume_number-order+1),ones(1,order)];
                end
                for rank_index=1:control_rank_number
                    curve=CurveBSpline([],[temp_control_X(rank_index,:)',temp_control_Y(rank_index,:)',temp_control_Z(rank_index,:)'],order,u_list);
                    control_X(rank_index,:)=curve.control_list(:,1)';
                    control_Y(rank_index,:)=curve.control_list(:,2)';
                    control_Z(rank_index,:)=curve.control_list(:,3)';
                end

                self.u_list=u_list;
                self.v_list=v_list;

                self.control_X=control_X;
                self.control_Y=control_Y;
                self.control_Z=control_Z;
                self.control_rank_number=control_rank_number;
                self.control_colume_number=control_colume_number;

            elseif ~isempty(control_X) && ~isempty(control_Y) && ~isempty(control_Z)
                % generate B spline by control point
                % input control_list, order, u_list(optional)
                %
                [control_rank_number,control_colume_number]=size(control_X);
                if any(size(control_Y) ~= [control_rank_number,control_colume_number]) ||...
                        any(size(control_Z) ~= [control_rank_number,control_colume_number])
                    error('SurfaceBSpline: size of control_X, control_Y, control_Z do not equal');
                end
                if control_rank_number < (order+1) || control_colume_number < (order+1)
                   error('SurfaceBSpline: control_number less than order+1'); 
                end
                self.control_X=control_X;
                self.control_Y=control_Y;
                self.control_Z=control_Z;
                self.control_rank_number=control_rank_number;
                self.control_colume_number=control_colume_number;
                self.order=order;
                
                if isempty(u_list)
                    u_list=[zeros(1,order),linspace(0,1,control_colume_number-order+1),ones(1,order)];
                end
                if isempty(v_list)
                    v_list=[zeros(1,order),linspace(0,1,control_rank_number-order+1),ones(1,order)];
                end

                if (length(u_list) ~= control_colume_number+order+1) ||...
                        (length(v_list) ~= control_rank_number+order+1)
                    error('CurveBSpline: length of do not equal to control_number+order+1');
                end

                self.u_list=u_list;
                self.v_list=v_list;

                % calculate node point list
                node_rank_number=control_rank_number-order+1;
                node_colume_number=control_colume_number-order+1;
                [U_x,V_x]=meshgrid(u_list(order+1:control_colume_number+1),v_list(order+1:control_rank_number+1));
                [self.node_X,self.node_Y,self.node_Z]=self.calPoint(U_x,V_x);
                self.node_rank_number=node_rank_number;
                self.node_colume_number=node_colume_number;

            else
                error('SurfaceBSpline: error input, lack control point or node point');
            end

        end

        function writeStep(self,step_filestr,start_index,head_flag)
            % write BSpline into step file
            %

            % cheak filename
            if length(step_filestr) > 4
                if ~strcmpi(step_filestr((end-3):end),'.inp')
                    step_filestr=[step_filestr,'.inp'];
                end
            else
                step_filestr=[step_filestr,'.inp'];
            end
        end

        function [step_str,object_index,OPEN_SHELL_index]=getStep(self,object_index)
            % write BSpline into step file
            %
            if nargin < 2
                object_index=1;
            end
            step_str=[];
            ctrl_rank_num=self.control_rank_number;
            ctrl_colume_num=self.control_colume_number;
            
            % generate CARTESIAN_POINT
            CARTESIAN_POINT_index=object_index;
            for rank_index=1:ctrl_rank_num
                for colume_index=1:ctrl_colume_num
                    str=[num2str(object_index,'#%d ='),' CARTESIAN_POINT ( ''NONE'', ',...
                        num2str(self.control_X(rank_index,colume_index),'( %.16f, '),...
                        num2str(self.control_Y(rank_index,colume_index),'%.16f, '),...
                        num2str(self.control_Z(rank_index,colume_index),'%.16f )'),...
                        ' );\n'];
                    step_str=[step_str,str];
                    object_index=object_index+1;
                end
            end

            step_str=[step_str,'\n'];

            % generate B_SPLINE_CURVE
            B_SPLINE_CURVE_index=object_index;
            str=getCurveStep(((1:ctrl_colume_num)-1)*ctrl_rank_num+CARTESIAN_POINT_index);
            step_str=[step_str,str];
            object_index=object_index+1;
            str=getCurveStep((((ctrl_colume_num-1)*ctrl_rank_num+1):ctrl_colume_num*ctrl_rank_num)+CARTESIAN_POINT_index-1);
            step_str=[step_str,str];
            object_index=object_index+1;
            str=getCurveStep(((ctrl_colume_num:-1:1))*ctrl_rank_num+CARTESIAN_POINT_index-1);
            step_str=[step_str,str];
            object_index=object_index+1;
            str=getCurveStep((ctrl_rank_num:-1:1)+CARTESIAN_POINT_index-1);
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];

            % generate VERTEX_POINT
            VERTEX_POINT_index=object_index;
            str=[num2str(object_index,'#%d'),' = VERTEX_POINT ',...
                num2str(CARTESIAN_POINT_index,'( ''NONE'', #%d );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = VERTEX_POINT ',...
                num2str((ctrl_colume_num-1)*ctrl_rank_num+CARTESIAN_POINT_index,'( ''NONE'', #%d );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = VERTEX_POINT ',...
                num2str(ctrl_colume_num*ctrl_rank_num+CARTESIAN_POINT_index-1,'( ''NONE'', #%d );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = VERTEX_POINT ',...
                num2str(ctrl_rank_num+CARTESIAN_POINT_index-1,'( ''NONE'', #%d );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];
            
            % generate EDGE_CURVE
            EDGE_CURVE_index=object_index;
            str=[num2str(object_index,'#%d'),' = EDGE_CURVE ',...
                num2str([VERTEX_POINT_index,VERTEX_POINT_index+1,B_SPLINE_CURVE_index],'( ''NONE'', #%d, #%d, #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = EDGE_CURVE ',...
                num2str([VERTEX_POINT_index+1,VERTEX_POINT_index+2,B_SPLINE_CURVE_index+1],'( ''NONE'', #%d, #%d, #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = EDGE_CURVE ',...
                num2str([VERTEX_POINT_index+2,VERTEX_POINT_index+3,B_SPLINE_CURVE_index+2],'( ''NONE'', #%d, #%d, #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = EDGE_CURVE ',...
                num2str([VERTEX_POINT_index+3,VERTEX_POINT_index,B_SPLINE_CURVE_index+3],'( ''NONE'', #%d, #%d, #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];

            % generate ORIENTED_EDGE
            ORIENTED_EDGE_index=object_index;
            str=[num2str(object_index,'#%d'),' = ORIENTED_EDGE ',...
                num2str(EDGE_CURVE_index,'( ''NONE'', *, *, #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = ORIENTED_EDGE ',...
                num2str(EDGE_CURVE_index+1,'( ''NONE'', *, *, #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = ORIENTED_EDGE ',...
                num2str(EDGE_CURVE_index+2,'( ''NONE'', *, *, #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = ORIENTED_EDGE ',...
                num2str(EDGE_CURVE_index+3,'( ''NONE'', *, *, #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];

            % generate EDGE_LOOP
            EDGE_LOOP_index=object_index;
            str=[num2str(object_index,'#%d'),' = EDGE_LOOP ',...
                num2str((0:3)+ORIENTED_EDGE_index,'( ''NONE'', ( #%d, #%d, #%d, #%d ) );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];

            % generate FACE_OUTER_BOUND
            FACE_OUTER_BOUND_index=object_index;
            str=[num2str(object_index,'#%d'),' = FACE_OUTER_BOUND ',...
                num2str(EDGE_LOOP_index,'( ''NONE'', #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];

            % generate B_SPLINE_SURFACE
            B_SPLINE_SURFACE_index=object_index;
            str=[num2str(object_index,'#%d'),' = B_SPLINE_SURFACE ( ''NONE'', 3, 3, (\n'];
            for rank_index=1:ctrl_rank_num-1
                str=[str,'( ',num2str((0:(ctrl_colume_num-2))*ctrl_rank_num+CARTESIAN_POINT_index+rank_index-1,'#%d, '),...
                    ' ',num2str((ctrl_colume_num-1)*ctrl_rank_num+CARTESIAN_POINT_index+rank_index-1,'#%d'),' ),\n'];
            end
            str=[str,'( ',num2str((0:(ctrl_colume_num-2))*ctrl_rank_num+CARTESIAN_POINT_index+ctrl_rank_num-1,'#%d, '),...
                ' ',num2str((ctrl_colume_num-1)*ctrl_rank_num+CARTESIAN_POINT_index+ctrl_rank_num-1,'#%d'),' )'];
            str=[str,'),\n.UNSPECIFIED., .F., .F., .F.);\n'];
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];

            % generate ADVANCED_FACE
            ADVANCED_FACE_index=object_index;
            str=[num2str(object_index,'#%d'),' = ADVANCED_FACE ',...
                num2str([FACE_OUTER_BOUND_index,B_SPLINE_SURFACE_index],'( ''NONE'', ( #%d ), #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];

            % generate OPEN_SHELL
            OPEN_SHELL_index=object_index;
            str=[num2str(object_index,'#%d'),' = OPEN_SHELL ',...
                num2str(ADVANCED_FACE_index,'( ''NONE'', ( #%d ) );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];

            function str=getCurveStep(index_list)
                str=[num2str(object_index,'#%d ='),' B_SPLINE_CURVE(''NONE'', ',...
                    num2str(self.order,'%d, '),...
                    '( ',num2str(index_list(1:end-1),'#%d, '),' ',num2str(index_list(end),'#%d'),' ),',...
                    '.UNSPECIFIED., .F., .F.',...
                    ');\n'];
            end
        end

        function [step_str,object_index,OPEN_SHELL_index]=getStepNode(self,object_index)
            % write BSpline into step file
            %
            if nargin < 2
                object_index=1;
            end
            step_str=[];
            node_rank_num=self.node_rank_number;
            node_colume_num=self.node_colume_number;
            u_node_list=self.u_list(self.order+1:self.control_colume_number+1);
            v_node_list=self.v_list(self.order+1:self.control_rank_number+1);
            
            % generate CARTESIAN_POINT
            CARTESIAN_POINT_index=object_index;

            for colume_index=1:node_colume_num
                for rank_index=1:node_rank_num
                    str=[num2str(object_index,'#%d ='),' CARTESIAN_POINT ( ''NONE'', ',...
                        num2str(self.node_X(rank_index,colume_index),'( %.16f, '),...
                        num2str(self.node_Y(rank_index,colume_index),'%.16f, '),...
                        num2str(self.node_Z(rank_index,colume_index),'%.16f )'),...
                        ' );\n'];
                    step_str=[step_str,str];
                    object_index=object_index+1;
                end
            end

            step_str=[step_str,'\n'];

            % generate B_SPLINE_CURVE_WITH_KNOTS
            B_SPLINE_CURVE_WITH_KNOTS_index=object_index;

            str=getCurveStep(((1:node_colume_num)-1)*node_rank_num+CARTESIAN_POINT_index,...
                node_colume_num,u_node_list);
            step_str=[step_str,str];
            object_index=object_index+1;
            str=getCurveStep((((node_colume_num-1)*node_rank_num+1):(node_colume_num*node_rank_num))+CARTESIAN_POINT_index-1,...
                node_rank_num,v_node_list);
            step_str=[step_str,str];
            object_index=object_index+1;
            str=getCurveStep((node_colume_num:-1:1)*node_rank_num+CARTESIAN_POINT_index-1,...
                node_colume_num,u_node_list);
            step_str=[step_str,str];
            object_index=object_index+1;

            str=getCurveStep((node_rank_num:-1:1)+CARTESIAN_POINT_index-1,...
                node_rank_num,v_node_list);
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];

            % generate VERTEX_POINT
            VERTEX_POINT_index=object_index;

            str=[num2str(object_index,'#%d'),' = VERTEX_POINT ',...
                num2str(CARTESIAN_POINT_index,'( ''NONE'', #%d );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = VERTEX_POINT ',...
                num2str((node_colume_num-1)*node_rank_num+CARTESIAN_POINT_index,'( ''NONE'', #%d );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = VERTEX_POINT ',...
                num2str(node_colume_num*node_rank_num+CARTESIAN_POINT_index-1,'( ''NONE'', #%d );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = VERTEX_POINT ',...
                num2str(node_rank_num+CARTESIAN_POINT_index-1,'( ''NONE'', #%d );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];
            
            % generate EDGE_CURVE
            EDGE_CURVE_index=object_index;

            str=[num2str(object_index,'#%d'),' = EDGE_CURVE ',...
                num2str([VERTEX_POINT_index,VERTEX_POINT_index+1,B_SPLINE_CURVE_WITH_KNOTS_index],'( ''NONE'', #%d, #%d, #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = EDGE_CURVE ',...
                num2str([VERTEX_POINT_index+1,VERTEX_POINT_index+2,B_SPLINE_CURVE_WITH_KNOTS_index+1],'( ''NONE'', #%d, #%d, #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = EDGE_CURVE ',...
                num2str([VERTEX_POINT_index+2,VERTEX_POINT_index+3,B_SPLINE_CURVE_WITH_KNOTS_index+2],'( ''NONE'', #%d, #%d, #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = EDGE_CURVE ',...
                num2str([VERTEX_POINT_index+3,VERTEX_POINT_index,B_SPLINE_CURVE_WITH_KNOTS_index+3],'( ''NONE'', #%d, #%d, #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];

            % generate ORIENTED_EDGE
            ORIENTED_EDGE_index=object_index;

            str=[num2str(object_index,'#%d'),' = ORIENTED_EDGE ',...
                num2str(EDGE_CURVE_index,'( ''NONE'', *, *, #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = ORIENTED_EDGE ',...
                num2str(EDGE_CURVE_index+1,'( ''NONE'', *, *, #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = ORIENTED_EDGE ',...
                num2str(EDGE_CURVE_index+2,'( ''NONE'', *, *, #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;
            str=[num2str(object_index,'#%d'),' = ORIENTED_EDGE ',...
                num2str(EDGE_CURVE_index+3,'( ''NONE'', *, *, #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];

            % generate EDGE_LOOP
            EDGE_LOOP_index=object_index;
            str=[num2str(object_index,'#%d'),' = EDGE_LOOP ',...
                num2str((0:3)+ORIENTED_EDGE_index,'( ''NONE'', ( #%d, #%d, #%d, #%d ) );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];

            % generate FACE_OUTER_BOUND
            FACE_OUTER_BOUND_index=object_index;
            str=[num2str(object_index,'#%d'),' = FACE_OUTER_BOUND ',...
                num2str(EDGE_LOOP_index,'( ''NONE'', #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];

            % generate B_SPLINE_SURFACE
            B_SPLINE_SURFACE_index=object_index;
            str=[num2str(object_index,'#%d'),' = B_SPLINE_SURFACE_WITH_KNOTS ( ''NONE'', 3, 3, (\n'];
            for rank_index=1:node_rank_num-1
                str=[str,'( ',num2str((0:(node_colume_num-2))*node_rank_num+CARTESIAN_POINT_index+rank_index-1,'#%d, '),...
                    ' ',num2str((node_colume_num-1)*node_rank_num+CARTESIAN_POINT_index+rank_index-1,'#%d'),' ),\n'];
            end
            str=[str,'( ',num2str((0:(node_colume_num-2))*node_rank_num+CARTESIAN_POINT_index+node_rank_num-1,'#%d, '),...
                ' ',num2str((node_colume_num-1)*node_rank_num+CARTESIAN_POINT_index+node_rank_num-1,'#%d'),' )'];

            v_node_list=linspace(0,1,node_rank_num-2);
            u_node_list=linspace(0,1,node_colume_num-2);

            str=[str,'),\n.UNSPECIFIED., .F., .F., .F.,\n',...
                '( 4, ',num2str(ones(1,node_rank_num-4),'%d, '),' 4 ), ',...
                '( 4, ',num2str(ones(1,node_colume_num-4),'%d, '),' 4 ), ',...
                '( ',num2str(v_node_list(1:end-1),'%.16f, '),num2str(v_node_list(end),'%.16f'),' ),',...
                '( ',num2str(u_node_list(1:end-1),'%.16f, '),num2str(u_node_list(end),'%.16f'),' ),',...
                ' \n.UNSPECIFIED.);\n'];
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];

            % generate ADVANCED_FACE
            ADVANCED_FACE_index=object_index;
            str=[num2str(object_index,'#%d'),' = ADVANCED_FACE ',...
                num2str([FACE_OUTER_BOUND_index,B_SPLINE_SURFACE_index],'( ''NONE'', ( #%d ), #%d, .T. );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];

            % generate OPEN_SHELL
            OPEN_SHELL_index=object_index;
            str=[num2str(object_index,'#%d'),' = OPEN_SHELL ',...
                num2str(ADVANCED_FACE_index,'( ''NONE'', ( #%d ) );'),'\n'];
            step_str=[step_str,str];
            object_index=object_index+1;

            step_str=[step_str,'\n'];

            function str=getCurveStep(index_list,node_number,u_node_list)
                u_node_list=linspace(0,1,node_number-2);
                str=[num2str(object_index,'#%d ='),' B_SPLINE_CURVE_WITH_KNOTS(''NONE'', ',...
                    num2str(self.order,'%d, '),...
                    '( ',num2str(index_list(1:end-1),'#%d, '),' ',num2str(index_list(end),'#%d'),' ),\n',...
                    '.UNSPECIFIED., .F., .F.,','( 4, ',num2str(ones(1,node_number-4),'%d, '),' 4 ), ',...
                    '( ',num2str(u_node_list(1:end-1),'%.16f, '),num2str(u_node_list(end),'%.16f'),' ), \n',...
                    '.UNSPECIFIED.);\n'];
            end
        end

    end

    % calculate point function
    methods
        function [X,Y,Z]=calSurface(self,varargin)
            % generate surface matrix by u_x_list or point_number
            %
            % default
            % xi_list=linspace(0,1,xi_gird_num(default is 50)+1)
            % psi_list=linspace(0,1,psi_gird_num(default is 50)+1)
            % [XI,PSI]=meshgrid(xi_list,psi_list), colume is LaWGS format line
            %
            % input:
            % U_x, V_x
            % u_x_number, v_x_number
            %
            % output:
            % X,Y,Z (colume is LaWGS format line)
            %

            % xi
            if nargin < 2 || isempty(varargin{1})
                % mean donot input xi_grid_number, use default number
                U_x=[];
                u_x_number=50;
            else
                if length(varargin{1}) == 1
                    U_x=[];
                    u_x_number=varargin{1};
                else
                    % mean input XI matrix
                    U_x=varargin{1};
                    u_x_number=size(U_x,2)-1;
                end
            end

            % psi
            if nargin < 3 || isempty(varargin{2})
                % mean donot input psi_grid_numebr, use default number
                V_x=[];
                v_x_number=50;
            else
                if length(varargin{2}) == 1
                    V_x=[];
                    v_x_number=varargin{2};
                else
                    % mean input PSI matrix
                    V_x=varargin{2};
                    v_x_number=size(V_x,1)-1;
                end
            end

            % calculate local coordinate matrix
            if isempty(U_x)
                U_x=repmat(linspace(0,1,u_x_number+1),v_x_number+1,1);
            end
            if isempty(V_x)
                V_x=repmat(linspace(0,1,v_x_number+1)',1,u_x_number+1);
            end

            [X,Y,Z]=calPoint(self,U_x,V_x);
        end

        function [X,Y,Z]=calPoint(self,U_x,V_x)
            % according u_x to calculate point
            % u_x_list is u_x_number x 1 matrix
            % point_list is point_number x dimension matrix
            %
            [rank_num,colume_num]=size(U_x);
            if any(size(V_x) ~= [rank_num,colume_num])
                error('SurfaceBSpline.calPoint: size of U_x do not equal to size of V_x');
            end

            X=zeros(rank_num,colume_num);
            Y=zeros(rank_num,colume_num);
            Z=zeros(rank_num,colume_num);

            N_u_list=zeros(self.order+1,1);
            N_v_list=zeros(1,self.order+1);
            for rank_idx=1:rank_num
                for colume_idx=1:colume_num
                    % local index of u_x in u_list, v_list
                    u_x=U_x(rank_idx,colume_idx);
                    v_x=V_x(rank_idx,colume_idx);
                    
                    %                     [index_end_v,index_end_u]=getIndex(); % y, x

                    index_end_u=self.control_colume_number; % is equal to the section index
                    while index_end_u > self.order+1 && u_x < self.u_list(index_end_u)
                        index_end_u=index_end_u-1;
                    end
                    index_end_v=self.control_rank_number; % is equal to the section index
                    while index_end_v > self.order+1 && v_x < self.v_list(index_end_v)
                        index_end_v=index_end_v-1;
                    end

                    index_start_u=index_end_u-self.order;
                    index_start_v=index_end_v-self.order;

                    % calculate base function
                    for N_index=1:self.order+1
                        N_u_list(N_index)=baseFunction(self.u_list,u_x,N_index+index_start_u-1,self.order);
                        N_v_list(N_index)=baseFunction(self.v_list,v_x,N_index+index_start_v-1,self.order);
                    end

                    X(rank_idx,colume_idx)=N_v_list*self.control_X(index_start_v:index_end_v,index_start_u:index_end_u)*N_u_list/sum(N_u_list)/sum(N_v_list);
                    Y(rank_idx,colume_idx)=N_v_list*self.control_Y(index_start_v:index_end_v,index_start_u:index_end_u)*N_u_list/sum(N_u_list)/sum(N_v_list);
                    Z(rank_idx,colume_idx)=N_v_list*self.control_Z(index_start_v:index_end_v,index_start_u:index_end_u)*N_u_list/sum(N_u_list)/sum(N_v_list);
                end
            end

%             function [index_end_v,index_end_u]=getIndex()
%                 for index_end_u=self.order+1:self.control_number % is equal to the section index
%                     for index_end_v=self.order+1:self.control_number % is equal to the section index
%                         if u_x >= self.U(index_end_v,index_end_u) && v_x >= self.V(index_end_v,index_end_u)
%                             return;
%                         end
%                     end
%                 end
%             end
        end

        function drawSurface(self,surface_option,control_option,figure_handle,U_x,V_x,node_option)
            % draw curve on figure handle
            %
            if nargin < 7
                node_option=[];
                if nargin < 5
                    V_x=[];
                    if nargin < 5
                        U_x=[];
                        if nargin < 4
                            figure_handle=[];
                            if nargin < 3
                                control_option=[];
                                if nargin < 2
                                    surface_option=[];
                                end
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
            if isempty(surface_option)
                surface_option=struct('LineStyle','none');
            end
            if isempty(node_option)
                node_option=struct('Marker','o','LineStyle','none','MarkerEdgeColor','b','FaceAlpha',0);
            end
            if isempty(control_option)
                control_option=struct('Marker','s','EdgeColor','r','LineStyle','--','MarkerEdgeColor','r','FaceAlpha',0);
            end

            % calculate point on curve
            [X,Y,Z]=calSurface(self,U_x,V_x);

            % plot surface
            surface(axes_handle,X,Y,Z,surface_option);
            surface(axes_handle,self.node_X,self.node_Y,self.node_Z,node_option);
            surface(axes_handle,self.control_X,self.control_Y,self.control_Z,control_option);
            view(3);

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

