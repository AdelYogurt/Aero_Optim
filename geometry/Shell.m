classdef Shell < handle
    % shell
    %
    properties
        name;
        face_list; % cell, can be surface CST3D or surface BSpline
    end

    % define Shell
    methods
        function self=Shell(name,face_list)
            if nargin < 2,face_list={};end
            self.name=name;
            self.face_list=face_list;
        end

        function fce_list=calShell(self,u_param,v_param)
            % calculate all face point
            %
            if nargin < 3
                v_param=[];
                if nargin < 2
                    u_param=[];
                end
            end

            fce_list=[];
            fce_num=length(self.face_list);
            for fce_idx=1:fce_num
                fce=self.face_list{fce_idx};
                if ~isempty(fce)
                    fce_list=[fce_list,{fce.calFace(u_param,v_param)}];
                end
            end
        end

        function [surf,surf_idx]=getFace(self,surf_name)
            % load surface from face_list base on input surface name
            %
            for surf_idx=1:length(self.face_list)
                surf=self.face_list{surf_idx};
                if strcmp(surf.name,surf_name)
                    return;
                end
            end
            surf=[];
            surf_idx=0;
        end
    end

    % visualizate function
    methods
        function drawShell(self,axe_hdl,u_param,v_param,crv_option,ctrl_option)
            % draw all surface of shell
            % wrapper of drawFace
            %
            if nargin < 6
                ctrl_option=[];
                if nargin < 5
                    crv_option=[];
                    if nargin < 4
                        v_param=[];
                        if nargin < 3
                            u_param=[];
                            if nargin < 2
                                axe_hdl=[];
                            end
                        end
                    end
                end
            end

            if isempty(axe_hdl),axe_hdl=axes(figure());end

            surf_num=length(self.face_list);
            for surf_idx=1:surf_num
                surf=self.face_list{surf_idx};
                if ~isempty(surf)
                    surf.drawFace(axe_hdl,u_param,v_param,crv_option,ctrl_option);
                end
            end
            
            xlabel('x');
            ylabel('y');
            zlabel('z');
            view(3);

            % axis equal;
            % x_range=xlim();
            % y_range=ylim();
            % z_range=zlim();
            % center=[mean(x_range),mean(y_range),mean(z_range)];
            % range=max([x_range(2)-x_range(1),y_range(2)-y_range(1),z_range(2)-z_range(1)])/2;
            % xlim([center(1)-range,center(1)+range]);
            % ylim([center(2)-range,center(2)+range]);
            % zlim([center(3)-range,center(3)+range]);
        end

        function obj_idx=writeStepHead(~,step_file,obj_idx,step_filename)
            % wite head of step
            %

            fprintf(step_file,'ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION (( ''STEP AP203'' ),''1'' );\nFILE_NAME (''%s'',''%s'',( '''' ),( '''' ),''Matlab step'',''Matlab'','''' );\nFILE_SCHEMA (( ''CONFIG_CONTROL_DESIGN'' ));\nENDSEC;\n',step_filename,date);
            fprintf(step_file,'\n');
            fprintf(step_file,'DATA;\n');

            fprintf(step_file,'#%d = CARTESIAN_POINT ( ''NONE'',  ( 0.0, 0.0, 0.0 ) ) ;\n',obj_idx);obj_idx=obj_idx+1;
            fprintf(step_file,'#%d = DIRECTION ( ''NONE'',  ( 0.0, 0.0, 1.0 ) ) ;\n',obj_idx);obj_idx=obj_idx+1;
            fprintf(step_file,'#%d = DIRECTION ( ''NONE'',  ( 1.0, 0.0, 0.0 ) ) ;\n',obj_idx);obj_idx=obj_idx+1;
            fprintf(step_file,'#%d = AXIS2_PLACEMENT_3D ( ''NONE'', #%d, #%d, #%d ) ;\n',obj_idx,1,2,3);obj_idx=obj_idx+1;
            fprintf(step_file,'\n');
            fprintf(step_file,'#%d =( LENGTH_UNIT ( ) NAMED_UNIT ( * ) SI_UNIT ( $., .METRE. ) );\n',obj_idx);obj_idx=obj_idx+1;
            fprintf(step_file,'#%d =( NAMED_UNIT ( * ) PLANE_ANGLE_UNIT ( ) SI_UNIT ( $, .RADIAN. ) );\n',obj_idx);obj_idx=obj_idx+1;
            fprintf(step_file,'#%d =( NAMED_UNIT ( * ) SI_UNIT ( $, .STERADIAN. ) SOLID_ANGLE_UNIT ( ) );\n',obj_idx);obj_idx=obj_idx+1;
            fprintf(step_file,'#%d = UNCERTAINTY_MEASURE_WITH_UNIT (LENGTH_MEASURE( 1.0E-05 ), #%d, ''distance_accuracy_value'', ''NONE'');\n',obj_idx,5);obj_idx=obj_idx+1;
            fprintf(step_file,'#%d =( GEOMETRIC_REPRESENTATION_CONTEXT ( 3 ) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT ( ( #%d ) ) GLOBAL_UNIT_ASSIGNED_CONTEXT ( ( #%d, #%d, #%d ) ) REPRESENTATION_CONTEXT ( ''NONE'', ''WORKASPACE'' ) );\n',obj_idx,8,5,6,7);obj_idx=obj_idx+1;
            fprintf(step_file,'\n');
        end

        function writeStepOpenShell(self,step_filestr)
            % write surface into step file
            %
            surf_num=length(self.face_list);
            [~,step_filename,~]=fileparts(step_filestr);

            % write head
            step_file=fopen(step_filestr,'w');
            obj_idx=1;
            obj_idx=writeStepHead(self,step_file,obj_idx,step_filename);

            % write surface
            ADVANCED_FACE_list=zeros(1,surf_num);
            
            for surf_idx=1:surf_num
                surf=self.face_list{surf_idx};

                [step_str,obj_idx,ADVANCED_FACE]=surf.getStep(obj_idx);
                fprintf(step_file,step_str);

                ADVANCED_FACE_list(surf_idx)=ADVANCED_FACE;
            end

            % generate OPEN_SHELL
            OPEN_SHELL=obj_idx;
            step_str=[num2str(obj_idx,'#%d ='),' OPEN_SHELL',...
                ' ( ',...
                '''NONE''',', ',...
                '( ',num2str(ADVANCED_FACE_list(1),'#%d'),num2str(ADVANCED_FACE_list(2:end),', #%d'),' )',...
                ' ) ;\n'];obj_idx=obj_idx+1;
            fprintf(step_file,step_str);
            fprintf(step_file,'\n');

            % write model
            SHELL_BASED_SURFACE_MODEL=obj_idx;
            step_str=[num2str(obj_idx,'#%d ='),' SHELL_BASED_SURFACE_MODEL',...
                ' ( ',...
                '''NONE''',', ',...
                '( ',num2str(OPEN_SHELL,'#%d'),' )',...
                ' ) ;\n'];obj_idx=obj_idx+1;
            fprintf(step_file,step_str);
            fprintf(step_file,'\n');

            % write end of step file
            writeStepEnd(self,step_file,obj_idx,step_filename,SHELL_BASED_SURFACE_MODEL);

            fclose(step_file);
            clear('step_file');
        end

        function obj_idx=writeStepEnd(~,step_file,obj_idx,step_filename,model_index)
            % write end of step file
            %

            % write product context 
            fprintf(step_file,'#%d = APPLICATION_CONTEXT ( ''configuration controlled 3d designs of mechanical parts and assemblies'' ) ;\n',obj_idx);obj_idx=obj_idx+1;
            fprintf(step_file,'#%d = MECHANICAL_CONTEXT ( ''NONE'', #%d, ''mechanical'' ) ;\n',obj_idx,obj_idx-1);obj_idx=obj_idx+1;
            fprintf(step_file,'#%d = PRODUCT ( ''%s'', ''%s'', '''', ( #%d ) ) ;\n',obj_idx,step_filename,step_filename,obj_idx-1);obj_idx=obj_idx+1;
            fprintf(step_file,'\n');
            fprintf(step_file,'#%d = APPLICATION_CONTEXT ( ''configuration controlled 3d designs of mechanical parts and assemblies'' ) ;\n',obj_idx);obj_idx=obj_idx+1;
            fprintf(step_file,'#%d = DESIGN_CONTEXT ( ''detailed design'', #%d, ''design'' ) ;\n',obj_idx,obj_idx-1);obj_idx=obj_idx+1;
            fprintf(step_file,'#%d = PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE ( ''NONE'', '''', #%d, .NOT_KNOWN. ) ;\n',obj_idx,obj_idx-3);obj_idx=obj_idx+1;
            fprintf(step_file,'#%d = PRODUCT_DEFINITION ( ''NONE'', '''', #%d, #%d ) ;\n',obj_idx,obj_idx-1,obj_idx-2);obj_idx=obj_idx+1;
            fprintf(step_file,'\n');
            fprintf(step_file,'#%d = PRODUCT_DEFINITION_SHAPE ( ''NONE'', ''NONE'',  #%d ) ;\n',obj_idx,obj_idx-1);obj_idx=obj_idx+1;
            fprintf(step_file,'#%d = MANIFOLD_SURFACE_SHAPE_REPRESENTATION ( ''test'', ( #%d, #%d ), #%d ) ;\n',obj_idx,model_index,4,9);obj_idx=obj_idx+1;
            fprintf(step_file,'#%d = SHAPE_DEFINITION_REPRESENTATION ( #%d, #%d ) ;\n',obj_idx,obj_idx-2,obj_idx-1);obj_idx=obj_idx+1;

            % write end
            fprintf(step_file,'\n');
            fprintf(step_file,'ENDSEC;\nEND-ISO-10303-21;\n');

        end
    end
end