classdef Body < handle
    properties
        surface_list; % cell, can be surface CST3D or surface BSpline
    end

    methods
        function self=Body(surface_list)
            if nargin < 1
                surface_list={};
            end

            self.surface_list=surface_list;
        end

        function [surf,surf_idx]=getSurface(self,surf_name)
            % load surface from surface_list base on input surface name
            %
            for surf_idx=1:length(self.surface_list)
                surf=self.surface_list{surf_idx};
                if strcmp(surf.name,surf_name)
                    return;
                end
            end
            surf=[];
            surf_idx=0;
        end

    end

    methods
        function object_index=writeStepHead(self,step_file,object_index,step_filename)
            % wite head of step
            %

            fprintf(step_file,'ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION (( ''STEP AP203'' ),''1'' );\nFILE_NAME (''%s'',''%s'',( '''' ),( '''' ),''Matlab step'',''Matlab'','''' );\nFILE_SCHEMA (( ''CONFIG_CONTROL_DESIGN'' ));\nENDSEC;\n',step_filename,date);
            fprintf(step_file,'\n');
            fprintf(step_file,'DATA;\n');
            
            fprintf(step_file,'#%d = CARTESIAN_POINT ( ''NONE'',  ( 0.0, 0.0, 0.0 ) ) ;\n',object_index);object_index=object_index+1;
            fprintf(step_file,'#%d = DIRECTION ( ''NONE'',  ( 0.0, 0.0, 1.0 ) ) ;\n',object_index);object_index=object_index+1;
            fprintf(step_file,'#%d = DIRECTION ( ''NONE'',  ( 1.0, 0.0, 0.0 ) ) ;\n',object_index);object_index=object_index+1;
            fprintf(step_file,'#%d = AXIS2_PLACEMENT_3D ( ''NONE'', #%d, #%d, #%d ) ;\n',object_index,1,2,3);object_index=object_index+1;
            fprintf(step_file,'\n');
            fprintf(step_file,'#%d =( LENGTH_UNIT ( ) NAMED_UNIT ( * ) SI_UNIT ( $., .METRE. ) );\n',object_index);object_index=object_index+1;
            fprintf(step_file,'#%d =( NAMED_UNIT ( * ) PLANE_ANGLE_UNIT ( ) SI_UNIT ( $, .RADIAN. ) );\n',object_index);object_index=object_index+1;
            fprintf(step_file,'#%d =( NAMED_UNIT ( * ) SI_UNIT ( $, .STERADIAN. ) SOLID_ANGLE_UNIT ( ) );\n',object_index);object_index=object_index+1;
            fprintf(step_file,'#%d = UNCERTAINTY_MEASURE_WITH_UNIT (LENGTH_MEASURE( 1.0E-05 ), #%d, ''distance_accuracy_value'', ''NONE'');\n',object_index,5);object_index=object_index+1;
            fprintf(step_file,'#%d =( GEOMETRIC_REPRESENTATION_CONTEXT ( 3 ) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT ( ( #%d ) ) GLOBAL_UNIT_ASSIGNED_CONTEXT ( ( #%d, #%d, #%d ) ) REPRESENTATION_CONTEXT ( ''NONE'', ''WORKASPACE'' ) );\n',object_index,8,5,6,7);object_index=object_index+1;
            fprintf(step_file,'\n');
        end

        function writeStepOpenShell(self,step_filestr)
            % write surface into step file
            %
            surf_num=length(self.surface_list);
            
            % check file name
            if length(step_filestr) > 4
                if ~strcmp(step_filestr((end-4):end),'.step')
                    step_filestr=[step_filestr,'.step'];
                end
            else
                step_filestr=[step_filestr,'.step'];
            end
            [~,step_filename,~]=fileparts(step_filestr);

            % write head
            step_file=fopen(step_filestr,'w');
            object_index=1;
            object_index=writeStepHead(self,step_file,object_index,step_filename);

            % write surface
            ADVANCED_FACE_index_list=zeros(1,surf_num);
            
            for surf_idx=1:surf_num
                surf=self.surface_list{surf_idx};

                [step_str,object_index,ADVANCED_FACE_index]=surf.getStepNode(object_index);
                fprintf(step_file,step_str);

                ADVANCED_FACE_index_list(surf_idx)=ADVANCED_FACE_index;
            end

            % generate OPEN_SHELL
            OPEN_SHELL_index=object_index;
            step_str=[num2str(object_index,'#%d'),' = OPEN_SHELL ',...
                '( ''NONE'', ',...
                '( ',num2str(ADVANCED_FACE_index_list(1:end-1),'#%d, '),' ',num2str(ADVANCED_FACE_index_list(end),'#%d'),' )',...
                ' );\n'];object_index=object_index+1;
            fprintf(step_file,step_str);
            fprintf(step_file,'\n');

            % write model
            SHELL_BASED_SURFACE_MODEL_index=object_index;
            step_str=[num2str(object_index,'#%d'),' = SHELL_BASED_SURFACE_MODEL ',...
                '( ''NONE'', ( ',...
                num2str(OPEN_SHELL_index,'#%d'),' )',...
                ' );\n'];object_index=object_index+1;
            fprintf(step_file,step_str);
            fprintf(step_file,'\n');

            % write end of step file
            writeStepEnd(self,step_file,object_index,step_filename,SHELL_BASED_SURFACE_MODEL_index);

            fclose(step_file);
            clear('step_file');
        end

        function writeStepEnd(self,step_file,object_index,step_filename,model_index)
            % write end of step file
            %

            % write product context 
            fprintf(step_file,'#%d = APPLICATION_CONTEXT ( ''configuration controlled 3d designs of mechanical parts and assemblies'' ) ;\n',object_index);object_index=object_index+1;
            fprintf(step_file,'#%d = MECHANICAL_CONTEXT ( ''NONE'', #%d, ''mechanical'' ) ;\n',object_index,object_index-1);object_index=object_index+1;
            fprintf(step_file,'#%d = PRODUCT ( ''%s'', ''%s'', '''', ( #%d ) ) ;\n',object_index,step_filename,step_filename,object_index-1);object_index=object_index+1;
            fprintf(step_file,'\n');
            fprintf(step_file,'#%d = APPLICATION_CONTEXT ( ''configuration controlled 3d designs of mechanical parts and assemblies'' ) ;\n',object_index);object_index=object_index+1;
            fprintf(step_file,'#%d = DESIGN_CONTEXT ( ''detailed design'', #%d, ''design'' ) ;\n',object_index,object_index-1);object_index=object_index+1;
            fprintf(step_file,'#%d = PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE ( ''NONE'', '''', #%d, .NOT_KNOWN. ) ;\n',object_index,object_index-3);object_index=object_index+1;
            fprintf(step_file,'#%d = PRODUCT_DEFINITION ( ''NONE'', '''', #%d, #%d ) ;\n',object_index,object_index-1,object_index-2);object_index=object_index+1;
            fprintf(step_file,'\n');
            fprintf(step_file,'#%d = PRODUCT_DEFINITION_SHAPE ( ''NONE'', ''NONE'',  #%d ) ;\n',object_index,object_index-1);object_index=object_index+1;
            fprintf(step_file,'#%d = MANIFOLD_SURFACE_SHAPE_REPRESENTATION ( ''test'', ( #%d, #%d ), #%d ) ;\n',object_index,model_index,4,9);object_index=object_index+1;
            fprintf(step_file,'#%d = SHAPE_DEFINITION_REPRESENTATION ( #%d, #%d ) ;\n',object_index,object_index-2,object_index-1);object_index=object_index+1;

            % write end
            fprintf(step_file,'\n');
            fprintf(step_file,'ENDSEC;\nEND-ISO-10303-21;\n');

        end
    end
end