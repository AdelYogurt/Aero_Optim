function writeGeomStep(geom_list,step_filestr)
% write step file
%
if nargin < 2, step_filestr=[];end
if isempty(step_filestr), step_filestr='geom.step';end

if ~iscell(geom_list), geom_list={geom_list};end

[~,step_filename,~]=fileparts(step_filestr);
obj_idx=1;geom_num=length(geom_list);

% write head of step file
[str_head,obj_idx]=writeStepHead(obj_idx,step_filename);

% write model
str_geom=[];SHELL_BASED_SURFACE_MODEL_list=zeros(1,geom_num);
for geom_idx=1:geom_num
    geom=geom_list{geom_idx};
    if isa(geom,'Shell')
        [str_new,SHELL_BASED_SURFACE_MODEL,obj_idx]=getOpenShell(geom,obj_idx);
    end
    str_geom=[str_geom,str_new,'\n'];
    SHELL_BASED_SURFACE_MODEL_list(geom_idx)=SHELL_BASED_SURFACE_MODEL;
end

% write end of step file
[str_end,obj_idx]=writeStepEnd(obj_idx,step_filename,SHELL_BASED_SURFACE_MODEL_list);

% write into file
str=[str_head,str_geom,str_end];
step_file=fopen(step_filestr,'w');
fprintf(step_file,str);
fclose(step_file);
clear('step_file');

    function [str,obj_idx]=writeStepHead(obj_idx,step_filename)
        % wite head of step
        %
        str=[];
        str=[str,sprintf('ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION (( ''STEP AP203'' ),''1'' );\nFILE_NAME (''%s'',''%s'',( '''' ),( '''' ),''Matlab step'',''Matlab'','''' );\nFILE_SCHEMA (( ''CONFIG_CONTROL_DESIGN'' ));\nENDSEC;\n',step_filename,date)];
        str=[str,newline];
        str=[str,sprintf('DATA;\n')];

        str=[str,sprintf('#%d = CARTESIAN_POINT ( ''NONE'',  ( 0.0, 0.0, 0.0 ) ) ;\n',obj_idx)]; obj_idx=obj_idx+1;
        str=[str,sprintf('#%d = DIRECTION ( ''NONE'',  ( 0.0, 0.0, 1.0 ) ) ;\n',obj_idx)]; obj_idx=obj_idx+1;
        str=[str,sprintf('#%d = DIRECTION ( ''NONE'',  ( 1.0, 0.0, 0.0 ) ) ;\n',obj_idx)]; obj_idx=obj_idx+1;
        str=[str,sprintf('#%d = AXIS2_PLACEMENT_3D ( ''NONE'', #%d, #%d, #%d ) ;\n',obj_idx,1,2,3)]; obj_idx=obj_idx+1;
        str=[str,newline];
        str=[str,sprintf('#%d =( LENGTH_UNIT ( ) NAMED_UNIT ( * ) SI_UNIT ( $., .METRE. ) );\n',obj_idx)]; obj_idx=obj_idx+1;
        str=[str,sprintf('#%d =( NAMED_UNIT ( * ) PLANE_ANGLE_UNIT ( ) SI_UNIT ( $, .RADIAN. ) );\n',obj_idx)]; obj_idx=obj_idx+1;
        str=[str,sprintf('#%d =( NAMED_UNIT ( * ) SI_UNIT ( $, .STERADIAN. ) SOLID_ANGLE_UNIT ( ) );\n',obj_idx)]; obj_idx=obj_idx+1;
        str=[str,sprintf('#%d = UNCERTAINTY_MEASURE_WITH_UNIT (LENGTH_MEASURE( 1.0E-05 ), #%d, ''distance_accuracy_value'', ''NONE'');\n',obj_idx,5)]; obj_idx=obj_idx+1;
        str=[str,sprintf('#%d =( GEOMETRIC_REPRESENTATION_CONTEXT ( 3 ) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT ( ( #%d ) ) GLOBAL_UNIT_ASSIGNED_CONTEXT ( ( #%d, #%d, #%d ) ) REPRESENTATION_CONTEXT ( ''NONE'', ''WORKASPACE'' ) );\n',obj_idx,8,5,6,7)]; obj_idx=obj_idx+1;
        str=[str,newline];
    end

    function [str,obj_idx]=writeStepEnd(obj_idx,step_filename,model_index)
        % write end of step file
        %

        % write product context
        str=[];
        str=[str,sprintf('#%d = APPLICATION_CONTEXT ( ''configuration controlled 3d designs of mechanical parts and assemblies'' ) ;\n',obj_idx)]; obj_idx=obj_idx+1;
        str=[str,sprintf('#%d = MECHANICAL_CONTEXT ( ''NONE'', #%d, ''mechanical'' ) ;\n',obj_idx,obj_idx-1)]; obj_idx=obj_idx+1;
        str=[str,sprintf('#%d = PRODUCT ( ''%s'', ''%s'', '''', ( #%d ) ) ;\n',obj_idx,step_filename,step_filename,obj_idx-1)]; obj_idx=obj_idx+1;
        str=[str,newline];
        str=[str,sprintf('#%d = APPLICATION_CONTEXT ( ''configuration controlled 3d designs of mechanical parts and assemblies'' ) ;\n',obj_idx)]; obj_idx=obj_idx+1;
        str=[str,sprintf('#%d = DESIGN_CONTEXT ( ''detailed design'', #%d, ''design'' ) ;\n',obj_idx,obj_idx-1)];obj_idx=obj_idx+1;
        str=[str,sprintf('#%d = PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE ( ''NONE'', '''', #%d, .NOT_KNOWN. ) ;\n',obj_idx,obj_idx-3)]; obj_idx=obj_idx+1;
        str=[str,sprintf('#%d = PRODUCT_DEFINITION ( ''NONE'', '''', #%d, #%d ) ;\n',obj_idx,obj_idx-1,obj_idx-2)]; obj_idx=obj_idx+1;
        str=[str,newline];
        str=[str,sprintf('#%d = PRODUCT_DEFINITION_SHAPE ( ''NONE'', ''NONE'',  #%d ) ;\n',obj_idx,obj_idx-1)]; obj_idx=obj_idx+1;
        str=[str,sprintf('#%d = MANIFOLD_SURFACE_SHAPE_REPRESENTATION ( ''test'', ( #%d, #%d ), #%d ) ;\n',obj_idx,model_index,4,9)]; obj_idx=obj_idx+1;
        str=[str,sprintf('#%d = SHAPE_DEFINITION_REPRESENTATION ( #%d, #%d ) ;\n',obj_idx,obj_idx-2,obj_idx-1)]; obj_idx=obj_idx+1;

        % write end
        str=[str,newline];
        str=[str,sprintf('ENDSEC;\nEND-ISO-10303-21;\n')];
        str=[str,newline];
    end
end

function [str,SHELL_BASED_SURFACE_MODEL,obj_idx]=getOpenShell(shell,obj_idx)
% write surface into step file
%
surf_num=length(shell.face_list);
if surf_num == 0
    str=[];SHELL_BASED_SURFACE_MODEL=[];
    return;
end

% write face
str_face=[];ADVANCED_FACE_list=zeros(1,surf_num);
for surf_idx=1:surf_num
    fce=shell.face_list{surf_idx};
    if ~isa(fce,'FaceCST')
        if all(fce.Weights == 1)
            [str_new,ADVANCED_FACE,obj_idx]=bsplineFaceStep(fce,obj_idx);
        else
            [str_new,ADVANCED_FACE,obj_idx]=NURBSFaceStep(fce,obj_idx);
        end
    else
        fce=fce.getNURBS();
        [str_new,ADVANCED_FACE,obj_idx]=bsplineFaceStep(fce,obj_idx);
    end
    str_face=[str_face,str_new,'\n'];
    ADVANCED_FACE_list(surf_idx)=ADVANCED_FACE;
end

% write OPEN_SHELL
OPEN_SHELL=obj_idx;
str_shell=[num2str(obj_idx,'#%d ='),' OPEN_SHELL',...
    ' ( ',...
    '''NONE''',', ',...
    '( ',num2str(ADVANCED_FACE_list(1),'#%d'),num2str(ADVANCED_FACE_list(2:end),', #%d'),' )',...
    ' ) ;\n']; obj_idx=obj_idx+1;

% write model
SHELL_BASED_SURFACE_MODEL=obj_idx;
str_model=[num2str(obj_idx,'#%d ='),' SHELL_BASED_SURFACE_MODEL',...
    ' ( ',...
    '''NONE''',', ',...
    '( ',num2str(OPEN_SHELL,'#%d'),' )',...
    ' ) ;\n']; obj_idx=obj_idx+1;

str=sprintf([str_face,str_shell,str_model]);

end

function [str,ADVANCED_FACE,obj_idx]=bsplineFaceStep(fce,obj_idx)
% write BSpline into step file
%
if nargin < 2,obj_idx=1;end
name=fce.name;
if isempty(name),name='NONE';end

[v_num,u_num,~]=size(fce.Poles);

% generate CARTESIAN_POINT
CARTESIAN_POINT=obj_idx;
Poles=fce.Poles;
Poles=reshape(Poles,v_num*u_num,3);
str_point=stepPoint(obj_idx,Poles); obj_idx=obj_idx+v_num*u_num;

% generate B_SPLINE_CURVE_WITH_KNOTS
B_SPLINE_CURVE_WITH_KNOTS=obj_idx;
str_curve=[];
str_curve=[str_curve,stepCurve(obj_idx,((0:1:(u_num-1))*v_num+1) +CARTESIAN_POINT-1,...
    fce.UDegree,fce.UMults,fce.UKnots)]; obj_idx=obj_idx+1;
str_curve=[str_curve,stepCurve(obj_idx,((v_num*u_num-v_num+1):1:(v_num*u_num)) +CARTESIAN_POINT-1,...
    fce.VDegree,fce.VMults,fce.VKnots)]; obj_idx=obj_idx+1;
str_curve=[str_curve,stepCurve(obj_idx,((u_num:-1:1)*v_num) +CARTESIAN_POINT-1,...
    fce.UDegree,fliplr(fce.UMults),fliplr(1-fce.UKnots))]; obj_idx=obj_idx+1;
str_curve=[str_curve,stepCurve(obj_idx,(v_num:-1:1) +CARTESIAN_POINT-1,...
    fce.VDegree,fliplr(fce.VMults),fliplr(1-fce.VKnots))]; obj_idx=obj_idx+1;

% generate VERTEX_POINT
VERTEX_POINT=obj_idx;
str_vertex=[];
str_vertex=[str_vertex,stepVertex(obj_idx,1 +CARTESIAN_POINT-1)]; obj_idx=obj_idx+1;
str_vertex=[str_vertex,stepVertex(obj_idx,v_num*u_num-v_num+1 +CARTESIAN_POINT-1)]; obj_idx=obj_idx+1;
str_vertex=[str_vertex,stepVertex(obj_idx,v_num*u_num +CARTESIAN_POINT-1)]; obj_idx=obj_idx+1;
str_vertex=[str_vertex,stepVertex(obj_idx,v_num +CARTESIAN_POINT-1)]; obj_idx=obj_idx+1;

% generate EDGE_CURVE
EDGE_CURVE=obj_idx;
str_edge=[];
str_edge=[str_edge,stepEdge(obj_idx,VERTEX_POINT,VERTEX_POINT+1,B_SPLINE_CURVE_WITH_KNOTS)]; obj_idx=obj_idx+1;
str_edge=[str_edge,stepEdge(obj_idx,VERTEX_POINT+1,VERTEX_POINT+2,B_SPLINE_CURVE_WITH_KNOTS+1)]; obj_idx=obj_idx+1;
str_edge=[str_edge,stepEdge(obj_idx,VERTEX_POINT+2,VERTEX_POINT+3,B_SPLINE_CURVE_WITH_KNOTS+2)]; obj_idx=obj_idx+1;
str_edge=[str_edge,stepEdge(obj_idx,VERTEX_POINT+3,VERTEX_POINT,B_SPLINE_CURVE_WITH_KNOTS+3)]; obj_idx=obj_idx+1;

% generate ORIENTED_EDGE
ORIENTED_EDGE=obj_idx;
str_oriented=[];
str_oriented=[str_oriented,stepOriented(obj_idx,EDGE_CURVE)]; obj_idx=obj_idx+1;
str_oriented=[str_oriented,stepOriented(obj_idx,EDGE_CURVE+1)]; obj_idx=obj_idx+1;
str_oriented=[str_oriented,stepOriented(obj_idx,EDGE_CURVE+2)]; obj_idx=obj_idx+1;
str_oriented=[str_oriented,stepOriented(obj_idx,EDGE_CURVE+3)]; obj_idx=obj_idx+1;

% generate EDGE_LOOP
EDGE_LOOP=obj_idx;
str_loop=[num2str(obj_idx,'#%d ='),' EDGE_LOOP',...
    ' ( ',...
    '''NONE''',', ',...
    '( ',num2str(ORIENTED_EDGE,'#%d'),', ',...
    num2str(ORIENTED_EDGE+1,'#%d'),', ',...
    num2str(ORIENTED_EDGE+2,'#%d'),', ',...
    num2str(ORIENTED_EDGE+3,'#%d'),') ',...
    ' ) ;\n']; obj_idx=obj_idx+1;

% generate FACE_OUTER_BOUND
FACE_OUTER_BOUND=obj_idx;
str_outer=[num2str(obj_idx,'#%d ='),' FACE_OUTER_BOUND',...
    ' ( ',...
    '''NONE''',', ',...
    num2str(EDGE_LOOP,'#%d'),', ',...
    '.T.',...
    ' ) ;\n']; obj_idx=obj_idx+1;

% generate B_SPLINE_SURFACE_WITH_KNOTS
B_SPLINE_SURFACE_WITH_KNOTS=obj_idx;
point_index=CARTESIAN_POINT-1+reshape((1:v_num*u_num),v_num,u_num);
str_surface=stepSurface(obj_idx,name,point_index,fce.UDegree,fce.VDegree,fce.UMults,fce.VMults,fce.UKnots,fce.VKnots); obj_idx=obj_idx+1;

% generate ADVANCED_FACE
ADVANCED_FACE=obj_idx;
str_face=[num2str(obj_idx,'#%d ='),' ADVANCED_FACE',...
    ' ( ',...
    '''',name,'''',', ',...
    '( ',num2str(FACE_OUTER_BOUND,'#%d'),' )',', '...
    num2str(B_SPLINE_SURFACE_WITH_KNOTS,'#%d'),', ',...
    '.T.',...
    ' ) ;\n']; obj_idx=obj_idx+1;

% generate all
str=sprintf([str_point,'\n',str_curve,'\n',str_vertex,'\n',...
    str_edge,'\n',str_oriented,'\n',str_loop,'\n',str_outer,'\n',...
    str_surface,'\n',str_face]);
end

function [str,ADVANCED_FACE,obj_idx]=NURBSFaceStep(fce,obj_idx)
% write BSpline into step file
%
if nargin < 2,obj_idx=1;end
name=fce.name;
if isempty(name),name='NONE';end

u_num=size(fce.Poles,2);
v_num=size(fce.Poles,1);

% generate CARTESIAN_POINT
CARTESIAN_POINT=obj_idx;
Poles=fce.Poles;
Poles=reshape(Poles,v_num*u_num,3);
str_point=stepPoint(obj_idx,Poles);obj_idx=obj_idx+v_num*u_num;

% generate B_SPLINE_CURVE_WITH_KNOTS
B_SPLINE_CURVE_WITH_KNOTS=obj_idx;
str_curve=[];
str_curve=[str_curve,stepCurveRational(obj_idx,((0:1:(u_num-1))*v_num+1) +CARTESIAN_POINT-1,...
    fce.UDegree,fce.UMults,fce.UKnots,fce.Weights(1,:))]; obj_idx=obj_idx+1;
str_curve=[str_curve,stepCurveRational(obj_idx,((v_num*u_num-v_num+1):1:(v_num*u_num)) +CARTESIAN_POINT-1,...
    fce.VDegree,fce.VMults,fce.VKnots,fce.Weights(:,end))]; obj_idx=obj_idx+1;
str_curve=[str_curve,stepCurveRational(obj_idx,((u_num:-1:1)*v_num) +CARTESIAN_POINT-1,...
    fce.UDegree,fliplr(fce.UMults),fliplr(1-fce.UKnots),fce.Weights(end,end:-1:1))]; obj_idx=obj_idx+1;
str_curve=[str_curve,stepCurveRational(obj_idx,(v_num:-1:1) +CARTESIAN_POINT-1,...
    fce.VDegree,fliplr(fce.VMults),fliplr(1-fce.VKnots),fce.Weights(end:-1:1,1))]; obj_idx=obj_idx+1;

% generate VERTEX_POINT
VERTEX_POINT=obj_idx;
str_vertex=[];
str_vertex=[str_vertex,stepVertex(obj_idx,1 +CARTESIAN_POINT-1)]; obj_idx=obj_idx+1;
str_vertex=[str_vertex,stepVertex(obj_idx,v_num*u_num-v_num+1 +CARTESIAN_POINT-1)]; obj_idx=obj_idx+1;
str_vertex=[str_vertex,stepVertex(obj_idx,v_num*u_num +CARTESIAN_POINT-1)]; obj_idx=obj_idx+1;
str_vertex=[str_vertex,stepVertex(obj_idx,v_num +CARTESIAN_POINT-1)]; obj_idx=obj_idx+1;

% generate EDGE_CURVE
EDGE_CURVE=obj_idx;
str_edge=[];
str_edge=[str_edge,stepEdge(obj_idx,VERTEX_POINT,VERTEX_POINT+1,B_SPLINE_CURVE_WITH_KNOTS)]; obj_idx=obj_idx+1;
str_edge=[str_edge,stepEdge(obj_idx,VERTEX_POINT+1,VERTEX_POINT+2,B_SPLINE_CURVE_WITH_KNOTS+1)]; obj_idx=obj_idx+1;
str_edge=[str_edge,stepEdge(obj_idx,VERTEX_POINT+2,VERTEX_POINT+3,B_SPLINE_CURVE_WITH_KNOTS+2)]; obj_idx=obj_idx+1;
str_edge=[str_edge,stepEdge(obj_idx,VERTEX_POINT+3,VERTEX_POINT,B_SPLINE_CURVE_WITH_KNOTS+3)]; obj_idx=obj_idx+1;

% generate ORIENTED_EDGE
ORIENTED_EDGE=obj_idx;
str_oriented=[];
str_oriented=[str_oriented,stepOriented(obj_idx,EDGE_CURVE)]; obj_idx=obj_idx+1;
str_oriented=[str_oriented,stepOriented(obj_idx,EDGE_CURVE+1)]; obj_idx=obj_idx+1;
str_oriented=[str_oriented,stepOriented(obj_idx,EDGE_CURVE+2)]; obj_idx=obj_idx+1;
str_oriented=[str_oriented,stepOriented(obj_idx,EDGE_CURVE+3)]; obj_idx=obj_idx+1;

% generate EDGE_LOOP
EDGE_LOOP=obj_idx;
str_loop=[num2str(obj_idx,'#%d ='),' EDGE_LOOP',...
    ' ( ',...
    '''NONE''',', ',...
    '( ',num2str(ORIENTED_EDGE,'#%d'),', ',...
    num2str(ORIENTED_EDGE+1,'#%d'),', ',...
    num2str(ORIENTED_EDGE+2,'#%d'),', ',...
    num2str(ORIENTED_EDGE+3,'#%d'),') ',...
    ' ) ;\n']; obj_idx=obj_idx+1;

% generate FACE_OUTER_BOUND
FACE_OUTER_BOUND=obj_idx;
str_outer=[num2str(obj_idx,'#%d ='),' FACE_OUTER_BOUND',...
    ' ( ',...
    '''NONE''',', ',...
    num2str(EDGE_LOOP,'#%d'),', ',...
    '.T.',...
    ' ) ;\n']; obj_idx=obj_idx+1;

% generate B_SPLINE_SURFACE_WITH_KNOTS
B_SPLINE_SURFACE_WITH_KNOTS=obj_idx;
point_index=CARTESIAN_POINT-1+reshape((1:v_num*u_num),v_num,u_num);
str_surface=stepSurfaceRational(obj_idx,name,point_index,fce.UDegree,fce.VDegree,fce.UMults,fce.VMults,fce.UKnots,fce.VKnots,fce.Weights); obj_idx=obj_idx+1;

% generate ADVANCED_FACE
ADVANCED_FACE=obj_idx;
str_face=[num2str(obj_idx,'#%d ='),' ADVANCED_FACE',...
    ' ( ',...
    '''',name,'''',', ',...
    '( ',num2str(FACE_OUTER_BOUND,'#%d'),' )',', '...
    num2str(B_SPLINE_SURFACE_WITH_KNOTS,'#%d'),', ',...
    '.T.',...
    ' ) ;\n']; obj_idx=obj_idx+1;

% generate all
str=sprintf([str_point,'\n',str_curve,'\n',str_vertex,'\n',...
    str_edge,'\n',str_oriented,'\n',str_loop,'\n',str_outer,'\n',...
    str_surface,'\n',str_face]);
end

function str_point=stepPoint(obj_idx,point)
% notice:
% point should be point_nunm x 3 matrix
%
obj_idx=obj_idx-1+(1:size(point,1));
data=[obj_idx;point'];
str_point=sprintf('#%d = CARTESIAN_POINT ( ''NONE'', ( %.16f, %.16f, %.16f ) ) ;\n',data);
end

function str_curve=stepCurve(obj_idx,point_index,Degree,Mults,Knots)
str_curve=[num2str(obj_idx,'#%d ='),' B_SPLINE_CURVE_WITH_KNOTS',...
    ' ( ',...
    '''NONE''',', ',...
    num2str(Degree,'%d'),', \n',...
    '( ',num2str(point_index(1),'#%d'),num2str(point_index(2:end),', #%d'),' ),\n',...
    '.UNSPECIFIED., .F., .F.,\n',...
    '( ',num2str(Mults(1),'%d'),num2str(Mults(2:end),', %d'),' ),\n',...
    '( ',num2str(Knots(1),'%.16f'),num2str(Knots(2:end),', %.16f'),' ),\n',...
    '.UNSPECIFIED.',...
    ' ) ;\n'];
end

function str_curve=stepCurveRational(obj_idx,point_index,Degree,Mults,Knots,Weights)
str_curve=[num2str(obj_idx,'#%d ='),' (\n'];Weights=Weights(:)';

% B_SPLINE_CURVE
str_curve=[str_curve,'B_SPLINE_CURVE ( ',...
    num2str(Degree,'%d'),', ',...
    '( ',num2str(point_index(1),'#%d'),num2str(point_index(2:end),', #%d'),' ),\n',...
    '.UNSPECIFIED., .F., .F. )\n'];

% B_SPLINE_CURVE_WITH_KNOTS
str_curve=[str_curve,'B_SPLINE_CURVE_WITH_KNOTS ( ( ',num2str(Mults(1),'%d'),num2str(Mults(2:end),', %d'),' ), ',...
    '( ',num2str(Knots(1),'%.16f'),num2str(Knots(2:end),', %.16f'),' ),\n',...
    '.UNSPECIFIED. )\n'];

% RATIONAL_B_SPLINE_CURVE
str_curve=[str_curve,'RATIONAL_B_SPLINE_CURVE ( (',num2str(Weights(1),' %.16f'),num2str(Weights(2:end),', %.16f'),' ) )\n'];

str_curve=[str_curve,') ;\n'];
end

function str_vertex=stepVertex(obj_idx,point_index)
str_vertex=[num2str(obj_idx,'#%d ='),' VERTEX_POINT',...
    ' ( ',...
    '''NONE''',', ',...
    num2str(point_index,'#%d'),...
    ' ) ;\n'];
end

function str_edge=stepEdge(obj_idx,start_vertex,end_vertex,edge_index)
str_edge=[num2str(obj_idx,'#%d ='),' EDGE_CURVE',...
    ' ( ',...
    '''NONE''',', ',...
    num2str(start_vertex,'#%d'),', ',...
    num2str(end_vertex,'#%d'),', ',...
    num2str(edge_index,'#%d'),', ',...
    '.T.',...
    ' ) ;\n'];

end

function str_oriented=stepOriented(obj_idx,edge_index)
str_oriented=[num2str(obj_idx,'#%d ='),' ORIENTED_EDGE',...
    ' ( ',...
    '''NONE''',', ',...
    '*',', ','*',', ',...
    num2str(edge_index,'#%d'),', ',...
    '.T.',...
    ' ) ;\n'];
end

function str_surface=stepSurface(obj_idx,name,point_index,UDegree,VDegree,UMults,VMults,UKnots,VKnots)
[v_num,~]=size(point_index);

str_surface=[num2str(obj_idx,'#%d ='),' B_SPLINE_SURFACE_WITH_KNOTS',...
    ' ( ',...
    '''',name,'''',', ',...
    num2str(UDegree,'%d'),', ',num2str(VDegree,'%d'),', \n'];
str_surface=[str_surface,'('];
out_format=['( ','#%d',repmat(', #%d',1,v_num-1),' ),\n'];
str_surface=[str_surface,sprintf(out_format,(point_index))];
str_surface=[str_surface(1:end-2),'),\n'];
str_surface=[str_surface,'.UNSPECIFIED.',', ','.F.',', ','.F.',', ','.F.',',\n',...
    '( ',num2str(UMults(1),'%d'),num2str(UMults(2:end),', %d'),' ),\n',...
    '( ',num2str(VMults(1),'%d'),num2str(VMults(2:end),', %d'),' ),\n',...
    '( ',num2str(UKnots(1),'%.16f'),num2str(UKnots(2:end),', %.16f'),' ),\n',...
    '( ',num2str(VKnots(1),'%.16f'),num2str(VKnots(2:end),', %.16f'),' ),\n',...
    '.UNSPECIFIED.',...
    ') ;\n'];
end

function str_surface=stepSurfaceRational(obj_idx,name,point_index,UDegree,VDegree,UMults,VMults,UKnots,VKnots,Weights)
[v_num,~]=size(point_index);
str_surface=[num2str(obj_idx,'#%d ='),' (\n'];

% B_SPLINE_SURFACE
str_surface=[str_surface,'B_SPLINE_SURFACE ( ',...
    '''',name,'''',', ',...
    num2str(UDegree,'%d'),', ',num2str(VDegree,'%d'),', \n'];
str_surface=[str_surface,'('];
out_format=['( ','#%d',repmat(', #%d',1,v_num-1),' ),\n'];
str_surface=[str_surface,sprintf(out_format,(point_index))];
str_surface=[str_surface(1:end-2),'),\n','.UNSPECIFIED.',', ','.F.',', ','.F.',', ','.F.',' )\n'];

% B_SPLINE_SURFACE_WITH_KNOTS
str_surface=[str_surface,'B_SPLINE_SURFACE_WITH_KNOTS ( \n',...
    '( ',num2str(UMults(1),'%d'),num2str(UMults(2:end),', %d'),' ),\n',...
    '( ',num2str(VMults(1),'%d'),num2str(VMults(2:end),', %d'),' ),\n',...
    '( ',num2str(UKnots(1),'%.16f'),num2str(UKnots(2:end),', %.16f'),' ),\n',...
    '( ',num2str(VKnots(1),'%.16f'),num2str(VKnots(2:end),', %.16f'),' ),\n',...
    '.UNSPECIFIED. )\n'];

% RATIONAL_B_SPLINE_SURFACE
str_surface=[str_surface,'RATIONAL_B_SPLINE_SURFACE ( \n'];
str_surface=[str_surface,'('];
out_format=['( ','%.16f',repmat(', %.16f',1,v_num-1),' ),\n'];
str_surface=[str_surface,sprintf(out_format,(Weights))];
str_surface=[str_surface(1:end-2),' ) )\n'];

str_surface=[str_surface,') ;\n'];
end
