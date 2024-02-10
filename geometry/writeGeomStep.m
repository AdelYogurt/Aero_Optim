function writeGeomStep(geom_list,step_filestr)
% write topology entity into step file
%
if nargin < 2, step_filestr=[];end
if isempty(step_filestr), step_filestr='geom.step';end

if ~iscell(geom_list), geom_list={geom_list};end

[~,step_filename,~]=fileparts(step_filestr);
obj_idx=1;geom_num=length(geom_list);

% write head of step file
[str_head,obj_idx]=writeStepHead(obj_idx,step_filename);

% write entity
str_entity='';WIRE=[];SHELL=[];
for geom_idx=1:geom_num
    geom=geom_list{geom_idx};

    % if isa(geom,'Edge')
    %     % GEOMETRIC_CURVE_SET
    %     [str_new,obj_idx,model_idx]=stepEdge(geom,obj_idx);
    %     WIRE=[WIRE,model_idx];
    % elseif isa(geom,'Wire')
    %     % EDGE_BASED_WIREFRAME_MODEL
    %     [str_new,obj_idx,model_idx]=stepLoop(geom,obj_idx);
    %     WIRE=[WIRE,model_idx];
    % elseif isa(geom,'Face')
    %     % 
    %     [str_new,obj_idx,model_idx]=stepFace(geom,obj_idx);
    %     SHELL=[SHELL,model_idx];
    % elseif isa(geom,'Shell')
    %     % SHELL_BASED_SURFACE_MODEL
    %     [str_new,obj_idx,model_idx]=stepShell(geom,obj_idx);
    %     SHELL=[SHELL,model_idx];
    % else
    %     error('writeGeomStep: unknown topology entity');
    % end

    if isa(geom,'Edge') || isa(geom,'Wire')
        % write wire
        if isa(geom,'Edge'), geom=Wire(geom.name,geom);end
        [str_new,obj_idx,model_idx]=writeLoop(geom,obj_idx);
        WIRE=[WIRE,model_idx];
    elseif isa(geom,'Face') || isa(geom,'Shell')
        % write shell
        if isa(geom,'Face'), geom=Shell(geom.name,geom);end
        [str_new,obj_idx,model_idx]=writeShell(geom,obj_idx);
        SHELL=[SHELL,model_idx];
    else
        error('writeGeomStep: unknown topology entity');
    end
    
    str_entity=[str_entity,str_new,'\n'];
end

% write model
if ~isempty(WIRE)
    WIRE_MODEL=obj_idx;
    format_out=['( ','#%d',repmat(', #%d',1,length(WIRE)-1),' )'];
    str_wire=[sprintf('#%d = EDGE_BASED_WIREFRAME_MODEL ( ''NONE'', ',WIRE_MODEL),sprintf(format_out,WIRE),sprintf(' ) ;\n')];
    obj_idx=obj_idx+1;
else
    WIRE_MODEL=[];
    str_wire='';
end

if ~isempty(SHELL)
    SHELL_MODEL=obj_idx;
    format_out=['( ','#%d',repmat(', #%d',1,length(SHELL)-1),' )'];
    str_shell=[sprintf('#%d = SHELL_BASED_SURFACE_MODEL ( ''NONE'', ',SHELL_MODEL),sprintf(format_out,SHELL),sprintf(' ) ;\n')];
    obj_idx=obj_idx+1;
else
    SHELL_MODEL=[];
    str_shell='';
end

% write end of step file
[str_end,obj_idx]=writeStepEnd(obj_idx,step_filename,WIRE_MODEL,SHELL_MODEL);

% write into file
str=[str_head,str_entity,str_wire,str_shell,str_end];
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

    function [str,obj_idx]=writeStepEnd(obj_idx,step_filename,WIRE_MODEL,SHELL_MODEL)
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
        if ~isempty(WIRE_MODEL)
            str=[str,sprintf('#%d = GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION ( ''geom'', ( #%d, #%d ), #%d ) ;\n',obj_idx,WIRE_MODEL,4,9)]; obj_idx=obj_idx+1;
        end
        if ~isempty(SHELL_MODEL)
            str=[str,sprintf('#%d = MANIFOLD_SURFACE_SHAPE_REPRESENTATION ( ''geom'', ( #%d, #%d ), #%d ) ;\n',obj_idx,SHELL_MODEL,4,9)]; obj_idx=obj_idx+1;
        end
        str=[str,sprintf('#%d = SHAPE_DEFINITION_REPRESENTATION ( #%d, #%d ) ;\n',obj_idx,obj_idx-2,obj_idx-1)]; obj_idx=obj_idx+1;

        % write end
        str=[str,newline];
        str=[str,sprintf('ENDSEC;\nEND-ISO-10303-21;\n')];
        str=[str,newline];
    end
end

%% topology step

function [str,obj_idx,VERTEX_POINT]=writeVertex(vtx,obj_idx)
% generate vertex step str
%
% notice:
% vertex should be point_num x 3 matrix
%

% write point
CARTESIAN_POINT=obj_idx;
[str_pnt,obj_idx]=writePoint(vtx,obj_idx);
Pnt_idx=CARTESIAN_POINT-1+(1:size(vtx,1));

VERTEX_POINT=obj_idx-1+(1:size(vtx,1));
out_data=[VERTEX_POINT;Pnt_idx];
str_vtx=sprintf('#%d = VERTEX_POINT ( ''NONE'', #%d ) ;\n',out_data);
obj_idx=obj_idx+size(vtx,1);

% generate all
str=[str_pnt,str_vtx];
end

function [str,obj_idx,ORIENTED_EDGE]=writeEdge(edg,obj_idx)
% generate edge step str
%

% load properties
name=edg.name;
crv=edg.curve;
topo=edg.topology;
vtx_list=edg.vertex_list;

% process properties
if isempty(name), name='NONE';end
vtx_list=vtx_list(topo.vertex,:);

% write curve
[str_crv,obj_idx]=writeCurve(crv,obj_idx);
CURVE=obj_idx-1;

% write vertex
[str_vtx,obj_idx,VERTEX_POINT]=writeVertex(vtx_list,obj_idx);

% write edge
EDGE_CURVE=obj_idx;
str_edg=sprintf('#%d = EDGE_CURVE ( ''NONE'', #%d, #%d, #%d, .T. ) ;\n',EDGE_CURVE,VERTEX_POINT(1),VERTEX_POINT(2),CURVE);
obj_idx=obj_idx+1;

% write oriented edg
ORIENTED_EDGE=obj_idx;
str_ori=sprintf('#%d = ORIENTED_EDGE ( ''%s'', *, *, #%d, .T. ) ;\n',ORIENTED_EDGE,name,EDGE_CURVE);
obj_idx=obj_idx+1;

% generate all
str=[str_crv,str_vtx,str_edg,str_ori];
end

function [str,obj_idx,EDGE_LOOP]=writeLoop(wir,obj_idx)
% generate loop step str
%

% load properties
name=wir.name;
edg_list=wir.edge_list;
topo=wir.topology;

% process properties
if isempty(name), name='NONE';end

% write edge
str_edg='';edg_num=length(edg_list);
ORIENTED_EDGE=zeros(1,edg_num);
for edg_idx=1:edg_num
    [str,obj_idx,ORIENTED_EDGE(edg_idx)]=writeEdge(edg_list(edg_idx),obj_idx);
    str_edg=[str_edg,str,'\n'];
end

% write loop
EDGE_LOOP=obj_idx;
format_out=['( ','#%d',repmat(', #%d',1,edg_num-1),' )'];
str_lop=[sprintf('#%d = EDGE_LOOP ( ''%s'', ',EDGE_LOOP,name),sprintf(format_out,ORIENTED_EDGE),' ) ;\n'];
obj_idx=obj_idx+1;

% generate all
str=[str_edg,str_lop];
end

function [str,obj_idx,ADVANCED_FACE]=writeFace(fce,obj_idx)
% generate face step str
%

% load properties
name=fce.name;
srf=fce.surface;
topo=fce.topology;
wire_list=fce.wire_list;

% process properties
if isempty(name), name='NONE';end

% write surface
[str_srf,obj_idx]=writeSurface(srf,obj_idx);
SURFACE=obj_idx-1;

% write edge loop
str_lop='';wir_num=length(wire_list);
EDGE_LOOP=zeros(1,wir_num);
for wir_idx=1:wir_num
    [str,obj_idx,EDGE_LOOP(wir_idx)]=writeLoop(wire_list(wir_idx),obj_idx);
    str_lop=[str_lop;str];
end

% generate FACE_OUTER_BOUND
FACE_OUTER_BOUND=obj_idx;
str_otr_bou=sprintf('#%d = FACE_OUTER_BOUND ( ''NONE'', #%d, .T. ) ;\n',FACE_OUTER_BOUND,EDGE_LOOP(1));
obj_idx=obj_idx+1;

% generate FACE_BOUND
if wir_num > 1
    FACE_BOUND=zeros(1,wir_num-1);
    str_bou='';
    for wir_idx=1:wir_num
        FACE_BOUND(wir_idx)=obj_idx;
        str_bou=[str_bou,sprintf('#%d = FACE_BOUND ( ''NONE'', #%d, .T. ) ;\n',FACE_BOUND(wir_idx),EDGE_LOOP(1+wir_idx))];
        obj_idx=obj_idx+1;
    end
else
    FACE_BOUND=[];
    str_bou='';
end

% generate ADVANCED_FACE
ADVANCED_FACE=obj_idx;
BOUND=[FACE_OUTER_BOUND,FACE_BOUND];
format_out=['( ','#%d',repmat(', #%d',1,wir_num-1),' )'];
str_fce=[sprintf('#%d = ADVANCED_FACE ( ''%s'', ',ADVANCED_FACE,name),sprintf(format_out,BOUND),sprintf(', #%d, .T. ) ;\n',SURFACE)];
obj_idx=obj_idx+1;

% generate all
str=[str_srf,str_lop,str_otr_bou,str_bou,str_fce];
end

function [str,obj_idx,SHELL]=writeShell(shl,obj_idx)
% write surface into step file
%

% load properties
name=shl.name;
fce_list=shl.face_list;
topo=shl.topology;
clsd=shl.closed;

% process properties
if isempty(name), name='NONE';end

% write face
str_fce='';fce_num=length(fce_list);
ADVANCED_FACE=zeros(1,fce_num);
for fce_idx=1:fce_num
    [str,obj_idx,ADVANCED_FACE(fce_idx)]=writeFace(fce_list(fce_idx),obj_idx);
    str_fce=[str_fce,str,'\n'];
end

% write shell
SHELL=obj_idx;
format_out=['( ','#%d',repmat(', #%d',1,fce_num-1),' )'];
if clsd
str_shell=[sprintf('#%d = CLOSED_SHELL ( ''%s'', ',SHELL,name),sprintf(format_out,ADVANCED_FACE),' ) ;\n'];
else
str_shell=[sprintf('#%d = OPEN_SHELL ( ''%s'', ',SHELL,name),sprintf(format_out,ADVANCED_FACE),' ) ;\n'];
end
obj_idx=obj_idx+1;

% generate all
str=[str_fce,str_shell];
end

%% geometry step

function [str,obj_idx]=writePoint(pnt,obj_idx)
% generate point step str
%
% notice:
% point should be point_num x 3 matrix
%
out_data=[obj_idx-1+(1:size(pnt,1));pnt'];
str=sprintf('#%d = CARTESIAN_POINT ( ''NONE'', ( %.16f, %.16f, %.16f ) ) ;\n',out_data);
obj_idx=obj_idx+size(pnt,1);
end

function [str,obj_idx]=writeCurve(crv,obj_idx)
% generate curve step str
%
if ~isa(crv,'Curve')
    error('stepCurve: only support BSpline curve');
end

% load properties
Poles=crv.Poles;
Degree=crv.Degree;
Mults=crv.Mults;
Knots=crv.Knots;
Weights=crv.Weights;

% write point
CARTESIAN_POINT=obj_idx;
[str,obj_idx]=writePoint(Poles,obj_idx);
Pnt_idx=CARTESIAN_POINT-1+(1:size(Poles,1));

if isempty(Weights)
    % unrational
    str=[str,num2str(obj_idx,'#%d ='),' B_SPLINE_CURVE_WITH_KNOTS',...
        ' ( ',...
        '''NONE''',', ',...
        num2str(Degree,'%d'),', \n',...
        '( ',num2str(Pnt_idx(1),'#%d'),num2str(Pnt_idx(2:end),', #%d'),' ),\n',...
        '.UNSPECIFIED., .F., .F.,\n',...
        '( ',num2str(Mults(1),'%d'),num2str(Mults(2:end),', %d'),' ),\n',...
        '( ',num2str(Knots(1),'%.16f'),num2str(Knots(2:end),', %.16f'),' ),\n',...
        '.UNSPECIFIED.',...
        ' ) ;\n']; obj_idx=obj_idx+1;
else
    % rational
    str=[str,num2str(obj_idx,'#%d ='),' (\n'];

    % B_SPLINE_CURVE
    str=[str,'B_SPLINE_CURVE ( ',...
        num2str(Degree,'%d'),', ',...
        '( ',num2str(Pnt_idx(1),'#%d'),num2str(Pnt_idx(2:end),', #%d'),' ),\n',...
        '.UNSPECIFIED., .F., .F. )\n'];

    % B_SPLINE_CURVE_WITH_KNOTS
    str=[str,'B_SPLINE_CURVE_WITH_KNOTS ( ( ',num2str(Mults(1),'%d'),num2str(Mults(2:end),', %d'),' ), ',...
        '( ',num2str(Knots(1),'%.16f'),num2str(Knots(2:end),', %.16f'),' ),\n',...
        '.UNSPECIFIED. )\n'];

    % RATIONAL_B_SPLINE_CURVE
    Weights=Weights(:)';
    str=[str,'RATIONAL_B_SPLINE_CURVE ( (',num2str(Weights(1),' %.16f'),num2str(Weights(2:end),', %.16f'),' ) )\n'];

    str=[str,') ;\n']; obj_idx=obj_idx+1;
end
str=sprintf(str);
end

function [str,obj_idx]=writeSurface(srf,obj_idx)
% generate surface step str
%
if ~isa(srf,'Surface')
    error('stepSurface: only support BSpline surface');
end

% load properties
Poles=srf.Poles;
UDegree=srf.UDegree;
VDegree=srf.VDegree;
UMults=srf.UMults;
VMults=srf.VMults;
UKnots=srf.UKnots;
VKnots=srf.VKnots;
Weights=srf.Weights;

% generate CARTESIAN_POINT
[v_num,u_num,~]=size(Poles);
CARTESIAN_POINT=obj_idx;
Poles=reshape(Poles,v_num*u_num,3);
[str,obj_idx]=writePoint(Poles,obj_idx);
Pnt_idx=CARTESIAN_POINT-1+(1:v_num*u_num);

if isempty(Weights)
    % unrational
    str=[str,num2str(obj_idx,'#%d ='),' B_SPLINE_SURFACE_WITH_KNOTS',...
        ' ( ',...
        '''NONE''',', ',...
        num2str(UDegree,'%d'),', ',num2str(VDegree,'%d'),', \n'];
    str=[str,'('];
    format_out=['( ','#%d',repmat(', #%d',1,v_num-1),' ),\n'];
    str=[str,sprintf(format_out,(Pnt_idx))];
    str=[str(1:end-2),'),\n'];
    str=[str,'.UNSPECIFIED.',', ','.F.',', ','.F.',', ','.F.',',\n',...
        '( ',num2str(UMults(1),'%d'),num2str(UMults(2:end),', %d'),' ),\n',...
        '( ',num2str(VMults(1),'%d'),num2str(VMults(2:end),', %d'),' ),\n',...
        '( ',num2str(UKnots(1),'%.16f'),num2str(UKnots(2:end),', %.16f'),' ),\n',...
        '( ',num2str(VKnots(1),'%.16f'),num2str(VKnots(2:end),', %.16f'),' ),\n',...
        '.UNSPECIFIED.',...
        ') ;\n']; obj_idx=obj_idx+1;
else
    % rational
    str=[str,num2str(obj_idx,'#%d ='),' (\n'];

    % B_SPLINE_SURFACE
    str=[str,'B_SPLINE_SURFACE ( ',...
        num2str(UDegree,'%d'),', ',num2str(VDegree,'%d'),', \n'];
    str=[str,'('];
    format_out=['( ','#%d',repmat(', #%d',1,v_num-1),' ),\n'];
    str=[str,sprintf(format_out,(Pnt_idx))];
    str=[str(1:end-2),'),\n','.UNSPECIFIED.',', ','.F.',', ','.F.',', ','.F.',' )\n'];

    % B_SPLINE_SURFACE_WITH_KNOTS
    str=[str,'B_SPLINE_SURFACE_WITH_KNOTS ( \n',...
        '( ',num2str(UMults(1),'%d'),num2str(UMults(2:end),', %d'),' ),\n',...
        '( ',num2str(VMults(1),'%d'),num2str(VMults(2:end),', %d'),' ),\n',...
        '( ',num2str(UKnots(1),'%.16f'),num2str(UKnots(2:end),', %.16f'),' ),\n',...
        '( ',num2str(VKnots(1),'%.16f'),num2str(VKnots(2:end),', %.16f'),' ),\n',...
        '.UNSPECIFIED. )\n'];

    % RATIONAL_B_SPLINE_SURFACE
    str=[str,'RATIONAL_B_SPLINE_SURFACE ( \n'];
    str=[str,'('];
    format_out=['( ','%.16f',repmat(', %.16f',1,v_num-1),' ),\n'];
    str=[str,sprintf(format_out,(Weights))];
    str=[str(1:end-2),' ) )\n'];

    str=[str,') ;\n']; obj_idx=obj_idx+1;
end
str=sprintf(str);
end
