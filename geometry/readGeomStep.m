function body_list=readGeomStep(step_filestr)
% read step file into
%

step_file=fopen(step_filestr,'r');

% read all data from step
Data={};Idx=[];
while ~feof(step_file)
    str=fgetl(step_file);

    % start read data
    if ~isempty(str) && str(1) == '#'
        str=strrep(str,' ','');
        while str(end) ~= ';'
            str_new=fgetl(step_file);
            str_new=strrep(str_new,' ','');
            str=[str,str_new];
        end

        data=cell(1,2);
        str_list=split(str,'=');
        Idx=[Idx;str2double(str_list{1}(2:end))];
        str_list=split(str_list{2},'(');
        data{1}=str_list{1};
        str_data=join(str_list(2:end),'(');
        data{2}=str_data{1}(1:end-2);
        Data=[Data;data];
    end
end

fclose(step_file);
clear('step_file');

[Idx,sort_idx]=sort(Idx);
Data=Data(sort_idx,:);

% add addition index
Data=[cell(Idx(1)-1,2);Data];

% create model
Idx_body=[find(strcmp(Data(:,1),'MANIFOLD_SOLID_BREP')),...
    find(strcmp(Data(:,1),'SHELL_BASED_SURFACE_MODEL'))];
body_list=cell(1,length(Idx_body));
for body_idx=1:length(Idx_body)
    idx_body=Idx_body(body_idx);
    body_type=Data{idx_body,1};

    switch body_type
        case 'SHELL_BASED_SURFACE_MODEL'
            str_body=Data{idx_body,2};
            str_list=split(str_body,',');
            shell_name=str_list{1};
            str_shell_list=join(str_list(2:end),',');
            Idx_shell=str2double(strsplit(strrep(str_shell_list{1}(2:end-1),'#',''),','));
            face_list=[];
            for shell_idx=1:length(Idx_shell)
                idx_shell=Idx_shell(shell_idx);
                % CLOSED/OPEN SHELL
                str_shell=Data{idx_shell,2};
                str_list=split(str_shell,',');str_shell=join(str_list(2:end),',');
                Idx_face=str2double(strsplit(strrep(str_shell{1}(2:end-1),'#',''),','));

                face_shell=cell(1,length(Idx_face));
                for face_idx=1:length(Idx_face)
                    idx_face=Idx_face(face_idx);

                    str_face=Data{idx_face,2};
                    str_list=split(str_face,',');
                    str_list=[str_list(1),join(str_list(2:end-2),','),str_list(end-1),str_list(end)];

                    % read surface and create BSpline surface
                    idx_surf=str2double(str_list{3}(2:end));
                    surf_type=Data{idx_surf,1};
                    switch surf_type
                        case 'B_SPLINE_SURFACE_WITH_KNOTS'
                            surf=getFaceBSpline(Data,idx_face,idx_surf);
                        case 'PLANE'
                            surf=getPLANE(Data,idx_face,idx_surf);
                        otherwise
                            continue;
                            % error('readGeomStep: unknown surface type');
                    end
                    face_shell{face_idx}=surf;
                end
                face_list=[face_list,face_shell];
            end

            body=Shell(shell_name,face_list);

        case 'MANIFOLD_SOLID_BREP'
            str_body=Data{idx_body,2};
            str_list=split(str_body,',');
            solid_name=str_list{1};
            Idx_solid=str2double(str_list{2}(2:end));

            body=Solid(solid_name);
        otherwise
            error('readGeomStep: unknown body type');
    end
    body_list{body_idx}=body;
end

end

function surf=getFaceBSpline(Data,idx_face,idx_surf)
% generate BSpline surface by BSpline surface string data
%

str_surf=Data{idx_surf,2};
% main properties
idx=find(str_surf == ',',1);
name=str_surf(1:idx-1);str_surf(1:idx)=[];
idx=find(str_surf == ',',1);
u_degree=str2double(str_surf(1:idx-1));str_surf(1:idx)=[];
idx=find(str_surf == ',',1);
v_degree=str2double(str_surf(1:idx-1));str_surf(1:idx)=[];

% read control point matrix
idx=2;par_sum=1;
while idx < length(str_surf)
    if str_surf(idx) == '(',par_sum=par_sum+1;elseif str_surf(idx) == ')',par_sum=par_sum-1;end
    if par_sum==0,break;end
    idx=idx+1;
end
str_control=strrep(str_surf(2:idx-1),'#','');
str_control=strsplit(str_control(2:end-1),'),(');
Idx_control=zeros(length(str_control),sum(str_control{1} == ',')+1);
for u_idx=1:length(str_control)
    Idx_control(u_idx,:)=str2double(strsplit(str_control{u_idx},','));
end
str_surf(1:idx+1)=[];

% properties
idx=find(str_surf == ',',1);
surface_form=(str_surf(1:idx-1));str_surf(1:idx)=[];
idx=find(str_surf == ',',1);
u_closed=(str_surf(1:idx-1));str_surf(1:idx)=[];
idx=find(str_surf == ',',1);
v_closed=(str_surf(1:idx-1));str_surf(1:idx)=[];
idx=find(str_surf == ',',1);
self_intersect=(str_surf(1:idx-1));str_surf(1:idx)=[];

% knot_multi
idx=find(str_surf == ')',1);
u_knot_multi=str2double(split(str_surf(2:idx-1),',',2));str_surf(1:idx+1)=[];
idx=find(str_surf == ')',1);
v_knot_multi=str2double(split(str_surf(2:idx-1),',',2));str_surf(1:idx+1)=[];

% knot_list
idx=find(str_surf == ')',1);
u_knot_list=str2double(split(str_surf(2:idx-1),',',2));str_surf(1:idx+1)=[];
idx=find(str_surf == ')',1);
v_knot_list=str2double(split(str_surf(2:idx-1),',',2));str_surf(1:idx+1)=[];

% knot_spec
knot_spec=str_surf;

% load control point data
[u_ctrl_num,v_ctrl_num]=size(Idx_control);
ctrl_X=zeros(v_ctrl_num,u_ctrl_num);
ctrl_Y=zeros(v_ctrl_num,u_ctrl_num);
ctrl_Z=zeros(v_ctrl_num,u_ctrl_num);
for u_idx=1:u_ctrl_num
    for v_idx=1:v_ctrl_num
        idx_point=Idx_control(u_idx,v_ctrl_num-v_idx+1);
        control_point=getPoint(Data,idx_point);
        ctrl_X(v_idx,u_idx)=control_point(1);
        ctrl_Y(v_idx,u_idx)=control_point(2);
        ctrl_Z(v_idx,u_idx)=control_point(3);
    end
end

surf=FaceBSpline(name,ctrl_X,ctrl_Y,ctrl_Z,...
    [],u_degree,v_degree,...
    u_knot_multi,v_knot_multi,u_knot_list,v_knot_list);

surf.surface_form=surface_form;
surf.u_closed=u_closed;
surf.v_closed=v_closed;
surf.self_intersect=self_intersect;
surf.knot_spec=knot_spec;
end

function surf=getPLANE(Data,idx_face,surf_idx)
% generate BSpline surface by plane string data
%

% ADVANCED_FACE
str_face=Data{idx_face,2};
str_list=split(str_face,',');
str_bound=join(str_list(2:end-2),',');str_bound=strrep(str_bound{1}(2:end-1),'#','');
Idx_bound=str2double(strsplit(str_bound,','));

% BOUND
edge_list=getBound(Data,Idx_bound);

str_surf=Data{surf_idx,2};
% main properties
idx=find(str_surf == ',',1);
name=str_surf(1:idx-1);str_surf(1:idx)=[];
idx_axis=str2double(str_surf(2:end));

% axis
str_axis=Data{idx_axis,2};

% process edge_list
if length(edge_list) > 4
    surf=[];
    return;
end
switch length(edge_list)
    case 2
        edge_1=edge_list{1};
        edge_3=edge_list{2};
        v_degree=max(edge_1.degree,edge_3.degree);
        u_degree=1;

        edge_1.changeDegree(v_degree);control_1=[edge_1.ctrl_X,edge_1.ctrl_Y,edge_1.ctrl_Z];
        edge_3.changeDegree(v_degree);control_3=[edge_3.ctrl_X,edge_3.ctrl_Y,edge_3.ctrl_Z];

        control_2=repmat(control_1(end,:),2,1);
        control_4=repmat(control_1(1,:),2,1);
    case 3
        edge_1=edge_list{1};
        edge_2=edge_list{2};
        edge_3=edge_list{3};
        u_degree=max(edge_2.degree);
        v_degree=max(edge_1.degree,edge_3.degree);

        edge_1.changeDegree(v_degree);control_1=[edge_1.ctrl_X,edge_1.ctrl_Y,edge_1.ctrl_Z];
        edge_2.changeDegree(u_degree);control_2=[edge_2.ctrl_X,edge_2.ctrl_Y,edge_2.ctrl_Z];
        edge_3.changeDegree(v_degree);control_3=[edge_3.ctrl_X,edge_3.ctrl_Y,edge_3.ctrl_Z];

        control_4=repmat(control_1(1,:),u_num,1);
    case 4
        edge_1=edge_list{1};
        edge_2=edge_list{2};
        edge_3=edge_list{3};
        edge_4=edge_list{4};
        u_degree=max(edge_2.degree,edge_4.degree);
        v_degree=max(edge_1.degree,edge_3.degree);

        edge_1.changeDegree(v_degree);control_1=[edge_1.ctrl_X,edge_1.ctrl_Y,edge_1.ctrl_Z];
        edge_2.changeDegree(u_degree);control_2=[edge_2.ctrl_X,edge_2.ctrl_Y,edge_2.ctrl_Z];
        edge_3.changeDegree(v_degree);control_3=[edge_3.ctrl_X,edge_3.ctrl_Y,edge_3.ctrl_Z];
        edge_4.changeDegree(u_degree);control_4=[edge_4.ctrl_X,edge_4.ctrl_Y,edge_4.ctrl_Z];
end

[ctrl_X,ctrl_Y,ctrl_Z]=geomMapGrid(control_2,control_3,control_4,control_1);

% process bound with axis
surf=FaceBSpline(name,ctrl_X,ctrl_Y,ctrl_Z,[],u_degree,v_degree);
end

function edge_list=getBound(Data,Idx_bound)
% get bound around face
%

edge_list=[];
% load bound
for bound_idx=1:length(Idx_bound)
    idx_bound=Idx_bound(bound_idx);
    str_bound=Data{idx_bound,2};

    % FACE_OUTER_BOUND
    idx=find(str_bound == ',',1);
    name=str_bound(1:idx-1);str_bound(1:idx)=[];
    idx=find(str_bound == ',',1);
    idx_loop=str2double(str_bound(2:idx-1));str_bound(1:idx)=[];
    orientation=str_bound;

    % EDGE_LOOP
    str_loop=Data{idx_loop,2};
    idx=find(str_loop == ',',1);
    name=str_loop(1:idx-1);str_loop(1:idx)=[];
    Idx_edge=str2double(strsplit(strrep(str_loop(2:end-1),'#',''),','));

    loop_edge_list=cell(1,length(Idx_edge));
    % ORIENTED_EDGE
    for edge_idx=1:length(Idx_edge)
        idx_edge=Idx_edge(edge_idx);

        % ORIENTED_EDGE
        str_oriented=Data{idx_edge,2};
        idx=find(str_oriented == ',',1);
        name=str_oriented(1:idx-1);str_oriented(1:idx)=[];
        idx=find(str_oriented == ',',1);
        idx_vertex_start=str2double(str_oriented(2:idx-1));str_oriented(1:idx)=[];
        idx=find(str_oriented == ',',1);
        idx_vertex_end=str2double(str_oriented(2:idx-1));str_oriented(1:idx)=[];
        idx=find(str_oriented == ',',1);
        idx_edge=str2double(str_oriented(2:idx-1));str_oriented(1:idx)=[];
        orientation=str_oriented;

        % EDGE_CURVE
        str_edge=Data{idx_edge,2};
        idx=find(str_edge == ',',1);
        name=str_edge(1:idx-1);str_edge(1:idx)=[];
        idx=find(str_edge == ',',1);
        idx_vertex_start=str2double(str_edge(2:idx-1));str_edge(1:idx)=[];
        idx=find(str_edge == ',',1);
        idx_vertex_end=str2double(str_edge(2:idx-1));str_edge(1:idx)=[];
        idx=find(str_edge == ',',1);
        idx_edge=str2double(str_edge(2:idx-1));str_edge(1:idx)=[];
        same_sense=str_edge;

        % CURVE
        edge_type=Data{idx_edge,1};
        switch edge_type
            case 'B_SPLINE_CURVE_WITH_KNOTS'
                edge=getEdgeBSpline(Data,idx_vertex_start,idx_vertex_end,idx_edge);
            case 'LINE'
                edge=getLINE(Data,idx_vertex_start,idx_vertex_end,idx_edge);
            otherwise
                error('getBound: unknown edge type');
        end

        % process orientation
        if contains(orientation,'F')
            % reverse
            edge.reverse();
        end

        loop_edge_list{edge_idx}=edge;
    end
    edge_list=[edge_list,loop_edge_list];
end
end

function edge=getEdgeBSpline(Data,idx_vertex_start,idx_vertex_end,idx_edge)
% generate BSpline edge by BSpline edge string data
%

str_edge=Data{idx_edge,2};

% main properties
idx=find(str_edge == ',',1);
name=str_edge(1:idx-1);str_edge(1:idx)=[];
idx=find(str_edge == ',',1);
degree=str2double(str_edge(1:idx-1));str_edge(1:idx)=[];

% read control point matrix
idx=find(str_edge == ')',1);
Idx_control=str2double(strsplit(strrep(str_edge(2:idx-1),'#',''),','));str_edge(1:idx+1)=[];

% properties
idx=find(str_edge == ',',1);
edge_form=(str_edge(1:idx-1));str_edge(1:idx)=[];
idx=find(str_edge == ',',1);
closed_edge=(str_edge(1:idx-1));str_edge(1:idx)=[];
idx=find(str_edge == ',',1);
self_intersect=(str_edge(1:idx-1));str_edge(1:idx)=[];

% knot_multi
idx=find(str_edge == ')',1);
knot_multi=str2double(split(str_edge(2:idx-1),',',2));str_edge(1:idx+1)=[];

% knot_list
idx=find(str_edge == ')',1);
knot_list=str2double(split(str_edge(2:idx-1),',',2));str_edge(1:idx+1)=[];

% knot_spec
knot_spec=str_edge;

% load control point data
ctrl_list=zeros(length(Idx_control),3);
for control_idx=1:length(Idx_control)
    idx_control=Idx_control(control_idx);
    ctrl_list(control_idx,:)=getPoint(Data,idx_control);
end

% VERTEX
str_vertex=Data{idx_vertex_start,2};
idx=find(str_vertex == ',',1);str_vertex(1:idx)=[];
idx_point_start=str2double(str_vertex(2:end));
point_start=getPoint(Data,idx_point_start);
str_vertex=Data{idx_vertex_end,2};
idx=find(str_vertex == ',',1);str_vertex(1:idx)=[];
idx_point_end=str2double(str_vertex(2:end));
point_end=getPoint(Data,idx_point_end);

% reverse
if norm(point_start-ctrl_list(1,:)) > 100*eps &&...
        norm(point_end-ctrl_list(end,:)) > 100*eps
    ctrl_list=flipud(ctrl_list);
    knot_multi=fliplr(knot_multi);
    knot_list=1-fliplr(knot_list);
end

edge=EdgeBSpline(name,ctrl_list(:,1),ctrl_list(:,2),ctrl_list(:,3),[],degree,...
    knot_multi,knot_list);

edge.edge_form=edge_form;
edge.closed_edge=closed_edge;
edge.self_intersect=self_intersect;
edge.knot_spec=knot_spec;
end

function edge=getLINE(Data,idx_vertex_start,idx_vertex_end,idx_edge)
% generate BSpline edge by line string data
%

str_edge=Data{idx_edge,2};

% name
idx=find(str_edge == ',',1);
name=str_edge(1:idx-1);str_edge(1:idx)=[];

% properties
idx=find(str_edge == ',',1);
idx_point=str2double(str_edge(2:idx-1));str_edge(1:idx)=[];
idx_vertor=str2double(str_edge(2:end));

% point start
point_start=getPoint(Data,idx_point);

% VECTOR
str_vector=Data{idx_vertor,2};
idx=find(str_vector == ',',1);str_vector(1:idx)=[];
idx=find(str_vector == ',',1);
idx_direction=str2double(str_vector(2:idx-1));str_vector(1:idx)=[];
magnitude=str2double(str_vector);

% DIRECTION
direction=getPoint(Data,idx_direction);

% VERTEX
str_vertex=Data{idx_vertex_start,2};
idx=find(str_vertex == ',',1);str_vertex(1:idx)=[];
idx_point_start=str2double(str_vertex(2:end));
point_start=getPoint(Data,idx_point_start);
str_vertex=Data{idx_vertex_end,2};
idx=find(str_vertex == ',',1);str_vertex(1:idx)=[];
idx_point_end=str2double(str_vertex(2:end));
point_end=getPoint(Data,idx_point_end);

point=[point_start;point_end];
edge=EdgeBSpline(name,point(:,1),point(:,2),point(:,3));

end

function point=getPoint(Data,idx_point)
str_point=Data{idx_point,2};
idx=find(str_point == ',',1);str_point(1:idx)=[];
point=str2double(strsplit(str_point(2:end-1),','));
end
