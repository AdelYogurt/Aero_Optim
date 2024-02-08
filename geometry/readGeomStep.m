function geom_list=readGeomStep(step_filestr)
% read step file
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
geom_list=cell(1,length(Idx_body));
for body_idx=1:length(Idx_body)
    idx_body=Idx_body(body_idx);
    body_type=Data{idx_body,1};

    switch body_type
        case 'SHELL_BASED_SURFACE_MODEL'
            str_body=Data{idx_body,2};
            str_list=split(str_body,',');
            shell_name=str_list{1}(2:end-1);
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
                            fce=getFaceNURBS(Data,idx_face,idx_surf);
                        case 'PLANE'
                            fce=getPLANE(Data,idx_face,idx_surf);
                        otherwise
                            continue;
                            % error('readGeomStep: unknown surface type');
                    end
                    face_shell{face_idx}=fce;
                end
                face_list=[face_list,face_shell];
            end

            body=Shell(shell_name,face_list);

        case 'MANIFOLD_SOLID_BREP'
            str_body=Data{idx_body,2};
            str_list=split(str_body,',');
            solid_name=str_list{1}(2:end-1);
            Idx_solid=str2double(str_list{2}(2:end));

            body=Solid(solid_name);
        otherwise
            error('readGeomStep: unknown body type');
    end
    geom_list{body_idx}=body;
end

end

function fce=getFaceNURBS(Data,idx_face,idx_surf)
% generate BSpline surface by BSpline surface string data
%

str_surf=Data{idx_surf,2};
% main properties
idx=find(str_surf == ',',1);
name=str_surf(2:idx-2);str_surf(1:idx)=[];
idx=find(str_surf == ',',1);
UDegree=str2double(str_surf(1:idx-1));str_surf(1:idx)=[];
idx=find(str_surf == ',',1);
VDegree=str2double(str_surf(1:idx-1));str_surf(1:idx)=[];

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
face_form=(str_surf(1:idx-1));str_surf(1:idx)=[];
idx=find(str_surf == ',',1);
UPeriodic=(str_surf(1:idx-1));str_surf(1:idx)=[];
idx=find(str_surf == ',',1);
VPeriodic=(str_surf(1:idx-1));str_surf(1:idx)=[];
idx=find(str_surf == ',',1);
self_intersect=(str_surf(1:idx-1));str_surf(1:idx)=[];

% Mults
idx=find(str_surf == ')',1);
UMults=str2double(split(str_surf(2:idx-1),',',2));str_surf(1:idx+1)=[];
idx=find(str_surf == ')',1);
VMults=str2double(split(str_surf(2:idx-1),',',2));str_surf(1:idx+1)=[];

% Knots
idx=find(str_surf == ')',1);
UKnots=str2double(split(str_surf(2:idx-1),',',2));str_surf(1:idx+1)=[];
idx=find(str_surf == ')',1);
VKnots=str2double(split(str_surf(2:idx-1),',',2));str_surf(1:idx+1)=[];

% knot_spec
knot_spec=str_surf;

% load control point data
[u_pole_num,v_pole_num]=size(Idx_control);
Poles=zeros(v_pole_num,u_pole_num,3);
for u_idx=1:u_pole_num
    for v_idx=1:v_pole_num
        idx_point=Idx_control(u_idx,v_idx);
        control_point=getPoint(Data,idx_point);
        Poles(v_idx,u_idx,:)=control_point;
    end
end

fce=Surface(name,Poles,UDegree,VDegree,UMults,VMults,UKnots,VKnots);

fce.face_form=face_form;
fce.UPeriodic=UPeriodic;
fce.VPeriodic=VPeriodic;
fce.self_intersect=self_intersect;
fce.knot_spec=knot_spec;
end

function fce=getPLANE(Data,idx_face,surf_idx)
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
name=str_surf(2:idx-2);str_surf(1:idx)=[];
idx_axis=str2double(str_surf(2:end));

% axis
str_axis=Data{idx_axis,2};

% process edge_list
if length(edge_list) > 4
    fce=[];
    return;
end

% if edge num equal to 4, edge order is 0v, u1, 1v, u0
switch length(edge_list)
    case 2
        edge_0v=edge_list{1};
        edge_1v=edge_list{2};
        VDegree=max(edge_0v.Degree,edge_1v.Degree);
        UDegree=1;

        edge_0v.addDegree(VDegree);Poles_0v=[edge_0v.Poles];
        edge_1v.addDegree(VDegree);Poles_1v=[edge_1v.Poles];

        Poles_u1=repmat(Poles_0v(end,:),2,1);
        Poles_u0=repmat(Poles_0v(1,:),2,1);
    case 3
        edge_0v=edge_list{1};
        edge_u1=edge_list{2};
        edge_1v=edge_list{3};
        UDegree=max(edge_u1.Degree);
        VDegree=max(edge_0v.Degree,edge_1v.Degree);

        edge_0v.addDegree(VDegree);Poles_0v=[edge_0v.Poles];
        edge_u1.addDegree(UDegree);Poles_u1=[edge_u1.Poles];
        edge_1v.addDegree(VDegree);Poles_1v=[edge_1v.Poles];

        Poles_u0=repmat(Poles_0v(1,:),u_num,1);
    case 4
        edge_0v=edge_list{1};
        edge_u1=edge_list{2};
        edge_1v=edge_list{3};
        edge_u0=edge_list{4};
        UDegree=max(edge_u0.Degree,edge_u1.Degree);
        VDegree=max(edge_0v.Degree,edge_1v.Degree);

        edge_0v.addDegree(VDegree);Poles_0v=[edge_0v.Poles];
        edge_u1.addDegree(UDegree);Poles_u1=[edge_u1.Poles];
        edge_1v.addDegree(VDegree);Poles_1v=[edge_1v.Poles];
        edge_u0.addDegree(UDegree);Poles_u0=[edge_u0.Poles];
end

Poles=GeomApp.MapGrid(Poles_u0,Poles_u1,Poles_0v,Poles_1v);

% process bound with axis
fce=Surface(name,Poles,UDegree,VDegree);
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
    name=str_bound(2:end-1);str_bound(1:idx)=[];
    idx=find(str_bound == ',',1);
    idx_loop=str2double(str_bound(2:idx-1));str_bound(1:idx)=[];
    orientation=str_bound;

    % EDGE_LOOP
    str_loop=Data{idx_loop,2};
    idx=find(str_loop == ',',1);
    name=str_loop(2:end-1);str_loop(1:idx)=[];
    Idx_edge=str2double(strsplit(strrep(str_loop(2:end-1),'#',''),','));

    loop_edge_list=cell(1,length(Idx_edge));
    % ORIENTED_EDGE
    for edge_idx=1:length(Idx_edge)
        idx_edge=Idx_edge(edge_idx);

        % ORIENTED_EDGE
        str_oriented=Data{idx_edge,2};
        idx=find(str_oriented == ',',1);
        name=str_oriented(2:end-1);str_oriented(1:idx)=[];
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
        name=str_edge(2:end-1);str_edge(1:idx)=[];
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
                edge=getEdgeNURBS(Data,idx_vertex_start,idx_vertex_end,idx_edge);
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

function edge=getEdgeNURBS(Data,idx_vertex_start,idx_vertex_end,idx_edge)
% generate BSpline edge by BSpline edge string data
%

str_edge=Data{idx_edge,2};

% main properties
idx=find(str_edge == ',',1);
name=str_edge(2:end-1);str_edge(1:idx)=[];
idx=find(str_edge == ',',1);
Degree=str2double(str_edge(1:idx-1));str_edge(1:idx)=[];

% read control point matrix
idx=find(str_edge == ')',1);
Idx_control=str2double(strsplit(strrep(str_edge(2:idx-1),'#',''),','));str_edge(1:idx+1)=[];

% properties
idx=find(str_edge == ',',1);
edge_form=(str_edge(1:idx-1));str_edge(1:idx)=[];
idx=find(str_edge == ',',1);
Periodic=(str_edge(1:idx-1));str_edge(1:idx)=[];
idx=find(str_edge == ',',1);
self_intersect=(str_edge(1:idx-1));str_edge(1:idx)=[];

% Mults
idx=find(str_edge == ')',1);
Mults=str2double(split(str_edge(2:idx-1),',',2));str_edge(1:idx+1)=[];

% Knots
idx=find(str_edge == ')',1);
Knots=str2double(split(str_edge(2:idx-1),',',2));str_edge(1:idx+1)=[];

% knot_spec
knot_spec=str_edge;

% load control point data
Poles=zeros(length(Idx_control),3);
for control_idx=1:length(Idx_control)
    idx_control=Idx_control(control_idx);
    Poles(control_idx,:)=getPoint(Data,idx_control);
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
if norm(point_start-Poles(1,:)) > 100*eps &&...
        norm(point_end-Poles(end,:)) > 100*eps
    Poles=flipud(Poles);
    Mults=fliplr(Mults);
    Knots=1-fliplr(Knots);
end

edge=EdgeNURBS(name,Poles,Degree,Mults,Knots);

edge.edge_form=edge_form;
edge.Periodic=Periodic;
edge.self_intersect=self_intersect;
edge.knot_spec=knot_spec;
end

function edge=getLINE(Data,idx_vertex_start,idx_vertex_end,idx_edge)
% generate BSpline edge by line string data
%

str_edge=Data{idx_edge,2};

% name
idx=find(str_edge == ',',1);
name=str_edge(2:end-1);str_edge(1:idx)=[];

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
edge=EdgeNURBS(name,point);

end

function point=getPoint(Data,idx_point)
str_point=Data{idx_point,2};
idx=find(str_point == ',',1);str_point(1:idx)=[];
point=str2double(strsplit(str_point(2:end-1),','));
end
