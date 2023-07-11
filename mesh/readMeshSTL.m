function part=readMeshSTL(mesh_filestr,scale,file_type,INFORMATION,geometry_torlance)
% read mash data from stl file
%
% input:
% filename_mesh(support .stl file), scale(geometry zoom scale), ...
% file_type(stl code format, binary/ASCII), INFORMATION(true or false), ...
% geometry_torlance(default is 1e-12)
%
% output:
% part(only one part)
%
% part(part.name, part.mesh_list{mesh.element_list, mesh.element_type, mesh.element_ID})
%
if nargin < 4
    INFORMATION=[];
    if nargin < 3
        file_type=[];
        if nargin < 2
            scale=[];
        end
    end
end

if nargin < 5
    geometry_torlance=1e-12;
end

if isempty(scale)
    scale=1;
end
if isempty(INFORMATION)
    INFORMATION=true(1);
end

% cheak file definition
if length(mesh_filestr) > 4
    if ~strcmpi(mesh_filestr((end-3):end),'.stl')
        mesh_filestr=[mesh_filestr,'.stl'];
    end
else
    mesh_filestr=[mesh_filestr,'.stl'];
end
if exist(mesh_filestr,'file')~=2
    error('readMeshSTL: mesh file do not exist')
end


% detech file type
if isempty(file_type)
    file_mesh=fopen(mesh_filestr,'rb');
    string_head=fread(file_mesh,80,'int8');

    if string_head(end)==32
        file_type='binary';
        fseek(file_mesh,0,'bof');
    else
        file_type='ASCII';
        fclose(file_mesh);
        clear('file_mesh');

        % reopen as ASCII file
        file_mesh=fopen(mesh_filestr,'r');
    end
end

if INFORMATION
    disp(['readMeshSTL: read mash data begin, file type read as ',file_type]);
end

if strcmp(file_type,'binary')
    % read head
    part_name=fread(file_mesh,80,'int8');
    part_name=char(part_name');
    part_name=strsplit(part_name);
    part_name=part_name{2};

    % read total element_number
    part_element_number=fread(file_mesh,1,'int32');
    mesh_element_list=zeros(3*part_element_number,3);

    % read element
    overlap_list=[];
    for element_index=1:part_element_number
        % read normal vector
        vector_normal=fread(file_mesh,3,'float32');

        % read point coordinate
        point_1=fread(file_mesh,3,'float32')*scale;
        point_2=fread(file_mesh,3,'float32')*scale;
        point_3=fread(file_mesh,3,'float32')*scale;

        mesh_element_list(3*element_index-2:3*element_index,:)=...
            [point_1,point_2,point_3]';

        attribute=fread(file_mesh,1,'int16');

        if norm(point_1-point_2) < geometry_torlance ||...
                norm(point_2-point_3) < geometry_torlance ||...
                norm(point_3-point_1) < geometry_torlance
            overlap_list=[overlap_list;element_index];
        end
    end

    mesh_element_list([3*overlap_list-2;3*overlap_list-1;overlap_list],:)=[];
    part_element_number=part_element_number-length(overlap_list);
else
    % read head
    part_name=fgetl(file_mesh);
    part_name=strsplit(part_name);
    part_name=part_name{2};

    % initial sort space
    mesh_element_list=zeros(99,3);
    element_number=0;

    % read normal vector
    vector_normal_string=strsplit(fgetl(file_mesh));
    while ~strcmp(vector_normal_string{1},'endsolid')
        fgetl(file_mesh);

        % read point1 coordinate
        point_1=getASCIIPoint(file_mesh)*scale;
        point_2=getASCIIPoint(file_mesh)*scale;
        point_3=getASCIIPoint(file_mesh)*scale;

        fgetl(file_mesh);fgetl(file_mesh);

        % read normal vector
        vector_normal_string=strsplit(fgetl(file_mesh));

        if norm(point_1-point_2) < geometry_torlance ||...
                norm(point_2-point_3) < geometry_torlance ||...
                norm(point_3-point_1) < geometry_torlance
        else
            element_number=element_number+1;
            mesh_element_list(3*element_number-2:3*element_number,:)=...
                [point_1,point_2,point_3]';
        end

        % add more sort space
        if element_number*3 == size(mesh_element_list,1)
            mesh_element_list=[mesh_element_list;zeros(99,3)];
        end
    end
    mesh_element_list=mesh_element_list(1:element_number*3,:);

    part_element_number=size(mesh_element_list,1)/3;
end

fclose(file_mesh);
clear('file_mesh');

mesh.element_list=mesh_element_list;
mesh.element_type='stl';
mesh.element_ID=int8(5);
mesh.element_number=part_element_number;

part.name=part_name;
part.mesh_list={mesh};
part.element_number=part_element_number;

if INFORMATION
    disp('readMeshSTL: read mash data done');
end

    function point=getASCIIPoint(file_mesh)
        % read ASCII point data
        %
        point_string=strsplit(fgetl(file_mesh));
        point=zeros(3,1);
        for node_index=1:3
            point(node_index)=str2double(point_string{2+node_index});
        end
    end
end