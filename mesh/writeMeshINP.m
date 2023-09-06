function writeMeshINP(mesh_filestr,mesh_data,marker_name_list)
% write mesh data into inp file(ABAQUS analysis format)
%
% input:
% mesh_filestr, mesh_data, marker_name_list(default all markers)
%
% notice:
% mesh_data(single zone): mesh_data.geometry, mesh_data.(marker)
% marker: marker.type, marker.ID, marker.element_list, marker.number_list
% geometry: point_list, dimension
%
if nargin < 3
    marker_name_list=[];
end
if isempty(marker_name_list)
    marker_name_list=fieldnames(mesh_data);
end
marker_index=1;
while marker_index <= length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    if strcmp(marker_name,'geometry')
        marker_name_list(marker_index)=[];
    else
        marker_index=marker_index+1;
    end
end

file_mesh=fopen(mesh_filestr,'w');

% write basic information of inp file
fprintf(file_mesh,'*Heading\n');
fprintf(file_mesh,'** Job name: Job-%s Model name: Model-%s\n',mesh_filestr,mesh_filestr);
fprintf(file_mesh,'** Generated by: matlab::writeMeshINP\n');
fprintf(file_mesh,'*Preprint, echo=NO, model=NO, history=NO, contact=NO\n');

% write start of part
fprintf(file_mesh,'**\n');
fprintf(file_mesh,'** PARTS\n');
fprintf(file_mesh,'**\n');

geometry=mesh_data.geometry;

% write all part
point_list=geometry.point_list;
for marker_index=1:length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    marker=mesh_data.(marker_name);

    % write part name
    fprintf(file_mesh,'*Part, name=%s\n',marker_name);

    % notice, element index in each part of inp file start from 1
    % search min and max element index of mesh
    element_list=marker.element_list;
    number_list=marker.number_list;
    if marker.ID == 20
        % for mixed element, separte it into element list and number list
        % than, element list will be sort by ID
        % the same ID element will be print in same *Element
        [element_list,ID_list]=getSeparteElement(element_list,number_list);
    end

    min_node_index=min(min(element_list));
    max_point_index=max(max(element_list));

    % write part point data
    fprintf(file_mesh,'*Node\n');
    for point_index=min_node_index:max_point_index
        fprintf(file_mesh,'% 7d,% 14f,% 14f,% 14f\n',point_index-min_node_index+1,point_list(point_index,1:3));
    end

    % write element data
    % notice, element index in each part of inp file start from 1
    if marker.ID == 20
        element_number=length(ID_list);
        [ID_list,index_list]=sort(ID_list);
        map_list=sort(index_list);
        node_index_list=cumsum(number_list);

        % sort element_list
        element_list_old=element_list;
        node_index=0;
        for element_index=1:element_number
            node_number=number_list(map_list(element_index));
            if map_list(element_index) ~= element_index
                element_list(node_index+(1:node_number))=element_list_old...
                    (node_index_list(map_list(element_index))-(node_number:-1:1)+1);
            end
            node_index=node_index+node_number;
        end

        id=ID_list(1);
        [type,node_number]=convertIDToType(id);
        fprintf(file_mesh,'*Element, type=%s\n',type);
        print_format=['%d',repmat(',%d',1,node_number),'\n'];

        % print element list
        node_index=0;
        for element_index=1:element_number
            % if element type change, state new type
            if ID_list(element_index) ~= id
                id=ID_list(element_index);
                [type,node_number]=convertIDToType(id);
                fprintf(file_mesh,'*Element, type=%s\n',type);
                print_format=['%d',repmat(',%d',1,node_number),'\n'];
            end

            fprintf(file_mesh,print_format,element_index,element_list(node_index+(1:node_number))-min_node_index+1);
            node_index=node_index+node_number;
        end
    else
        fprintf(file_mesh,'*Element, type=%s\n',convertTypeToType(marker.type));

        [element_number,node_number]=size(element_list);
        print_format=['%d',repmat(',%d',1,node_number),'\n'];
        for element_index=1:element_number
            fprintf(file_mesh,print_format,element_index,element_list(element_index,1:node_number)-min_node_index+1);
        end
    end

    % write node setting
    if isfield(marker,'Nset') && ~isempty(marker.Nset)
        Nset=marker.Nset;
        Nset_name_list=fieldnames(Nset);
        for Nset_index=1:length(Nset_name_list)
            Nset_name=Nset_name_list{Nset_index};
            node_index=Nset.(Nset_name);

            if length(node_index) == node_index(end)
                fprintf(file_mesh,'*Nset, nset=%s, generate\n',Nset_name);
                fprintf(file_mesh,'%d,%d,%d\n',node_index(1),node_index(end),1);
            else
                fprintf(file_mesh,'*Nset, nset=%s\n',Nset_name);
                node_number=length(node_index);
                if node_number > 10
                    number=floor(node_number/10)*10;
                    fprintf(file_mesh,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n',node_index(1:number));
                    node_index=node_index(number+1:end);
                end
                fprintf(file_mesh,'%d,',node_index(1:end-1));fprintf(file_mesh,'%d\n',node_index(end));
            end
        end
    end

    % write element setting
    if isfield(marker,'Elset') && ~isempty(marker.Elset)
        Elset=marker.Elset;
        Elset_name_list=fieldnames(Elset);
        for Elset_index=1:length(Elset_name_list)
            Elset_name=Elset_name_list{Elset_index};
            node_index=Elset.(Elset_name);

            if length(node_index) == node_index(end)
                fprintf(file_mesh,'*Elset, elset=%s, generate\n',Elset_name);
                fprintf(file_mesh,'%d,%d,%d\n',node_index(1),node_index(end),1);
            else
                fprintf(file_mesh,'*Elset, elset=%s\n',Elset_name);
                node_number=length(node_index);
                if node_number > 10
                    number=floor(node_number/10)*10;
                    fprintf(file_mesh,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n',node_index(1:number));
                    node_index=node_index(number+1:end);
                end
                fprintf(file_mesh,'%d,',node_index(1:end-1));fprintf(file_mesh,'%d\n',node_index(end));
            end
        end
    end

    % write material setting
    if isfield(marker,'material') && ~isempty(marker.material)
        material=marker.material;
        Elset_name_list=fieldnames(material);
        for material_index=1:length(Elset_name_list)
            Elset_name=Elset_name_list{material_index};
            fprintf(file_mesh,'*Solid Section, elset=%s, material=%s\n',Elset_name,material.(Elset_name));
        end
        fprintf(file_mesh,',\n');
    end

    % end part
    fprintf(file_mesh,'*End Part\n');
    fprintf(file_mesh,'**\n');
end

if isfield(geometry,'residual_data')
    fprintf(file_mesh,'%s\n',geometry.residual_data);
else
    % write start of assembly
    fprintf(file_mesh,'**\n');
    fprintf(file_mesh,'** ASSEMBLY\n');
    fprintf(file_mesh,'**\n');
    fprintf(file_mesh,'*Assembly, name=Assembly\n');

    % import part
    for marker_index=1:length(marker_name_list)
        marker_name=marker_name_list{marker_index};
        if strcmp(marker_name,'point_list') || strcmp(marker_name,'name')
            continue;
        end

        fprintf(file_mesh,'**  \n');
        fprintf(file_mesh,'*Instance, name=%s, part=%s\n',marker_name,marker_name);
        fprintf(file_mesh,'*End Instance\n');
        fprintf(file_mesh,'**\n');
    end

    % write end of assembly
    fprintf(file_mesh,'*End Assembly\n');
end

fclose(file_mesh);
clear('file_mesh');

end

function [element_list,ID_list]=getSeparteElement(element_list,number_list)
% separate mixed into element_list, number_list, ID_list
%

data_number=length(element_list);
% calculate element number
ID_list=zeros(length(number_list),1,'uint32');
index_list=zeros(length(number_list),1,'uint32');

element_index=1;
index=1;
while index < data_number
    ID_list(element_index)=element_list(index);
    index_list(element_index)=index;

    index=index+number_list(element_index)+1;
    element_index=element_index+1;
end

element_list(index_list)=[];
end

function [type,node_number]=convertIDToType(id)
% inp version of type and ID converter
%
switch id
    case 3
        type='B21';
        node_number=2;
    case 5
        type='S3';
        node_number=3;
    case 7
        type='S4R';
        node_number=4;
    case 8
        type='S8R';
        node_number=8;
    case 9
        type='S9R';
        node_number=9;
    case 10
        type='C3D4';
        node_number=4;
    case 17
        type='C3D8R';
        node_number=8;
    otherwise
        error('idType: unknown identifier')
end
end

function type=convertTypeToType(type)
% inp version of type and ID converter
%
switch type
    case 'BAR_2'
        type='B21';
    case 'TRI_3'
        type='S3';
    case 'QUAD_4'
        type='S4R';
    case 'QUAD_8'
        type='S8R';
    case 'QUAD_9'
        type='S9R';
    case 'TETRA_4'
        type='C3D4';
    case 'HEXA_8'
        type='C3D8R';
    otherwise
end
end
