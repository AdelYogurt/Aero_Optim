function [part_list,point_list,geometry]=convertWGSToMesh...
    (part_list,remove_redundance,geometry_torlance)
% convert LaWGS format part into mesh format part
%
% if remove_redundance, function will call ADTreePoint to remove repeat
% point, change element_idx to same coordinate point
%
if nargin < 2
    remove_redundance=true(1);
end

if ~iscell(part_list)
    part_list={part_list};
end

if nargin < 3 || isempty(geometry_torlance)
    geometry_torlance=1e-12;
end

% add all surface to point list
point_list=[];
for part_idx=1:length(part_list)
    mesh_list=part_list{part_idx}.mesh_list;

    % process each mesh
    for mesh_idx=1:length(mesh_list)
        mesh=mesh_list{mesh_idx};
        if ~strcmp(mesh.element_type,'wgs')
            error('convertWGSToSTL: element_type of mesh is not LaWGS format');
        end

        X=mesh.X;
        Y=mesh.Y;
        Z=mesh.Z;

        if any(size(X) ~= size(Y)) || any(size(X) ~= size(Z))
            error('convertWGSToSTL: size of X,Y,Z in mesh are not equal');
        end

        point_mesh_list=[X(:),Y(:),Z(:)];

        % add this mesh point into point list
        point_list=[point_list;point_mesh_list];
    end
end

if remove_redundance
    % delete same coordination point
    % need ADTreePoint
    [ADT,idx_list]=ADTreePoint(point_list,[],[],geometry_torlance);
    
    exist_point_number=int32(0); % offset idx of point idx list
    for part_idx=1:length(part_list)
        mesh_list=part_list{part_idx}.mesh_list;
        part_element_number=0;
        
        % process each mesh
        for mesh_idx=1:length(mesh_list)
            mesh=mesh_list{mesh_idx};

            [point_number,line_number]=size(mesh.X);
            mesh_point_number=point_number*line_number;

            % initialize element sort space
            calculate_element=(point_number-1)*(line_number-1)*2;
            element_number=0;
            element_list=zeros(calculate_element,3);
            
            % sort element(S3)
            for line_idx=1:line_number-1
                base_offset=point_number*(line_idx-1);
                for point_idx=1:point_number-1
                    point1_idx=idx_list(base_offset+point_idx+exist_point_number);
                    point2_idx=idx_list(base_offset+point_idx+point_number+exist_point_number);
                    point3_idx=idx_list(base_offset+point_idx+1+point_number+exist_point_number);
                    point4_idx=idx_list(base_offset+point_idx+1+exist_point_number);

                    if point1_idx == point2_idx || point2_idx == point3_idx || point1_idx == point3_idx
                        % small element degeneration to line, discard it
                    else
                        element_number=element_number+1;
                        element_list(element_number,:)=...
                            [int32(point1_idx),int32(point2_idx),int32(point3_idx)];
                    end

                    if point3_idx == point4_idx || point4_idx == point1_idx || point1_idx == point3_idx
                        % small element degeneration to line, discard it
                    else
                        element_number=element_number+1;
                        element_list(element_number,:)=...
                            [int32(point3_idx),int32(point4_idx),int32(point1_idx)];
                    end
                end
            end

            mesh=rmfield(mesh,'X');
            mesh=rmfield(mesh,'Y');
            mesh=rmfield(mesh,'Z');

            % sort element
            mesh.element_list=element_list(1:element_number,:);
            mesh.element_type='S3';
            mesh.element_ID=int8(5);

            mesh_list{mesh_idx}=mesh;

            exist_point_number=exist_point_number+int32(mesh_point_number);
            part_element_number=part_element_number+element_number;
        end

        part_list{part_idx}.mesh_list=mesh_list;
        part_list{part_idx}.element_number=part_element_number;
    end

    geometry.min_bou=ADT.min_bou;
    geometry.max_bou=ADT.max_bou;
    geometry.dimension=3;
else
    % do not delete same coordination point
    exist_point_number=int32(0); % offset idx of point idx list
    for part_idx=1:length(part_list)
        mesh_list=part_list{part_idx}.mesh_list;

        % process each mesh
        delete_idx=[];
        for mesh_idx=1:length(mesh_list)
            mesh=mesh_list{mesh_idx};

            [point_number,line_number]=size(mesh.X);
            mesh_point_number=point_number*line_number;

            % initialize element sort space
            calculate_element=(point_number-1)*(line_number-1)*2;
            element_number=0;
            element_list=zeros(calculate_element,3);
            
            % sort element(S3)
            for line_idx=1:line_number-1
                base_offset=point_number*(line_idx-1);
                for point_idx=1:point_number-1
                    point1_idx=base_offset+point_idx+exist_point_number;
                    point2_idx=base_offset+point_idx+point_number+exist_point_number;
                    point3_idx=base_offset+point_idx+1+point_number+exist_point_number;
                    point4_idx=base_offset+point_idx+1+exist_point_number;

                    d12=point_list(point2_idx,:)-point_list(point1_idx,:);
                    d23=point_list(point3_idx,:)-point_list(point2_idx,:);

                    if norm(d12) < eps || norm(d23) < eps
                        % small element degeneration to line, discard it
                    else
                        element_number=element_number+1;
                        element_list(element_number,:)=...
                            [int32(point1_idx),int32(point2_idx),int32(point3_idx)];
                    end

                    d34=point_list(point4_idx,:)-point_list(point3_idx,:);
                    d41=point_list(point1_idx,:)-point_list(point4_idx,:);

                    if norm(d34) < eps || norm(d41) < eps
                        % small element degeneration to line, discard it
                    else
                        element_number=element_number+1;
                        element_list(element_number,:)=...
                            [int32(point3_idx),int32(point4_idx),int32(point1_idx)];
                    end
                end
            end

            if element_number == 0
                delete_idx=[delete_idx;mesh_idx];
                continue;
            end

            mesh=rmfield(mesh,'X');
            mesh=rmfield(mesh,'Y');
            mesh=rmfield(mesh,'Z');

            % sort element
            mesh.element_list=element_list(1:element_number,:);
            mesh.element_type='S3';
            mesh.element_ID=int8(5);

            mesh_list{mesh_idx}=mesh;

            exist_point_number=exist_point_number+int32(mesh_point_number);
        end

        mesh_list(delete_idx)=[];

        part_list{part_idx}.mesh_list=mesh_list;
    end

    geometry.min_bou=min(point_list);
    geometry.max_bou=max(point_list);
    geometry.dimension=3;
end

if length(part_list) == 1
    part_list=part_list{1};
end

end
