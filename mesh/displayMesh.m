function displayMesh(mesh_data,marker_name_list,fig_hdl,draw_option)
% draw mesh in mesh_data
%
if nargin < 4
    draw_option=struct();
    if nargin < 3
        fig_hdl=[];
        if nargin < 2
            marker_name_list=[];
        end
    end

end

if isempty(marker_name_list)
    marker_name_list=fieldnames(mesh_data);
end
if ischar(marker_name_list)
    marker_name_list={marker_name_list};
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

if isfield(mesh_data,'geometry') && isfield(mesh_data.geometry,'point_list')
    point_list=mesh_data.geometry.point_list;
    dimension=size(point_list,2);
else
    point_list=[];
    dimension=3;
end

for marker_index=1:length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    if strcmp(marker_name,'point_list') || strcmp(marker_name,'name')
        continue;
    end

    % check if part.name exist in part_name_list
    if ~isfield(mesh_data,marker_name)
        continue;
    end

    marker=mesh_data.(marker_name);

    if strcmp(marker.type,'scatter')
        % scatter format data
        hold on;
        if isfield(marker,'Z')
            scatter3(marker.X,marker.Y,marker.Z);
        else
            scatter(marker.X,marker.Y);
            dimension=2;
        end
        hold off;
    elseif strcmp(marker.type,'wgs')
        % LaWGS format data
        hold on;
        if isfield(marker,'Z')
            surface(marker.X,marker.Y,marker.Z,'FaceColor','none');
        else
            surface(marker.X,marker.Y,'FaceColor','none');
            dimension=2;
        end
        hold off;
    elseif strcmp(marker.type,'stl')
        % stl format data
        element_list=marker.element_list;
        element_number=size(element_list,1)/3;

        for element_index=1:element_number
            point_index_list=[3*element_index-2:3*element_index,3*element_index-2];
            line(element_list(point_index_list,1), ...
                element_list(point_index_list,2), ...
                element_list(point_index_list,3))
        end
    elseif strcmp(marker.type,'MIXED2') || strcmp(marker.type,'MIXED3')
        element_list=marker.element_list;
        number_list=marker.number_list;
        element_number=length(number_list);

        node_idx=0;
        for element_index=1:element_number
            node_num=number_list(element_index);

            % element point index
            point_index_list=[element_list(node_idx+(1:node_num)+1);element_list(node_idx+2)];

            if dimension == 2
                line(point_list(point_index_list,1), ...
                    point_list(point_index_list,2));
            else
                line(point_list(point_index_list,1), ...
                    point_list(point_index_list,2), ...
                    point_list(point_index_list,3));
            end

            node_idx=node_idx+node_num+1;
        end
    else
        element_list=marker.element_list;
        element_number=size(element_list,1);

        for element_index=1:element_number
            % element point index
            point_index_list=[element_list(element_index,:),element_list(element_index,1)];

            if dimension == 2
                line(point_list(point_index_list,1), ...
                    point_list(point_index_list,2));
            else
                line(point_list(point_index_list,1), ...
                    point_list(point_index_list,2), ...
                    point_list(point_index_list,3));
            end
        end
    end
end

axis equal;
xlabel('x');
ylabel('y');

x_range=xlim();
y_range=ylim();
if dimension ~= 2
    zlabel('z');
    view(3);
    z_range=zlim();
    center=[mean(x_range),mean(y_range),mean(z_range)];
    range=max([x_range(2)-x_range(1),y_range(2)-y_range(1),z_range(2)-z_range(1)])/2;
    xlim([center(1)-range,center(1)+range]);
    ylim([center(2)-range,center(2)+range]);
    zlim([center(3)-range,center(3)+range]);
else
    center=[mean(x_range),mean(y_range)];
    range=max([x_range(2)-x_range(1),y_range(2)-y_range(1)])/2;
    xlim([center(1)-range,center(1)+range]);
    ylim([center(2)-range,center(2)+range]);
end

end


