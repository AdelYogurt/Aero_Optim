function mesh_data=calMesh(point_list,marker_index_list)
% calculate mesh data by point list and marker index list
%
marker_name_list=fieldnames(marker_index_list);
for marker_idx=1:length(marker_name_list)
    marker_name=marker_name_list{marker_idx};
    index=marker_index_list.(marker_name);

    mesh_data.(marker_name).index=index;
    mesh_data.(marker_name).X=point_list(index,1);
    mesh_data.(marker_name).Y=point_list(index,2);
    if size(point_list,2) == 3
        mesh_data.(marker_name).Z=point_list(index,3);
    end
     mesh_data.(marker_name).type='scatter';
end
end