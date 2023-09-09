clc;
clear;
% close all hidden;

%% read CGNS file

% mesh_filestr='hermes.cgns';
% marker_name_list={'HERMES'};

mesh_filestr='../W/W.cgns';
marker_name_list={'head_up','head_low','head_side','back'};

[point_list,marker_index_list]=readBCCGNS(mesh_filestr,marker_name_list);

%% generate coord

% mesh_point is struct
% mesh_point.marker which contain index, X, Y, Z
marker_name_list=fieldnames(marker_index_list);
for marker_idx=1:length(marker_name_list)
    marker_name=marker_name_list{marker_idx};
    point_idx=marker_index_list.(marker_name);
    mesh_point.(marker_name).index=point_idx;
    mesh_point.(marker_name).X=point_list(point_idx,1);
    mesh_point.(marker_name).Y=point_list(point_idx,2);
    mesh_point.(marker_name).Z=point_list(point_idx,3);
    mesh_point.(marker_name).type='scatter';
end

% % write mesh point
% writePoint(mesh_point,'hermes_deform.dat',3)

% draw mesh point
drawPoint(mesh_point,3)
