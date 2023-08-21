clc;
clear;
close all hidden;


%% read CGNS file

cgns_filestr='hermes.cgns';

marker_name_list={'HERMES'};

[point_list,marker_index_list]=readMeshCGNSMarker(cgns_filestr,marker_name_list);

%% generate coord

% mesh_point is struct
% mesh_point.marker which contain index, X, Y, Z
hermes_idx=marker_index_list.HERMES;
mesh_point.HERMES.index=hermes_idx;
mesh_point.HERMES.X=point_list(hermes_idx,1);
mesh_point.HERMES.Y=point_list(hermes_idx,2);
mesh_point.HERMES.Z=point_list(hermes_idx,3);

% % write mesh point
% writePoint(mesh_point,'hermes_deform.dat',3)

% draw mesh point
drawPoint(mesh_point,3)


