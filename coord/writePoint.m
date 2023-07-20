function writePoint(mesh_point,point_filestr,dimension)
% write mesh point into file
%
[point_filedir,~,~]=fileparts(point_filestr);
if ~isempty(point_filedir) && ~exist(point_filedir,"dir")
    mkdir(point_filedir);
end

switch dimension
    case 2
        curve_name_list=fieldnames(mesh_point);

        point_file=fopen(point_filestr,'w');
        for curve_idx=1:length(curve_name_list)
            curve_name=curve_name_list{curve_idx};
            index=mesh_point.(curve_name).index;
            X=mesh_point.(curve_name).X;
            Y=mesh_point.(curve_name).Y;

            for point_index=1:size(index,1)
                fprintf(point_file,'%d %f %f\n',index(point_index)-1,X(point_index),Y(point_index));
            end
        end
        fclose(point_file);

    case 3

        surf_name_list=fieldnames(mesh_point);

        point_file=fopen(point_filestr,'w');
        for surf_idx=1:length(surf_name_list)
            surf_name=surf_name_list{surf_idx};
            index=mesh_point.(surf_name).index;
            X=mesh_point.(surf_name).X;
            Y=mesh_point.(surf_name).Y;
            Z=mesh_point.(surf_name).Z;

            for point_index=1:size(index,1)
                fprintf(point_file,'%d %f %f %f\n',index(point_index)-1,X(point_index),Y(point_index),Z(point_index));
            end
        end
        fclose(point_file);

    otherwise
        error('pointWrite: unknow dimension');
end

end