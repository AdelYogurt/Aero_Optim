function writeCoord(coord_data,dimension,coor_filedir)
% writer local coordinate into txt file
%
if ~exist(coor_filedir,"dir")
    mkdir(coor_filedir);
end

switch dimension
    case 1
        curve_name_list=fieldnames(coord_data);

        for curve_idx=1:length(curve_name_list)
            curve_name=curve_name_list{curve_idx};
            index=coord_data.(curve_name).index;
            U=coord_data.(curve_name).U;
            
            % write coord
            coord_file=fopen(fullfile(coor_filedir,[curve_name,'_local_coord.txt']),'w');
            for point_index=1:length(index)
                fprintf(coord_file,'%d %.12f\n',index(point_index),U(point_index));
            end
            fclose(coord_file);
        end

    case 2
        surface_name_list=fieldnames(coord_data);

        for surf_idx=1:length(surface_name_list)
            surf_name=surface_name_list{surf_idx};
            index=coord_data.(surf_name).index;
            U=coord_data.(surf_name).U;
            V=coord_data.(surf_name).V;

            % write coord
            coord_file=fopen(fullfile(coor_filedir,[surf_name,'_local_coord.txt']),'w');
            for point_index=1:length(index)
                fprintf(coord_file,'%d %.12f %.12f\n',index(point_index),U(point_index),V(point_index));
            end
            fclose(coord_file);
        end
    otherwise
        error('coordWrite: unknow dimension');
end

end