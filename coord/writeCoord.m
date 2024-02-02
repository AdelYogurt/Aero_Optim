function writeCoord(coord_data,dimension,coor_filedir)
% writer local coordinate into txt file
%
if ~exist(coor_filedir,"dir")
    mkdir(coor_filedir);
end

switch dimension
    case 1
        edge_name_list=fieldnames(coord_data);

        for edge_idx=1:length(edge_name_list)
            edge_name=edge_name_list{edge_idx};
            index=coord_data.(edge_name).index;
            U=coord_data.(edge_name).U;
            
            % write coord
            coord_file=fopen(fullfile(coor_filedir,[edge_name,'_local_coord.txt']),'w');
            fprintf(coord_file,'%d %.12f\n',[double(index),U]');
            fclose(coord_file);
        end

    case 2
        face_name_list=fieldnames(coord_data);

        for face_idx=1:length(face_name_list)
            face_name=face_name_list{face_idx};
            index=coord_data.(face_name).index;
            U=coord_data.(face_name).U;
            V=coord_data.(face_name).V;

            % write coord
            coord_file=fopen(fullfile(coor_filedir,[face_name,'_local_coord.txt']),'w');
            fprintf(coord_file,'%d %.12f %.12f\n',[double(index),U,V]');
            fclose(coord_file);
        end
    otherwise
        error('coordWrite: unknow dimension');
end

end