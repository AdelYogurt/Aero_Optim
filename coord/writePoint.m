function writePoint(mesh_point,point_filestr,dimension)
% write mesh point into file
%
[point_filedir,~,~]=fileparts(point_filestr);
if ~isempty(point_filedir) && ~exist(point_filedir,"dir")
    mkdir(point_filedir);
end

switch dimension
    case 2
        edge_name_list=fieldnames(mesh_point);

        Idx_origin=[];X_origin=[];Y_origin=[];
        for edge_idx=1:length(edge_name_list)
            edge_name=edge_name_list{edge_idx};

            Idx_origin=[Idx_origin;mesh_point.(edge_name).index];
            X_origin=[X_origin;mesh_point.(edge_name).X];
            Y_origin=[Y_origin;mesh_point.(edge_name).Y];
        end

        [Idx,~,Map]=unique(Idx_origin);
        Repeat=zeros(length(Idx),1);X=zeros(length(Idx),1);Y=zeros(length(Idx),1);
        for point_index=1:length(Idx_origin)
            idx=Map(point_index);
            Repeat(idx)=Repeat(idx)+1;
            X(idx)=X(idx)+X_origin(point_index);
            Y(idx)=Y(idx)+Y_origin(point_index);
        end
        X=X./Repeat;Y=Y./Repeat;

        point_file=fopen(point_filestr,'w');
        for point_index=1:size(Idx,1)
            fprintf(point_file,'%d %f %f\n',Idx(point_index)-1,X(point_index),Y(point_index));
        end
        fclose(point_file);
    case 3
        face_name_list=fieldnames(mesh_point);

        Idx_origin=[];X_origin=[];Y_origin=[];Z_origin=[];
        for face_idx=1:length(face_name_list)
            face_name=face_name_list{face_idx};

            Idx_origin=[Idx_origin;mesh_point.(face_name).index];
            X_origin=[X_origin;mesh_point.(face_name).X];
            Y_origin=[Y_origin;mesh_point.(face_name).Y];
            Z_origin=[Z_origin;mesh_point.(face_name).Z];
        end

        [Idx,~,Map]=unique(Idx_origin);
        Repeat=zeros(length(Idx),1);X=zeros(length(Idx),1);Y=zeros(length(Idx),1);Z=zeros(length(Idx),1);
        for point_index=1:length(Idx_origin)
            idx=Map(point_index);
            Repeat(idx)=Repeat(idx)+1;
            X(idx)=X(idx)+X_origin(point_index);
            Y(idx)=Y(idx)+Y_origin(point_index);
            Z(idx)=Z(idx)+Z_origin(point_index);
        end
        X=X./Repeat;Y=Y./Repeat;Z=Z./Repeat;

        point_file=fopen(point_filestr,'w');
        for point_index=1:size(Idx,1)
            fprintf(point_file,'%d %f %f %f\n',Idx(point_index)-1,X(point_index),Y(point_index),Z(point_index));
        end
        fclose(point_file);
    otherwise
        error('pointWrite: unknow dimension');
end

end