function drawPoint(mesh_point,dimension)
% write mesh point into file
%

switch dimension
    case 2
        curve_name_list=fieldnames(mesh_point);

        hold on;
        for curve_idx=1:length(curve_name_list)
            curve_name=curve_name_list{curve_idx};
            X=mesh_point.(curve_name).X;
            Y=mesh_point.(curve_name).Y;

            scatter(X,Y)
        end
        hold off;

    case 3

        surface_name_list=fieldnames(mesh_point);

        hold on;
        for surface_idx=1:length(surface_name_list)
            surface_name=surface_name_list{surface_idx};
            X=mesh_point.(surface_name).X;
            Y=mesh_point.(surface_name).Y;
            Z=mesh_point.(surface_name).Z;

            scatter3(X,Y,Z)
        end
        hold off;

        zlabel('z');view(3);axis equal;
    otherwise
        error('pointWrite: unknow dimension');
end
xlabel('x');ylabel('y');

end