classdef WaveriderBase
    properties
        % CST parameter
        total_length;
        par_width;
        par_hight_up;
        par_hight_low;
        par_T;
        par_M_up;
        par_M_low;
        par_N_up;
        par_N_low;
        par_R;
            
        blunt_xi;
        blunt_eta;

        % face object list
        surface_list;

        % function
        shape_head_side_x
        shape_head_side_y
        shape_back_side_y
    end

    methods
        function self=WaveriderBase...
                (total_length,par_width,par_hight_up,par_hight_low,...
                par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R)
            % generate waverider wing wgs data
            % u, x=X(u)
            % v, y=Y(u)
            % w, z=Z(u,v)
            %
            % notice:
            % colume vector in matrix is wgs format
            % par_rho1 should between 0 and 1
            %
            global_symmetry_y=true(1);

            self.total_length=total_length;
            self.par_width=par_width;
            self.par_hight_up=par_hight_up;
            self.par_hight_low=par_hight_low;
            self.par_T=par_T;
            self.par_M_up=par_M_up;
            self.par_M_low=par_M_low;
            self.par_N_up=par_N_up;
            self.par_N_low=par_N_low;
            self.par_R=par_R;

            self.blunt_xi=0.002; % is local length parameter
            self.blunt_eta=0.005;

            % calculate blunt radius parameter
            % local coordinate, local x is global x, local y is global z
            % head radius
            % calculate gradient of up and low discrete surface to local radius center
            if par_M_up > 1
                radius_center_up=0;
            else
                dz_dx=(par_hight_up*(self.blunt_xi)^par_M_up)/(total_length*self.blunt_xi);
                radius_center_up=par_R*dz_dx;
            end
            if par_M_low > 1
                radius_center_low=0;
            else
                dz_dx=(par_hight_low*(self.blunt_xi)^par_M_low)/(total_length*self.blunt_xi);
                radius_center_low=par_R*dz_dx;
            end
            
            radius_center_head_sied_x=max(radius_center_up,radius_center_low);
            radius_head_side_x_sq=(radius_center_head_sied_x*radius_center_head_sied_x+par_R*par_R);
            radius_head_side_x=sqrt(radius_head_side_x_sq);
            shape_head_side_x=@(Y) (sqrt(radius_head_side_x_sq-Y.^2))-radius_center_head_sied_x;
            self.shape_head_side_x=shape_head_side_x;

            % calculate head side blunt radius of yoz plane
            dz_dx=(self.blunt_eta)^par_N_up*(1-self.blunt_eta)^par_N_up/(0.5)^(2*par_N_up)*par_hight_up/(par_width*self.blunt_eta);
            radius_center_up=par_R*dz_dx;
            dz_dx=(self.blunt_eta)^par_N_low*(1-self.blunt_eta)^par_N_low/(0.5)^(2*par_N_low)*par_hight_low/(par_width*self.blunt_eta);
            radius_center_low=par_R*dz_dx;
            
            radius_center_head_side_y=max(radius_center_low,radius_center_up);
            radius_head_side_y_sq=(radius_center_head_side_y*radius_center_head_side_y+par_R*par_R);
            radius_head_side_y=sqrt(radius_head_side_y_sq);
            shape_head_side_y=@(Y) (sqrt(radius_head_side_y_sq-Y.^2))-radius_center_head_side_y;
            self.shape_head_side_y=shape_head_side_y;

            shape_back_side_y=@(X) (sqrt(radius_head_side_y_sq-(X+radius_center_head_side_y).^2));
            self.shape_back_side_y=shape_back_side_y;

            side_high=radius_head_side_y-radius_center_head_side_y;

            %% define head
            head_up=SurfaceCST3D(total_length,par_width,par_hight_up,[],[par_T,0],[par_M_up,0],{[par_N_up,par_N_up]},global_symmetry_y);
            head_low=SurfaceCST3D(total_length,par_width,-par_hight_low,[],[par_T,0],[par_M_low,0],{[par_N_low,par_N_low]},global_symmetry_y);
            
            % blunt translation
            head_up.addTranslation(0,0,par_R);
            head_low.addTranslation(0,0,-par_R);

            %% define blunt head side
            % local coordinate, local x is global -x, local y is global z, local z is global y
            shape_fcn_X=@(V) 1+shape_head_side_x((V-0.5)*par_R*2)/total_length;
            shape_fcn_Z=@(U,V) shape_head_side_y((V-0.5)*par_R*2);
            head_side=SurfaceCST3D(total_length,2*par_R,1,shape_fcn_X,[],[0.0,0.01],shape_fcn_Z);

            Y_edge=@(U) U.^par_T*par_width/2;
            tran_fun_Z=@(U,V) Y_edge(1-U);
            % deform surface
            head_side.addDeform([],[],tran_fun_Z);

            % rotation surface
            head_side.addRotation(90,180,0);

            % translation surface
            head_side.addTranslation(total_length,0,-par_R);

            %% define back
            shape_fcn_Y_up=@(U) par_hight_up*(0.5+U/2).^par_N_up.*(0.5-U/2).^par_N_up/(0.5)^(2*par_N_up)+par_R;
            shape_fcn_Y_low=@(U) par_hight_low*(0.5+U/2).^par_N_low.*(0.5-U/2).^par_N_low/(0.5)^(2*par_N_low)+par_R;

            % local coordinate, local x is global y, local y is global z
            back_up=SurfaceCST3D(par_width/2,1,0,[],shape_fcn_Y_up);
            back_low=SurfaceCST3D(par_width/2,-1,0,[],shape_fcn_Y_low);

            % rotation to global coordinate
            back_up.addRotation(90,90,0);
            back_low.addRotation(90,90,0);

            % translation to global coordination
            back_up.addTranslation(total_length,0,0);
            back_low.addTranslation(total_length,0,0);


            %% define back side
            % local coordinate, local x is global -x, local y is global z
            shape_fcn_Y_up=@(U) shape_back_side_y(U*side_high);
            shape_fcn_Y_low=@(U) shape_back_side_y(U*side_high);

            back_side_up=SurfaceCST3D(side_high,1,0,[],shape_fcn_Y_up);
            back_side_low=SurfaceCST3D(side_high,-1,0,[],shape_fcn_Y_low);

            back_side_up.addRotation(90,90,0);
            back_side_low.addRotation(90,90,0);

            % translation to global coordinate
            back_side_up.addTranslation(total_length,par_width/2,0);
            back_side_low.addTranslation(total_length,par_width/2,0);

            %% sort data
            self.surface_list=struct('head_up',head_up,'head_low',head_low,'head_side',head_side,...
                'back_up',back_up,'back_low',back_low,...
                'back_side_up',back_side_up,'back_side_low',back_side_low);

        end

        function [X,Y,Z]=calSurfaceMatrix...
                (self,xi_grid_num_head,psi_grid_num_head,edge_gird_num)
            % generate waverider wing wgs data
            % u, x=X(u)
            % v, y=Y(u)
            % w, z=Z(u,v)
            %
            % notice colume vector in matrix is wgs format
            %

            % auto allocation u
            head_xi_list=linspace(0,1,xi_grid_num_head+1);
            head_xi_list=head_xi_list.^(1/self.par_T);

            % calculate head
            [U,V]=meshgrid(head_xi_list,linspace(0.5,1,psi_grid_num_head+1));
            [X_head,Y_head,Z_head_up]=self.surface_list.head_up.calSurface(U,V);
            [X_head,Y_head,Z_head_low]=self.surface_list.head_low.calSurface(U,V);

            % calculate blunt head side
            head_side_xi_list=1-head_xi_list;
            [U,W]=meshgrid(head_side_xi_list,linspace(0,1,edge_gird_num*2+1));
            [X_head_side,Y_head_side,Z_head_side]=self.surface_list.head_side.calSurface(U,W);

            % calculate back
            [X_back,Y_back,Z_back_up]=self.surface_list.back_up.calSurface(psi_grid_num_head,edge_gird_num);
            [X_back,Y_back,Z_back_low]=self.surface_list.back_low.calSurface(psi_grid_num_head,edge_gird_num);

            % calculate back side
            back_side_eta_list=linspace(1,0,edge_gird_num+1);
            back_side_eta_list=self.shape_head_side_y(back_side_eta_list*self.par_R)/self.surface_list.back_side_up.LX;
            [X_back_side,Y_back_side,Z_back_side_up]=self.surface_list.back_side_up.calSurface(back_side_eta_list,edge_gird_num);
            [X_back_side,Y_back_side,Z_back_side_low]=self.surface_list.back_side_low.calSurface(back_side_eta_list,edge_gird_num);

            %% sort data
            X=struct();Y=struct();Z=struct();
            [X,Y,Z]=self.sortSurfaceUpLow(X,Y,Z,'head',X_head,Y_head,Z_head_up,Z_head_low);
            [X,Y,Z]=self.sortSurface(X,Y,Z,'head_side',X_head_side,Y_head_side,Z_head_side);
            [X,Y,Z]=self.sortSurfaceUpLow(X,Y,Z,'back',X_back,Y_back,Z_back_up,Z_back_low);
            [X,Y,Z]=self.sortSurfaceUpLow(X,Y,Z,'back_side',X_back_side,Y_back_side,Z_back_side_up,Z_back_side_low);

        end
    end

    methods
        function [X,Y,Z]=sortSurface(self,X,Y,Z,name,X_surf,Y_surf,Z_surf)
            X.(name)=X_surf;
            Y.(name)=Y_surf;
            Z.(name)=Z_surf;
        end

        function [X,Y,Z]=sortSurfaceUpLow(self,X,Y,Z,name,X_surf,Y_surf,Z_surf_up,Z_surf_low)
            X.([name,'_up'])=X_surf;
            Y.([name,'_up'])=Y_surf;
            Z.([name,'_up'])=Z_surf_up;
            X.([name,'_low'])=X_surf;
            Y.([name,'_low'])=Y_surf;
            Z.([name,'_low'])=Z_surf_low;
        end
    
        function part=getWGSMesh(self,part_name,...
                xi_grid_num_head,psi_grid_num_head,edge_gird_num)
            surface_name_list=fieldnames(self.surface_list);
            surface_number=length(surface_name_list);

            [X_total,Y_total,Z_total]=self.calSurfaceMatrix...
                (xi_grid_num_head,psi_grid_num_head,edge_gird_num);

            mesh_list=cell(surface_number,1);
            for surface_index=1:surface_number
                surface_name=surface_name_list{surface_index};
                if contains(surface_name,'up') || contains(surface_name,'front')
                    mesh.X=X_total.(surface_name);
                    mesh.Y=Y_total.(surface_name);
                    mesh.Z=Z_total.(surface_name);
                else
                    mesh.X=fliplr(X_total.(surface_name));
                    mesh.Y=fliplr(Y_total.(surface_name));
                    mesh.Z=fliplr(Z_total.(surface_name));
                end
                mesh.element_type='wgs';
                mesh_list{surface_index,1}=mesh;
            end

            part.name=part_name;
            part.mesh_list=mesh_list;
        end

        function drawBody(self,...
                xi_grid_num_head,psi_grid_num_head,edge_gird_num)
            surface_name_list=fieldnames(self.surface_list);
            surface_number=length(surface_name_list);

            [X_total,Y_total,Z_total]=self.calSurfaceMatrix...
                (xi_grid_num_head,psi_grid_num_head,edge_gird_num);

            hold on;
            for surface_index=1:surface_number
                surface_name=surface_name_list{surface_index};
                surf(X_total.(surface_name),Y_total.(surface_name),Z_total.(surface_name))
            end
            hold off;

            xlabel('x');
            ylabel('y');
            zlabel('z');
            axis equal;
            view(3);
        end

    end
end