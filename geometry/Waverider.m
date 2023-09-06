classdef Waverider < Body
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
        par_G;
        par_F;
    end

    properties
        % calculate parameter
        head_length;
        body_length;

        % local coordinate parameter
        blunt_u=0.002; % is local length parameter
        blunt_eta=0.005; % is local length parameter

        % function
        body_v_max_fun;
        shape_head_edge_x;
        shape_tri_wing_edge_x;
        shape_wing_edge_x;
        shape_wing_fcn;
    end

    % define function
    methods
        function self=Waverider(total_length,varargin)
            % generate parameterized waverider
            % u, x=X(u)
            % v, y=Y(u)
            % w, z=Z(u,v)
            %
            [par_width,par_hight_up,par_hight_low,...
                par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
                par_G,par_F]=Waverider.decode(varargin{1});

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
            self.par_G=par_G;
            self.par_F=par_F;

            % calculate blunt radius parameter
            % calculate gradient of up and low discrete surface to local radius center

            % calculate head radius of zox plane of head
            if par_M_up > 1
                radius_center_up=0;
            else
                dz_dx=(par_hight_up*(self.blunt_u)^par_M_up)/(total_length*self.blunt_u);
                radius_center_up=par_R*dz_dx;
            end
            if par_M_low > 1
                radius_center_low=0;
            else
                dz_dx=(par_hight_low*(self.blunt_u)^par_M_low)/(total_length*self.blunt_u);
                radius_center_low=par_R*dz_dx;
            end
            radius_center_head=max(radius_center_up,radius_center_low);
            radius_head_sq=(radius_center_head*radius_center_head+par_R*par_R);
            shape_head_edge_x=@(Y) (sqrt(radius_head_sq-Y.^2))-radius_center_head;
            self.shape_head_edge_x=shape_head_edge_x;

            % calculate head side blunt radius of yoz plane
            dz_dx=(self.blunt_eta)^par_N_up*(1-self.blunt_eta)^par_N_up/(0.5)^(2*par_N_up);
            radius_center_up=par_R*dz_dx;
            dz_dx=(self.blunt_eta)^par_N_low*(1-self.blunt_eta)^par_N_low/(0.5)^(2*par_N_low);
            radius_center_low=par_R*dz_dx;
            radius_center_head_side=max(radius_center_low,radius_center_up);
            radius_head_side_sq=(radius_center_head_side*radius_center_head_side+par_R*par_R);
            radius_head_side=sqrt(radius_head_side_sq);
            shape_head_side_x=@(Y) (sqrt(radius_head_side_sq-Y.^2))-radius_center_head_side;
            side_length=radius_head_side-radius_center_head_side;

            %% define head
            head_up=SurfaceCST3D('head_up',total_length,par_width,par_hight_up,[],[par_T,0],[par_M_up,0],[par_N_up,par_N_up],global_symmetry_y);
            head_low=SurfaceCST3D('head_low',total_length,par_width,-par_hight_low,[],[par_T,0],[par_M_low,0],[par_N_low,par_N_low],global_symmetry_y);
            
            tran_fun_Z_up=@(U,V) par_G*U.^par_F;
            tran_fun_Z_low=@(U,V) par_G*U.^par_F;
            % deform surface
            head_up.addDeform([],[],tran_fun_Z_up)
            head_low.addDeform([],[],tran_fun_Z_low)

            % blunt translation
            head_up.addTranslation(0,0,par_R);
            head_low.addTranslation(0,0,-par_R);

            %% define blunt head side
            % local coordinate, local x is global -x, local y is global z, local z is global y
            shape_fcn_X=@(V) 1+shape_head_side_x((V-0.5)*par_R*2)/total_length;
            shape_fcn_Z=@(U,V) (sqrt(radius_head_side_sq-(((V-0.5)*2)*par_R).^2)-radius_center_head_side).*(1-U).^0.3;
            head_side=SurfaceCST3D('head_side',total_length,2*par_R,1,shape_fcn_X,[],shape_fcn_Z,[0.0,0.0]);

            Y_edge=@(U) U.^par_T*par_width/2;
            tran_fun_Y_up=@(U) par_G*U.^par_F;
            tran_fun_Z=@(U,V) Y_edge(1-U);
            % deform surface
            head_side.addDeform([],tran_fun_Y_up,tran_fun_Z);

            % rotation surface
            head_side.addRotation(90,180,0);

            % translation surface
            head_side.addTranslation(total_length,0,-par_R);

            %% define back
            shape_fcn_Y_up=@(U) par_hight_up*(0.5+U/2).^par_N_up.*(0.5-U/2).^par_N_up/(0.5)^(2*par_N_up)+par_R;
            shape_fcn_Y_low=@(U) par_hight_low*(0.5+U/2).^par_N_low.*(0.5-U/2).^par_N_low/(0.5)^(2*par_N_low)+par_R;

            % local coordinate, local x is global y, local y is global z
            back_up=SurfaceCST3D('body_back_up',par_width/2,1,0,[],shape_fcn_Y_up);
            back_low=SurfaceCST3D('body_back_low',par_width/2,-1,0,[],shape_fcn_Y_low);

            % rotation to global coordinate
            back_up.addRotation(90,90,0);
            back_low.addRotation(90,90,0);

            % translation to global coordination
            back_up.addTranslation(total_length,0,par_G);
            back_low.addTranslation(total_length,0,par_G);

            %% define back side
            % local coordinate, local x is global -x, local y is global z
            shape_fcn_X=@(V) shape_head_side_x((V-0.5)*2*par_R);

            back_side=SurfaceCST3D('back_side',side_length,2*par_R,0,shape_fcn_X,[]);

            back_side.addRotation(90,90,0);

            % translation to global coordinate
            back_side.addTranslation(total_length,par_width/2,-par_R+par_G);

            %% sort data
            self.surface_list={head_up,head_low,head_side,back_up,back_low,back_side};
        end
    
    end

    % visualization function
    methods
        function surf_total=calSurfaceMatrix(self,varargin)
            % generate waverider wing wgs data
            % u, x=X(u)
            % v, y=Y(u)
            % w, z=Z(u,v)
            %
            % notice colume vector in matrix is wgs format
            %

            surf_num=length(self.surface_list);
            surf_total=cell(surf_num,1);
            if nargin <= 2
                if nargin == 1
                    value_torl=1e-3;
                else
                    value_torl=varargin{2};
                end

                % calculate all surface matrix
                for surf_idx=1:surf_num
                    surf.name=self.surface_list{surf_idx}.name;
                    [surf.X,surf.Y,surf.Z,surf.U,surf.V]=self.surface_list{surf_idx}.calSurface(value_torl);
                    surf.element_type='wgs';
                    surf_total{surf_idx}=surf;
                end
            else
                u_grid_num_head=varargin{1};
                v_grid_num_head=varargin{2};
                edge_gird_num=varargin{3};

                num_list=[
                    u_grid_num_head,v_grid_num_head;
                    u_grid_num_head,v_grid_num_head;
                    u_grid_num_head,edge_gird_num*2;
                    v_grid_num_head,edge_gird_num;
                    v_grid_num_head,edge_gird_num;
                    edge_gird_num,edge_gird_num*2];
                name_list={'head_up';'head_low';'head_side';'back_up';'back_low';'back_side'};

                for surf_idx=1:length(name_list)
                    surf_CST=self.getSurface(name_list{surf_idx});
                    surf.name=surf_CST.name;
                    [surf.X,surf.Y,surf.Z,surf.U,surf.V]=surf_CST.calSurface(num_list(surf_idx,1),num_list(surf_idx,2));
                    surf.element_type='wgs';
                    surf_total{surf_idx}=surf;
                end
            end

            % fix normal vector
            for surf_idx=1:surf_num
                surf=surf_total{surf_idx};
                if contains(surf.name,'low')
                    surf.X=flipud(surf.X);
                    surf.Y=flipud(surf.Y);
                    surf.Z=flipud(surf.Z);
                    surf.U=flipud(surf.U);
                    surf.V=flipud(surf.V);
                end
                surf_total{surf_idx}=surf;
            end
        end

        function mesh_data=getWGSMesh(self,varargin)
            % interface of calSurfaceMatrix to get WGS part
            %
            if length(varargin) > 1
                surf_total=self.calSurfaceMatrix...
                    (varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
            elseif length(varargin) == 1
                surf_total=self.calSurfaceMatrix(varargin{1});
            else
                surf_total=self.calSurfaceMatrix();
            end

            mesh_data=struct();
            for surf_index=1:length(surf_total)
                surf=surf_total{surf_index};
                marker_name=surf.name;

                mesh_data.(marker_name).X=surf.X;
                mesh_data.(marker_name).Y=surf.Y;
                mesh_data.(marker_name).Z=surf.Z;
                mesh_data.(marker_name).type='wgs';
            end
        end

        function writeStepOpenShell(self,step_filestr,varargin)
            % write surface into step file
            %
            [~,step_filename,~]=fileparts(step_filestr);

            % write head
            step_file=fopen(step_filestr,'w');
            object_index=1;
            object_index=writeStepHead(self,step_file,object_index,step_filename);

            % write surface
            ADVANCED_FACE_index_list=zeros(1,length(self.surface_list));
            if length(varargin) > 1
                surf_total=self.calSurfaceMatrix...
                    (varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
            elseif length(varargin) == 1
                surf_total=self.calSurfaceMatrix(varargin{1});
            else
                surf_total=self.calSurfaceMatrix();
            end

            for surf_idx=1:length(surf_total)
                surf=surf_total{surf_idx};
                if size(surf.X,2) == 2
                    u_list=[0,0,1,1];
                else
                    u_list=[0,0,0,surf.U(1,:),1,1,1];
                end
                if size(surf.X,1) == 2
                    v_list=[0,0,1,1];
                else
                    v_list=[0,0,0,surf.V(:,1)',1,1,1];
                end
                surf=SurfaceBSpline(surf_name,[],[],[],surf.X,surf.Y,surf.Z,[],[],u_list,v_list);
                [step_str,object_index,ADVANCED_FACE_index_list(surf_idx)]=surf.getStepNode(object_index);
                fprintf(step_file,step_str);
                fprintf(step_file,'\n');
            end

            % generate OPEN_SHELL
            OPEN_SHELL_index=object_index;
            step_str=[num2str(object_index,'#%d'),' = OPEN_SHELL ',...
                '( ''NONE'', ',...
                '( ',num2str(ADVANCED_FACE_index_list(1:end-1),'#%d, '),' ',num2str(ADVANCED_FACE_index_list(end),'#%d'),' )',...
                ' );\n'];object_index=object_index+1;
            fprintf(step_file,step_str);
            fprintf(step_file,'\n');

            % write model
            SHELL_BASED_SURFACE_MODEL_index=object_index;
            step_str=[num2str(object_index,'#%d'),' = SHELL_BASED_SURFACE_MODEL ',...
                '( ''NONE'', ( ',...
                num2str(OPEN_SHELL_index,'#%d'),' )',...
                ' );\n'];object_index=object_index+1;
            fprintf(step_file,step_str);
            fprintf(step_file,'\n');

            % write end of step file
            writeStepEnd(self,step_file,object_index,step_filename,SHELL_BASED_SURFACE_MODEL_index);

            fclose(step_file);
            clear('step_file');

        end

        function drawBody(self,varargin)
            % show all surface of body
            %
            if length(varargin) > 1
                surf_total=self.calSurfaceMatrix...
                    (varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
            elseif length(varargin) == 1
                surf_total=self.calSurfaceMatrix(varargin{1});
            else
                surf_total=self.calSurfaceMatrix();
            end

            for surf_index=1:length(surf_total)
                surf=surf_total{surf_index};
                surface(surf.X,surf.Y,surf.Z);
            end

            xlabel('x');
            ylabel('y');
            zlabel('z');
            axis equal;
            view(3);
        end

    end

    % parameterized function
    methods
        function WWD_coord=calCoord(self,point_list,surf_index_list)
            % calculate all input point local coordinate
            %
            WWD_coord=struct();

            surf_name_list=fieldnames(surf_index_list);
            for surf_idx=1:length(surf_name_list)
                % get surf
                surf_name=surf_name_list{surf_idx};
                surf=self.getSurface(surf_name);
                if isempty(surf)
                    continue;
                end
                point_idx=surf_index_list.(surf_name);

                % calculate coordinate
                point=point_list(point_idx,1:3);
                [U,V,~,~,~]=surf.calCoordinate(point(:,1),point(:,2),point(:,3));
                WWD_coord.(surf_name).index=point_idx;
                WWD_coord.(surf_name).U=U;
                WWD_coord.(surf_name).V=V;
            end

        end

        function mesh_data=calMeshPoint(self,WWD_coord)
            % calculate all mesh point by input WWD coord
            %
            mesh_data=struct();

            surf_name_list=fieldnames(WWD_coord);
            for surf_idx=1:length(surf_name_list)
                % get surface
                surf_name=surf_name_list{surf_idx};
                surf=self.getSurface(surf_name);
                if isempty(surf)
                    continue;
                end
                point_idx=WWD_coord.(surf_name).index;
                U=WWD_coord.(surf_name).U;
                V=WWD_coord.(surf_name).V;

                % calculate coordinate
                [X,Y,Z]=surf.calPoint(U,V);
                mesh_data.(surf_name).type='scatter';
                mesh_data.(surf_name).index=point_idx;
                mesh_data.(surf_name).X=X;
                mesh_data.(surf_name).Y=Y;
                mesh_data.(surf_name).Z=Z;
            end

        end
    end

    % decode x function
    methods(Static)
        function [par_width,par_hight_up,par_hight_low,...
                par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
                par_G,par_F]=decode(param)
            % decode x into parameter
            %
            if isnumeric(param)
                param=num2cell(param,[1,11]);
            end

            % length parameter
            par_M_up=param{1};
            par_M_low=param{2};
            % width parameter
            par_width=param{3};
            par_T=param{4};
            par_N_up=param{5};
            par_N_low=param{6};
            % height parameter
            par_hight_up=param{7};
            par_hight_low=param{8};
            par_R=param{9};
            % wing parameter
            par_G=param{10};
            par_F=param{11};
        end
    end
end