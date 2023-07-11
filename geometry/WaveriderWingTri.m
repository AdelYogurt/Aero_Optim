classdef WaveriderWingTri
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
        par_rho1;
        par_rho12;
        par_WS1;
        par_WS2;
        par_WT_up;
        par_WT_low;
        global_symmetry_y;

        % face object list
        surface_list;

        % calculate parameter
        head_length;
        y_cut;
        body_psi_max_fun;

        % local coordinate parameter
        first_xi;
        second_xi;
    end

    methods
        function self=WaveriderWingTri...
                (total_length,par_width,par_hight_up,par_hight_low,...
                par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
                par_rho1,par_rho12,par_WS1,par_WS2,par_WT_up,par_WT_low)
            % generate waverider wing wgs data
            % xi, x=X(xi)
            % psi, y=Y(xi)
            % zeta, z=Z(xi,psi)
            %
            % notice colume vector in matrix is wgs format
            %
            if par_rho1 >= 1 || par_rho1 <= 0
                error('genWaveriderWing: do not support pure waverider');
            end
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
            self.par_rho1=par_rho1;
            self.par_rho12=par_rho12;
            self.par_WS1=par_WS1;
            self.par_WS2=par_WS2;
            self.par_WT_up=par_WT_up;
            self.par_WT_low=par_WT_low;

            self.global_symmetry_y=global_symmetry_y;

            self.first_xi=0.002;
            self.second_xi=0.010;

            head_length=(1-par_rho1)*total_length;
            self.head_length=head_length;

            y_cut=(1-par_rho1).^par_T*par_width/2;
            self.y_cut=y_cut;
            body_psi_max_fun=@(xi) y_cut./(xi.^par_T*par_width)-eps;
            self.body_psi_max_fun=body_psi_max_fun;

            % wing base parameter
            psi_end_cut=self.body_psi_max_fun(1);

            tri_wing_length=par_rho1*total_length;
            tri_wing_width=par_WS1;
            tri_wing_height_up=(psi_end_cut+0.5)^par_N_up*(0.5-psi_end_cut)^par_N_up/(0.5)^(2*par_N_up)*par_hight_up;
            tri_wing_height_low=(psi_end_cut+0.5)^par_N_low*(0.5-psi_end_cut)^par_N_low/(0.5)^(2*par_N_low)*par_hight_low;

            wing_length=par_rho1*par_rho12*total_length;
            wing_width=par_WS2;
            wing_height_up=par_WT_up;
            wing_height_low=par_WT_low;

            % calculate blunt radius parameter
            % local coordinate, local x is global x, local y is global z
            % head radius
            % calculate gradient of up and low discrete surface to local radius center
            if par_M_up > 1
                radius_center_up=0;
            else
                dz_dx=(par_hight_up*(self.first_xi)^par_M_up)/(total_length*self.first_xi);
                radius_center_up=par_R*dz_dx;
            end
            if par_M_low > 1
                radius_center_low=0;
            else
                dz_dx=(par_hight_low*(self.first_xi)^par_M_low)/(total_length*self.first_xi);
                radius_center_low=par_R*dz_dx;
            end
            radius_center_head=(radius_center_up+radius_center_low)/2;
            radius_head_sq=(sum(radius_center_head*radius_center_head+par_R*par_R));
            radius_head=sqrt(radius_head_sq);
            shape_head_edge_x=@(Y) real(sqrt(radius_head_sq-Y.^2))-radius_center_head;

            % calculate blunt radius of zox plane of tri wing
            tan_theta=((1-par_rho12)*par_rho1*total_length)/tri_wing_width;
            cos_theta=1/sqrt(tan_theta^2+1);
            radius_center_tri_wing=0;
            radius_tri_wing_sq=(par_R/cos_theta).^2;
            radius_tri_wing=sqrt(radius_tri_wing_sq);
            shape_tri_wing_edge_x=@(Y) real(sqrt(1-Y.^2/par_R/par_R))*radius_tri_wing-radius_center_tri_wing;

            % calculate blunt radius and center of zox plane of wing
            dz_dx=par_WT_up/wing_length;
            radius_center_up=par_R*dz_dx;
            dz_dx=par_WT_low/wing_length;
            radius_center_low=par_R*dz_dx;
            radius_center_wing=(radius_center_up+radius_center_low)/2;
            radius_wing_sq=(sum(radius_center_wing*radius_center_wing+par_R*par_R));
            radius_wing=sqrt(radius_wing_sq);
            shape_wing_edge_x=@(Y) real(sqrt(radius_wing_sq-Y.^2))-radius_center_wing;

            %% calculate head
            head_up=CST3DSurface(total_length,par_width,par_hight_up,[],[par_T,0],[par_M_up,0],{[par_N_up,par_N_up]},global_symmetry_y);
            head_low=CST3DSurface(total_length,par_width,-par_hight_low,[],[par_T,0],[par_M_low,0],{[par_N_low,par_N_low]},global_symmetry_y);
            
            % blunt translation
            head_up.addTranslation(0,0,par_R);
            head_low.addTranslation(0,0,-par_R);

            %% calculate blunt head side
            % local coordinate, local x is global -x, local y is global z, local z is global y
            shape_fcn_X=@(PSI) (shape_head_edge_x((PSI-0.5)*par_R*2)/head_length)-(shape_tri_wing_edge_x((PSI-0.5)*par_R*2)/head_length)+1;
            head_side=CST3DSurface(head_length,2*par_R,par_R,shape_fcn_X,[],[0.0,0.1],{[0.5,0.5]});

            tran_fun_X=@(PSI) shape_tri_wing_edge_x((PSI-0.5)*par_R*2);
            Y_edge=@(XI) XI.^par_T*par_width/2;
            tran_fun_Z=@(XI,PSI) Y_edge((1-XI)*(1-par_rho1));
            % deform surface
            head_side.addDeform(tran_fun_X,[],tran_fun_Z);

            % rotation surface
            head_side.addRotation(90,180,0);

            % translation surface
            head_side.addTranslation(head_length,0,-par_R);

            %% calculate body
            body_up=CST3DSurface(total_length,par_width,par_hight_up,[],[par_T,0],[par_M_up,0],{[par_N_up,par_N_up]},global_symmetry_y);
            body_low=CST3DSurface(total_length,par_width,-par_hight_low,[],[par_T,0],[par_M_low,0],{[par_N_low,par_N_low]},global_symmetry_y);

            % blunt translation
            body_up.addTranslation(0,0,par_R);
            body_low.addTranslation(0,0,-par_R);

            %% calculate body back
            shape_fcn_Y_up=@(XI) par_hight_up*(0.5+XI).^par_N_up.*(0.5-XI).^par_N_up/(0.5)^(2*par_N_up)+par_R;
            shape_fcn_Y_low=@(XI) par_hight_low*(0.5+XI).^par_N_low.*(0.5-XI).^par_N_low/(0.5)^(2*par_N_low)+par_R;

            % local coordinate, local x is global y, local y is global z
            body_back_up=CST3DSurface(par_width,1,0,[],shape_fcn_Y_up);
            body_back_low=CST3DSurface(par_width,-1,0,[],shape_fcn_Y_low);

            % rotation to global coordinate
            body_back_up.addRotation(90,90,0);
            body_back_low.addRotation(90,90,0);

            % translation to global coordination
            body_back_up.addTranslation(total_length,0,0);
            body_back_low.addTranslation(total_length,0,0);

            %% calculate transition wing
            % local coordination, local x is global y, local y is global -x
            shape_fcn_Y=@(XI) (1-XI*(1-par_rho12));
            shape_body_size_up=@(PSI) body_up.LZ.*...
                body_up.shape_fcn_Z(1-PSI*par_rho1,0.5+body_psi_max_fun(1-PSI*par_rho1)).*...
                body_up.class_fcn_Z(1-PSI*par_rho1,0.5+body_psi_max_fun(1-PSI*par_rho1));
            shape_body_size_low=@(PSI) -body_low.LZ.*...
                body_low.shape_fcn_Z(1-PSI*par_rho1,0.5+body_psi_max_fun(1-PSI*par_rho1)).*...
                body_low.class_fcn_Z(1-PSI*par_rho1,0.5+body_psi_max_fun(1-PSI*par_rho1));
            shape_fcn_Z_up=@(XI,PSI) (1-XI).*shape_body_size_up(PSI)+XI.*self.shapeWing(1-PSI)*(wing_height_up);
            shape_fcn_Z_low=@(XI,PSI) (1-XI).*shape_body_size_low(PSI)+XI.*self.shapeWing(1-PSI)*(wing_height_low);

            tri_wing_up=CST3DSurface(tri_wing_width,tri_wing_length,1,[],shape_fcn_Y,shape_fcn_Z_up,{[0,0]});
            tri_wing_low=CST3DSurface(tri_wing_width,tri_wing_length,-1,[],shape_fcn_Y,shape_fcn_Z_low,{[0,0]});

            % add rotation
            tri_wing_up.addRotation(0,0,90);
            tri_wing_low.addRotation(0,0,90);

            % add translation
            tri_wing_up.addTranslation(total_length,y_cut,par_R);
            tri_wing_low.addTranslation(total_length,y_cut,-par_R);

            %% calculate blunt tir wing front
            % local coordinate, local x is global -y, local y is global z, local z is global -x
            shape_fcn_X=@(PSI) 1-par_R*sqrt(1-(2*(PSI-0.5)).^2)/tri_wing_width; % if do not want blunt intersection, remove this and change blunt head side class N1
            tri_wing_front=CST3DSurface(tri_wing_width,2*par_R,radius_tri_wing,shape_fcn_X,[],[],{[0.5,0.5]});

            % deform surface
            tran_fun_Z=@(XI,PSI) XI*(par_rho1*(1-par_rho12)*total_length);
            tri_wing_front.addDeform([],[],tran_fun_Z);

            % rotation surface
            tri_wing_front.addRotation(90,-90,0);

            % translation surface
            tri_wing_front.addTranslation(total_length-wing_length,y_cut+par_WS1,-par_R);

            %% calculate tri wing back
            % local coordination, local x is global y, local y is global z
            shape_fcn_Y_up=@(XI) (1-XI)*(tri_wing_height_up+par_R)+XI*(wing_height_up+par_R);
            shape_fcn_Y_low=@(XI) (1-XI)*(tri_wing_height_low+par_R)+XI*(wing_height_low+par_R);

            tri_wing_back_up=CST3DSurface(tri_wing_width,1,0,[],shape_fcn_Y_up);
            tri_wing_back_low=CST3DSurface(tri_wing_width,-1,0,[],shape_fcn_Y_low);

            % rotation to global coordinate
            tri_wing_back_up.addRotation(90,90,0);
            tri_wing_back_low.addRotation(90,90,0);

            % translation to global coordination
            tri_wing_back_up.addTranslation(total_length,y_cut,0);
            tri_wing_back_low.addTranslation(total_length,y_cut,0);

            %% calculate wing
            % local coordination, local x is global y, local y is global -x
            shape_fcn_Z=@(XI,PSI) self.shapeWing(1-PSI);

            wing_up=CST3DSurface(wing_width,wing_length,wing_height_up,[],[],shape_fcn_Z,{[0,0],[0,0]});
            wing_low=CST3DSurface(wing_width,wing_length,-wing_height_low,[],[],shape_fcn_Z,{[0,0],[0,0]});

            % add rotation
            wing_up.addRotation(0,0,90);
            wing_low.addRotation(0,0,90);

            % add translation
            wing_up.addTranslation(total_length,y_cut+par_WS1,par_R);
            wing_low.addTranslation(total_length,y_cut+par_WS1,-par_R);

            %% calculate blunt wing front
            % local coordinate, local x is global -y, local y is global z, local z is global -x
            shape_fcn_Z=@(XI,PSI) (1-XI).*shape_wing_edge_x((PSI-0.5)*2*par_R)+XI.*shape_tri_wing_edge_x((PSI-0.5)*2*par_R);
            wing_front=CST3DSurface(wing_width,2*par_R,1,[],[],shape_fcn_Z,{[0,0]});

            % rotation surface
            wing_front.addRotation(90,-90,0);

            % translation surface
            wing_front.addTranslation(total_length-wing_length,y_cut+par_WS1+par_WS2,-par_R);

            %% calculate wing back
            % local coordination, local x is global y, local y is global z
            wing_back_up=CST3DSurface(wing_width,wing_height_up+par_R,0);
            wing_back_low=CST3DSurface(wing_width,-wing_height_low-par_R,0);

            % rotation to global coordinate
            wing_back_up.addRotation(90,90,0);
            wing_back_low.addRotation(90,90,0);

            % translation to global coordination
            wing_back_up.addTranslation(total_length,y_cut+par_WS1,0);
            wing_back_low.addTranslation(total_length,y_cut+par_WS1,0);

            %% calculate wing side
            % local coordinate, local x is global -x, local y is global z
            shape_fcn_Y_up=@(PSI) self.shapeWing(1-PSI)*par_WT_up+par_R;
            shape_fcn_Y_low=@(PSI) self.shapeWing(1-PSI)*par_WT_low+par_R;

            wing_side_up=CST3DSurface(wing_length,1,0,[],shape_fcn_Y_up);
            wing_side_low=CST3DSurface(wing_length,-1,0,[],shape_fcn_Y_low);

            % wing side blunt front
            shape_fcn_X=@(PSI) shape_wing_edge_x((PSI-0.5)*2*par_R);
            wing_side_front=CST3DSurface(1,2*par_R,0,shape_fcn_X);

            wing_side_up.addRotation(90,180,0);
            wing_side_low.addRotation(90,180,0);
            wing_side_front.addRotation(90,180,0);

            % translation to global coordinate
            wing_side_up.addTranslation(total_length,y_cut+par_WS1+par_WS2,0);
            wing_side_low.addTranslation(total_length,y_cut+par_WS1+par_WS2,0);
            wing_side_front.addTranslation(total_length-wing_length,y_cut+par_WS1+par_WS2,-par_R);

            %% sort data
            self.surface_list=struct('head_up',head_up,'head_low',head_low,'head_side',head_side,...
                'body_up',body_up,'body_low',body_low,'body_back_up',body_back_up,'body_back_low',body_back_low,...
                'tri_wing_up',tri_wing_up,'tri_wing_low',tri_wing_low,'tri_wing_front',tri_wing_front,'tri_wing_back_up',tri_wing_back_up,'tri_wing_back_low',tri_wing_back_low,...
                'wing_up',wing_up,'wing_low',wing_low,'wing_front',wing_front,'wing_back_up',wing_back_up,'wing_back_low',wing_back_low,...
                'wing_side_up',wing_side_up,'wing_side_low',wing_side_low,'wing_side_front',wing_side_front);

        end

        function [X,Y,Z]=calSurfaceMatrix...
                (self,xi_grid_num_head,psi_grid_num_head,xi_grid_num_body,psi_grid_num_wing,edge_gird_num)
            % generate waverider wing wgs data
            % xi, x=X(xi)
            % psi, y=Y(xi)
            % zeta, z=Z(xi,psi)
            %
            % notice colume vector in matrix is wgs format
            %

            % calculate head
            head_xi_list=[0,self.first_xi,self.second_xi,linspace(0.015,(1-self.par_rho1),xi_grid_num_head-2)];
            [XI,PSI]=meshgrid(head_xi_list,linspace(0.5,1,psi_grid_num_head+1));
            [X_head,Y_head,Z_head_up]=self.surface_list.head_up.calSurface(XI,PSI);
            [X_head,Y_head,Z_head_low]=self.surface_list.head_low.calSurface(XI,PSI);

            % calculate blunt head side
            head_xi_list=fliplr(1-head_xi_list/(1-self.par_rho1));
            [XI,ZETA]=meshgrid(head_xi_list,linspace(0,1,edge_gird_num*2+1));
            [X_head_side,Y_head_side,Z_head_side]=self.surface_list.head_side.calSurface(XI,ZETA);

            % calculate body
            xi_list=linspace((1-self.par_rho1),1,xi_grid_num_body+1);
            XI=repmat(xi_list,psi_grid_num_head+1,1);
            PSI=zeros(psi_grid_num_head+1,xi_grid_num_body+1);
            for xi_index=1:xi_grid_num_body+1
                psi_max=self.body_psi_max_fun(xi_list(xi_index));
                PSI(:,xi_index)=linspace(0.5,psi_max+0.5,psi_grid_num_head+1);
            end

            [X_body,Y_body,Z_body_up]=self.surface_list.body_up.calSurface(XI,PSI);
            [X_body,Y_body,Z_body_low]=self.surface_list.body_low.calSurface(XI,PSI);

            % calculate body back
            psi_max=self.body_psi_max_fun(1);
            psi_list=linspace(0,psi_max,psi_grid_num_head+1);
            [XI,ZETA]=meshgrid(psi_list,linspace(0,1,edge_gird_num+1));

            [X_body_back,Y_body_back,Z_body_back_up]=self.surface_list.body_back_up.calSurface(XI,ZETA);
            [X_body_back,Y_body_back,Z_body_back_low]=self.surface_list.body_back_low.calSurface(XI,ZETA);

            % calculate transition wing
            [X_tri_wing,Y_tri_wing,Z_tri_wing_up]=self.surface_list.tri_wing_up.calSurface(psi_grid_num_wing,xi_grid_num_body);
            [X_tri_wing,Y_tri_wing,Z_tri_wing_low]=self.surface_list.tri_wing_low.calSurface(psi_grid_num_wing,xi_grid_num_body);

            % calculate blunt tir wing front
            [X_tri_wing_front,Y_tri_wing_front,Z_tri_wing_front]=self.surface_list.tri_wing_front.calSurface(psi_grid_num_wing,2*edge_gird_num);

            % calculate tri wing back
            [Y_tri_wing_back,Y_tri_wing_back,Z_tri_wing_back_up]=self.surface_list.tri_wing_back_up.calSurface(psi_grid_num_wing,edge_gird_num);
            [X_tri_wing_back,Y_tri_wing_back,Z_tri_wing_back_low]=self.surface_list.tri_wing_back_low.calSurface(psi_grid_num_wing,edge_gird_num);

            % calculate wing
            [X_wing,Y_wing,Z_wing_up]=self.surface_list.wing_up.calSurface(psi_grid_num_wing,xi_grid_num_body);
            [X_wing,Y_wing,Z_wing_low]=self.surface_list.wing_low.calSurface(psi_grid_num_wing,xi_grid_num_body);

            % calculate blunt wing front
            [X_wing_front,Y_wing_front,Z_wing_front]=self.surface_list.wing_front.calSurface(psi_grid_num_wing,2*edge_gird_num);

            % calculate wing back
            [X_wing_back,Y_wing_back,Z_wing_back_up]=self.surface_list.wing_back_up.calSurface(psi_grid_num_wing,edge_gird_num);
            [X_wing_back,Y_wing_back,Z_wing_back_low]=self.surface_list.wing_back_low.calSurface(psi_grid_num_wing,edge_gird_num);

            % calculate wing side
            [X_wing_side,Y_wing_side,Z_wing_side_up]=self.surface_list.wing_side_up.calSurface(xi_grid_num_body,edge_gird_num);
            [X_wing_side,Y_wing_side,Z_wing_side_low]=self.surface_list.wing_side_low.calSurface(xi_grid_num_body,edge_gird_num);
            [X_wing_side_front,Y_wing_side_front,Z_wing_side_front]=self.surface_list.wing_side_front.calSurface(edge_gird_num,2*edge_gird_num);

            %% sort data
            X=struct();Y=struct();Z=struct();
            [X,Y,Z]=self.sortSurfaceUpLow(X,Y,Z,'head',X_head,Y_head,Z_head_up,Z_head_low);
            [X,Y,Z]=self.sortSurface(X,Y,Z,'head_side',X_head_side,Y_head_side,Z_head_side);
            [X,Y,Z]=self.sortSurfaceUpLow(X,Y,Z,'body',X_body,Y_body,Z_body_up,Z_body_low);
            [X,Y,Z]=self.sortSurfaceUpLow(X,Y,Z,'body_back',X_body_back,Y_body_back,Z_body_back_up,Z_body_back_low);
            [X,Y,Z]=self.sortSurfaceUpLow(X,Y,Z,'tri_wing',X_tri_wing,Y_tri_wing,Z_tri_wing_up,Z_tri_wing_low);
            [X,Y,Z]=self.sortSurface(X,Y,Z,'tri_wing_front',X_tri_wing_front,Y_tri_wing_front,Z_tri_wing_front);
            [X,Y,Z]=self.sortSurfaceUpLow(X,Y,Z,'tri_wing_back',X_tri_wing_back,Y_tri_wing_back,Z_tri_wing_back_up,Z_tri_wing_back_low);
            [X,Y,Z]=self.sortSurfaceUpLow(X,Y,Z,'wing',X_wing,Y_wing,Z_wing_up,Z_wing_low);
            [X,Y,Z]=self.sortSurface(X,Y,Z,'wing_front',X_wing_front,Y_wing_front,Z_wing_front);
            [X,Y,Z]=self.sortSurfaceUpLow(X,Y,Z,'wing_back',X_wing_back,Y_wing_back,Z_wing_back_up,Z_wing_back_low);
            [X,Y,Z]=self.sortSurfaceUpLow(X,Y,Z,'wing_side',X_wing_side,Y_wing_side,Z_wing_side_up,Z_wing_side_low);
            [X,Y,Z]=self.sortSurface(X,Y,Z,'wing_side_front',X_wing_side_front,Y_wing_side_front,Z_wing_side_front);

        end

        function PSI=shapeWing(self,XI)
            PSI=XI;
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

        function part=getWGSPart(self,part_name,...
                xi_grid_num_head,psi_grid_num_head,xi_grid_num_body,psi_grid_num_wing,edge_gird_num)
            surface_name_list=fieldnames(self.surface_list);
            surface_number=length(surface_name_list);

            [X_total,Y_total,Z_total]=self.calSurfaceMatrix...
                (xi_grid_num_head,psi_grid_num_head,xi_grid_num_body,psi_grid_num_wing,edge_gird_num);

            mesh_list=cell(surface_number,1);
            for surface_index=1:surface_number
                surface_name=surface_name_list{surface_index};
                if contains(surface_name,'up') || contains(surface_name,'front') || strcmp(surface_name,'head_side')
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
                xi_grid_num_head,psi_grid_num_head,xi_grid_num_body,psi_grid_num_wing,edge_gird_num)
            surface_name_list=fieldnames(self.surface_list);
            surface_number=length(surface_name_list);

            [X_total,Y_total,Z_total]=self.calSurfaceMatrix...
                (xi_grid_num_head,psi_grid_num_head,xi_grid_num_body,psi_grid_num_wing,edge_gird_num);

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