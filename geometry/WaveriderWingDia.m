classdef WaveriderWingDia < Body
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
        par_rho23;
        par_WS1;
        par_WS2;
        par_WT_up;
        par_WT_low;

        % calculate parameter
        head_length;
        head_side_length;
        stag_length;
        body_length
        y_cut;
        body_psi_max_fun;

        % local coordinate parameter
        blunt_xi;
        blunt_eta;
        stag_xi;

        % function
        shape_head_edge_x
        shape_tri_wing_edge_x
        shape_wing_edge_x

    end

    methods
        function self=WaveriderWingDia...
                (total_length,par_width,par_hight_up,par_hight_low,...
                par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
                par_rho1,par_rho12,par_rho23,par_WS1,par_WS2,par_WT_up,par_WT_low)
            % generate waverider wing wgs data
            % xi, x=X(xi)
            % psi, y=Y(xi)
            % zeta, z=Z(xi,psi)
            %
            % notice:
            % colume vector in matrix is wgs format
            % par_rho1 should between 0 and 1
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
            self.par_rho23=par_rho23;
            self.par_WS1=par_WS1;
            self.par_WS2=par_WS2;
            self.par_WT_up=par_WT_up;
            self.par_WT_low=par_WT_low;

            self.blunt_xi=0.002; % is local length parameter
            self.blunt_eta=0.005; % is local length parameter
            self.stag_xi=0.010; % is local length parameter
            
            stag_length=(self.stag_xi*(1-par_rho1))^par_T*par_width/2;
            self.stag_length=stag_length;

            head_length=(1-par_rho1)*total_length;
            self.head_length=head_length;

            head_side_length=head_length*(1-self.stag_xi);
            self.head_side_length=head_side_length;

            body_length=par_rho1*total_length;
            self.body_length=body_length;

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
                dz_dx=(par_hight_up*(self.blunt_xi)^par_M_up)/(total_length*self.blunt_xi);
                radius_center_up=par_R*dz_dx;
            end
            if par_M_low > 1
                radius_center_low=0;
            else
                dz_dx=(par_hight_low*(self.blunt_xi)^par_M_low)/(total_length*self.blunt_xi);
                radius_center_low=par_R*dz_dx;
            end
            radius_center_head=max(radius_center_up,radius_center_low);
            radius_head_sq=(radius_center_head*radius_center_head+par_R*par_R);
            radius_head=sqrt(radius_head_sq);
            shape_head_edge_x=@(Y) (sqrt(radius_head_sq-Y.^2))-radius_center_head;
            self.shape_head_edge_x=shape_head_edge_x;

            % calculate blunt radius of zox plane of tri wing
            dz_dx=(par_M_up*(1-par_rho1)^(par_M_up-1)/total_length);
            radius_center_up=par_R*dz_dx;
            dz_dx=(par_M_low*(1-par_rho1)^(par_M_low-1)/total_length);
            radius_center_low=par_R*dz_dx;
            radius_center_tri_wing=max(radius_center_up,radius_center_low);
            radius_tri_wing_sq=(radius_center_tri_wing*radius_center_tri_wing+par_R*par_R);
            shape_tri_wing_edge_x=@(Y) sqrt(radius_tri_wing_sq-Y.^2)-radius_center_tri_wing;
            self.shape_tri_wing_edge_x=shape_tri_wing_edge_x;

            % calculate blunt radius and center of zox plane of wing
            dz_dx=par_WT_up/wing_length*2;
            radius_center_up=par_R*dz_dx;
            dz_dx=par_WT_low/wing_length*2;
            radius_center_low=par_R*dz_dx;
            radius_center_wing=max(radius_center_up,radius_center_low);
            radius_wing_sq=(radius_center_wing*radius_center_wing+par_R*par_R);
            radius_wing=sqrt(radius_wing_sq);
            shape_wing_edge_x=@(Y) (sqrt(radius_wing_sq-Y.^2))-radius_center_wing;
            self.shape_wing_edge_x=shape_wing_edge_x;

            % calculate head side blunt radius of yoz plane
            dz_dx=(self.blunt_eta)^par_N_up*(1-self.blunt_eta)^par_N_up/(0.5)^(2*par_N_up)*par_hight_up*(1-par_rho1)^par_M_up/(2*self.blunt_eta*y_cut);
            radius_center_up=par_R*dz_dx;
            dz_dx=(self.blunt_eta)^par_N_low*(1-self.blunt_eta)^par_N_low/(0.5)^(2*par_N_low)*par_hight_low*(1-par_rho1)^par_M_low/(2*self.blunt_eta*y_cut);
            radius_center_low=par_R*dz_dx;
            radius_center_head_side=max(radius_center_low,radius_center_up);
            radius_head_side_sq=(radius_center_head_side*radius_center_head_side+par_R*par_R);
            radius_head_side=sqrt(radius_head_side_sq);

            %% define head
            shape_fcn_Y=@(XI) (XI*(1-par_rho1)).^par_T;
            shape_fcn_Z_up=@(XI,PSI) (XI*(1-par_rho1)).^par_M_up;
            shape_fcn_Z_low=@(XI,PSI) (XI*(1-par_rho1)).^par_M_low;
            head_up=SurfaceCST3D(head_length,par_width,par_hight_up,[],shape_fcn_Y,shape_fcn_Z_up,{[par_N_up,par_N_up]},global_symmetry_y);
            head_low=SurfaceCST3D(head_length,par_width,-par_hight_low,[],shape_fcn_Y,shape_fcn_Z_low,{[par_N_low,par_N_low]},global_symmetry_y);
            
            % blunt translation
            head_up.addTranslation(0,0,par_R);
            head_low.addTranslation(0,0,-par_R);

            %% define blunt head side
            % local coordinate, local x is global -x, local y is global z, local z is global y
            shape_fcn_X=@(PSI) 1-shape_tri_wing_edge_x((PSI-0.5)*par_R*2)/head_length+shape_head_edge_x((PSI-0.5)*par_R*2)/head_length;
            shape_fcn_Z=@(XI,PSI) (sqrt(radius_head_side_sq-(((PSI-0.5)*2)*par_R).^2)-radius_center_head_side).*(1-XI).^(0.01);
            head_side=SurfaceCST3D(head_length,2*par_R,1,shape_fcn_X,[],shape_fcn_Z,{[0.0,0.0],[0.0,0.0]});

            tran_fun_X=@(PSI) shape_tri_wing_edge_x((PSI-0.5)*par_R*2);
            Y_edge=@(XI) XI.^par_T*par_width/2;
            tran_fun_Z=@(XI,PSI) Y_edge((1-XI)*(1-par_rho1));
            % deform surface
            head_side.addDeform(tran_fun_X,[],tran_fun_Z);

            % rotation surface
            head_side.addRotation(90,180,0);

            % translation surface
            head_side.addTranslation(head_length,0,-par_R);

            %% define body
            shape_fcn_Z_up=@(XI,PSI) (XI*par_rho1+(1-par_rho1)).^par_M_up;
            class_fcn_Z_up=@(XI,PSI) ((PSI-0.5)*2.*body_psi_max_fun(XI*par_rho1+(1-par_rho1))+0.5).^par_N_up...
                .*(1-((PSI-0.5)*2.*body_psi_max_fun(XI*par_rho1+(1-par_rho1))+0.5)).^par_N_up/(0.5)^(2*par_N_up);
            shape_fcn_Z_low=@(XI,PSI) (XI*par_rho1+(1-par_rho1)).^par_M_low;
            class_fcn_Z_low=@(XI,PSI) ((PSI-0.5)*2.*body_psi_max_fun(XI*par_rho1+(1-par_rho1))+0.5).^par_N_low...
                .*(1-((PSI-0.5)*2.*body_psi_max_fun(XI*par_rho1+(1-par_rho1))+0.5)).^par_N_low/(0.5)^(2*par_N_low);
            body_up=SurfaceCST3D(body_length,y_cut*2,par_hight_up,[],[],shape_fcn_Z_up,class_fcn_Z_up,global_symmetry_y);
            body_low=SurfaceCST3D(body_length,y_cut*2,-par_hight_low,[],[],shape_fcn_Z_low,class_fcn_Z_low,global_symmetry_y);

            % blunt translation
            body_up.addTranslation(head_length,0,par_R);
            body_low.addTranslation(head_length,0,-par_R);

            %% define body back
            shape_fcn_Y_up=@(XI) par_hight_up*(0.5+XI*(psi_end_cut)).^par_N_up.*(0.5-XI*(psi_end_cut)).^par_N_up/(0.5)^(2*par_N_up)+par_R;
            shape_fcn_Y_low=@(XI) par_hight_low*(0.5+XI*(psi_end_cut)).^par_N_low.*(0.5-XI*(psi_end_cut)).^par_N_low/(0.5)^(2*par_N_low)+par_R;

            % local coordinate, local x is global y, local y is global z
            body_back_up=SurfaceCST3D(y_cut,1,0,[],shape_fcn_Y_up);
            body_back_low=SurfaceCST3D(y_cut,-1,0,[],shape_fcn_Y_low);

            % rotation to global coordinate
            body_back_up.addRotation(90,90,0);
            body_back_low.addRotation(90,90,0);

            % translation to global coordination
            body_back_up.addTranslation(total_length,0,0);
            body_back_low.addTranslation(total_length,0,0);

            %% define transition wing
            % local coordination, local x is global y, local y is global -x, local z is global z
            shape_fcn_Y=@(XI) (1-XI*(1-par_rho12));
            shape_body_size_up=@(PSI) body_up.LZ.*...
                body_up.shape_fcn_Z(1-PSI,1).*...
                body_up.class_fcn_Z(1-PSI,1);
            shape_body_size_low=@(PSI) -body_low.LZ.*...
                body_low.shape_fcn_Z(1-PSI,1).*...
                body_low.class_fcn_Z(1-PSI,1);
            shape_fcn_Z_up=@(XI,PSI) (1-XI).*shape_body_size_up(PSI)+XI.*self.shapeWing(1-PSI)*(wing_height_up);
            shape_fcn_Z_low=@(XI,PSI) (1-XI).*shape_body_size_low(PSI)+XI.*self.shapeWing(1-PSI)*(wing_height_low);

            tri_wing_up=SurfaceCST3D(tri_wing_width,tri_wing_length,1,[],shape_fcn_Y,shape_fcn_Z_up,{[0,0]});
            tri_wing_low=SurfaceCST3D(tri_wing_width,tri_wing_length,-1,[],shape_fcn_Y,shape_fcn_Z_low,{[0,0]});

            % add rotation
            tri_wing_up.addRotation(0,0,90);
            tri_wing_low.addRotation(0,0,90);

            % add translation
            tri_wing_up.addTranslation(total_length,y_cut,par_R);
            tri_wing_low.addTranslation(total_length,y_cut,-par_R);

            %% define blunt tir wing front
            % local coordinate, local x is global -y, local y is global z, local z is global -x
            shape_fcn_X=@(PSI) 1-(sqrt(radius_head_side_sq-(((PSI-0.5)*2)*par_R).^2)-radius_center_head_side)/tri_wing_width; % if do not want blunt intersection, remove this and change blunt head side class N1
            shape_fcn_Z=@(XI,PSI) (1-XI).*shape_wing_edge_x((PSI-0.5)*2*par_R)+XI.*shape_tri_wing_edge_x((PSI-0.5)*2*par_R);
            tri_wing_front=SurfaceCST3D(tri_wing_width,2*par_R,1,shape_fcn_X,[],[],shape_fcn_Z);

            % deform surface
            tran_fun_Z=@(XI,PSI) XI*(par_rho1*(1-par_rho12)*total_length);
            tri_wing_front.addDeform([],[],tran_fun_Z);

            % rotation surface
            tri_wing_front.addRotation(90,-90,0);

            % translation surface
            tri_wing_front.addTranslation(total_length-wing_length,y_cut+par_WS1,-par_R);

            %% define tri wing back
            % local coordination, local x is global y, local y is global z
            shape_fcn_Y_up=@(XI) (1-XI)*(tri_wing_height_up+par_R)+XI*(par_R);
            shape_fcn_Y_low=@(XI) (1-XI)*(tri_wing_height_low+par_R)+XI*(par_R);

            tri_wing_back_up=SurfaceCST3D(tri_wing_width,1,0,[],shape_fcn_Y_up);
            tri_wing_back_low=SurfaceCST3D(tri_wing_width,-1,0,[],shape_fcn_Y_low);

            % rotation to global coordinate
            tri_wing_back_up.addRotation(90,90,0);
            tri_wing_back_low.addRotation(90,90,0);

            % translation to global coordination
            tri_wing_back_up.addTranslation(total_length,y_cut,0);
            tri_wing_back_low.addTranslation(total_length,y_cut,0);

            %% define wing
            % local coordination, local x is global y, local y is global -x
            shape_fcn_Y=@(XI) (1-XI*(1-par_rho23));
            shape_fcn_Z=@(XI,PSI) self.shapeWing(1-PSI).*(1-XI*(1-par_rho23));

            wing_up=SurfaceCST3D(wing_width,wing_length,wing_height_up,[],shape_fcn_Y,shape_fcn_Z,{[0,0],[0,0]});
            wing_low=SurfaceCST3D(wing_width,wing_length,-wing_height_low,[],shape_fcn_Y,shape_fcn_Z,{[0,0],[0,0]});

            % add rotation
            wing_up.addRotation(0,0,90);
            wing_low.addRotation(0,0,90);

            % add translation
            wing_up.addTranslation(total_length,y_cut+par_WS1,par_R);
            wing_low.addTranslation(total_length,y_cut+par_WS1,-par_R);

            %% define blunt wing front
            % local coordinate, local x is global -y, local y is global z, local z is global -x
            shape_fcn_Z=@(XI,PSI) shape_wing_edge_x((PSI-0.5)*2*par_R);
            wing_front=SurfaceCST3D(wing_width,2*par_R,1,[],[],shape_fcn_Z,{[0,0]});

            % deform surface
            tran_fun_Z=@(XI,PSI) XI*(par_rho1*par_rho12*(1-par_rho23)*total_length);
            wing_front.addDeform([],[],tran_fun_Z);

            % rotation surface
            wing_front.addRotation(90,-90,0);

            % translation surface
            wing_front.addTranslation(total_length-wing_length*par_rho23,y_cut+par_WS1+par_WS2,-par_R);

            %% define wing back
            % local coordination, local x is global y, local y is global z
            wing_back_up=SurfaceCST3D(wing_width,par_R,0);
            wing_back_low=SurfaceCST3D(wing_width,-par_R,0);

            % rotation to global coordinate
            wing_back_up.addRotation(90,90,0);
            wing_back_low.addRotation(90,90,0);

            % translation to global coordination
            wing_back_up.addTranslation(total_length,y_cut+par_WS1,0);
            wing_back_low.addTranslation(total_length,y_cut+par_WS1,0);

            %% define wing side
            % local coordinate, local x is global -x, local y is global z
            shape_fcn_Y_up=@(PSI) self.shapeWing(1-PSI)*par_WT_up*par_rho23+par_R;
            shape_fcn_Y_low=@(PSI) self.shapeWing(1-PSI)*par_WT_low*par_rho23+par_R;

            wing_side_up=SurfaceCST3D(wing_length*par_rho23,1,0,[],shape_fcn_Y_up);
            wing_side_low=SurfaceCST3D(wing_length*par_rho23,-1,0,[],shape_fcn_Y_low);

            % wing side blunt front
            shape_fcn_X=@(PSI) shape_wing_edge_x((PSI-0.5)*2*par_R);
            wing_side_front=SurfaceCST3D(1,2*par_R,0,shape_fcn_X);

            wing_side_up.addRotation(90,180,0);
            wing_side_low.addRotation(90,180,0);
            wing_side_front.addRotation(90,180,0);

            % translation to global coordinate
            wing_side_up.addTranslation(total_length,y_cut+par_WS1+par_WS2,0);
            wing_side_low.addTranslation(total_length,y_cut+par_WS1+par_WS2,0);
            wing_side_front.addTranslation(total_length-wing_length*par_rho23,y_cut+par_WS1+par_WS2,-par_R);

            %% sort data
            self.surface_list=struct('head_up',head_up,'head_low',head_low,'head_side',head_side,...
                'body_up',body_up,'body_low',body_low,'body_back_up',body_back_up,'body_back_low',body_back_low,...
                'tri_wing_up',tri_wing_up,'tri_wing_low',tri_wing_low,'tri_wing_front',tri_wing_front,'tri_wing_back_up',tri_wing_back_up,'tri_wing_back_low',tri_wing_back_low,...
                'wing_up',wing_up,'wing_low',wing_low,'wing_front',wing_front,'wing_back_up',wing_back_up,'wing_back_low',wing_back_low,...
                'wing_side_up',wing_side_up,'wing_side_low',wing_side_low,'wing_side_front',wing_side_front);

        end
        
        function PSI=shapeWing(self,XI)
            % height normalize to 1
            %
            PSI=1-abs(1-2*XI);
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
            [X_head,Y_head,Z_head_up]=self.surface_list.head_up.calSurface(xi_grid_num_head,psi_grid_num_head);
            [X_head,Y_head,Z_head_low]=self.surface_list.head_low.calSurface(xi_grid_num_head,psi_grid_num_head);

            % calculate blunt head side
            [X_head_side,Y_head_side,Z_head_side]=self.surface_list.head_side.calSurface(xi_grid_num_head,edge_gird_num*2);

            % calculate body
            [X_body,Y_body,Z_body_up]=self.surface_list.body_up.calSurface(xi_grid_num_body,psi_grid_num_head);
            [X_body,Y_body,Z_body_low]=self.surface_list.body_low.calSurface(xi_grid_num_body,psi_grid_num_head);

            % calculate body back
            [X_body_back,Y_body_back,Z_body_back_up]=self.surface_list.body_back_up.calSurface(psi_grid_num_head,edge_gird_num);
            [X_body_back,Y_body_back,Z_body_back_low]=self.surface_list.body_back_low.calSurface(psi_grid_num_head,edge_gird_num);

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
            [X,Y,Z]=sortSurfaceUpLow(X,Y,Z,'head',X_head,Y_head,Z_head_up,Z_head_low);
            [X,Y,Z]=sortSurface(X,Y,Z,'head_side',X_head_side,Y_head_side,Z_head_side);
            [X,Y,Z]=sortSurfaceUpLow(X,Y,Z,'body',X_body,Y_body,Z_body_up,Z_body_low);
            [X,Y,Z]=sortSurfaceUpLow(X,Y,Z,'body_back',X_body_back,Y_body_back,Z_body_back_up,Z_body_back_low);
            [X,Y,Z]=sortSurfaceUpLow(X,Y,Z,'tri_wing',X_tri_wing,Y_tri_wing,Z_tri_wing_up,Z_tri_wing_low);
            [X,Y,Z]=sortSurface(X,Y,Z,'tri_wing_front',X_tri_wing_front,Y_tri_wing_front,Z_tri_wing_front);
            [X,Y,Z]=sortSurfaceUpLow(X,Y,Z,'tri_wing_back',X_tri_wing_back,Y_tri_wing_back,Z_tri_wing_back_up,Z_tri_wing_back_low);
            [X,Y,Z]=sortSurfaceUpLow(X,Y,Z,'wing',X_wing,Y_wing,Z_wing_up,Z_wing_low);
            [X,Y,Z]=sortSurface(X,Y,Z,'wing_front',X_wing_front,Y_wing_front,Z_wing_front);
            [X,Y,Z]=sortSurfaceUpLow(X,Y,Z,'wing_back',X_wing_back,Y_wing_back,Z_wing_back_up,Z_wing_back_low);
            [X,Y,Z]=sortSurfaceUpLow(X,Y,Z,'wing_side',X_wing_side,Y_wing_side,Z_wing_side_up,Z_wing_side_low);
            [X,Y,Z]=sortSurface(X,Y,Z,'wing_side_front',X_wing_side_front,Y_wing_side_front,Z_wing_side_front);

        end

    end

    methods    
        function part=getWGSMesh(self,part_name,...
                xi_grid_num_head,psi_grid_num_head,xi_grid_num_body,psi_grid_num_wing,edge_gird_num)
            surface_name_list=fieldnames(self.surface_list);
            surface_number=length(surface_name_list);

            [X_total,Y_total,Z_total]=self.calSurfaceMatrix...
                (xi_grid_num_head,psi_grid_num_head,xi_grid_num_body,psi_grid_num_wing,edge_gird_num);

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

        function writeStepShell(self,step_filestr)
            % write surface into step file
            %
            write_surface_list=self.surface_list;
            surface_name_list=fieldnames(write_surface_list);
            surf_num=length(surface_name_list);
            
            % check file name
            if length(step_filestr) > 4
                if ~strcmpi(step_filestr((end-3):end),'.step')
                    step_filestr=[step_filestr,'.step'];
                end
            else
                step_filestr=[step_filestr,'.step'];
            end
            [~,step_filename,~]=fileparts(step_filestr);

            % write head
            step_file=fopen(step_filestr,'w');
            fprintf(step_file,'ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION (( ''STEP AP203'' ),''1'' );\nFILE_NAME (''%s'',''%s'',( '''' ),( '''' ),''Matlab step'',''Matlab'','''' );\nFILE_SCHEMA (( ''CONFIG_CONTROL_DESIGN'' ));\nENDSEC;\n',step_filestr,date);
            fprintf(step_file,'\n');
            fprintf(step_file,'DATA;\n');
            object_index=1;
            fprintf(step_file,'#%d = CARTESIAN_POINT ( ''NONE'',  ( 0.000000000000000000, 0.000000000000000000, 0.000000000000000000 ) ) ;\n',object_index);object_index=object_index+1;
            fprintf(step_file,'#%d = DIRECTION ( ''NONE'',  ( 0.000000000000000000, 0.000000000000000000, 1.000000000000000000 ) ) ;\n',object_index);object_index=object_index+1;
            fprintf(step_file,'#%d = DIRECTION ( ''NONE'',  ( 1.000000000000000000, 0.000000000000000000, 0.000000000000000000 ) ) ;\n',object_index);object_index=object_index+1;
            fprintf(step_file,'#%d = AXIS2_PLACEMENT_3D ( ''NONE'', #%d, #%d, #%d ) ;\n',object_index,1,2,3);object_index=object_index+1;
            fprintf(step_file,'\n');
            fprintf(step_file,'#%d =( LENGTH_UNIT ( ) NAMED_UNIT ( * ) SI_UNIT ( $., .METRE. ) );\n',object_index);object_index=object_index+1;
            fprintf(step_file,'#%d =( NAMED_UNIT ( * ) PLANE_ANGLE_UNIT ( ) SI_UNIT ( $, .RADIAN. ) );\n',object_index);object_index=object_index+1;
            fprintf(step_file,'#%d =( NAMED_UNIT ( * ) SI_UNIT ( $, .STERADIAN. ) SOLID_ANGLE_UNIT ( ) );\n',object_index);object_index=object_index+1;
            fprintf(step_file,'#%d = UNCERTAINTY_MEASURE_WITH_UNIT (LENGTH_MEASURE( 1.00000000000000000E-05 ), #%d, ''distance_accuracy_value'', ''NONE'');\n',object_index,5);object_index=object_index+1;
            fprintf(step_file,'#%d =( GEOMETRIC_REPRESENTATION_CONTEXT ( 3 ) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT ( ( #%d ) ) GLOBAL_UNIT_ASSIGNED_CONTEXT ( ( #%d, #%d, #%d ) ) REPRESENTATION_CONTEXT ( ''NONE'', ''WORKASPACE'' ) );\n',object_index,8,5,6,7);object_index=object_index+1;
            fprintf(step_file,'\n');

            % write surface
            OPEN_SHELL_index_list=zeros(1,surf_num);
            surf_idx=1;

            % write head_up
            surf=write_surface_list.('head_up');
            [X,Y,Z,XI,PSI]=surf.calSurface(1e-5);
            u_list=[0,0,0,XI(1,:),0,0,0];
            v_list=[0,0,0,PSI(:,1)',0,0,0];
            surf=SurfaceBSpline([],[],[],X,Y,Z,3,u_list,v_list);

            [step_str,object_index,OPEN_SHELL_index]=surf.getStepNode(object_index);
            fprintf(step_file,step_str);

            OPEN_SHELL_index_list(surf_idx)=OPEN_SHELL_index;
            surf_idx=surf_idx+1;
            write_surface_list=rmfield(write_surface_list,'head_up');

            % write head_low
            surf=write_surface_list.('head_low');
            [X,Y,Z,XI,PSI]=surf.calSurface(XI,PSI);
            u_list=[0,0,0,XI(1,:),0,0,0];
            v_list=[0,0,0,PSI(:,1)',0,0,0];
            surf=SurfaceBSpline([],[],[],X,Y,Z,3,u_list,v_list);

            [step_str,object_index,OPEN_SHELL_index]=surf.getStepNode(object_index);
            fprintf(step_file,step_str);

            OPEN_SHELL_index_list(surf_idx)=OPEN_SHELL_index;
            surf_idx=surf_idx+1;
            write_surface_list=rmfield(write_surface_list,'head_low');

            % write head_side
            surf=write_surface_list.('head_side');
            [X,Y,Z,XI,PSI]=surf.calSurface(1-XI,20);
            u_list=[0,0,0,XI(1,:),0,0,0];
            v_list=[0,0,0,PSI(:,1)',0,0,0];
            surf=SurfaceBSpline([],[],[],X,Y,Z,3,u_list,v_list);

            [step_str,object_index,OPEN_SHELL_index]=surf.getStepNode(object_index);
            fprintf(step_file,step_str);

            OPEN_SHELL_index_list(surf_idx)=OPEN_SHELL_index;
            surf_idx=surf_idx+1;
            write_surface_list=rmfield(write_surface_list,'head_side');

            write_surface_list=struct2cell(write_surface_list);

            for surf_idx=1:length(write_surface_list)
                surf=write_surface_list{surf_idx};
                surf=surf.getBSpline(20,20);

                [step_str,object_index,OPEN_SHELL_index]=surf.getStepNode(object_index);
                fprintf(step_file,step_str);

                OPEN_SHELL_index_list(surf_idx+3)=OPEN_SHELL_index;
            end

            % write model
            SHELL_BASED_SURFACE_MODEL_index=object_index;
            step_str=[num2str(object_index,'#%d'),' = SHELL_BASED_SURFACE_MODEL ',...
                '( ''NONE'', ( ',num2str(OPEN_SHELL_index_list(1:end-1),'#%d, '),' ',num2str(OPEN_SHELL_index_list(end),'#%d'),' ) );\n'];object_index=object_index+1;
            fprintf(step_file,step_str);
            fprintf(step_file,'\n');

            % write product context 
            fprintf(step_file,'#%d = APPLICATION_CONTEXT ( ''configuration controlled 3d designs of mechanical parts and assemblies'' ) ;\n',object_index);object_index=object_index+1;
            fprintf(step_file,'#%d = MECHANICAL_CONTEXT ( ''NONE'', #%d, ''mechanical'' ) ;\n',object_index,object_index-1);object_index=object_index+1;
            fprintf(step_file,'#%d = PRODUCT ( ''%s'', ''%s'', '''', ( #%d ) ) ;\n',object_index,step_filename,step_filename,object_index-1);object_index=object_index+1;
            fprintf(step_file,'\n');
            fprintf(step_file,'#%d = APPLICATION_CONTEXT ( ''configuration controlled 3d designs of mechanical parts and assemblies'' ) ;\n',object_index);object_index=object_index+1;
            fprintf(step_file,'#%d = DESIGN_CONTEXT ( ''detailed design'', #%d, ''design'' ) ;\n',object_index,object_index-1);object_index=object_index+1;
            fprintf(step_file,'#%d = PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE ( ''NONE'', '''', #%d, .NOT_KNOWN. ) ;\n',object_index,object_index-3);object_index=object_index+1;
            fprintf(step_file,'#%d = PRODUCT_DEFINITION ( ''NONE'', '''', #%d, #%d ) ;\n',object_index,object_index-1,object_index-2);object_index=object_index+1;
            fprintf(step_file,'\n');
            fprintf(step_file,'#%d = PRODUCT_DEFINITION_SHAPE ( ''NONE'', ''NONE'',  #%d ) ;\n',object_index,object_index-1);object_index=object_index+1;
            fprintf(step_file,'#%d = MANIFOLD_SURFACE_SHAPE_REPRESENTATION ( ''test'', ( #%d, #%d ), #%d ) ;\n',object_index,SHELL_BASED_SURFACE_MODEL_index,4,9);object_index=object_index+1;
            fprintf(step_file,'#%d = SHAPE_DEFINITION_REPRESENTATION ( #%d, #%d ) ;\n',object_index,object_index-2,object_index-1);object_index=object_index+1;

            % write end
            fprintf(step_file,'\n');
            fprintf(step_file,'ENDSEC;\nEND-ISO-10303-21;\n');
            fclose(step_file);
            clear('step_file');

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

function [X,Y,Z]=sortSurface(X,Y,Z,name,X_surf,Y_surf,Z_surf)
X.(name)=X_surf;
Y.(name)=Y_surf;
Z.(name)=Z_surf;
end

function [X,Y,Z]=sortSurfaceUpLow(X,Y,Z,name,X_surf,Y_surf,Z_surf_up,Z_surf_low)
X.([name,'_up'])=X_surf;
Y.([name,'_up'])=Y_surf;
Z.([name,'_up'])=Z_surf_up;
X.([name,'_low'])=X_surf;
Y.([name,'_low'])=Y_surf;
Z.([name,'_low'])=Z_surf_low;
end

