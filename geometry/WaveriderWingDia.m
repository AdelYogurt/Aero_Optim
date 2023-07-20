classdef WaveriderWingDia < Body
    % parameterized body of waverider with dia cross section wing
    %
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
    end

    properties
        % calculate parameter
        head_length;
        head_side_length;
        stag_length;
        body_length
        y_cut;
        body_psi_max_fun;

        % local coordinate parameter
        blunt_xi=0.002; % is local length parameter
        blunt_eta=0.005; % is local length parameter

        % function
        shape_head_edge_x
        shape_tri_wing_edge_x
        shape_wing_edge_x
    end

    % define function
    methods
        function self=WaveriderWingDia...
                (total_length,varargin)
            % generate waverider wing wgs data
            % xi, x=X(xi)
            % psi, y=Y(xi)
            % zeta, z=Z(xi,psi)
            %
            % notice:
            % colume vector in matrix is wgs format
            % par_rho1 should between 0 and 1
            %
            if nargin == 2
                par_input=varargin{1};
            else
                par_input=varargin;
            end

            [par_width,par_hight_up,par_hight_low,...
                par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
                par_rho1,par_rho12,par_rho23,par_WS1,par_WS2,par_WT_up,par_WT_low]=WaveriderWingDia.decode(par_input);

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

            head_length=(1-par_rho1)*total_length;
            self.head_length=head_length;

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
            head_up=SurfaceCST3D('head_up',head_length,par_width,par_hight_up,[],shape_fcn_Y,shape_fcn_Z_up,[par_N_up,par_N_up],global_symmetry_y);
            head_low=SurfaceCST3D('head_low',head_length,par_width,-par_hight_low,[],shape_fcn_Y,shape_fcn_Z_low,[par_N_low,par_N_low],global_symmetry_y);
            
            % blunt translation
            head_up.addTranslation(0,0,par_R);
            head_low.addTranslation(0,0,-par_R);

            %% define blunt head side
            % local coordinate, local x is global -x, local y is global z, local z is global y
            shape_fcn_X=@(PSI) 1-shape_tri_wing_edge_x((PSI-0.5)*par_R*2)/head_length+shape_head_edge_x((PSI-0.5)*par_R*2)/head_length;
            shape_fcn_Z=@(XI,PSI) (sqrt(radius_head_side_sq-(((PSI-0.5)*2)*par_R).^2)-radius_center_head_side).*(1-XI);
            head_side=SurfaceCST3D('head_side',head_length,2*par_R,1,shape_fcn_X,[],shape_fcn_Z,{[0.0,0.0],[0.0,0.0]});

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
            body_up=SurfaceCST3D('body_up',body_length,y_cut*2,par_hight_up,[],[],shape_fcn_Z_up,class_fcn_Z_up,global_symmetry_y);
            body_low=SurfaceCST3D('body_low',body_length,y_cut*2,-par_hight_low,[],[],shape_fcn_Z_low,class_fcn_Z_low,global_symmetry_y);

            % blunt translation
            body_up.addTranslation(head_length,0,par_R);
            body_low.addTranslation(head_length,0,-par_R);

            %% define body back
            shape_fcn_Y_up=@(XI) par_hight_up*(0.5+XI*(psi_end_cut)).^par_N_up.*(0.5-XI*(psi_end_cut)).^par_N_up/(0.5)^(2*par_N_up)+par_R;
            shape_fcn_Y_low=@(XI) par_hight_low*(0.5+XI*(psi_end_cut)).^par_N_low.*(0.5-XI*(psi_end_cut)).^par_N_low/(0.5)^(2*par_N_low)+par_R;

            % local coordinate, local x is global y, local y is global z
            body_back_up=SurfaceCST3D('body_back_up',y_cut,1,0,[],shape_fcn_Y_up);
            body_back_low=SurfaceCST3D('body_back_low',y_cut,-1,0,[],shape_fcn_Y_low);

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

            tri_wing_up=SurfaceCST3D('tri_wing_up',tri_wing_width,tri_wing_length,1,[],shape_fcn_Y,shape_fcn_Z_up,{[0,0]});
            tri_wing_low=SurfaceCST3D('tri_wing_low',tri_wing_width,tri_wing_length,-1,[],shape_fcn_Y,shape_fcn_Z_low,{[0,0]});

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
            tri_wing_front=SurfaceCST3D('tri_wing_front',tri_wing_width,2*par_R,1,shape_fcn_X,[],[],shape_fcn_Z);

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

            tri_wing_back_up=SurfaceCST3D('tri_wing_back_up',tri_wing_width,1,0,[],shape_fcn_Y_up);
            tri_wing_back_low=SurfaceCST3D('tri_wing_back_low',tri_wing_width,-1,0,[],shape_fcn_Y_low);

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

            wing_up=SurfaceCST3D('wing_up',wing_width,wing_length,wing_height_up,[],shape_fcn_Y,shape_fcn_Z,{[0,0],[0,0]});
            wing_low=SurfaceCST3D('wing_low',wing_width,wing_length,-wing_height_low,[],shape_fcn_Y,shape_fcn_Z,{[0,0],[0,0]});

            % add rotation
            wing_up.addRotation(0,0,90);
            wing_low.addRotation(0,0,90);

            % add translation
            wing_up.addTranslation(total_length,y_cut+par_WS1,par_R);
            wing_low.addTranslation(total_length,y_cut+par_WS1,-par_R);

            %% define blunt wing front
            % local coordinate, local x is global -y, local y is global z, local z is global -x
            shape_fcn_Z=@(XI,PSI) shape_wing_edge_x((PSI-0.5)*2*par_R);
            wing_front=SurfaceCST3D('wing_front',wing_width,2*par_R,1,[],[],shape_fcn_Z,{[0,0]});

            % deform surface
            tran_fun_Z=@(XI,PSI) XI*(par_rho1*par_rho12*(1-par_rho23)*total_length);
            wing_front.addDeform([],[],tran_fun_Z);

            % rotation surface
            wing_front.addRotation(90,-90,0);

            % translation surface
            wing_front.addTranslation(total_length-wing_length*par_rho23,y_cut+par_WS1+par_WS2,-par_R);

            %% define wing back
            % local coordination, local x is global y, local y is global z
            wing_back_up=SurfaceCST3D('wing_back_up',wing_width,par_R,0);
            wing_back_low=SurfaceCST3D('wing_back_low',wing_width,-par_R,0);

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

            wing_side_up=SurfaceCST3D('wing_side_up',wing_length*par_rho23,1,0,[],shape_fcn_Y_up);
            wing_side_low=SurfaceCST3D('wing_side_low',wing_length*par_rho23,-1,0,[],shape_fcn_Y_low);

            % wing side blunt front
            shape_fcn_X=@(PSI) shape_wing_edge_x((PSI-0.5)*2*par_R);
            wing_side_front=SurfaceCST3D('wing_side_front',1,2*par_R,0,shape_fcn_X);

            wing_side_up.addRotation(90,180,0);
            wing_side_low.addRotation(90,180,0);
            wing_side_front.addRotation(90,180,0);

            % translation to global coordinate
            wing_side_up.addTranslation(total_length,y_cut+par_WS1+par_WS2,0);
            wing_side_low.addTranslation(total_length,y_cut+par_WS1+par_WS2,0);
            wing_side_front.addTranslation(total_length-wing_length*par_rho23,y_cut+par_WS1+par_WS2,-par_R);

            %% sort data
            %             self.surface_list=struct('head_up',head_up,'head_low',head_low,'head_side',head_side,...
            %                 'body_up',body_up,'body_low',body_low,'body_back_up',body_back_up,'body_back_low',body_back_low,...
            %                 'tri_wing_up',tri_wing_up,'tri_wing_low',tri_wing_low,'tri_wing_front',tri_wing_front,'tri_wing_back_up',tri_wing_back_up,'tri_wing_back_low',tri_wing_back_low,...
            %                 'wing_up',wing_up,'wing_low',wing_low,'wing_front',wing_front,'wing_back_up',wing_back_up,'wing_back_low',wing_back_low,...
            %                 'wing_side_up',wing_side_up,'wing_side_low',wing_side_low,'wing_side_front',wing_side_front);

            self.surface_list={head_up,head_low,head_side,...
                body_up,body_low,body_back_up,body_back_low,...
                tri_wing_up,tri_wing_low,tri_wing_front,tri_wing_back_up,tri_wing_back_low,...
                wing_up,wing_low,wing_front,wing_back_up,wing_back_low,...
                wing_side_up,wing_side_low,wing_side_front};

        end
        
        function PSI=shapeWing(self,XI)
            % height normalize to 1
            %
            PSI=1-abs(1-2*XI);
        end
    
    end

    % output function
    methods
        function [X_total,Y_total,Z_total,XI_total,PSI_total]=calSurfaceMatrix(self,varargin)
            % generate waverider wing wgs data
            % xi, x=X(xi)
            % psi, y=Y(xi)
            % zeta, z=Z(xi,psi)
            %
            % notice colume vector in matrix is wgs format
            %
            X_total=struct();Y_total=struct();Z_total=struct();XI_total=struct();PSI_total=struct();

            if nargin <= 2
                if nargin == 1
                    value_torl=1e-3;
                else
                    value_torl=varargin{1};
                end

                % calculate all surface matrix
                for surf_idx=1:length(self.surface_list)
                    surf=self.surface_list{surf_idx};
                    [X_surf,Y_surf,Z_surf,XI_surf,PSI_surf]=surf.calSurface(value_torl);
                    X_total.(surf.name)=X_surf;Y_total.(surf.name)=Y_surf;Z_total.(surf.name)=Z_surf;XI_total.(surf.name)=XI_surf;PSI_total.(surf.name)=PSI_surf;
                end

            else
                xi_grid_num_head=varargin{1};
                psi_grid_num_head=varargin{2};
                xi_grid_num_body=varargin{3};
                psi_grid_num_wing=varargin{4};
                edge_gird_num=varargin{5};

                % calculate head
                [X_total.('head_up'),Y_total.('head_up'),Z_total.('head_up'),XI_total.('head_up'),PSI_total.('head_up')]=self.getSurface('head_up').calSurface(xi_grid_num_head,psi_grid_num_head);
                [X_total.('head_low'),Y_total.('head_low'),Z_total.('head_low'),XI_total.('head_low'),PSI_total.('head_low')]=self.getSurface('head_low').calSurface(xi_grid_num_head,psi_grid_num_head);

                % calculate blunt head side
                [X_total.('head_side'),Y_total.('head_side'),Z_total.('head_side'),XI_total.('head_side'),PSI_total.('head_side')]=self.getSurface('head_side').calSurface(xi_grid_num_head,edge_gird_num*2);

                % calculate body
                [X_total.('body_up'),Y_total.('body_up'),Z_total.('body_up'),XI_total.('body_up'),PSI_total.('body_up')]=self.getSurface('body_up').calSurface(xi_grid_num_body,psi_grid_num_head);
                [X_total.('body_low'),Y_total.('body_low'),Z_total.('body_low'),XI_total.('body_low'),PSI_total.('body_low')]=self.getSurface('body_low').calSurface(xi_grid_num_body,psi_grid_num_head);

                % calculate body back
                [X_total.('body_back_up'),Y_total.('body_back_up'),Z_total.('body_back_up'),XI_total.('body_back_up'),PSI_total.('body_back_up')]=self.getSurface('body_back_up').calSurface(psi_grid_num_head,edge_gird_num);
                [X_total.('body_back_low'),Y_total.('body_back_low'),Z_total.('body_back_low'),XI_total.('body_back_low'),PSI_total.('body_back_low')]=self.getSurface('body_back_low').calSurface(psi_grid_num_head,edge_gird_num);

                % calculate transition wing
                [X_total.('tri_wing_up'),Y_total.('tri_wing_up'),Z_total.('tri_wing_up'),XI_total.('tri_wing_up'),PSI_total.('tri_wing_up')]=self.getSurface('tri_wing_up').calSurface(psi_grid_num_wing,xi_grid_num_body);
                [X_total.('tri_wing_low'),Y_total.('tri_wing_low'),Z_total.('tri_wing_low'),XI_total.('tri_wing_low'),PSI_total.('tri_wing_low')]=self.getSurface('tri_wing_low').calSurface(psi_grid_num_wing,xi_grid_num_body);

                % calculate blunt tir wing front
                [X_total.('tri_wing_front'),Y_total.('tri_wing_front'),Z_total.('tri_wing_front'),XI_total.('tri_wing_front'),PSI_total.('tri_wing_front')]=self.getSurface('tri_wing_front').calSurface(psi_grid_num_wing,2*edge_gird_num);

                % calculate tri wing back
                [X_total.('tri_wing_back_up'),Y_total.('tri_wing_back_up'),Z_total.('tri_wing_back_up'),XI_total.('tri_wing_back_up'),PSI_total.('tri_wing_back_up')]=self.getSurface('tri_wing_back_up').calSurface(psi_grid_num_wing,edge_gird_num);
                [X_total.('tri_wing_back_low'),Y_total.('tri_wing_back_low'),Z_total.('tri_wing_back_low'),XI_total.('tri_wing_back_low'),PSI_total.('tri_wing_back_low')]=self.getSurface('tri_wing_back_low').calSurface(psi_grid_num_wing,edge_gird_num);

                % calculate wing
                [X_total.('wing_up'),Y_total.('wing_up'),Z_total.('wing_up'),XI_total.('wing_up'),PSI_total.('wing_up')]=self.getSurface('wing_up').calSurface(psi_grid_num_wing,xi_grid_num_body);
                [X_total.('wing_low'),Y_total.('wing_low'),Z_total.('wing_low'),XI_total.('wing_low'),PSI_total.('wing_low')]=self.getSurface('wing_low').calSurface(psi_grid_num_wing,xi_grid_num_body);

                % calculate blunt wing front
                [X_total.('wing_front'),Y_total.('wing_front'),Z_total.('wing_front'),XI_total.('wing_front'),PSI_total.('wing_front')]=self.getSurface('wing_front').calSurface(psi_grid_num_wing,2*edge_gird_num);

                % calculate wing back
                [X_total.('wing_back_up'),Y_total.('wing_back_up'),Z_total.('wing_back_up'),XI_total.('wing_back_up'),PSI_total.('wing_back_up')]=self.getSurface('wing_back_up').calSurface(psi_grid_num_wing,edge_gird_num);
                [X_total.('wing_back_low'),Y_total.('wing_back_low'),Z_total.('wing_back_low'),XI_total.('wing_back_low'),PSI_total.('wing_back_low')]=self.getSurface('wing_back_low').calSurface(psi_grid_num_wing,edge_gird_num);

                % calculate wing side
                [X_total.('wing_side_up'),Y_total.('wing_side_up'),Z_total.('wing_side_up'),XI_total.('wing_side_up'),PSI_total.('wing_side_up')]=self.getSurface('wing_side_up').calSurface(xi_grid_num_body,edge_gird_num);
                [X_total.('wing_side_low'),Y_total.('wing_side_low'),Z_total.('wing_side_low'),XI_total.('wing_side_low'),PSI_total.('wing_side_low')]=self.getSurface('wing_side_low').calSurface(xi_grid_num_body,edge_gird_num);
                [X_total.('wing_side_front'),Y_total.('wing_side_front'),Z_total.('wing_side_front'),XI_total.('wing_side_front'),PSI_total.('wing_side_front')]=self.getSurface('wing_side_front').calSurface(edge_gird_num,2*edge_gird_num);
            end

            % fix normal vector
            for surf_idx=1:length(self.surface_list)
                surf_name=self.surface_list{surf_idx}.name;
                if contains(surf_name,'low') || contains(surf_name,'side')
                    X_total.(surf_name)=fliplr(X_total.(surf_name));
                    Y_total.(surf_name)=fliplr(Y_total.(surf_name));
                    Z_total.(surf_name)=fliplr(Z_total.(surf_name));
                    XI_total.(surf_name)=fliplr(XI_total.(surf_name));
                    PSI_total.(surf_name)=fliplr(PSI_total.(surf_name));
                end
            end

        end

        function part=getWGSMesh(self,part_name,varargin)
            % get mesh data in LaWGS format
            %

            if length(varargin) > 1
                [X_total,Y_total,Z_total]=self.calSurfaceMatrix...
                    (varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
            elseif length(varargin) == 1
                [X_total,Y_total,Z_total]=self.calSurfaceMatrix(varargin{1});
            else
                [X_total,Y_total,Z_total]=self.calSurfaceMatrix();
            end

            mesh_list=cell(surface_number,1);
            for surf_index=1:length(self.surface_list)
                surf_name=self.surface_list{surf_index}.name;
                mesh.X=X_total.(surf_name);
                mesh.Y=Y_total.(surf_name);
                mesh.Z=Z_total.(surf_name);
                mesh.element_type='wgs';
                mesh_list{surf_index,1}=mesh;
            end

            part.name=part_name;
            part.mesh_list=mesh_list;
        end

        function writeStepOpenShell(self,step_filestr,varargin)
            % write surface into step file
            %

            % check file name
            if length(step_filestr) > 4
                if ~strcmp(step_filestr((end-4):end),'.step')
                    step_filestr=[step_filestr,'.step'];
                end
            else
                step_filestr=[step_filestr,'.step'];
            end
            [~,step_filename,~]=fileparts(step_filestr);

            % write head
            step_file=fopen(step_filestr,'w');
            object_index=1;
            object_index=writeStepHead(self,step_file,object_index,step_filename);

            % write surface
            ADVANCED_FACE_index_list=zeros(1,length(self.surface_list));
            if length(varargin) > 1
                [X_total,Y_total,Z_total,XI_total,PSI_total]=self.calSurfaceMatrix...
                    (varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
            elseif length(varargin) == 1
                [X_total,Y_total,Z_total,XI_total,PSI_total]=self.calSurfaceMatrix(varargin{1});
            else
                [X_total,Y_total,Z_total,XI_total,PSI_total]=self.calSurfaceMatrix();
            end

            for surf_idx=1:length(self.surface_list)
                surf_name=self.surface_list{surf_idx}.name;
                if contains(surf_name,'wing') && ~contains(surf_name,'front')
                    [X_total.(surf_name),Y_total.(surf_name),Z_total.(surf_name),XI_total.(surf_name),PSI_total.(surf_name)]=self.surface_list{surf_idx}.calSurface(21,21);
                end
            end

            for surf_idx=1:length(self.surface_list)
                surf_name=self.surface_list{surf_idx}.name;
                if size(X_total.(surf_name),2) == 2
                    u_list=[0,0,1,1];
                else
                    u_list=[0,0,0,XI_total.(surf_name)(1,:),1,1,1];
                end
                if size(X_total.(surf_name),1) == 2
                    v_list=[0,0,1,1];
                else
                    v_list=[0,0,0,PSI_total.(surf_name)(:,1)',1,1,1];
                end
                surf=SurfaceBSpline(surf_name,[],[],[],X_total.(surf_name),Y_total.(surf_name),Z_total.(surf_name),[],[],u_list,v_list);
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
                [X_total,Y_total,Z_total]=self.calSurfaceMatrix...
                    (varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
            elseif length(varargin) == 1
                [X_total,Y_total,Z_total]=self.calSurfaceMatrix(varargin{1});
            else
                [X_total,Y_total,Z_total]=self.calSurfaceMatrix();
            end

            for surf_index=1:length(self.surface_list)
                surf_name=self.surface_list{surf_index}.name;
                surface(X_total.(surf_name),Y_total.(surf_name),Z_total.(surf_name))
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
                [U,V]=surf.calCoordinate(point(:,1),point(:,2),point(:,3));
                WWD_coord.(surf_name).index=point_idx;
                WWD_coord.(surf_name).U=U;
                WWD_coord.(surf_name).V=V;
            end

        end

        function mesh_point=calMeshPoint(self,WWD_coord)
            % calculate all mesh point by input WWD coord
            %
            mesh_point=struct();

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
                mesh_point.(surf_name).index=point_idx;
                mesh_point.(surf_name).X=X;
                mesh_point.(surf_name).Y=Y;
                mesh_point.(surf_name).Z=Z;
            end

        end
      end

    % decode x function
    methods(Static)
        function [par_width,par_hight_up,par_hight_low,...
                par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
                par_rho1,par_rho12,par_rho23,par_WS1,par_WS2,par_WT_up,par_WT_low]=decode(par_input)
            % decode x into parameter
            %

            if isnumeric(par_input)
                par_input=num2cell(par_input,[1,16]);
            end

            % length parameter
            par_M_up=par_input{1};
            par_M_low=par_input{2};
            % width parameter
            par_width=par_input{3};
            par_T=par_input{4};
            par_N_up=par_input{5};
            par_N_low=par_input{6};
            % height parameter
            par_hight_up=par_input{7};
            par_hight_low=par_input{8};
            par_R=par_input{9};
            % wing parameter
            par_rho1=par_input{10};
            par_rho12=par_input{11};
            par_rho23=par_input{12};
            par_WS1=par_input{13};
            par_WS2=par_input{14};
            par_WT_up=par_input{15};
            par_WT_low=par_input{16};
        end
    end
end
