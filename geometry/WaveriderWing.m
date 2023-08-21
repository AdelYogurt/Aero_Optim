classdef WaveriderWing < Body
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
        body_v_max_fun;

        % local coordinate parameter
        blunt_u=0.002; % is local length parameter
        blunt_eta=0.005; % is local length parameter

        % function
        shape_head_edge_x;
        shape_tri_wing_edge_x;
        shape_wing_edge_x;
        shape_wing_fcn;
    end

    % define function
    methods
        function self=WaveriderWing(total_length,varargin)
            % generate parameterized waverider with wing
            % u, x=X(u)
            % v, y=Y(u)
            % w, z=Z(u,v)
            %
            % notice:
            % par_rho1 should between 0 and 1
            %
            if nargin == 3
                par_input=varargin{1};
                waverider_type=varargin{2};
            else
                par_input=varargin(1:end-1);
                waverider_type=varargin{end};
            end

            [par_width,par_hight_up,par_hight_low,...
                par_T,par_M_up,par_M_low,par_N_up,par_N_low,par_R,...
                par_rho1,par_rho12,par_rho23,par_WS1,par_WS2,par_WT_up,par_WT_low]=WaveriderWing.decode(par_input);

            if par_rho1 >= 1 || par_rho1 <= 0
                error('WaveriderWing: do not support pure waverider');
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
            body_v_max_fun=@(u) y_cut./(u.^par_T*par_width)-eps;
            self.body_v_max_fun=body_v_max_fun;

            % wing base parameter
            v_end_cut=self.body_v_max_fun(1);

            tri_wing_length=par_rho1*total_length;
            tri_wing_width=par_WS1;
            tri_wing_height_up=(v_end_cut+0.5)^par_N_up*(0.5-v_end_cut)^par_N_up/(0.5)^(2*par_N_up)*par_hight_up;
            tri_wing_height_low=(v_end_cut+0.5)^par_N_low*(0.5-v_end_cut)^par_N_low/(0.5)^(2*par_N_low)*par_hight_low;

            wing_length=par_rho1*par_rho12*total_length;
            wing_width=par_WS2;
            wing_height_up=par_WT_up;
            wing_height_low=par_WT_low;

            switch waverider_type
                case 'Fei'
                    shape_wing_fcn=@self.shapeWingFei;
                case 'Dia'
                    shape_wing_fcn=@self.shapeWingDia;
                case 'Tri'
                    shape_wing_fcn=@self.shapeWingTri;
                otherwise
                    error('WaveriderWing: unknown waverider type');
            end
            self.shape_wing_fcn=shape_wing_fcn;

            if ~strcmp(waverider_type,'Fei')
                % calculate blunt radius parameter
                % local coordinate, local x is global x, local y is global z
                % head radius
                % calculate gradient of up and low discrete surface to local radius center
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
                shape_wing_edge_x=@(Y) (sqrt(radius_wing_sq-Y.^2))-radius_center_wing;
                self.shape_wing_edge_x=shape_wing_edge_x;

                % calculate head side blunt radius of yoz plane
                dz_dx=(self.blunt_eta)^par_N_up*(1-self.blunt_eta)^par_N_up/(0.5)^(2*par_N_up)*par_hight_up*(1-par_rho1)^par_M_up/(2*self.blunt_eta*y_cut);
                radius_center_up=par_R*dz_dx;
                dz_dx=(self.blunt_eta)^par_N_low*(1-self.blunt_eta)^par_N_low/(0.5)^(2*par_N_low)*par_hight_low*(1-par_rho1)^par_M_low/(2*self.blunt_eta*y_cut);
                radius_center_low=par_R*dz_dx;
                radius_center_head_side=max(radius_center_low,radius_center_up);
                radius_head_side_sq=(radius_center_head_side*radius_center_head_side+par_R*par_R);
            end

            %% define head
            shape_fcn_Y=@(U) (U*(1-par_rho1)).^par_T;
            shape_fcn_Z_up=@(U,V) (U*(1-par_rho1)).^par_M_up;
            shape_fcn_Z_low=@(U,V) (U*(1-par_rho1)).^par_M_low;
            head_up=SurfaceCST3D('head_up',head_length,par_width,par_hight_up,[],shape_fcn_Y,shape_fcn_Z_up,[par_N_up,par_N_up],global_symmetry_y);
            head_low=SurfaceCST3D('head_low',head_length,par_width,-par_hight_low,[],shape_fcn_Y,shape_fcn_Z_low,[par_N_low,par_N_low],global_symmetry_y);
            
            if strcmp(waverider_type,'Fei')
                tran_fun_Z_up=@(U,V) -par_R*(1-U).^100;
                tran_fun_Z_low=@(U,V) par_R*(1-U).^100;
                % blunt deform
                head_up.addDeform([],[],tran_fun_Z_up)
                head_low.addDeform([],[],tran_fun_Z_low)
            end

            % blunt translation
            head_up.addTranslation(0,0,par_R);
            head_low.addTranslation(0,0,-par_R);

            %% define blunt head side
            % local coordinate, local x is global -x, local y is global z, local z is global y
            if strcmp(waverider_type,'Fei')
                shape_fcn_Y=@(U) 1-U.^100;
                head_side=SurfaceCST3D('head_side',head_length,2*par_R,0,[],shape_fcn_Y,[0,0],[0,0]);
            else
                shape_fcn_X=@(V) 1-shape_tri_wing_edge_x((V-0.5)*par_R*2)/head_length+shape_head_edge_x((V-0.5)*par_R*2)/head_length;
                shape_fcn_Z=@(U,V) (sqrt(radius_head_side_sq-(((V-0.5)*2)*par_R).^2)-radius_center_head_side).*(1-U);
                % shape_fcn_Z=@(U,V) (sqrt(radius_head_side_sq-(((V-0.5)*2)*par_R).^2)-radius_center_head_side).*(1-U).^0.3;
                head_side=SurfaceCST3D('head_side',head_length,2*par_R,1,shape_fcn_X,[],shape_fcn_Z,{[0.0,0.0],[0.0,0.0]});
            end

            if strcmp(waverider_type,'Fei')
                tran_fun_X=[];
                tran_fun_Y=@(U) -par_R*(1-U.^100);
            else
                tran_fun_X=@(V) shape_tri_wing_edge_x((V-0.5)*par_R*2);
                tran_fun_Y=[];
            end
            Y_edge=@(U) U.^par_T*par_width/2;
            tran_fun_Z=@(U,V) Y_edge((1-U)*(1-par_rho1));
            % deform surface
            head_side.addDeform(tran_fun_X,tran_fun_Y,tran_fun_Z);

            % rotation surface
            head_side.addRotation(90,180,0);

            % translation surface
            if strcmp(waverider_type,'Fei')
                head_side.addTranslation(head_length,0,0);
            else
                head_side.addTranslation(head_length,0,-par_R);
            end

            %% define body
            shape_fcn_Z_up=@(U,V) (U*par_rho1+(1-par_rho1)).^par_M_up;
            class_fcn_Z_up=@(U,V) ((V-0.5)*2.*body_v_max_fun(U*par_rho1+(1-par_rho1))+0.5).^par_N_up...
                .*(1-((V-0.5)*2.*body_v_max_fun(U*par_rho1+(1-par_rho1))+0.5)).^par_N_up/(0.5)^(2*par_N_up);
            shape_fcn_Z_low=@(U,V) (U*par_rho1+(1-par_rho1)).^par_M_low;
            class_fcn_Z_low=@(U,V) ((V-0.5)*2.*body_v_max_fun(U*par_rho1+(1-par_rho1))+0.5).^par_N_low...
                .*(1-((V-0.5)*2.*body_v_max_fun(U*par_rho1+(1-par_rho1))+0.5)).^par_N_low/(0.5)^(2*par_N_low);
            body_up=SurfaceCST3D('body_up',body_length,y_cut*2,par_hight_up,[],[],shape_fcn_Z_up,class_fcn_Z_up,global_symmetry_y);
            body_low=SurfaceCST3D('body_low',body_length,y_cut*2,-par_hight_low,[],[],shape_fcn_Z_low,class_fcn_Z_low,global_symmetry_y);

            % blunt translation
            body_up.addTranslation(head_length,0,par_R);
            body_low.addTranslation(head_length,0,-par_R);

            %% define body back
            shape_fcn_Y_up=@(U) par_hight_up*(0.5+U*(v_end_cut)).^par_N_up.*(0.5-U*(v_end_cut)).^par_N_up/(0.5)^(2*par_N_up)+par_R;
            shape_fcn_Y_low=@(U) par_hight_low*(0.5+U*(v_end_cut)).^par_N_low.*(0.5-U*(v_end_cut)).^par_N_low/(0.5)^(2*par_N_low)+par_R;

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
            shape_fcn_Y=@(U) (1-U*(1-par_rho12));
            shape_body_size_up=@(V) body_up.LZ.*...
                body_up.shape_fcn_Z(1-V,1).*...
                body_up.class_fcn_Z(1-V,1);
            shape_body_size_low=@(V) -body_low.LZ.*...
                body_low.shape_fcn_Z(1-V,1).*...
                body_low.class_fcn_Z(1-V,1);

            if strcmp(waverider_type,'Fei')
                shape_fcn_Z_up=@(U,V) (1-U).*shape_body_size_up(V)+U.*(self.shape_wing_fcn(1-V)*(wing_height_up)-par_R);
                shape_fcn_Z_low=@(U,V) (1-U).*shape_body_size_low(V)+U.*(wing_height_low-par_R);
            else
                shape_fcn_Z_up=@(U,V) (1-U).*shape_body_size_up(V)+U.*self.shape_wing_fcn(1-V)*(wing_height_up);
                shape_fcn_Z_low=@(U,V) (1-U).*shape_body_size_low(V)+U.*self.shape_wing_fcn(1-V)*(wing_height_low);
            end

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
            if strcmp(waverider_type,'Fei')
                shape_fcn_Y=@(U) (1-U)*wing_height_low*2+U*par_R*2;
                tri_wing_front=SurfaceCST3D('tri_wing_front',tri_wing_width,1,0,[],shape_fcn_Y,[],{[0,0]});
            else
                shape_fcn_X=@(V) 1-(sqrt(radius_head_side_sq-(((V-0.5)*2)*par_R).^2)-radius_center_head_side)/tri_wing_width; % if do not want blunt intersection, remove this and change blunt head side class N1
                shape_fcn_Z=@(U,V) (1-U).*shape_wing_edge_x((V-0.5)*2*par_R)+U.*shape_tri_wing_edge_x((V-0.5)*2*par_R);
                tri_wing_front=SurfaceCST3D('tri_wing_front',tri_wing_width,2*par_R,1,shape_fcn_X,[],[],shape_fcn_Z);
            end

            % deform surface
            if strcmp(waverider_type,'Fei')
                tran_fun_Y=@(U) U*(wing_height_low-par_R);
            else
                tran_fun_Y=[];
            end
            tran_fun_Z=@(U,V) U*(par_rho1*(1-par_rho12)*total_length);
            tri_wing_front.addDeform([],tran_fun_Y,tran_fun_Z);

            % rotation surface
            tri_wing_front.addRotation(90,-90,0);

            % translation surface
            if strcmp(waverider_type,'Fei')
                tri_wing_front.addTranslation(total_length-wing_length,y_cut+par_WS1,-wing_height_low);
            else
                tri_wing_front.addTranslation(total_length-wing_length,y_cut+par_WS1,-par_R);
            end

            %% define tri wing back
            % local coordination, local x is global y, local y is global z
            if strcmp(waverider_type,'Fei')
                shape_fcn_Y_up=@(U) (1-U)*(tri_wing_height_up+par_R)+U*wing_height_low;
                shape_fcn_Y_low=@(U) (1-U)*(tri_wing_height_low+par_R)+U*wing_height_low;
            elseif strcmp(waverider_type,'Dia')
                shape_fcn_Y_up=@(U) (1-U)*(tri_wing_height_up+par_R)+U*(par_R);
                shape_fcn_Y_low=@(U) (1-U)*(tri_wing_height_low+par_R)+U*(par_R);
            elseif strcmp(waverider_type,'Tri')
                shape_fcn_Y_up=@(U) (1-U)*(tri_wing_height_up+par_R)+U*(wing_height_up+par_R);
                shape_fcn_Y_low=@(U) (1-U)*(tri_wing_height_low+par_R)+U*(wing_height_low+par_R);
            end

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
            shape_fcn_Y=@(U) (1-U*(1-par_rho23));
            if strcmp(waverider_type,'Fei')
                shape_fcn_Z_up=@(U,V) self.shape_wing_fcn(1-V).*(1-U*(1-par_rho23));
                shape_fcn_Z_low=@(U,V) (1-U*(1-par_rho23));
            else
                shape_fcn_Z_up=@(U,V) self.shape_wing_fcn(1-V).*(1-U*(1-par_rho23));
                shape_fcn_Z_low=@(U,V) self.shape_wing_fcn(1-V).*(1-U*(1-par_rho23));
            end

            wing_up=SurfaceCST3D('wing_up',wing_width,wing_length,wing_height_up,[],shape_fcn_Y,shape_fcn_Z_up,{[0,0],[0,0]});
            wing_low=SurfaceCST3D('wing_low',wing_width,wing_length,-wing_height_low,[],shape_fcn_Y,shape_fcn_Z_low,{[0,0],[0,0]});

            % add rotation
            wing_up.addRotation(0,0,90);
            wing_low.addRotation(0,0,90);

            % add translation
            if strcmp(waverider_type,'Fei')
                wing_up.addTranslation(total_length,y_cut+par_WS1,0);
                wing_low.addTranslation(total_length,y_cut+par_WS1,0);
            else
                wing_up.addTranslation(total_length,y_cut+par_WS1,par_R);
                wing_low.addTranslation(total_length,y_cut+par_WS1,-par_R);
            end

            %% define blunt wing front
            % local coordinate, local x is global -y, local y is global z, local z is global -x
            if strcmp(waverider_type,'Fei')
                wing_front=SurfaceCST3D('wing_front',wing_width,2*wing_height_low,0,[],[],[],{[0,0]});
            else
                shape_fcn_Z=@(U,V) shape_wing_edge_x((V-0.5)*2*par_R);
                wing_front=SurfaceCST3D('wing_front',wing_width,2*par_R,1,[],[],shape_fcn_Z,{[0,0]});

                % deform surface
                tran_fun_Z=@(U,V) U*(par_rho1*par_rho12*(1-par_rho23)*total_length);
                wing_front.addDeform([],[],tran_fun_Z);
            end

            % rotation surface
            wing_front.addRotation(90,-90,0);

            % translation surface
            if strcmp(waverider_type,'Fei')
                wing_front.addTranslation(total_length-wing_length*par_rho23,y_cut+par_WS1+par_WS2,-wing_height_low);
            else
                wing_front.addTranslation(total_length-wing_length*par_rho23,y_cut+par_WS1+par_WS2,-par_R);
            end

            %% define wing back
            % local coordination, local x is global y, local y is global z
            if strcmp(waverider_type,'Fei')
                wing_back_up=SurfaceCST3D('wing_back_up',wing_width,wing_height_low,0);
                wing_back_low=SurfaceCST3D('wing_back_up',wing_width,-wing_height_low,0);
            elseif strcmp(waverider_type,'Dia')
                wing_back_up=SurfaceCST3D('wing_back_up',wing_width,par_R,0);
                wing_back_low=SurfaceCST3D('wing_back_low',wing_width,-par_R,0);
            elseif strcmp(waverider_type,'Tri')
                wing_back_up=SurfaceCST3D('wing_back_up',wing_width,(wing_height_up+par_R),0);
                wing_back_low=SurfaceCST3D('wing_back_low',wing_width,-(wing_height_low+par_R),0);
            end

            % rotation to global coordinate
            wing_back_up.addRotation(90,90,0);
            wing_back_low.addRotation(90,90,0);

            % translation to global coordination
            wing_back_up.addTranslation(total_length,y_cut+par_WS1,0);
            wing_back_low.addTranslation(total_length,y_cut+par_WS1,0);

            %% define wing side
            % local coordinate, local x is global -x, local y is global z
            if strcmp(waverider_type,'Fei')
                shape_fcn_Y_up=@(V) self.shape_wing_fcn(1-V)*(wing_height_up);
                shape_fcn_Y_low=@(V) wing_height_low;
            else
                shape_fcn_Y_up=@(V) self.shape_wing_fcn(1-V)*par_WT_up*par_rho23+par_R;
                shape_fcn_Y_low=@(V) self.shape_wing_fcn(1-V)*par_WT_low*par_rho23+par_R;
            end

            wing_side_up=SurfaceCST3D('wing_side_up',wing_length*par_rho23,1,0,[],shape_fcn_Y_up);
            wing_side_low=SurfaceCST3D('wing_side_low',wing_length*par_rho23,-1,0,[],shape_fcn_Y_low);

            % rotation to global coordinate
            wing_side_up.addRotation(90,180,0);
            wing_side_low.addRotation(90,180,0);

            % translation to global coordinate
            wing_side_up.addTranslation(total_length,y_cut+par_WS1+par_WS2,0);
            wing_side_low.addTranslation(total_length,y_cut+par_WS1+par_WS2,0);
            
            %% sort data
            if strcmp(waverider_type,'Fei')
                self.surface_list={head_up,head_low,head_side,...
                    body_up,body_low,body_back_up,body_back_low,...
                    tri_wing_up,tri_wing_low,tri_wing_front,tri_wing_back_up,tri_wing_back_low,...
                    wing_up,wing_low,wing_front,wing_back_up,wing_back_low,...
                    wing_side_up,wing_side_low};
            else
                % wing side blunt front
                shape_fcn_X=@(V) shape_wing_edge_x((V-0.5)*2*par_R);
                wing_side_front=SurfaceCST3D('wing_side_front',1,2*par_R,0,shape_fcn_X);

                % rotation to global coordinate
                wing_side_front.addRotation(90,180,0);

                % translation to global coordinate
                wing_side_front.addTranslation(total_length-wing_length*par_rho23,y_cut+par_WS1+par_WS2,-par_R);

                self.surface_list={head_up,head_low,head_side,...
                body_up,body_low,body_back_up,body_back_low,...
                tri_wing_up,tri_wing_low,tri_wing_front,tri_wing_back_up,tri_wing_back_low,...
                wing_up,wing_low,wing_front,wing_back_up,wing_back_low,...
                wing_side_up,wing_side_low,wing_side_front};
            end

        end
        
        function V=shapeWingFei(~,U)
            % height normalize to 1
            %
            V=0.25-0.75*(abs(1-2*U)-1);
        end

        function V=shapeWingDia(~,U)
            % height normalize to 1
            %
            V=1-abs(1-2*U);
        end

        function V=shapeWingTri(~,U)
            % height normalize to 1
            %
            V=U;
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

            surf_total=cell(length(self.surface_list),1);
            if nargin <= 2
                if nargin == 1
                    value_torl=1e-3;
                else
                    value_torl=varargin{2};
                end

                % calculate all surface matrix
                for surf_idx=1:length(self.surface_list)
                    surf.name=self.surface_list{surf_idx}.name;
                    [surf.X,surf.Y,surf.Z,surf.U,surf.V]=self.surface_list{surf_idx}.calSurface(value_torl);
                    surf.element_type='wgs';
                    surf_total{surf_idx}=surf;
                end
            else
                u_grid_num_head=varargin{1};
                v_grid_num_head=varargin{2};
                u_grid_num_body=varargin{3};
                v_grid_num_wing=varargin{4};
                edge_gird_num=varargin{5};

                num_list=[
                    u_grid_num_head,v_grid_num_head;
                    u_grid_num_head,v_grid_num_head;
                    u_grid_num_head,edge_gird_num*2;
                    u_grid_num_body,v_grid_num_head;
                    u_grid_num_body,v_grid_num_head;
                    v_grid_num_head,edge_gird_num;
                    v_grid_num_head,edge_gird_num;
                    v_grid_num_wing,u_grid_num_body;
                    v_grid_num_wing,u_grid_num_body;
                    v_grid_num_wing,2*edge_gird_num;
                    v_grid_num_wing,edge_gird_num;
                    v_grid_num_wing,edge_gird_num;
                    v_grid_num_wing,u_grid_num_body;
                    v_grid_num_wing,u_grid_num_body;
                    v_grid_num_wing,2*edge_gird_num;
                    v_grid_num_wing,edge_gird_num;
                    v_grid_num_wing,edge_gird_num;
                    u_grid_num_body,edge_gird_num;
                    u_grid_num_body,edge_gird_num;
                    edge_gird_num,2*edge_gird_num];
                name_list={'head_up';'head_low';'head_side';
                    'body_up';'body_low';'body_back_up';'body_back_low';
                    'tri_wing_up';'tri_wing_low';'tri_wing_front';
                    'tri_wing_back_up';'tri_wing_back_low';
                    'wing_up';'wing_low';'wing_front';
                    'wing_back_up';'wing_back_low';
                    'wing_side_up';'wing_side_low';'wing_side_front'};

                for surf_idx=1:length(name_list)
                    surf_CST=self.getSurface(name_list{surf_idx});
                    surf.name=surf_CST.name;
                    [surf.X,surf.Y,surf.Z,surf.U,surf.V]=surf_CST.calSurface(num_list(surf_idx,1),num_list(surf_idx,2));
                    surf.element_type='wgs';
                    surf_total{surf_idx}=surf;
                end
            end

            % fix normal vector
            for surf_idx=1:length(surf_total)
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
                surf_total=self.calSurfaceMatrix...
                    (varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
            elseif length(varargin) == 1
                surf_total=self.calSurfaceMatrix(varargin{1});
            else
                surf_total=self.calSurfaceMatrix();
            end

            for surf_idx=1:length(surf_total)
                surf_name=surf_total{surf_idx}.name;
                if contains(surf_name,'wing') && ~contains(surf_name,'front')
                    [surf_total{surf_idx}.X,surf_total{surf_idx}.Y,surf_total{surf_idx}.Z,...
                        surf_total{surf_idx}.U,surf_total{surf_idx}.V]=self.surface_list{surf_idx}.calSurface(21,21);
                end
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
                par_rho1,par_rho12,par_rho23,par_WS1,par_WS2,par_WT_up,par_WT_low]=decode(param)
            % decode x into parameter
            %

            if isnumeric(param)
                param=num2cell(param,[1,17]);
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
            par_rho1=param{10};
            par_rho12=param{11};
            par_rho23=param{12};
            par_WS1=param{13};
            par_WS2=param{14};
            par_WT_up=param{15};
            par_WT_low=param{16};
        end
    end
end
