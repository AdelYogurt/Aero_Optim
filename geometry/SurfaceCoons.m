classdef SurfaceCoons < SurfaceBSpline
    % Coons surface
    %
    properties
        point_00;
        point_10;
        point_01;
        point_11;

        curve_u0;
        curve_u1;
        curve_0v;
        curve_1v;
    end
    methods
        function self=SurfaceCoons(name,curve_u0,curve_u1,curve_0v,curve_1v)
            % generate Coons surface
            %
            self=self@SurfaceBSpline(name);

            geom_torl=100*eps;
            % get vertex
            [x0,y0,z0]=curve_u0(0);
            [x1,y1,z1]=curve_0v(0);
            if norm([x0,y0,z0]-[x1,y1,z1]) > geom_torl
                error('SurfaceCoons: curve do not connect');
            else
                self.point_00=([x0,y0,z0]+[x1,y1,z1])/2;
            end

            [x0,y0,z0]=curve_u0(1);
            [x1,y1,z1]=curve_1v(0);
            if norm([x0,y0,z0]-[x1,y1,z1]) > geom_torl
                error('SurfaceCoons: curve do not connect');
            else
                self.point_10=([x0,y0,z0]+[x1,y1,z1])/2;
            end

            [x0,y0,z0]=curve_u1(0);
            [x1,y1,z1]=curve_0v(1);
            if norm([x0,y0,z0]-[x1,y1,z1]) > geom_torl
                error('SurfaceCoons: curve do not connect');
            else
                self.point_01=([x0,y0,z0]+[x1,y1,z1])/2;
            end

            [x0,y0,z0]=curve_u1(1);
            [x1,y1,z1]=curve_1v(1);
            if norm([x0,y0,z0]-[x1,y1,z1]) > geom_torl
                error('SurfaceCoons: curve do not connect');
            else
                self.point_11=([x0,y0,z0]+[x1,y1,z1])/2;
            end

            self.curve_u0=curve_u0;
            self.curve_u1=curve_u1;
            self.curve_0v=curve_0v;
            self.curve_1v=curve_1v;
        end

        function [X,Y,Z]=calPoint(self,U,V)
            [rank_num,colume_num]=size(U);
            X=zeros(rank_num,colume_num);
            Y=zeros(rank_num,colume_num);
            Z=zeros(rank_num,colume_num);

            x00=self.point_00(1);
            y00=self.point_00(2);
            z00=self.point_00(3);

            x10=self.point_10(1);
            y10=self.point_10(2);
            z10=self.point_10(3);

            x01=self.point_01(1);
            y01=self.point_01(2);
            z01=self.point_01(3);

            x11=self.point_11(1);
            y11=self.point_11(2);
            z11=self.point_11(3);

            for rank_idx=1:rank_num
                for colume_idx=1:colume_num
                    u_x=U(rank_idx,colume_idx);
                    v_x=V(rank_idx,colume_idx);

                    [xu0,yu0,zu0]=self.curve_u0(u_x);
                    [xu1,yu1,zu1]=self.curve_u1(u_x);
                    [x0v,y0v,z0v]=self.curve_0v(v_x);
                    [x1v,y1v,z1v]=self.curve_1v(v_x);

                    blend_u=[-1,1-u_x,u_x];
                    blend_v=[-1;1-v_x;v_x];

                    X(rank_idx,colume_idx)=-blend_u*[
                        0.0,xu0,xu1;
                        x0v,x00,x01;
                        x1v,x10,x11]*blend_v;
                    Y(rank_idx,colume_idx)=-blend_u*[
                        0.0,yu0,yu1;
                        y0v,y00,y01;
                        y1v,y10,y11]*blend_v;
                    Z(rank_idx,colume_idx)=-blend_u*[
                        0.0,zu0,zu1;
                        z0v,z00,z01;
                        z1v,z10,z11]*blend_v;
                end
            end
        end
    end
end