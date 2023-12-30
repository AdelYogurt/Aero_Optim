classdef FaceCoons < FaceNURBS
    % Coons surface
    %
    properties
        vertex_00;
        vertex_10;
        vertex_01;
        vertex_11;

        edge_u0;
        edge_u1;
        edge_0v;
        edge_1v;
    end

    methods
        function self=FaceCoons(name,edge_u0,edge_u1,edge_0v,edge_1v)
            % generate Coons surface
            %
            self=self@FaceNURBS(name);

            geom_torl=100*eps;

            % get vertex
            vtx_0=edge_u0(0);
            vtx_1=edge_0v(0);
            if norm(vtx_0-vtx_1) > geom_torl
                error('FaceCoons: edge do not connect');
            else
                self.vertex_00=(vtx_0+vtx_1)/2;
            end

            vtx_0=edge_u0(1);
            vtx_1=edge_1v(0);
            if norm(vtx_0-vtx_1) > geom_torl
                error('FaceCoons: edge do not connect');
            else
                self.vertex_10=(vtx_0+vtx_1)/2;
            end

            vtx_0=edge_u1(0);
            vtx_1=edge_0v(1);
            if norm(vtx_0-vtx_1) > geom_torl
                error('FaceCoons: edge do not connect');
            else
                self.vertex_01=(vtx_0+vtx_1)/2;
            end

            vtx_0=edge_u1(1);
            vtx_1=edge_1v(1);
            if norm(vtx_0-vtx_1) > geom_torl
                error('FaceCoons: edge do not connect');
            else
                self.vertex_11=(vtx_0+vtx_1)/2;
            end

            self.edge_u0=edge_u0;
            self.edge_u1=edge_u1;
            self.edge_0v=edge_0v;
            self.edge_1v=edge_1v;

            % get dimension
            self.dimension=length(vtx_0);
        end
    end

    % calculate point
    methods
        function Point=calPoint(self,U,V)
            [rank_num,colume_num]=size(U);
            Point=zeros(rank_num,colume_num,self.dimension);

            zo_dim=zeros(1,1,self.dimension);

            vtx_00=reshape(self.vertex_00,1,1,3);
            vtx_10=reshape(self.vertex_10,1,1,3);
            vtx_01=reshape(self.vertex_01,1,1,3);
            vtx_11=reshape(self.vertex_11,1,1,3);

            Pnt_u0=self.edge_u0(U);Pnt_u0=reshape(Pnt_u0,[rank_num,colume_num,self.dimension]);
            Pnt_u1=self.edge_u1(U);Pnt_u1=reshape(Pnt_u1,[rank_num,colume_num,self.dimension]);
            Pnt_0v=self.edge_0v(V);Pnt_0v=reshape(Pnt_0v,[rank_num,colume_num,self.dimension]);
            Pnt_1v=self.edge_1v(V);Pnt_1v=reshape(Pnt_1v,[rank_num,colume_num,self.dimension]);

            for rdx=1:rank_num
                for cdx=1:colume_num
                    u_x=U(rdx,cdx);
                    v_x=V(rdx,cdx);

                    blend_u=[-1,1-u_x,u_x];
                    blend_v=[-1;1-v_x;v_x];

                    Point(rdx,cdx,:)=-pagemtimes(pagemtimes(blend_u,[
                        zo_dim,Pnt_u0(rdx,cdx,:),Pnt_u1(rdx,cdx,:);
                        Pnt_0v(rdx,cdx,:),vtx_00,vtx_01;
                        Pnt_1v(rdx,cdx,:),vtx_10,vtx_11]),blend_v);
                end
            end
        end
    end

    % visualizate function
    methods
        function fce=getNURBS(self,u_param,v_param)
            % convert Coons Face into NURBS Face
            %
            % input:
            % u_param, v_param
            % or:
            % value_torl, []
            %
            if nargin < 3
                v_param=[];
                if nargin < 2
                    u_param=[];
                end
            end

            if ~isempty(u_param) && length(u_param) > 1 && u_param(1,1) > u_param(1,2)
                u_param=fliplr(u_param);
            end
            if ~isempty(v_param) && length(v_param) > 1 && v_param(1,1) > v_param(2,1)
                v_param=flipud(v_param);
            end
            [Nodes,U_node,V_node]=calFace(self,u_param,v_param);

            UDegree=min(size(U_node,2)-1,3);VDegree=min(size(V_node,1)-1,3);
            fce=GeomApp.VertexToFace(self.name,Nodes,UDegree,VDegree,[],[],U_node(1,:),V_node(:,1));
        end

        function [step_str,obj_idx,ADVANCED_FACE]=getStep(self,obj_idx)
            % interface of BSpline surface getStep function
            %
            surf_BSpline=self.getNURBS(1e-2);
            [step_str,obj_idx,ADVANCED_FACE]=surf_BSpline.getStep(obj_idx);
        end
    end
end