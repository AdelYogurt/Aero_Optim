classdef Surface < handle
    % Non-Uniform Rational B-Splines Face
    % define reference to step standard
    %
    properties % Explicit Attributes
        name='';
        UDegree; % u_degree
        VDegree; % v_degree
        Poles; % control_points_list
        face_form='.UNSPECIFIED.'; % face_form
        UPeriodic='.F.'; % closed_surface (boolean)
        VPeriodic='.F.'; % closed_surface (boolean)
        self_intersect='.F.'; % self_intersect (boolean)
        Weights; % weights_data

        UMults; % u_knot_multiplicities
        VMults; % v_knot_multiplicities
        UKnots; % UKnots
        VKnots; % VKnots
        knot_spec='.UNSPECIFIED.';

        u_list; % u_knot_vector
        v_list; % v_knot_vector
    end

    % define surface
    methods
        function self=Surface(Poles,UDegree,VDegree,UMults,VMults,UKnots,VKnots,Weights)
            % generate Non-Uniform Rational B-Splines Face
            %
            % input:
            % Poles (matrix): control point, v_pole_num x u_pole_num x dimension matrix
            % UDegree (matrix): optional input
            % VDegree (matrix): optional input
            % UMults (matrix): optional input
            % VMults (matrix): optional input
            % UKnots (matrix): optional input
            % VKnots (matrix): optional input
            % Weights (matrix): optional input
            %
            % output:
            % Surface
            %
            % notice:
            % colume of X, Y, Z is LaWGS format
            % colume of node_X, node_Y, node_Z will by reseve calculate
            % colume direction is u(x), rank direction is v(y)
            %
            if nargin < 8
                Weights=[];
                if nargin < 7
                    VKnots=[];
                    if nargin < 6
                        VMults=[];
                        if nargin < 5
                            UKnots=[];
                            if nargin < 4
                                UMults=[];
                                if nargin < 3
                                    VDegree=[];
                                    if nargin < 2
                                        UDegree=[];
                                    end
                                end
                            end
                        end
                    end
                end
            end

            [v_pole_num,u_pole_num,~]=size(Poles);

            % default value
            if isempty(UDegree),UDegree=u_pole_num-1;end
            if isempty(VDegree),VDegree=v_pole_num-1;end
            if isempty(UMults),UMults=[UDegree+1,ones(1,u_pole_num-UDegree-1),UDegree+1];end;UMults=UMults(:)';
            if isempty(VMults),VMults=[VDegree+1,ones(1,v_pole_num-VDegree-1),VDegree+1];end;VMults=VMults(:)';
            if isempty(UKnots),UKnots=linspace(0,1,u_pole_num-UDegree+1);end;UKnots=UKnots(:)';
            if isempty(VKnots),VKnots=linspace(0,1,v_pole_num-VDegree+1);end;VKnots=VKnots(:)';
            if isempty(Weights), Weights=ones(v_pole_num,u_pole_num);end
            u_list=baseKnotVec(UMults,UKnots);
            v_list=baseKnotVec(VMults,VKnots);

            if u_pole_num < (UDegree+1) || v_pole_num < (VDegree+1)
                error('Surface: pole_num less than Degree+1');
            end

            if length(u_list) ~= u_pole_num+UDegree+1 || length(v_list) ~= v_pole_num+VDegree+1
                error('Surface: knot_num is not equal to pole_num+Degree+1');
            end

            % Explicit Attributes
            self.UDegree=UDegree;
            self.VDegree=VDegree;
            self.Poles=Poles;
            self.UMults=UMults;
            self.VMults=VMults;
            self.UKnots=UKnots;
            self.VKnots=VKnots;
            self.Weights=Weights;

            % Derived Attributes
            self.u_list=u_list;
            self.v_list=v_list;
        end

        function Points=calPoint(self,U_x,V_x)
            % calculate point on Non-Uniform Rational B-Splines Surface
            %
            [v_num,u_num]=size(V_x);
            uv_num=v_num*u_num;
            if any(size(U_x) ~= size(V_x))
                error('Surface.calPoint: U_x and V_x have different size')
            end
            [v_pole_num,u_pole_num,dim]=size(self.Poles);
            Ctrls=cat(3,self.Poles.*self.Weights,self.Weights);

            % evaluate along the v direction
            Ctrls_V=reshape(Ctrls,v_pole_num,(dim+1)*u_pole_num);
            Points_V=calBSpline(Ctrls_V,self.VDegree,self.v_list,V_x(:));
            Points_V=reshape(Points_V,[uv_num,u_pole_num,(dim+1)]);

            % evaluate along the u direction
            Points=zeros(uv_num,(dim+1));
            for uv_idx=1:uv_num
                Ctrls_U=reshape(Points_V(uv_idx,:,:),[u_pole_num,(dim+1)]);
                Points(uv_idx,:)=calBSpline(Ctrls_U,self.UDegree,self.u_list,U_x(uv_idx));
            end
            Points=reshape(Points,[v_num,u_num,(dim+1)]);
            
            Points=Points(:,:,1:end-1)./Points(:,:,end);

            function Points=calBSpline(Ctrls,Degree,u_list,U_x)
                [N_list,idx_srt,idx_end]=baseFcnN(U_x,Degree,u_list);
                Points=zeros(length(U_x),size(Ctrls,2));
                for deg_idx=1:Degree+1
                    Points=Points+N_list(:,deg_idx).*Ctrls(idx_srt+(deg_idx-1),:);
                end
            end

            % [rank_num,colume_num]=size(U);
            % Point=zeros(rank_num,colume_num,size(self.Poles,3));
            % 
            % N_u_list=zeros(1,self.UDegree+1);
            % N_v_list=zeros(self.VDegree+1,1);
            % for rank_idx=1:rank_num
            %     for colume_idx=1:colume_num
            %         % local index of u_x in u_list, v_list
            %         u_x=U(rank_idx,colume_idx);
            %         v_x=V(rank_idx,colume_idx);
            % 
            %         % [index_end_v,index_end_u]=getIndex(); % y, x
            % 
            %         idx_end_u=size(self.Poles,2); % is equal to the section index
            %         while idx_end_u > self.UDegree+1 && u_x < self.u_list(idx_end_u)
            %             idx_end_u=idx_end_u-1;
            %         end
            %         idx_end_v=size(self.Poles,1); % is equal to the section index
            %         while idx_end_v > self.VDegree+1 && v_x < self.v_list(idx_end_v)
            %             idx_end_v=idx_end_v-1;
            %         end
            % 
            %         idx_start_u=idx_end_u-self.UDegree;
            %         idx_start_v=idx_end_v-self.VDegree;
            % 
            %         % calculate base function
            %         for N_idx=1:self.UDegree+1
            %             N_u_list(N_idx)=baseFcnN(self.u_list,u_x,N_idx+idx_start_u-1,self.UDegree);
            %         end
            % 
            %         for N_idx=1:self.VDegree+1
            %             N_v_list(N_idx)=baseFcnN(self.v_list,v_x,N_idx+idx_start_v-1,self.VDegree);
            %         end
            %         NW_list=N_u_list.*self.Weights(idx_start_v:idx_end_v,idx_start_u:idx_end_u).*N_v_list;
            % 
            %         Point(rank_idx,colume_idx,:)=sum(NW_list.*self.Poles(idx_start_v:idx_end_v,idx_start_u:idx_end_u,:),[1,2])/sum(NW_list,'all');
            %     end
            % end
        end

        function srf_hdl=displayGeom(self,axe_hdl,srf_option,u_param,v_param)
            % draw surface on axes handle
            %
            if nargin < 5
                v_param=[];
                if nargin < 4
                    u_param=[];
                    if nargin < 3
                        srf_option=[];
                        if nargin < 2
                            axe_hdl=[];
                        end
                    end
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end

            % default draw option
            if isempty(srf_option), srf_option=struct('LineStyle','none');end

            % calculate point on surface
            if isempty(u_param), u_param=101;end
            if length(u_param) == 1, u_param=linspace(min(self.UKnots),max(self.UKnots),u_param);end
            if isempty(v_param), v_param=101;end
            if length(v_param) == 1, v_param=linspace(min(self.VKnots),max(self.VKnots),v_param);end
            [u_param,v_param]=meshgrid(u_param,v_param);
            Points=self.calPoint(u_param,v_param);

            % draw Points on axe_hdl
            if size(self.Poles,3) == 2
                srf_hdl=surface(axe_hdl,Points(:,:,1),Points(:,:,2),srf_option);
            else
                srf_hdl=surface(axe_hdl,Points(:,:,1),Points(:,:,2),Points(:,:,3),srf_option);
                zlabel('z');
                view(3);
            end
            xlabel('x');
            ylabel('y');
        end
        
        function srf_hdl=displayPoles(self,axe_hdl,pole_option)
            % draw curve on figure handle
            %
            if nargin < 3
                pole_option=[];
                if nargin < 2
                    axe_hdl=[];
                end
            end
            if isempty(axe_hdl),axe_hdl=gca();end

            % default draw option
            if isempty(pole_option)
                pole_option=struct('Marker','s','MarkerEdgeColor','r','EdgeColor','r','LineStyle','--','FaceAlpha',0);
            end

            % draw Poles on axe_hdl
            if size(self.Poles,2) == 2
                srf_hdl=surface(axe_hdl,self.Poles(:,:,1),self.Poles(:,:,2),pole_option);
            else
                srf_hdl=surface(axe_hdl,self.Poles(:,:,1),self.Poles(:,:,2),self.Poles(:,:,3),pole_option);
                zlabel('z');
                view(3);
            end
            xlabel('x');
            ylabel('y');
        end

    end

    % control surface
    methods
        function self=reverseU(self)
            % revese U direction of surface
            %
            self.Poles=self.Poles(:,end:-1:1,:);
            self.UMults=self.UMults(end:-1:1);
            self.UKnots=max(self.UKnots)-self.UKnots(end:-1:1);
            self.Weights=self.Weights(:,end:-1:1,:);
            self.u_list=max(self.u_list)-self.u_list(end:-1:1);
        end

        function self=reverseV(self)
            % revese V direction of surface
            %
            self.Poles=self.Poles(end:-1:1,:,:);
            self.VMults=self.VMults(end:-1:1);
            self.VKnots=max(self.VKnots)-self.VKnots(end:-1:1);
            self.Weights=self.Weights(end:-1:1,:,:);
            self.v_list=max(self.v_list)-self.v_list(end:-1:1);
        end

        function addDegree(self,Udeg_tar,Vdeg_tar)
            % increase surface degree
            %
            if Vdeg_tar <= self.VDegree && Udeg_tar <= self.UDegree, return;end

            [v_pole_num,u_pole_num,dim]=size(self.Poles);
            Ctrls=cat(3,self.Poles.*self.Weights,self.Weights);

            % modify VDegree, VMults and Ctrls
            if self.VDegree < Vdeg_tar
                Ctrls=permute(Ctrls,[1,2,3]);
                Ctrls=reshape(Ctrls,v_pole_num,(dim+1)*u_pole_num);
                [Ctrls,self.VMults]=GeomApp.addDegree(self.VDegree,Ctrls,self.VMults,self.VKnots,Vdeg_tar);
                self.VDegree=Vdeg_tar;
                self.v_list=baseKnotVec(self.VMults,self.VKnots);

                v_pole_num=size(Ctrls,1);
                Ctrls=reshape(Ctrls,[v_pole_num,u_pole_num,(dim+1)]);
                Ctrls=permute(Ctrls,[1,2,3]);
            end

            % modify UDegree, UMults and Ctrls
            if self.UDegree < Udeg_tar
                Ctrls=permute(Ctrls,[2,1,3]);
                Ctrls=reshape(Ctrls,u_pole_num,(dim+1)*v_pole_num);
                [Ctrls,self.UMults]=GeomApp.addDegree(self.UDegree,Ctrls,self.UMults,self.UKnots,Udeg_tar);
                self.VDegree=Vdeg_tar;
                self.u_list=baseKnotVec(self.UMults,self.UKnots);

                u_pole_num=size(Ctrls,1);
                Ctrls=reshape(Ctrls,[u_pole_num,v_pole_num,(dim+1)]);
                Ctrls=permute(Ctrls,[2,1,3]);
            end

            % modify Poles and Weights
            self.Poles=Ctrls(:,:,1:end-1)./Ctrls(:,end);
            self.Weights=Ctrls(:,:,end);
        end

        function insertKnot(self,U_ins,V_ins)
            % insert knot to surface
            %
            if isempty(U_ins) && isempty(V_ins), return;end
            [v_pole_num,u_pole_num,dim]=size(self.Poles);
            Ctrls=[self.Poles.*self.Weights,self.Weights];

            % insert knots along the v direction
            if ~isempty(V_ins)
                Ctrls=permute(Ctrls,[1,2,3]);
                Ctrls=reshape(Ctrls,v_pole_num,(dim+1)*u_pole_num);
                [Ctrls,self.v_list]=GeomApp.insertKnot(self.VDegree,Ctrls,self.v_list,V_ins);
                self.VKnots=unique(self.v_list);
                self.VMults=ones(size(self.VKnots));
                for k_idx=1:length(self.VKnots)
                    self.VMults(k_idx)=sum(self.v_list == self.VKnots(k_idx));
                end

                v_pole_num=size(Ctrls,1);
                Ctrls=reshape(Ctrls,[v_pole_num,u_pole_num,(dim+1)]);
                Ctrls=permute(Ctrls,[1,2,3]);
            end

            % insert knots along the u direction
            if ~isempty(U_ins)
                Ctrls=permute(Ctrls,[2,1,3]);
                Ctrls=reshape(Ctrls,u_pole_num,(dim+1)*v_pole_num);
                [Ctrls,self.u_list]=GeomApp.insertKnot(Udeg,Ctrls,self.u_list,U_ins);
                self.UKnots=unique(self.u_list);
                self.UMults=ones(size(self.UKnots));
                for k_idx=1:length(self.UKnots)
                    self.UMults(k_idx)=sum(self.u_list == self.UKnots(k_idx));
                end

                u_pole_num=size(Ctrls,1);
                Ctrls=reshape(Ctrls,[u_pole_num,v_pole_num,(dim+1)]);
                Ctrls=permute(Ctrls,[2,1,3]);
            end

            % modify Poles and Weights
            self.Poles=Ctrls(:,:,1:end-1)./Ctrls(:,:,end);
            self.Weights=Ctrls(:,:,end);
        end

        function self=translate(self,tran_vctr)
            self.Poles=self.Poles+reshape(tran_vctr,1,1,size(self.Poles,3));
        end

        function self=rotate(self,rotate_matrix)
            [v_num,u_num,dim]=size(self.Poles);
            self.Poles=reshape(reshape(self.Poles,[],size(self.Poles,3))*rotate_matrix',v_num,u_num,dim);
        end

        function crv=getBoundCurve(self,type)
            switch type
                case 'u0' % 1
                    crv=Curve(reshape(self.Poles(1,:,:),size(self.Poles,2),[]),self.UDegree,self.UMults,self.UKnots,self.Weights(1,:));
                case 'u1' % 3.reverse
                    crv=Curve(reshape(self.Poles(end,:,:),size(self.Poles,2),[]),self.UDegree,self.UMults,self.UKnots,self.Weights(end,:));
                case '0v' % 4.reverse
                    crv=Curve(reshape(self.Poles(:,1,:),size(self.Poles,1),[]),self.VDegree,self.VMults,self.VKnots,self.Weights(:,1));
                case '1v' % 2
                    crv=Curve(reshape(self.Poles(:,end,:),size(self.Poles,1),[]),self.VDegree,self.VMults,self.VKnots,self.Weights(:,end));
            end
        end
    end

end