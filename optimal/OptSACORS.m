classdef OptSACORS < handle
    % SACO-RS optimization algorithm
    % only support constraints problem
    %
    % Copyright 2023 Adel
    %
    % abbreviation:
    % obj: objective, con: constraint, iter: iteration, torl: torlance
    % fcn: function, surr: surrogate
    % lib: library, init: initial, rst: restart, potl: potential
    %
    properties
        % basic parameter
        NFE_max;
        iter_max;
        obj_torl;
        con_torl;
        data_lib;

        DRAW_FIGURE_FLAG=0; % whether draw data
        INFORMATION_FLAG=1; % whether print data
        CONVERGENCE_JUDGMENT_FLAG=0; % whether judgment convergence
        WRIRE_FILE_FLAG=0; % whether write to file
        save_description=[]; % iter save mat name

        nomlz_value=10; % max obj when normalize obj,con,coneq
        protect_range=1e-14; % surrogate add point protect range
        identiy_torl=1e-2; % if inf norm of point less than identiy_torlance, point will be consider as same local best
        min_r_interest=1e-3;
        max_r_interest=1e-1;

        str_data_file='result_total.txt'
    end

    properties
        % problem parameter
        GPC_simplify_flag=true(1);
        expensive_con_flag;

        X_local_best=[];
        Obj_local_best=[];

        X_potl=[];
        Obj_potl=[];
        Vio_potl=[];

        detect_local_flag=true(1);

        model_GPC=struct('hyp',struct('mean',0,'cov',[0,0]));

        obj_fcn_surr;
        con_fcn_surr;
        ks_fcn_surr;

        obj_surr;
        con_surr_list;
        coneq_surr_list;
        ks_surr;
    end

    % main function
    methods
        function self=OptSACORS(NFE_max,iter_max,obj_torl,con_torl,data_lib,save_description)
            % initialize optimization
            %
            if nargin < 6
                save_description=[];
                if nargin < 5
                    data_lib=[];
                    if nargin < 4
                        con_torl=[];
                        if nargin < 3
                            obj_torl=[];
                            if nargin < 2
                                iter_max=[];
                                if nargin < 1
                                    NFE_max=[];
                                end
                            end
                        end
                    end
                end
            end

            if isempty(con_torl)
                con_torl=1e-3;
            end
            if isempty(obj_torl)
                obj_torl=1e-6;
            end

            self.NFE_max=NFE_max;
            self.iter_max=iter_max;
            self.obj_torl=obj_torl;
            self.con_torl=con_torl;
            self.data_lib=data_lib;
            self.save_description=save_description;
        end

        function [x_best,obj_best,NFE,output]=optimize(self,varargin)
            % main optimize function
            %
            if length(varargin) == 1
                % input is struct or object
                problem=varargin{1};
                if isstruct(problem)
                    fields=fieldnames(problem);
                    if all(~contains(fields,'objcon_fcn')), error('optSACORS.optimize: input problem lack objcon_fcn'); end
                    objcon_fcn=problem.objcon_fcn;
                    if all(~contains(fields,'vari_num')), error('optSACORS.optimize: input problem lack vari_num'); end
                    if all(~contains(fields,'low_bou')), error('optSACORS.optimize: input problem lack low_bou'); end
                    if all(~contains(fields,'up_bou')), error('optSACORS.optimize: input problem lack up_bou'); end
                    if all(~contains(fields,'con_fcn_cheap')), problem.con_fcn_cheap=[]; end
                    con_fcn_cheap=problem.con_fcn_cheap;
                else
                    m=methods(problem);p=properties(problem);
                    if all(~contains(m,'objcon_fcn')) && all(~contains(p,'objcon_fcn')), error('optSACORS.optimize: input problem lack objcon_fcn'); end
                    objcon_fcn=@problem.objcon_fcn;
                    if all(~contains(p,'vari_num')), error('optSACORS.optimize: input problem lack vari_num'); end
                    if all(~contains(p,'low_bou')), error('optSACORS.optimize: input problem lack low_bou'); end
                    if all(~contains(p,'up_bou')), error('optSACORS.optimize: input problem lack up_bou'); end
                    if all(~contains(m,'con_fcn_cheap')) && all(~contains(p,'con_fcn_cheap'))
                        con_fcn_cheap=[];
                    else
                        con_fcn_cheap=@problem.con_fcn_cheap;
                    end
                end
                vari_num=problem.vari_num;
                low_bou=problem.low_bou;
                up_bou=problem.up_bou;
            else
                % multi input
                varargin=[varargin,repmat({[]},1,5-length(varargin))];
                [objcon_fcn,vari_num,low_bou,up_bou,con_fcn_cheap]=varargin{:};
            end

            % NFE and iteration setting
            if isempty(self.NFE_max)
                self.NFE_max=10+10*vari_num;
            end

            if isempty(self.iter_max)
                self.iter_max=20+20*vari_num;
            end

            % hyper parameter
            sample_num_init=6+3*vari_num;

            sample_num_rst=sample_num_init;
            sample_num_add=ceil(log(sample_num_init)/2);

            trial_num=min(100*vari_num,100);

            done=0;NFE=0;iter=0;

            % step 1
            % generate initial data library
            if isempty(self.data_lib)
                self.data_lib=struct('objcon_fcn',objcon_fcn,'vari_num',vari_num,'low_bou',low_bou,'up_bou',up_bou,'con_torl',self.con_torl,...
                    'result_best_idx',[],'X',[],'Obj',[],'Con',[],'Coneq',[],'Vio',[],'Ks',[]);
            else
                % load exist data
                self.data_lib=struct('objcon_fcn',objcon_fcn,'vari_num',vari_num,'low_bou',low_bou,'up_bou',up_bou,'con_torl',self.con_torl,...
                    'result_best_idx',[],'X',self.data_lib.X,'Obj',self.data_lib.Obj,'Con',self.data_lib.Con,'Coneq',self.data_lib.Coneq,'Vio',self.data_lib.Vio,'Ks',self.data_lib.Ks);
            end

            if size(self.data_lib.X,1) < sample_num_init
                % use latin hypercube method to get enough initial sample x_list
                X_updata=lhsdesign(sample_num_init-size(self.data_lib.X,1),vari_num,'iterations',50,'criterion','maximin').*(up_bou-low_bou)+low_bou;

                % updata data lib by x_list
                [self.data_lib,~,~,~,~,~,~,~,NFE_updata]=dataUpdata(self.data_lib,X_updata,0,self.WRIRE_FILE_FLAG,self.str_data_file);
                NFE=NFE+NFE_updata;
            end

            % detech expensive constraints
            if ~isempty(self.data_lib.Vio)
                self.expensive_con_flag=true(1);
            else
                self.expensive_con_flag=false(1);
            end

            % find fesiable data in current data lib
            if self.expensive_con_flag
                Bool_feas=self.data_lib.Vio == 0;
            end
            Bool_conv=false(size(self.data_lib.X,1),1);

            fmincon_options=optimoptions('fmincon','Display','none','Algorithm','sqp','ConstraintTolerance',0);

            result_x_best=zeros(self.iter_max,vari_num);
            result_obj_best=zeros(self.iter_max,1);

            iter=iter+1;

            while ~done
                % step 2
                % nomalization all data by max obj and to create surrogate objcon_fcn
                [X,Obj,Con,Coneq,Vio,Ks]=dataLoad(self.data_lib);

                [X_surr,Obj_surr,Con_surr,Coneq_surr,Vio_surr,Ks_surr,...
                    obj_max,con_max_list,coneq_max_list,vio_max_list,ks_max_list]=self.getModelData...
                    (X,Obj,Con,Coneq,Vio,Ks,self.nomlz_value);

                % construct surrogate objcon_fcn
                [self.obj_fcn_surr,self.con_fcn_surr,output_model]=self.getSurrFcnRBF...
                    (X_surr,Obj_surr,Con_surr,Coneq_surr);
                self.obj_surr=output_model.obj_surr;
                self.con_surr_list=output_model.con_surr_list;
                self.coneq_surr_list=output_model.coneq_surr_list;
                [self.ks_fcn_surr,self.ks_surr]=self.interpRadialBasisPreModel(X_surr,Ks_surr);

                % step 3
                % updata or identify potential local best
                self.updataLocalPotential...
                    (X_surr,Vio_surr,vari_num,low_bou,up_bou,con_fcn_cheap,fmincon_options,...
                    obj_max,con_max_list,coneq_max_list);

                % step 4
                % select best potential point as x_infill
                x_infill=self.X_potl(1,:);
                obj_infill_pred=self.Obj_potl(1,:);

                % updata infill point
                [self.data_lib,x_infill,obj_infill,con_infill,coneq_infill,vio_infill,ks_infill,repeat_idx,NFE_updata]=...
                    dataUpdata(self.data_lib,x_infill,self.protect_range,self.WRIRE_FILE_FLAG,self.str_data_file);
                NFE=NFE+NFE_updata;

                if isempty(x_infill)
                    % process error
                    x_infill=self.data_lib.X(repeat_idx,:);
                    obj_infill=self.data_lib.Obj(repeat_idx,:);
                    if ~isempty(Con)
                        con_infill=self.data_lib.Con(repeat_idx,:);
                    end
                    if ~isempty(Coneq)
                        coneq_infill=self.data_lib.Coneq(repeat_idx,:);
                    end
                    if ~isempty(Vio)
                        vio_infill=self.data_lib.Vio(repeat_idx,:);
                    end
                else
                    if ~isempty(vio_infill) && vio_infill > 0
                        Bool_feas=[Bool_feas;false(1)];
                    else
                        Bool_feas=[Bool_feas;true(1)];
                    end
                    Bool_conv=[Bool_conv;false(1)];
                end
                self.Obj_potl(1,:)=obj_infill;
                self.Vio_potl(1,:)=vio_infill;

                if self.DRAW_FIGURE_FLAG && vari_num < 3
                    interpVisualize(self.obj_surr,low_bou,up_bou);
                    line(x_infill(1),x_infill(2),obj_infill./obj_max*self.nomlz_value,'Marker','o','color','r','LineStyle','none');
                    view(2);
                end

                % find best result to record
                [X,Obj,~,~,Vio,~]=dataLoad(self.data_lib);
                X_unconv=X(~Bool_conv,:);
                Obj_unconv=Obj(~Bool_conv,:);
                if ~isempty(Vio)
                    Vio_unconv=Vio(~Bool_conv,:);
                else
                    Vio_unconv=[];
                end
                idx=find(Vio_unconv == 0);
                if isempty(idx)
                    [vio_best,min_idx]=min(Vio_unconv);
                    obj_best=Obj_unconv(min_idx);
                    x_best=X_unconv(min_idx,:);
                else
                    [obj_best,min_idx]=min(Obj_unconv(idx));
                    vio_best=0;
                    x_best=X_unconv(idx(min_idx),:);
                end

                if self.INFORMATION_FLAG
                    fprintf('obj:    %f    violation:    %f    NFE:    %-3d\n',obj_best,vio_best,NFE);
                    %         fprintf('iteration:          %-3d    NFE:    %-3d\n',iteration,NFE);
                    %         fprintf('x:          %s\n',num2str(x_infill));
                    %         fprintf('value:      %f\n',obj_infill);
                    %         fprintf('violation:  %s  %s\n',num2str(con_infill),num2str(coneq_infill));
                    %         fprintf('\n');
                end

                result_x_best(iter,:)=x_best;
                result_obj_best(iter,:)=obj_best;
                iter=iter+1;

                % forced interrupt
                if iter > self.iter_max || NFE >= self.NFE_max
                    done=1;
                end

                % convergence judgment
                if self.CONVERGENCE_JUDGMENT_FLAG
                    if ( ((iter > 2) && (abs((obj_infill-obj_infill_old)/obj_infill_old) < self.torl)) && ...
                            ((~isempty(vio_infill) && vio_infill == 0) || isempty(vio_infill)) )
                        done=1;
                    end
                end

                if ~done
                    [X,Obj,~,~,Vio,~]=dataLoad(self.data_lib);
                    X_add=[];

                    % check if converage
                    if  ( ((iter > 2) && (abs((obj_infill-obj_infill_old)/obj_infill_old) < self.obj_torl)) ...
                            && ((~isempty(vio_infill) && vio_infill == 0) || isempty(vio_infill)))
                        % resample LHD
                        % step 6.1
                        self.X_local_best=[self.X_local_best;x_infill];
                        self.Obj_local_best=[self.Obj_local_best;obj_infill];
                        self.X_potl=self.X_potl(2:end,:);
                        self.Obj_potl=self.Obj_potl(2:end,:);
                        self.Vio_potl=self.Vio_potl(2:end,:);

                        if isempty(self.X_potl)
                            self.detect_local_flag=true(1);
                        end

                        % step 6.2
                        % detech converage
                        for x_idx=1:size(X,1)
                            if ~Bool_conv(x_idx)
                                x_single_pred=fmincon(self.obj_fcn_surr,X(x_idx,:),[],[],[],[],low_bou,up_bou,self.con_fcn_surr,fmincon_options);

                                converage_flag=false(1);

                                for x_check_idx=1:size(self.X_local_best,1)
                                    if sum(abs(self.X_local_best(x_check_idx,:)-x_single_pred),2)/vari_num < self.identiy_torl
                                        converage_flag=true(1);
                                        break;
                                    end
                                end

                                if converage_flag
                                    % if converage to local minimum, set to infeasible
                                    Bool_conv(x_idx)=true(1);
                                end
                            end
                        end

                        % step 6.3
                        % use GPC to limit do not converage to exist local best
                        if ~all(Bool_conv)
                            Class=-1*ones(size(X,1),1);
                            Class(Bool_conv)=1; % cannot go into converage area

                            [pred_fcn_GPC,self.model_GPC]=self.classifyGaussProcess(X,Class,self.model_GPC.hyp,self.GPC_simplify_flag);
                            con_GPC_fcn=@(x) self.conFcnGPC(x,pred_fcn_GPC);
                        else
                            con_GPC_fcn=[];
                        end

                        % step 6.4
                        % resample latin hypercubic and updata into data lib
                        try
                            X_add=self.lhdESLHS(min(floor(sample_num_rst),self.NFE_max-NFE),vari_num,...
                                low_bou,up_bou,X,con_GPC_fcn);
                        catch
                            X_add=lhsdesign(min(floor(sample_num_rst),self.NFE_max-NFE),vari_num).*(up_bou-low_bou)+low_bou;
                        end

                        if self.DRAW_FIGURE_FLAG && vari_num < 3
                            classifyVisualization(self.model_GPC,low_bou,up_bou);
                            line(X(:,1),X(:,2),'Marker','o','color','k','LineStyle','none');
                        end
                    else
                        % step 5.1
                        % check if improve
                        improve=0;
                        if isempty(repeat_idx)
                            Bool_comp=(~Bool_conv)&Bool_feas;
                            Bool_comp(end)=false(1);
                            if self.expensive_con_flag
                                min_vio=min(Vio(~Bool_conv(1:end-1)));
                                min_obj=min(Obj(Bool_comp));

                                % if all point is infeasible,violation of point infilled is
                                % less than min violation of all point means improve.if
                                % feasible point exist,obj of point infilled is less than min
                                % obj means improve
                                if vio_infill == 0 || vio_infill < min_vio
                                    if ~isempty(min_obj)
                                        if obj_infill < min_obj
                                            % improve, continue local search
                                            improve=1;
                                        end
                                    else
                                        % improve, continue local search
                                        improve=1;
                                    end
                                end
                            else
                                min_obj=min(Obj(~Bool_conv));

                                % obj of point infilled is less than min obj means improve
                                if obj_infill < min_obj
                                    % imporve, continue local search
                                    improve=1;
                                end
                            end
                        end

                        % step 5.2
                        % if obj no improve, use GPC to identify interest area
                        % than, imporve interest area surrogate quality
                        if ~improve
                            % construct GPC
                            train_num=min(size(self.data_lib.X,1),11*vari_num-1+25);
                            [pred_fcn_GPC,x_pareto_center,pareto_idx_list]=self.trainFilter(x_infill,train_num,Bool_conv);
                            con_GPC_fcn=@(x) self.conFcnGPC(x,pred_fcn_GPC);

                            % step 5.3
                            % identify interest area
                            center_point=fmincon(con_GPC_fcn,x_pareto_center,[],[],[],[],low_bou,up_bou,con_fcn_cheap,fmincon_options);

                            r_interest=sum(abs(center_point-x_infill)./(up_bou-low_bou))/vari_num;
                            r_interest=max(self.min_r_interest,r_interest);
                            r_interest=min(self.max_r_interest,r_interest);

                            % generate trial point
                            trial_point=repmat(x_infill,trial_num,1);
                            for variable_idx=1:vari_num
                                trial_point(:,variable_idx)=trial_point(:,variable_idx)+...
                                    normrnd(0,r_interest,[trial_num,1]);
                            end
                            trial_point=max(trial_point,low_bou);
                            trial_point=min(trial_point,up_bou);

                            Bool_negetive=pred_fcn_GPC(trial_point) == -1;
                            if sum(Bool_negetive) < sample_num_add
                                value=con_GPC_fcn(trial_point);
                                thres=quantile(value,0.25);
                                Bool_negetive=value<thres;
                            end
                            trial_point_filter=trial_point(Bool_negetive,:);

                            % step 5.4
                            % select point
                            if size(trial_point_filter,1) <= min(sample_num_add,self.NFE_max-NFE)
                                X_add=trial_point_filter;
                            else
                                max_dist=0;
                                iteration_select=1;
                                while iteration_select < 100
                                    select_idx=randperm(size(trial_point_filter,1),min(sample_num_add,self.NFE_max-NFE));
                                    dist=self.calMinDistanceIter(trial_point_filter(select_idx,:),X);
                                    if max_dist < dist
                                        X_add=trial_point_filter(select_idx,:);
                                        max_dist=dist;
                                    end

                                    iteration_select=iteration_select+1;
                                end
                            end

                            if self.DRAW_FIGURE_FLAG && vari_num < 3
                                classifyVisualization(self.model_GPC,low_bou,up_bou);
                                line(trial_point(:,1),trial_point(:,2),'Marker','o','color','k','LineStyle','none');
                                line(X_add(:,1),X_add(:,2),'Marker','o','color','g','LineStyle','none');
                            end
                        end
                    end

                    % step 7
                    % updata data lib
                    [self.data_lib,X_add,~,~,~,Vio_add,~,~,NFE_updata]=...
                        dataUpdata(self.data_lib,X_add,self.protect_range,self.WRIRE_FILE_FLAG,self.str_data_file);
                    NFE=NFE+NFE_updata;
                    Bool_feas=[Bool_feas;Vio_add == 0];
                    Bool_conv=[Bool_conv;false(size(X_add,1),1)];

                    % forced interrupt
                    if iter > self.iter_max || NFE >= self.NFE_max
                        done=1;
                    end
                end

                obj_infill_old=obj_infill;

                % save iteration
                if ~isempty(self.save_description)
                    data_lib=self.data_lib;
                    save(self.save_description,'data_lib');
                end
            end

            % find best result to record
            x_best=self.data_lib.X(self.data_lib.result_best_idx(end),:);
            obj_best=self.data_lib.Obj(self.data_lib.result_best_idx(end),:);

            result_x_best=result_x_best(1:iter-1,:);
            result_obj_best=result_obj_best(1:iter-1);

            output.result_x_best=result_x_best;
            output.result_obj_best=result_obj_best;

            output.x_local_best=self.X_local_best;
            output.obj_local_best=self.Obj_local_best;

            output.data_lib=self.data_lib;
        end

        function updataLocalPotential...
                (self,X_surr,Vio_model,vari_num,low_bou,up_bou,con_fcn_cheap,fmincon_options,...
                obj_max,con_max_list,coneq_max_list)
            if ~isempty(self.con_fcn_surr) || ~isempty(con_fcn_cheap)
                con_fcn=@(x) self.conFcnTotal...
                    (x,self.con_fcn_surr,con_fcn_cheap);
            else
                con_fcn=[];
            end

            if self.detect_local_flag
                % detech potential local best point
                for x_idx=1:size(X_surr,1)

                    x_initial=X_surr(x_idx,:);
                    [x_potl,obj_potl_pred,exit_flag,output_fmincon]=fmincon(self.obj_fcn_surr,x_initial,[],[],[],[],...
                        low_bou,up_bou,con_fcn,fmincon_options);

                    if exit_flag == 1 || exit_flag == 2
                        % check if x_potl have existed
                        add_flag=true(1);
                        for x_check_idx=1:size(self.X_potl,1)
                            if sum(abs(self.X_potl(x_check_idx,:)-x_potl),2)/vari_num < self.identiy_torl
                                add_flag=false(1);
                                break;
                            end
                        end
                        for x_check_idx=1:size(self.X_local_best,1)
                            if sum(abs(self.X_local_best(x_check_idx,:)-x_potl),2)/vari_num < self.identiy_torl
                                add_flag=false(1);
                                break;
                            end
                        end

                        % updata into self.X_potl
                        if add_flag
                            self.X_potl=[self.X_potl;x_potl];
                            self.Obj_potl=[self.Obj_potl;obj_potl_pred/self.nomlz_value.*obj_max];
                            [con_potential,coneq_potential]=self.con_fcn_surr(x_potl);
                            if ~isempty(con_potential)
                                con_potential=con_potential/self.nomlz_value.*con_max_list;
                            end
                            if ~isempty(coneq_potential)
                                coneq_potential=coneq_potential/self.nomlz_value.*coneq_max_list;
                            end
                            self.Vio_potl=[self.Vio_potl;calViolation(con_potential,coneq_potential,self.con_torl)];
                        end
                    end
                end

                % if self.X_potl is empty, try to use KS surrogate as x potential
                if isempty(self.X_potl)
                    [~,x_idx]=min(Vio_model);
                    x_initial=X_surr(x_idx,:);
                    [x_potl,~,exit_flag,output_fmincon]=fmincon(self.ks_fcn_surr,x_initial,[],[],[],[],...
                        low_bou,up_bou,con_fcn_cheap,fmincon_options);
                    obj_potl_pred=self.obj_fcn_surr(x_potl);

                    self.X_potl=[self.X_potl;x_potl];
                    self.Obj_potl=[self.Obj_potl;obj_potl_pred/self.nomlz_value*obj_max];
                    [con_potential,coneq_potential]=self.con_fcn_surr(x_potl);
                    if ~isempty(con_potential)
                        con_potential=con_potential/self.nomlz_value.*con_max_list;
                    end
                    if ~isempty(coneq_potential)
                        coneq_potential=coneq_potential/self.nomlz_value.*coneq_max_list;
                    end
                    self.Vio_potl=[self.Vio_potl;calViolation(con_potential,coneq_potential,self.con_torl)];
                end

                % sort X potential by Vio
                [self.Vio_potl,idx]=sort(self.Vio_potl);
                self.Obj_potl=self.Obj_potl(idx,:);
                self.X_potl=self.X_potl(idx,:);

                self.detect_local_flag=false(1);
            else
                % updata X potential
                for x_idx=1:size(self.X_potl,1)
                    x_potl=self.X_potl(x_idx,:);

                    [x_potl,obj_potl_pred,exit_flag,output_fmincon]=fmincon(self.obj_fcn_surr,x_potl,[],[],[],[],...
                        low_bou,up_bou,con_fcn,fmincon_options);

                    self.X_potl(x_idx,:)=x_potl;
                    self.Obj_potl(x_idx,:)=obj_potl_pred/self.nomlz_value*obj_max;
                    [con_potential,coneq_potential]=self.con_fcn_surr(x_potl);
                    if ~isempty(con_potential)
                        con_potential=con_potential/self.nomlz_value.*con_max_list;
                    end
                    if ~isempty(coneq_potential)
                        coneq_potential=coneq_potential/self.nomlz_value.*coneq_max_list;
                    end
                    self.Vio_potl(x_idx,:)=calViolation(con_potential,coneq_potential,self.con_torl);
                end

                % merge X potential
                % Upward merge
                for x_idx=size(self.X_potl,1):-1:1
                    x_potl=self.X_potl(x_idx,:);

                    % check if x_potl have existed
                    merge_flag=false(1);
                    for x_check_idx=1:x_idx-1
                        if sum(abs(self.X_potl(x_check_idx,:)-x_potl),2)/vari_num < self.identiy_torl
                            merge_flag=true(1);
                            break;
                        end
                    end

                    % updata into self.X_potl
                    if merge_flag
                        self.X_potl(x_idx,:)=[];
                        self.Obj_potl(x_idx,:)=[];
                        self.Vio_potl(x_idx,:)=[];
                    end
                end
            end

            % sort X potential by Obj
            flag=find(self.Vio_potl == 0, 1, 'last' );
            if isempty(flag) % mean do not have fesaible point
                [self.Obj_potl,idx]=sort(self.Obj_potl);
                self.Vio_potl=self.Vio_potl(idx,:);
                self.X_potl=self.X_potl(idx,:);
            else
                [self.Obj_potl(1:flag,:),idx_feas]=sort(self.Obj_potl(1:flag,:));
                [self.Obj_potl(flag+1:end,:),idx_infeas]=sort(self.Obj_potl(flag+1:end,:));
                idx=[idx_feas;idx_infeas+flag];

                self.Vio_potl=self.Vio_potl(idx,:);
                self.X_potl=self.X_potl(idx,:);
            end

        end

        function [pred_func_GPC,x_pareto_center,pareto_idx_list]=trainFilter(self,x_infill,train_num,Bool_conv)
            % train filter of gaussian process classifer
            %
            [X,Obj,~,~,~,Ks]=dataLoad(self.data_lib);
            % base on distance select usable point
            x_dist=sum(abs(X-x_infill),2);
            [~,idx]=sort(x_dist);
            Obj=Obj(idx(1:train_num),:);
            Ks=Ks(idx(1:train_num),:);
            X=X(idx(1:train_num),:);
            Bool_conv=Bool_conv(idx(1:train_num),:);

            if self.expensive_con_flag
                % base on filter to decide which x should be choose
                %     pareto_idx_list=self.getParetoFront([Obj(~Bool_feas),Ks(~Bool_feas)]);
                pareto_idx_list=self.getParetoFront(Obj,Ks);

                Class=ones(size(X,1),1);
                Class(pareto_idx_list)=-1;
                %     Class(Bool_feas)=-1; % can go into feasiable area
                Class(Bool_conv)=1; % cannot go into convarage area

                x_pareto_center=sum(X(pareto_idx_list,:),1)/length(pareto_idx_list);

                [pred_func_GPC,self.model_GPC]=self.classifyGaussProcess(X,Class,self.model_GPC.hyp,self.GPC_simplify_flag);
            else
                obj_threshold=prctile(Obj,50-40*sqrt(NFE/self.NFE_max));

                Class=ones(size(X,1),1);
                Class(Obj < obj_threshold)=-1;
                Class(Bool_conv)=1; % cannot go into convarage area

                x_pareto_center=sum(X(Obj < obj_threshold,:),1)/sum(Obj < obj_threshold);

                [pred_func_GPC,self.model_GPC]=self.classifyGaussProcess(X,Class,self.model_GPC.hyp,self.GPC_simplify_flag);
            end

            pareto_idx_list=idx(pareto_idx_list);
        end

        function [obj_fcn_surr,con_fcn_surr,output]=getSurrFcnRBF...
                (self,x_list,obj_list,con_list,coneq_list)
            % generate surrogate function of objective and constraints
            %
            % output:
            % obj_fcn_surr(output is obj_pred),...
            % con_fcn_surr(output is con_pred, coneq_pred)
            %

            % generate obj surrogate
            [pred_fcn_obj,obj_RBF]=self.interpRadialBasisPreModel(x_list,obj_list);

            % generate con surrogate
            if ~isempty(con_list)
                pred_fcn_con_list=cell(size(con_list,2),1);
                con_RBF_list=cell(size(con_list,2),1);

                for con_idx=1:size(con_list,2)
                    [pred_fcn_con_list{con_idx},con_RBF_list{con_idx}]=self.interpRadialBasisPreModel...
                        (x_list,con_list(:,con_idx));
                end
            else
                pred_fcn_con_list=[];
                con_RBF_list=[];
            end

            % generate coneq surrogate
            if ~isempty(coneq_list)
                pred_fcn_coneq_list=cell(size(coneq_list,2),1);
                coneq_RBF_list=cell(size(coneq_list,2),1);

                for coneq_idx=1:size(coneq_list,2)
                    [pred_fcn_coneq_list{coneq_idx},coneq_RBF_list{coneq_idx}]=self.interpRadialBasisPreModel...
                        (x_list,coneq_list(:,coneq_idx));
                end
            else
                pred_fcn_coneq_list=[];
                coneq_RBF_list=[];
            end

            obj_fcn_surr=@(X_pred) objFcnSurr(X_pred,pred_fcn_obj);
            if isempty(pred_fcn_con_list) && isempty(pred_fcn_coneq_list)
                con_fcn_surr=[];
            else
                con_fcn_surr=@(X_pred) conFcnSurr(X_pred,pred_fcn_con_list,pred_fcn_coneq_list);
            end

            output.obj_surr=obj_RBF;
            output.con_surr_list=con_RBF_list;
            output.coneq_surr_list=coneq_RBF_list;

            function obj_pred=objFcnSurr...
                    (X_pred,pred_fcn_obj)
                % connect all predict favl
                %
                obj_pred=pred_fcn_obj(X_pred);
            end

            function [con_pred,coneq_pred]=conFcnSurr...
                    (X_pred,pred_fcn_con,pred_fcn_coneq)
                % connect all predict con and coneq
                %
                if isempty(pred_fcn_con)
                    con_pred=[];
                else
                    con_pred=zeros(size(X_pred,1),length(pred_fcn_con));
                    for con_idx__=1:length(pred_fcn_con)
                        con_pred(:,con_idx__)=....
                            pred_fcn_con{con_idx__}(X_pred);
                    end
                end
                if isempty(pred_fcn_coneq)
                    coneq_pred=[];
                else
                    coneq_pred=zeros(size(X_pred,1),length(pred_fcn_coneq));
                    for coneq_idx__=1:length(pred_fcn_coneq)
                        coneq_pred(:,coneq_idx__)=...
                            pred_fcn_coneq{coneq_idx__}(X_pred);
                    end
                end
            end
        end

    end

    % auxiliary function
    methods(Static)
        function [X_surr,Obj_surr,Con_surr,Coneq_surr,Vio_model,Ks_model,...
                obj_max,con_max_list,coneq_max_list,vio_max_list,ks_max_list]=getModelData...
                (X,Obj,Con,Coneq,Vio,Ks,nomlz_obj)
            % normalize data to construct surrogate model
            %
            X_surr=X;
            obj_max=max(abs(Obj),[],1);
            Obj_surr=Obj/obj_max*nomlz_obj;
            if ~isempty(Con)
                con_max_list=max(abs(Con),[],1);
                Con_surr=Con./con_max_list*nomlz_obj;
            else
                con_max_list=[];
                Con_surr=[];
            end
            if ~isempty(Coneq)
                coneq_max_list=max(abs(Coneq),[],1);
                Coneq_surr=Coneq./coneq_max_list*nomlz_obj;
            else
                coneq_max_list=[];
                Coneq_surr=[];
            end
            if ~isempty(Vio)
                vio_max_list=max(abs(Vio),[],1);
                Vio_model=Vio./vio_max_list*nomlz_obj;
            else
                vio_max_list=[];
                Vio_model=[];
            end
            if ~isempty(Ks)
                ks_max_list=max(abs(Ks),[],1);
                Ks_model=Ks./ks_max_list*nomlz_obj;
            else
                ks_max_list=[];
                Ks_model=[];
            end

        end

        function ND_idx_list=getNoDominate(data_list)
            % distinguish no dominate front of data list
            % data_list is x_number x data_number matrix
            % notice if all data of x1 is less than x2,x1 domain x2
            %
            % notice: there are any way to find pareto base on followed
            % rule: point of no dominate will no be dominate by all point
            %
            x_number=size(data_list,1);
            ND_idx_list=[]; % sort all idx of filter point list

            % select no domain filter
            for x_idx=1:x_number
                data=data_list(x_idx,:);
                pareto_idx=1;
                add_filter_flag=1;
                while pareto_idx <= length(ND_idx_list)
                    % compare x with exit pareto front point
                    x_pareto_idx=ND_idx_list(pareto_idx,:);

                    % contain constraint of x_filter
                    data_pareto=data_list(x_pareto_idx,:);

                    % compare x with x_pareto
                    judge=data > data_pareto;
                    if all(judge)
                        add_filter_flag=0;
                        break;
                    end

                    % if better than exit pareto point,reject pareto point
                    judge=data < data_pareto;
                    if all(judge)
                        ND_idx_list(pareto_idx)=[];
                        pareto_idx=pareto_idx-1;
                    end

                    pareto_idx=pareto_idx+1;
                end

                % add into pareto list if possible
                if add_filter_flag
                    ND_idx_list=[ND_idx_list;x_idx];
                end
            end
        end

        function pareto_idx_list=getParetoFront(obj_list,ks_list)
            % distinguish pareto front of data list
            % dominate define as followed
            % Solution i is feasible and solution j is not.
            % Solutions i and j are both infeasible,...
            % but solution i has a smaller overall constraint violation.
            % Solutions i and j are feasible and solution i dominates solution j
            %
            x_number=size(obj_list,1);
            pareto_idx_list=[]; % sort all idx of filter point list

            % select no domain filter
            for x_idx=1:x_number
                obj=obj_list(x_idx,:);
                ks=ks_list(x_idx,:);

                pareto_idx=1;
                add_filter_flag=true(1);
                while pareto_idx <= length(pareto_idx_list)
                    % compare x with exit pareto front point
                    x_pareto_idx=pareto_idx_list(pareto_idx,:);

                    % contain constraint of x_filter
                    obj_pareto=obj_list(x_pareto_idx,:);
                    ks_pareto=ks_list(x_pareto_idx,:);

                    % compare x with x_pareto
                    if ks_pareto <= 0
                        if obj > obj_pareto || ks > 0
                            add_filter_flag=false(1);
                            break;
                        end
                    else
                        if obj > obj_pareto && ks > ks_pareto
                            add_filter_flag=false(1);
                            break;
                        end
                    end

                    % if better than exit pareto point,reject pareto point
                    delete_filter_flag=false(1);
                    if ks <= 0
                        if obj_pareto > obj || ks_pareto > 0
                            delete_filter_flag=true(1);
                        end
                    else
                        if obj_pareto > obj && ks_pareto > ks
                            delete_filter_flag=true(1);
                        end
                    end
                    if delete_filter_flag
                        pareto_idx_list(pareto_idx)=[];
                        pareto_idx=pareto_idx-1;
                    end

                    pareto_idx=pareto_idx+1;
                end

                % add into pareto list if possible
                if add_filter_flag
                    pareto_idx_list=[pareto_idx_list;x_idx];
                end
            end

        end

        function [con,coneq]=conFcnTotal...
                (x,con_fcn,con_fcn_cheap)
            con=[];
            coneq=[];

            if ~isempty(con_fcn)
                [expencon,expenconeq]=con_fcn(x);
                con=[con,expencon];
                coneq=[coneq,expenconeq];
            end

            if ~isempty(con_fcn_cheap)
                [expencon,expenconeq]=con_fcn_cheap(x);
                con=[con,expencon];
                coneq=[coneq,expenconeq];
            end

        end

        function [con,coneq]=conFcnGPC(x,pred_fcn_GPC)
            % function to obtain probability predict function
            %
            [~,~,con]=pred_fcn_GPC(x);
            coneq=[];
        end

        function dist_min=calMinDistanceIter(X,X_exist)
            % get distance min from x_list
            % this version do not consider distance between x exist
            %

            % sort x_supply_list_initial to decrese distance calculate times
            X=sortrows(X,1);
            [sample_num,vari_num]=size(X);
            dist_min=vari_num;
            for x_idx=1:sample_num
                x_curr=X(x_idx,:);
                x_next_idx=x_idx + 1;
                % only search in min_distance(x_list had been sort)
                search_range=vari_num;
                while x_next_idx <= sample_num &&...
                        (X(x_next_idx,1)-X(x_idx,1))^2 ...
                        < search_range
                    x_next=X(x_next_idx,:);
                    distance_temp=sum((x_next-x_curr).^2);
                    if distance_temp < dist_min
                        dist_min=distance_temp;
                    end
                    if distance_temp < search_range
                        search_range=distance_temp;
                    end
                    x_next_idx=x_next_idx+1;
                end
                for x_exist_idx=1:size(X_exist,1)
                    x_next=X_exist(x_exist_idx,:);
                    distance_temp=sum((x_next-x_curr).^2);
                    if distance_temp < dist_min
                        dist_min=distance_temp;
                    end
                end
            end
            dist_min=sqrt(dist_min);
        end

    end

    % machine learning
    methods(Static)
        function [pred_fcn,model_GPC]=classifyGaussProcess...
                (X,Class,hyp,simplify_flag)
            % generate gaussian process classifier model
            % version 7, this version is assembly of gpml-3.6 EP method
            % only support binary classification, -1 and 1
            %
            % input:
            % X(x_num x vari_num matrix), Class(x_num x 1 matrix),...
            % hyp(mean, cov(len, eta)), simplify_flag(means whether use one hyerparameter of all variable)
            %
            % output:
            % pred_fcn, model_GPC
            %
            % abbreviation:
            % pred: predicted, nomlz: normalization,num: number
            % var: variance, fcn: function
            %
            [~,variable_number]=size(X);
            if nargin < 4
                simplify_flag=false(0);
            end

            if nargin < 3 || isempty(hyp)
                hyp.mean=0;
                if simplify_flag
                    hyp.cov=zeros(1,2);
                else
                    hyp.cov=zeros(1,variable_number+1);
                end
            end

            % normalization data
            aver_X=mean(X);
            stdD_X=std(X);
            index__=find(stdD_X == 0);
            if  ~isempty(index__),stdD_X(index__)=1; end
            X_nomlz=(X-aver_X)./stdD_X;

            object_function=@(x) objectNLLGPC(x,{@ OptSACORS.infEP},{@ OptSACORS.meanConst},{@calCov},{@ OptSACORS.likErf},X_nomlz,Class);
            hyp_x=[hyp.mean,hyp.cov];

            % [fval,gradient]=object_function(hyp_x)
            % [fval_differ,gradient_differ]=differ(object_function,hyp_x)

            if simplify_flag
                low_bou_hyp=-3*ones(1,3);
                up_bou_hyp=3*ones(1,3);
            else
                low_bou_hyp=-3*ones(1,variable_number+2);
                up_bou_hyp=3*ones(1,variable_number+2);
            end

            hyp_x=fmincon(object_function,hyp_x,[],[],[],[],low_bou_hyp,up_bou_hyp,[],...
                optimoptions('fmincon','Display','none','SpecifyObjectiveGradient',true,...
                'MaxFunctionEvaluations',20,'OptimalityTolerance',1e-6));

            %             hyp.mean=hyp_x(1);
            hyp.mean=0;
            hyp.cov=hyp_x(2:end);
            hyp.lik=[];
            post=OptSACORS.infEP(hyp,{@ OptSACORS.meanConst},{@calCov},{@ OptSACORS.likErf},X_nomlz,Class);
            pred_fcn=@(x_pred) classifyGaussPredictor...
                (x_pred,hyp,{@ OptSACORS.meanConst},{@calCov},{@ OptSACORS.likErf},post,X_nomlz,aver_X,stdD_X);

            % output model
            model_GPC.X=X;
            model_GPC.Class=Class;
            model_GPC.X_nomlz=X_nomlz;
            model_GPC.aver_X=aver_X;
            model_GPC.stdD_X=stdD_X;
            model_GPC.predict_function=pred_fcn;
            model_GPC.hyp=hyp;
            model_GPC.post=post;
            model_GPC.simplify_flag=simplify_flag;

            function [fval,gradient]=objectNLLGPC(x,inf,mean,cov,lik,X,Y)
                hyp_iter.mean=x(1);
                hyp_iter.cov=x(2:end);
                hyp_iter.lik=[];

                if nargout < 2
                    [~,nlZ]=feval(inf{:},hyp_iter,mean,cov,lik,X,Y);
                    fval=nlZ;
                elseif nargout < 3
                    [~,nlZ,dnlZ]=feval(inf{:},hyp_iter,mean,cov,lik,X,Y);
                    fval=nlZ;
                    gradient=[dnlZ.mean,dnlZ.cov];
                end
            end

            function [class,possibility,miu_pre,var_pre]=classifyGaussPredictor...
                    (X_pred,hyp,mean,cov,lik,post,X,aver_X,stdD_X)
                % predict function
                %
                X_pred_nomlz=(X_pred-aver_X)./stdD_X;
                pred_num=size(X_pred_nomlz,1);
                ys=ones(pred_num,1);

                alpha=post.alpha; L=post.L; sW=post.sW;
                nz=true(size(alpha,1),1);               % non-sparse representation
                %verify whether L contains valid Cholesky decomposition or something different
                Lchol=isnumeric(L) && all(all(tril(L,-1)==0)&diag(L)'>0&isreal(diag(L))');
                ns=size(X_pred_nomlz,1);                                       % number of data points
                nperbatch=1000;                       % number of data points per mini batch
                nact=0;                       % number of already processed test data points
                ymu=zeros(ns,1); ys2=ymu; miu_pre=ymu; var_pre=ymu; possibility=ymu;   % allocate mem
                while nact<ns               % process minibatches of test cases to save memory
                    id=(nact+1):min(nact+nperbatch,ns);               % data points to process
                    kss=feval(cov{:},hyp.cov,X_pred_nomlz(id,:),'diag');              % self-variance
                    Ks=feval(cov{:},hyp.cov,X(nz,:),X_pred_nomlz(id,:));        % avoid computation
                    ms=feval(mean{:},hyp.mean,X_pred_nomlz(id,:));
                    N=size(alpha,2);  % number of alphas (usually 1; more in case of sampling)
                    Fmu=repmat(ms,1,N) + Ks'*full(alpha(nz,:));        % conditional mean fs|f
                    miu_pre(id)=sum(Fmu,2)/N;                                   % predictive means
                    if Lchol    % L contains chol decomp => use Cholesky parameters (alpha,sW,L)
                        V =L'\(repmat(sW,1,length(id)).*Ks);
                        var_pre(id)=kss - sum(V.*V,1)';                       % predictive variances
                    else                % L is not triangular => use alternative parametrisation
                        if isnumeric(L),LKs=L*Ks; else LKs=L(Ks); end    % matrix or callback
                        var_pre(id)=kss + sum(Ks.*LKs,1)';                    % predictive variances
                    end
                    var_pre(id)=max(var_pre(id),0);   % remove numerical noise i.e. negative variances
                    Fs2=repmat(var_pre(id),1,N);     % we have multiple values in case of sampling
                    if nargin<9
                        [Lp,Ymu,Ys2]=feval(lik{:},hyp.lik,[],Fmu(:),Fs2(:));
                    else
                        Ys=repmat(ys(id),1,N);
                        [Lp,Ymu,Ys2]=feval(lik{:},hyp.lik,Ys(:),Fmu(:),Fs2(:));
                    end
                    possibility(id) =sum(reshape(Lp,[],N),2)/N;    % log probability; sample averaging
                    ymu(id)=sum(reshape(Ymu,[],N),2)/N;          % predictive mean ys|y and ..
                    ys2(id)=sum(reshape(Ys2,[],N),2)/N;                          % .. variance
                    nact=id(end);          % set counter to index of last processed data point
                end

                possibility=exp(possibility);
                class=ones(pred_num,1);
                class(possibility < 0.5)=-1;
            end

            function [K,dK_dcov]=calCov(cov,X,Z)
                % obtain covariance of x
                % cov: eta,len(equal to 1/len.^2)
                %
                % k=eta*exp(-sum(x_dis*len)/vari_num);
                %
                [x_num,vari_num]=size(X);

                len=exp(cov(1:end-1));
                eta=exp(cov(end));

                if nargin > 2 && nargout < 2 && ~isempty(Z)
                    % predict
                    if strcmp(Z,'diag')
                        K=eta;
                    else
                        [z_number,vari_num]=size(Z);
                        % initializate square of X inner distance/ vari_num
                        K=zeros(x_num,z_number);

                        if simplify_flag
                            for len_index=1:vari_num
                                K=K+(X(:,len_index)-Z(:,len_index)').^2*len/vari_num;
                            end
                        else
                            for len_index=1:vari_num
                                K=K+(X(:,len_index)-Z(:,len_index)').^2*len(len_index)/vari_num;
                            end
                        end

                        K=eta*exp(-K);
                    end
                else
                    % train
                    % initializate square of X inner distance sq
                    sq_dis_v=zeros(x_num,x_num,vari_num);
                    for len_index=1:vari_num
                        sq_dis_v(:,:,len_index)=(X(:,len_index)-X(:,len_index)').^2/vari_num;
                    end

                    % exp of x__x with theta
                    exp_dis=zeros(x_num);

                    if simplify_flag
                        for len_index=1:vari_num
                            exp_dis=exp_dis+sq_dis_v(:,:,len_index)*len;
                        end
                    else
                        for len_index=1:vari_num
                            exp_dis=exp_dis+sq_dis_v(:,:,len_index)*len(len_index);
                        end
                    end

                    exp_dis=exp(-exp_dis);
                    K=exp_dis*eta;

                    if nargout >= 2
                        if simplify_flag
                            dK_dcov=cell(1,2);
                        else
                            dK_dcov=cell(1,vari_num+1);
                        end

                        % len
                        if simplify_flag
                            dK_dlen=zeros(x_num,x_num);
                            for len_index=1:vari_num
                                dK_dlen=dK_dlen + sq_dis_v(:,:,len_index);
                            end
                            dK_dlen=-dK_dlen.*K*len;
                            dK_dcov{1}=dK_dlen;
                        else
                            for len_index=1:vari_num
                                dK_dcov{len_index}=-K.*sq_dis_v(:,:,len_index)*len(len_index);
                            end
                        end

                        % eta
                        dK_dcov{end}=K;
                    end
                end

            end

        end

        function [post nlZ dnlZ]=infEP(hyp,mean,cov,lik,x,y)
            % Expectation Propagation approximation to the posterior Gaussian Process.
            % The function takes a specified covariance function (see covFunctions.m) and
            % likelihood function (see likFunctions.m),and is designed to be used with
            % gp.m. See also infMethods.m. In the EP algorithm,the sites are
            % updated in random order,for better performance when cases are ordered
            % according to the targets.
            %
            % Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2013-09-13.
            %
            % See also INFMETHODS.M.
            %
            persistent last_ttau last_tnu              % keep tilde parameters between calls
            tol=1e-4; max_sweep=10; min_sweep=2;     % tolerance to stop EP iterations

            inf='infEP';
            n=size(x,1);
            if isnumeric(cov),K=cov;                    % use provided covariance matrix
            else K=feval(cov{:},hyp.cov,x); end       % evaluate the covariance matrix
            if isnumeric(mean),m=mean;                         % use provided mean vector
            else m=feval(mean{:},hyp.mean,x); end             % evaluate the mean vector

            % A note on naming: variables are given short but descriptive names in
            % accordance with Rasmussen & Williams "GPs for Machine Learning" (2006): mu
            % and s2 are mean and variance,nu and tau are natural parameters. A leading t
            % means tilde,a subscript _ni means "not i" (for cavity parameters),or _n
            % for a vector of cavity parameters. N(f|mu,Sigma) is the posterior.

            % marginal likelihood for ttau=tnu=zeros(n,1); equals n*log(2) for likCum*
            nlZ0=-sum(feval(lik{:},hyp.lik,y,m,diag(K),inf));
            if any(size(last_ttau) ~= [n 1])      % find starting point for tilde parameters
                ttau=zeros(n,1); tnu =zeros(n,1);        % init to zero if no better guess
                Sigma=K;                     % initialize Sigma and mu,the parameters of ..
                mu=m; nlZ=nlZ0;                  % .. the Gaussian posterior approximation
            else
                ttau=last_ttau; tnu =last_tnu;   % try the tilde values from previous call
                [Sigma,mu,L,alpha,nlZ]=epComputeParams(K,y,ttau,tnu,lik,hyp,m,inf);
                if nlZ > nlZ0                                           % if zero is better ..
                    ttau=zeros(n,1); tnu =zeros(n,1);       % .. then init with zero instead
                    Sigma=K;                   % initialize Sigma and mu,the parameters of ..
                    mu=m; nlZ=nlZ0;                % .. the Gaussian posterior approximation
                end
            end

            nlZ_old=Inf; sweep=0;               % converged,max. sweeps or min. sweeps?
            while (abs(nlZ-nlZ_old) > tol && sweep < max_sweep) || sweep<min_sweep
                nlZ_old=nlZ; sweep=sweep+1;
                for i=randperm(n)       % iterate EP updates (in random order) over examples
                    tau_ni=1/Sigma(i,i)-ttau(i);      %  first find the cavity distribution ..
                    nu_ni=mu(i)/Sigma(i,i)-tnu(i);                % .. params tau_ni and nu_ni

                    % compute the desired derivatives of the indivdual log partition function
                    [lZ,dlZ,d2lZ]=feval(lik{:},hyp.lik,y(i),nu_ni/tau_ni,1/tau_ni,inf);
                    ttau_old=ttau(i); tnu_old=tnu(i);  % find the new tilde params,keep old
                    ttau(i)=                    -d2lZ  /(1+d2lZ/tau_ni);
                    ttau(i)=max(ttau(i),0); % enforce positivity i.e. lower bound ttau by zero
                    tnu(i) =( dlZ - nu_ni/tau_ni*d2lZ )/(1+d2lZ/tau_ni);

                    dtt=ttau(i)-ttau_old; dtn=tnu(i)-tnu_old;      % rank-1 update Sigma ..
                    si=Sigma(:,i); ci=dtt/(1+dtt*si(i));
                    Sigma=Sigma - ci*si*si';                         % takes 70% of total time
                    mu=mu - (ci*(mu(i)+si(i)*dtn)-dtn)*si;               % .. and recompute mu
                end
                % recompute since repeated rank-one updates can destroy numerical precision
                [Sigma,mu,L,alpha,nlZ]=epComputeParams(K,y,ttau,tnu,lik,hyp,m,inf);
            end

            %             if sweep == max_sweep && abs(nlZ-nlZ_old) > tol
            %                 error('maximum number of sweeps exceeded in function infEP')
            %             end

            last_ttau=ttau; last_tnu=tnu;                       % remember for next call
            post.alpha=alpha; post.sW=sqrt(ttau); post.L=L;  % return posterior params

            if nargout>2                                           % do we want derivatives?
                dnlZ=hyp;                                   % allocate space for derivatives
                tau_n=1./diag(Sigma)-ttau;             % compute the log marginal likelihood
                nu_n =mu./diag(Sigma)-tnu;                    % vectors of cavity parameters
                sW=sqrt(ttau);
                F=alpha*alpha'-repmat(sW,1,n).*(L\(L'\diag(sW)));   % covariance hypers
                [K,dK]=feval(cov{:},hyp.cov,x,[]);
                for i=1:length(hyp.cov)
                    dnlZ.cov(i)=-sum(sum(F.*dK{i}))/2;
                end
                for i=1:numel(hyp.lik)                                   % likelihood hypers
                    dlik=feval(lik{:},hyp.lik,y,nu_n./tau_n,1./tau_n,inf,i);
                    dnlZ.lik(i)=-sum(dlik);
                end
                [junk,dlZ]=feval(lik{:},hyp.lik,y,nu_n./tau_n,1./tau_n,inf);% mean hyps
                for i=1:numel(hyp.mean)
                    dm=feval(mean{:},hyp.mean,x,i);
                    dnlZ.mean(i)=-dlZ'*dm;
                end
            end

            function [Sigma,mu,L,alpha,nlZ]=epComputeParams(K,y,ttau,tnu,lik,hyp,m,inf)
                % function to compute the parameters of the Gaussian approximation,Sigma and
                % mu,and the negative log marginal likelihood,nlZ,from the current site
                % parameters,ttau and tnu. Also returns L (useful for predictions).
                %
                n=length(y);                                      % number of training cases
                sW=sqrt(ttau);                                        % compute Sigma and mu
                L=chol(eye(n)+sW*sW'.*K);                            % L'*L=B=eye(n)+sW*K*sW
                V=L'\(repmat(sW,1,n).*K);
                Sigma=K - V'*V;
                alpha=tnu-sW.*(L\(L'\(sW.*(K*tnu+m))));
                mu=K*alpha+m; v=diag(Sigma);

                tau_n=1./diag(Sigma)-ttau;             % compute the log marginal likelihood
                nu_n =mu./diag(Sigma)-tnu;                    % vectors of cavity parameters
                lZ=feval(lik{:},hyp.lik,y,nu_n./tau_n,1./tau_n,inf);
                p=tnu-m.*ttau; q=nu_n-m.*tau_n;                        % auxiliary vectors
                nlZ=sum(log(diag(L))) - sum(lZ) - p'*Sigma*p/2 + (v'*p.^2)/2 ...
                    - q'*((ttau./tau_n.*q-2*p).*v)/2 - sum(log(1+ttau./tau_n))/2;
            end

        end

        function A=meanConst(hyp,x,i)

            % Constant mean function. The mean function is parameterized as:
            %
            % m(x)=c
            %
            % The hyperparameter is:
            %
            % hyp=[ c ]
            %
            % Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch,2010-08-04.
            %
            % See also MEANFUNCTIONS.M.

            if nargin<2,A='1'; return; end             % report number of hyperparameters
            if numel(hyp)~=1,error('Exactly one hyperparameter needed.'),end
            c=hyp;
            if nargin==2
                A=c*ones(size(x,1),1);                                       % evaluate mean
            else
                if i==1
                    A=ones(size(x,1),1);                                          % derivative
                else
                    A=zeros(size(x,1),1);
                end
            end
        end

        function [varargout]=likErf(hyp,y,mu,s2,inf,i)
            % likErf - Error function or cumulative Gaussian likelihood function for binary
            % classification or probit regression. The expression for the likelihood is
            %   likErf(t)=(1+erf(t/sqrt(2)))/2=normcdf(t).
            %
            % Several modes are provided,for computing likelihoods,derivatives and moments
            % respectively,see likFunctions.m for the details. In general,care is taken
            % to avoid numerical issues when the arguments are extreme.
            %
            % Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch,2014-03-19.
            %
            % See also LIKFUNCTIONS.M.
            %
            if nargin<3,varargout={'0'}; return; end   % report number of hyperparameters
            if nargin>1,y=sign(y); y(y==0)=1; else y=1; end % allow only +/- 1 values
            if numel(y)==0,y=1; end

            if nargin<5                              % prediction mode if inf is not present
                y=y.*ones(size(mu));                                       % make y a vector
                s2zero=1; if nargin>3&&numel(s2)>0&&norm(s2)>eps,s2zero=0; end  % s2==0 ?
                if s2zero                                         % log probability evaluation
                    lp=logphi(y.*mu);
                else                                                              % prediction
                    lp=OptSACORS.likErf(hyp,y,mu,s2,'infEP');
                end
                p=exp(lp); ymu={}; ys2={};
                if nargout>1
                    ymu=2*p-1;                                                % first y moment
                    if nargout>2
                        ys2=4*p.*(1-p);                                        % second y moment
                    end
                end
                varargout={lp,ymu,ys2};
            else                                                            % inference mode
                switch inf
                    case 'infLaplace'
                        if nargin<6                                             % no derivative mode
                            f=mu; yf=y.*f;                            % product latents and labels
                            varargout=cell(nargout,1); [varargout{:}]=logphi(yf);   % query logphi
                            if nargout>1
                                varargout{2}=y.*varargout{2};
                                if nargout>3,varargout{4}=y.*varargout{4}; end
                            end
                        else                                                       % derivative mode
                            varargout={[],[],[]};                         % derivative w.r.t. hypers
                        end

                    case 'infEP'
                        if nargin<6                                             % no derivative mode
                            z=mu./sqrt(1+s2); dlZ={}; d2lZ={};
                            if numel(y)>0,z=z.*y; end
                            if nargout <= 1,lZ=logphi(z);                         % log part function
                            else          [lZ,n_p]=logphi(z); end
                            if nargout > 1
                                if numel(y)==0,y=1; end
                                dlZ=y.*n_p./sqrt(1+s2);                      % 1st derivative wrt mean
                                if nargout>2,d2lZ=-n_p.*(z+n_p)./(1+s2); end         % 2nd derivative
                            end
                            varargout={lZ,dlZ,d2lZ};
                        else                                                       % derivative mode
                            varargout={[]};                                     % deriv. wrt hyp.lik
                        end
                end
            end
            function [lp,dlp,d2lp,d3lp]=logphi(z)
                % Safe computation of logphi(z)=log(normcdf(z)) and its derivatives
                %                    dlogphi(z)=normpdf(x)/normcdf(x).
                % The function is based on index 5725 in Hart et al. and gsl_sf_log_erfc_e.
                %
                % Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch,2013-11-13.
                %
                z=real(z);                                 % support for real arguments only
                lp=zeros(size(z));                                         % allocate memory
                id1=z.*z<0.0492;                                 % first case: close to zero
                lp0=-z(id1)/sqrt(2*pi);
                c=[ 0.00048204; -0.00142906; 0.0013200243174; 0.0009461589032;
                    -0.0045563339802; 0.00556964649138; 0.00125993961762116;
                    -0.01621575378835404; 0.02629651521057465; -0.001829764677455021;
                    2*(1-pi/3); (4-pi)/3; 1; 1];
                f=0; for i=1:14,f=lp0.*(c(i)+f); end,lp(id1)=-2*f-log(2);
                id2=z<-11.3137;                                    % second case: very small
                r=[ 1.2753666447299659525; 5.019049726784267463450;
                    6.1602098531096305441; 7.409740605964741794425;
                    2.9788656263939928886 ];
                q=[ 2.260528520767326969592;  9.3960340162350541504;
                    12.048951927855129036034; 17.081440747466004316;
                    9.608965327192787870698;  3.3690752069827527677 ];
                num=0.5641895835477550741; for i=1:5,num=-z(id2).*num/sqrt(2) + r(i); end
                den=1.0;                   for i=1:6,den=-z(id2).*den/sqrt(2) + q(i); end
                e=num./den; lp(id2)=log(e/2) - z(id2).^2/2;
                id3=~id2 & ~id1; lp(id3)=log(erfc(-z(id3)/sqrt(2))/2);  % third case: rest
                if nargout>1                                        % compute first derivative
                    dlp=zeros(size(z));                                      % allocate memory
                    dlp( id2)=abs(den./num) * sqrt(2/pi); % strictly positive first derivative
                    dlp(~id2)=exp(-z(~id2).*z(~id2)/2-lp(~id2))/sqrt(2*pi); % safe computation
                    if nargout>2                                     % compute second derivative
                        d2lp=-dlp.*abs(z+dlp);             % strictly negative second derivative
                        if nargout>3                                    % compute third derivative
                            d3lp=-d2lp.*abs(z+2*dlp)-dlp;     % strictly positive third derivative
                        end
                    end
                end
            end

        end

    end

    % surrogate function
    methods(Static)
        function [pred_fcn,radialbasis_model]=interpRadialBasisPreModel...
                (X,Y,basis_function)
            % radial basis function interp pre objcon_fcn function
            % input initial data X,Y,which are real data
            % X,Y are x_number x variable_number matrix
            % aver_X,stdD_X is 1 x x_number matrix
            % output is a radial basis objcon_fcn,include X,Y,base_function
            % and predict_function
            %
            % Copyright 2023 Adel
            %
            if nargin < 3
                basis_function=[];
            end

            [x_number,variable_number]=size(X);

            % normalize data
            aver_X=mean(X);
            stdD_X=std(X);
            aver_Y=mean(Y);
            stdD_Y=std(Y);
            idx__=find(stdD_X == 0);
            if ~isempty(idx__),stdD_X(idx__)=1;end
            idx__=find(stdD_Y == 0);
            if ~isempty(idx__),stdD_Y(idx__)=1;end
            X_nomlz=(X-aver_X)./stdD_X;
            Y_nomlz=(Y-aver_Y)./stdD_Y;

            if isempty(basis_function)
                %     c=(prod(max(X_nomlz)-min(X_nomlz))/x_number)^(1/variable_number);
                %     basis_function=@(r) exp(-(r.^2)/c);
                basis_function=@(r) r.^3;
            end

            % initialization distance of all X
            X_dis=zeros(x_number,x_number);
            for variable_idx=1:variable_number
                X_dis=X_dis+(X_nomlz(:,variable_idx)-X_nomlz(:,variable_idx)').^2;
            end
            X_dis=sqrt(X_dis);

            [beta,rdibas_matrix,inv_rdibas_matrix]=interpRadialBasis...
                (X_dis,Y_nomlz,basis_function,x_number);

            % initialization predict function
            pred_fcn=@(X_predict) interpRadialBasisPredictor...
                (X_predict,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
                x_number,variable_number,beta,basis_function);

            radialbasis_model.X=X;
            radialbasis_model.Y=Y;
            radialbasis_model.radialbasis_matrix=rdibas_matrix;
            radialbasis_model.inv_radialbasis_matrix=inv_rdibas_matrix;
            radialbasis_model.beta=beta;

            radialbasis_model.aver_X=aver_X;
            radialbasis_model.stdD_X=stdD_X;
            radialbasis_model.aver_Y=aver_Y;
            radialbasis_model.stdD_Y=stdD_Y;
            radialbasis_model.basis_function=basis_function;

            radialbasis_model.predict_function=pred_fcn;

            % abbreviation:
            % num: number,pred: predict,vari: variable
            function [beta,rdibas_matrix,inv_rdibas_matrix]=interpRadialBasis...
                    (X_dis,Y,basis_function,x_number)
                % interp polynomial responed surface core function
                % calculation beta
                %
                % Copyright 2022 Adel
                %
                rdibas_matrix=basis_function(X_dis);

                % stabilize matrix
                rdibas_matrix=rdibas_matrix+eye(x_number)*1e-9;

                % get inverse matrix
                inv_rdibas_matrix=rdibas_matrix\eye(x_number);

                % solve beta
                beta=inv_rdibas_matrix*Y;
            end

            function [Y_pred]=interpRadialBasisPredictor...
                    (X_pred,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
                    x_num,vari_num,beta,basis_function)
                % radial basis function interpolation predict function
                %
                [x_pred_num,~]=size(X_pred);

                % normalize data
                X_pred_nomlz=(X_pred-aver_X)./stdD_X;

                % calculate distance
                X_dis_pred=zeros(x_pred_num,x_num);
                for vari_idx=1:vari_num
                    X_dis_pred=X_dis_pred+...
                        (X_pred_nomlz(:,vari_idx)-X_nomlz(:,vari_idx)').^2;
                end
                X_dis_pred=sqrt(X_dis_pred);

                % predict variance
                Y_pred=basis_function(X_dis_pred)*beta;

                % normalize data
                Y_pred=Y_pred*stdD_Y+aver_Y;
            end

        end

    end

    % LHD
    methods(Static)
        function [X_sample,dist_min_nomlz,X_total]=lhdESLHS...
                (sample_number,variable_number,...
                low_bou,up_bou,X_exist,cheapcon_function)
            % generate latin hypercube design
            % ESLHS method is used(sample and iteration)
            % election combination mode of point and find best combination
            %
            % input:
            % sample_number(new point to sample), variable_number, ...
            % low_bou, up_bou, X_exist(exist point), cheapcon_function
            %
            % output:
            % X_sample, dist_min_nomlz(min distance of normalize data), ...
            % X_total(include all data in area)
            %
            % reference: [1] LONG T, LI X, SHI R, et al., Gradient-Free
            % Trust-Region-Based Adaptive Response Surface Method for Expensive
            % Aircraft Optimization[J]. AIAA Journal, 2018, 56(2): 862-73.
            %
            % Copyright 2023 03 Adel
            %
            if nargin < 6
                cheapcon_function=[];
                if nargin < 5
                    X_exist=[];
                    if nargin < 4
                        up_bou=ones(1,variable_number);
                        if nargin < 3
                            low_bou=zeros(1,variable_number);
                            if nargin < 2
                                error('getLatinHypercube: lack variable_number');
                            end
                        end
                    end
                end
            end

            iteration_max=min(100*sample_number,100);

            % check x_exist_list if meet boundary
            if ~isempty(X_exist)
                if size(X_exist,2) ~= variable_number
                    error('getLatinHypercube: x_exist_list variable_number error');
                end
                X_exist_nomlz=(X_exist-low_bou)./(up_bou-low_bou);
            else
                X_exist_nomlz=[];
            end

            exist_number=size(X_exist,1);
            total_number=sample_number+exist_number;
            if sample_number <= 0
                X_total=X_exist;
                X_sample=[];
                dist_min_nomlz=getMinDistance(X_exist_nomlz);
                return;
            end

            % get quasi-feasible point
            x_initial_number=10*sample_number;
            x_quasi_number=3*sample_number;
            if ~isempty(cheapcon_function)
                X_supply_quasi_nomlz=[];

                % check if have enough X_supply_nomlz
                iteration=0;
                while size(X_supply_quasi_nomlz,1) < x_quasi_number && iteration < 100
                    X_supply_initial_nomlz=lhsdesign(x_initial_number,variable_number);

                    qusai_index=[];
                    for x_index=1:size(X_supply_initial_nomlz,1)
                        if cheapcon_function(X_supply_initial_nomlz(x_index,:).*(up_bou-low_bou)+low_bou) <= 0
                            qusai_index=[qusai_index,x_index];
                        end
                    end
                    X_supply_quasi_nomlz=[X_supply_quasi_nomlz;X_supply_initial_nomlz(qusai_index,:)];

                    iteration=iteration+1;
                end

                if iteration == 100 && isempty(X_supply_quasi_nomlz)
                    error('getLatinHypercube: feasible quasi point cannot be found');
                end
            else
                X_supply_quasi_nomlz=lhsdesign(x_quasi_number,variable_number);
            end

            % iterate and get final x_supply_list
            iteration=0;
            x_supply_quasi_number=size(X_supply_quasi_nomlz,1);
            dist_min_nomlz=0;
            X_sample_nomlz=[];

            % dist_min_nomlz_result=zeros(1,iteration);
            while iteration <= iteration_max
                % random select x_new_number X to X_trial_nomlz
                x_select_index=randperm(x_supply_quasi_number,sample_number);

                % get distance min itertion X_
                distance_min_iteration=getMinDistanceIter...
                    (X_supply_quasi_nomlz(x_select_index,:),X_exist_nomlz);

                % if distance_min_iteration is large than last time
                if distance_min_iteration > dist_min_nomlz
                    dist_min_nomlz=distance_min_iteration;
                    X_sample_nomlz=X_supply_quasi_nomlz(x_select_index,:);
                end

                iteration=iteration+1;
                %     dist_min_nomlz_result(iteration)=dist_min_nomlz;
            end
            dist_min_nomlz=getMinDistance([X_sample_nomlz;X_exist_nomlz]);
            X_sample=X_sample_nomlz.*(up_bou-low_bou)+low_bou;
            X_total=[X_exist;X_sample];

            function distance_min__=getMinDistance(x_list__)
                % get distance min from x_list
                %
                if isempty(x_list__)
                    distance_min__=[];
                    return;
                end

                % sort x_supply_list_initial to decrese distance calculate times
                x_list__=sortrows(x_list__,1);
                [sample_number__,variable_number__]=size(x_list__);
                distance_min__=variable_number__;
                for x_index__=1:sample_number__
                    x_curr__=x_list__(x_index__,:);
                    x_next_index__=x_index__ + 1;
                    % only search in min_distance(x_list had been sort)
                    search_range__=variable_number__;
                    while x_next_index__ <= sample_number__ &&...
                            (x_list__(x_next_index__,1)-x_list__(x_index__,1))^2 ...
                            < search_range__
                        x_next__=x_list__(x_next_index__,:);
                        distance_temp__=sum((x_next__-x_curr__).^2);
                        if distance_temp__ < distance_min__
                            distance_min__=distance_temp__;
                        end
                        if distance_temp__ < search_range__
                            search_range__=distance_temp__;
                        end
                        x_next_index__=x_next_index__+1;
                    end
                end
                distance_min__=sqrt(distance_min__);
            end
            function distance_min__=getMinDistanceIter...
                    (x_list__,x_exist_list__)
                % get distance min from x_list
                % this version do not consider distance between x exist
                %

                % sort x_supply_list_initial to decrese distance calculate times
                x_list__=sortrows(x_list__,1);
                [sample_number__,variable_number__]=size(x_list__);
                distance_min__=variable_number__;
                for x_index__=1:sample_number__
                    x_curr__=x_list__(x_index__,:);
                    x_next_index__=x_index__ + 1;
                    % only search in min_distance(x_list had been sort)
                    search_range__=variable_number__;
                    while x_next_index__ <= sample_number__ &&...
                            (x_list__(x_next_index__,1)-x_list__(x_index__,1))^2 ...
                            < search_range__
                        x_next__=x_list__(x_next_index__,:);
                        distance_temp__=sum((x_next__-x_curr__).^2);
                        if distance_temp__ < distance_min__
                            distance_min__=distance_temp__;
                        end
                        if distance_temp__ < search_range__
                            search_range__=distance_temp__;
                        end
                        x_next_index__=x_next_index__+1;
                    end
                    for x_exist_index=1:size(x_exist_list__,1)
                        x_next__=x_exist_list__(x_exist_index,:);
                        distance_temp__=sum((x_next__-x_curr__).^2);
                        if distance_temp__ < distance_min__
                            distance_min__=distance_temp__;
                        end
                    end
                end
                distance_min__=sqrt(distance_min__);
            end
        end

    end
end

%% data library function
function [data_lib,x_new,obj_new,con_new,coneq_new,vio_new,ks_new,repeat_idx,NFE]=dataUpdata...
    (data_lib,x_origin_new,protect_range,write_file_flag,str_data_file)
% updata data lib
% updata format:
% variable_number,obj_number,con_number,coneq_number
% x,obj,con,coneq
%
if nargin < 3
    protect_range=0;
end

[x_new_num,~]=size(x_origin_new);
x_new=[];
obj_new=[];
con_new=[];
coneq_new=[];
vio_new=[];
ks_new=[];
repeat_idx=[];
NFE=0;

if write_file_flag
    file_data=fopen(str_data_file,'a');
end

% updata format:
% variable_number,obj_number,con_number,coneq_number
% x,obj,con,coneq
for x_idx=1:x_new_num
    x=x_origin_new(x_idx,:);

    if protect_range ~= 0
        % updata data with same_point_avoid protect
        % check x_potl if exist in data lib
        % if exist, jump updata
        distance=sum((abs(x-data_lib.X)./(data_lib.up_bou-data_lib.low_bou)),2);
        [distance_min,min_idx]=min(distance);
        if distance_min < data_lib.vari_num*protect_range
            % distance to exist point of point to add is small than protect_range
            repeat_idx=[repeat_idx;min_idx];
            continue;
        end
    end

    [obj,con,coneq]=data_lib.objcon_fcn(x); % eval value
    NFE=NFE+1;

    con=con(:)';
    coneq=coneq(:)';
    % calculate vio
    if isempty(con) && isempty(coneq)
        vio=[];
        ks=[];
    else
        vio=calViolation(con,coneq,data_lib.con_torl);
        ks=max([con,coneq]);
    end


    x_new=[x_new;x];
    obj_new=[obj_new;obj];
    if ~isempty(con)
        con_new=[con_new;con];
    end
    if ~isempty(coneq)
        coneq_new=[coneq_new;coneq];
    end
    if ~isempty(vio)
        vio_new=[vio_new;vio];
    end
    if ~isempty(ks)
        ks_new=[ks_new;ks];
    end

    if write_file_flag
        % write data to txt_result
        fprintf(file_data,'%d ',repmat('%.8e ',1,data_lib.vari_num));
        fprintf(file_data,'%d ',length(obj));
        fprintf(file_data,'%d ',length(con));
        fprintf(file_data,'%d ',length(coneq));

        fprintf(file_data,data_lib.x_format,x);
        fprintf(file_data,repmat('%.8e ',1,length(obj)),obj);
        fprintf(file_data,repmat('%.8e ',1,length(con)),con);
        fprintf(file_data,repmat('%.8e ',1,length(coneq)),coneq);
        fprintf(file_data,'\n');
    end

    data_lib=dataJoin(data_lib,x,obj,con,coneq,vio,ks);

    % record best
    if isempty(data_lib.result_best_idx)
        data_lib.result_best_idx=1;
    else
        if isempty(vio) || vio == 0
            if obj <= data_lib.Obj(data_lib.result_best_idx(end))
                data_lib.result_best_idx=[data_lib.result_best_idx;size(data_lib.X,1)];
            else
                data_lib.result_best_idx=[data_lib.result_best_idx;data_lib.result_best_idx(end)];
            end
        else
            if vio <= data_lib.Vio(data_lib.result_best_idx(end))
                data_lib.result_best_idx=[data_lib.result_best_idx;size(data_lib.X,1)];
            else
                data_lib.result_best_idx=[data_lib.result_best_idx;data_lib.result_best_idx(end)];
            end
        end
    end
end

if write_file_flag
    fclose(file_data);
    clear('file_data');
end
end

function [x_list,obj_list,con_list,coneq_list,vio_list,ks_list]=dataLoad...
    (data_lib,low_bou,up_bou)
% updata data to exist data lib
%
if nargin < 3
    up_bou=realmax;
    if nargin < 2
        low_bou=-realmax;
    end
end

idx=[];
for x_idx=1:size(data_lib.X,1)
    x=data_lib.X(x_idx,:);
    if all(x > low_bou) && all(x < up_bou)
        idx=[idx;x_idx];
    end
end

x_list=data_lib.X(idx,:);
obj_list=data_lib.Obj(idx,:);
if ~isempty(data_lib.Con)
    con_list=data_lib.Con(idx,:);
else
    con_list=[];
end
if ~isempty(data_lib.Coneq)
    coneq_list=data_lib.Coneq(idx,:);
else
    coneq_list=[];
end
if ~isempty(data_lib.Vio)
    vio_list=data_lib.Vio(idx,:);
else
    vio_list=[];
end
if ~isempty(data_lib.Ks)
    ks_list=data_lib.Ks(idx);
else
    ks_list=[];
end
end

function data_lib=dataJoin(data_lib,x,obj,con,coneq,vio,ks)
% updata data to exist data lib
%
data_lib.X=[data_lib.X;x];
data_lib.Obj=[data_lib.Obj;obj];
if ~isempty(data_lib.Con) || ~isempty(con)
    data_lib.Con=[data_lib.Con;con];
end
if ~isempty(data_lib.Coneq) || ~isempty(coneq)
    data_lib.Coneq=[data_lib.Coneq;coneq];
end
if ~isempty(data_lib.Vio) || ~isempty(vio)
    data_lib.Vio=[data_lib.Vio;vio];
end
if ~isempty(data_lib.Ks) || ~isempty(ks)
    data_lib.Ks=[data_lib.Ks;ks];
end
end

function vio_list=calViolation(con_list,coneq_list,con_torl)
% calculate violation of data
%
if isempty(con_list) && isempty(coneq_list)
    vio_list=[];
else
    vio_list=zeros(max(size(con_list,1),size(coneq_list,1)),1);
    if ~isempty(con_list)
        vio_list=vio_list+sum(max(con_list-con_torl,0),2);
    end
    if ~isempty(coneq_list)
        vio_list=vio_list+sum((abs(coneq_list)-con_torl),2);
    end
end
end

function ks=calKs(objcon_fcn,x)
[obj,con,coneq]=objcon_fcn(x);
ks=max([con,coneq]);
end

function [con,coneq]=getCon(objcon_fcn,x)
[obj,con,coneq]=objcon_fcn(x);
end
