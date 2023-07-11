classdef OptSADERBF < handle
    % RBF-CDE optimization algorithm
    % use RBF model instead of origin kriging model
    %
    % referance: [1] 叶年辉,龙腾,武宇飞,et al.
    % 基于Kriging代理模型的约束差分进化算法 [J]. 航空学报,2021,42(6): 13.
    %
    % Copyright 2023 Adel
    %
    % abbreviation:
    % obj: objective, con: constraint, iter: iteration, tor: torlance
    % fcn: function
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

        nomlz_value=10; % max obj when normalize obj,con,coneq
        protect_range=1e-16; % surrogate add point protect range

        % differ evoluation parameter
        scaling_factor = 0.8; % F
        cross_rate = 0.8;

        str_data_file='result_total.txt'

    end

    properties
        % problem parameter
        expensive_con_flag;

        obj_fcn_surr;
        con_fcn_surr;
        ks_fcn_surr;

        obj_surr;
        con_surr_list;
        coneq_surr_list;

    end

    methods % main
        function self=OptSADERBF(NFE_max,iter_max,obj_torl,con_torl,data_lib)
            % initialize optimization
            %
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

        end

        function [x_best,obj_best,NFE,output]=optimize(self,varargin)
            % main optimize function
            %
            if length(varargin) == 1
                % input is struct or object
                problem=varargin{1};
                if isstruct(problem)
                    fields=fieldnames(problem);
                    if ~contains(fields,'objcon_fcn'), error('optSACORS.optimize: input problem lack objcon_fcn'); end
                    objcon_fcn=problem.objcon_fcn;
                    if ~contains(fields,'vari_num'), error('optSACORS.optimize: input problem lack vari_num'); end
                    if ~contains(fields,'low_bou'), error('optSACORS.optimize: input problem lack low_bou'); end
                    if ~contains(fields,'up_bou'), error('optSACORS.optimize: input problem lack up_bou'); end
                    if ~contains(fields,'con_fcn_cheap'), problem.con_fcn_cheap=[]; end
                    con_fcn_cheap=problem.con_fcn_cheap;
                else
                    m=methods(problem);
                    if ~contains(m,'objcon_fcn'), error('optSACORS.optimize: input problem lack objcon_fcn'); end
                    objcon_fcn=@(x) problem.objcon_fcn(x);
                    p=properties(problem);
                    if ~contains(p,'vari_num'), error('optSACORS.optimize: input problem lack vari_num'); end
                    if ~contains(p,'low_bou'), error('optSACORS.optimize: input problem lack low_bou'); end
                    if ~contains(p,'up_bou'), error('optSACORS.optimize: input problem lack up_bou'); end
                    if ~contains(m,'con_fcn_cheap')
                        con_fcn_cheap=[];
                    else
                        con_fcn_cheap=@(x) problem.con_fcn_cheap(x);
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

            % hyper parameter
            pop_num = min(100,10*vari_num);
            RBF_number = max(100,(vari_num+1)*(vari_num+2)/2);

            % NFE and iteration setting
            if isempty(self.NFE_max)
                self.NFE_max = 10+10*vari_num;
            end

            if isempty(self.iter_max)
                self.iter_max = 20+20*vari_num;
            end

            done = 0;NFE = 0;iter = 0;

            % step 2
            % generate initial sample X
            if isempty(self.data_lib)
                self.data_lib=struct('objcon_fcn',objcon_fcn,'vari_num',vari_num,'low_bou',low_bou,'up_bou',up_bou,'con_torl',self.con_torl,...
                    'result_best_idx',[],'X',[],'Obj',[],'Con',[],'Coneq',[],'Vio',[],'Ks',[]);

                X_updata=lhsdesign(pop_num,vari_num,'iterations',50,'criterion','maximin').*(up_bou-low_bou)+low_bou;

                % updata data lib by x_list
                [self.data_lib,~,~,~,~,~,~,~,NFE_updata]=dataUpdata(self.data_lib,X_updata,0,self.WRIRE_FILE_FLAG,self.str_data_file);
                NFE=NFE+NFE_updata;
            else
                % load exist data
                self.data_lib=struct('objcon_fcn',objcon_fcn,'vari_num',vari_num,'low_bou',low_bou,'up_bou',up_bou,'con_torl',self.con_torl,...
                    'result_best_idx',[],'X',self.data_lib.X,'Obj',self.data_lib.Obj,'Con',self.data_lib.Con,'Coneq',self.data_lib.Coneq,'Vio',self.data_lib.Vio,'Ks',self.data_lib.Ks);
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

            fmincon_options=optimoptions('fmincon','Display','none','Algorithm','sqp','ConstraintTolerance',0);

            result_x_best = zeros(self.iter_max,vari_num);
            result_obj_best = zeros(self.iter_max,1);

            iter = iter+1;

            next_search_mode = 'G'; % 'G' is global search,'l' is local search
            while ~done
                search_mode = next_search_mode;
                % nomalization all data by max obj and to create surrogate model
                [X,Obj,Con,Coneq,Vio,Ks]=dataLoad(self.data_lib);

                [X_surr,Obj_surr,Con_surr,Coneq_surr,~,~,...
                    obj_max,con_max_list,coneq_max_list,vio_max_list,ks_max_list]=self.getModelData...
                    (X,Obj,Con,Coneq,Vio,Ks,self.nomlz_value);

                if search_mode == 'G'
                    % global search
                    x_infill = self.searchGlobal...
                        (X_surr,Obj_surr,Con_surr,Coneq_surr,...
                        vari_num,low_bou,up_bou,con_fcn_cheap,...
                        pop_num);
                else
                    % local search
                    x_infill = self.searchLocal...
                        (X_surr,Obj_surr,Con_surr,Coneq_surr,...
                        vari_num,low_bou,up_bou,con_fcn_cheap,...
                        pop_num,RBF_number,fmincon_options);
                end

                % updata infill point
                [self.data_lib,x_infill,obj_infill,~,~,vio_infill,~,repeat_idx,NFE_updata] = ...
                    dataUpdata(self.data_lib,x_infill,self.protect_range,self.WRIRE_FILE_FLAG,self.str_data_file);
                NFE = NFE+NFE_updata;

                % process error
                if isempty(x_infill)
                    % continue;
                    x_infill = X(repeat_idx,:);
                    obj_infill = Obj(repeat_idx,:);
                    if ~isempty(Con)
                        con_infill = Con(repeat_idx,:);
                    end
                    if ~isempty(Coneq)
                        coneq_infill = Coneq(repeat_idx,:);
                    end
                    if ~isempty(Vio)
                        vio_local_infill = Vio(repeat_idx,:);
                    end
                else
                    if self.expensive_con_flag
                        Bool_feas = [Bool_feas;vio_infill == 0];
                    end
                end

                improve_flag = false(1);
                if self.expensive_con_flag
                    min_vio = min(Vio);
                    min_obj = min(Obj([Bool_feas(1:end-1);false(1)]),[],1);

                    % if all point is infeasible,violation of point infilled is
                    % less than min violation of all point means improve.if
                    % feasible point exist,obj of point infilled is less than min
                    % obj means improve
                    if vio_infill < min_vio
                        if ~isempty(min_obj)
                            if obj_infill < min_obj
                                % improve, continue local search
                                improve_flag = true(1);
                            end
                        else
                            % improve, continue local search
                            improve_flag = true(1);
                        end
                    end
                else
                    min_obj = min(Obj(1:end-1));

                    % obj of point infilled is less than min obj means improve
                    if obj_infill < min_obj
                        % imporve, continue local search
                        improve_flag = true(1);
                    end
                end

                % if no imporve begin local search or global search
                if ~improve_flag
                    next_search_mode = flipSearchMode(next_search_mode);
                end

                if self.DRAW_FIGURE_FLAG && vari_num < 3
                    interpVisualize(self.obj_surr,low_bou,up_bou);
                    line(x_infill(1),x_infill(2),obj_infill./obj_max*self.nomlz_value,'Marker','o','color','r');
                end

                % find best result to record
                [X,Obj,~,~,Vio,~]=dataLoad(self.data_lib);
                idx=find(Vio == 0);
                if isempty(idx)
                    [vio_best,min_idx]=min(Vio);
                    obj_best=Obj(min_idx);
                    x_best=X(min_idx,:);
                else
                    [obj_best,min_idx]=min(Obj(idx));
                    vio_best=0;
                    x_best=X(idx(min_idx),:);
                end

                if self.INFORMATION_FLAG
                    fprintf('model: %s obj:    %f    violation:    %f    NFE:    %-3d\n',search_mode,obj_best,vio_best,NFE);
                    %         fprintf('iteration:          %-3d    NFE:    %-3d\n',iteration,NFE);
                    %         fprintf('x:          %s\n',num2str(x_infill));
                    %         fprintf('value:      %f\n',obj_infill);
                    %         fprintf('violation:  %s  %s\n',num2str(con_infill),num2str(coneq_infill));
                    %         fprintf('\n');
                end

                result_x_best(iter,:) = x_best;
                result_obj_best(iter,:) = obj_best;
                iter = iter+1;

                % forced interrupt
                if iter > self.iter_max || NFE >= self.NFE_max
                    done = 1;
                end

                % convergence judgment
                if self.CONVERGENCE_JUDGMENT_FLAG
                    if ( iter > 2 && ...
                            abs((obj_best-obj_best_old)/obj_best_old) < self.obj_torl && ...
                            ((~isempty(vio_best) && vio_best == 0) || isempty(vio_best)) )
                        done = 1;
                    end
                end

                obj_best_old = obj_best;
            end

            result_x_best = result_x_best(1:iter-1,:);
            result_obj_best = result_obj_best(1:iter-1);

            output.result_x_best = result_x_best;
            output.result_obj_best = result_obj_best;
            output.data_lib = self.data_lib;

            function next_search_flag=flipSearchMode(next_search_flag)
                switch next_search_flag
                    case 'L'
                        next_search_flag = 'G';
                    case 'G'
                        next_search_flag = 'L';
                end
            end

        end

        function x_infill = searchGlobal...
                (self,X_surr,Obj_surr,Con_surr,Coneq_surr,...
                vari_num,low_bou,up_bou,con_fcn_cheap,...
                pop_num)
            % find global infill point function
            %

            % step 5
            % rank x_list data
            [x_rank_list,~,~,~] = self.rankData...
                (X_surr,Obj_surr,Con_surr,Coneq_surr,...
                con_fcn_cheap,self.con_torl);

            % step 6
            % only the first population_number will be use
            X_best_pop = x_rank_list(1:pop_num,:);

            % differ evolution mutations
            X_new_R1 = DERand...
                (low_bou,up_bou,X_best_pop,self.scaling_factor,pop_num,1);
            X_new_R2 = DERand...
                (low_bou,up_bou,X_best_pop,self.scaling_factor,pop_num,2);
            X_new_CR = DECurrentRand...
                (low_bou,up_bou,X_best_pop,self.scaling_factor);
            X_new_CB = DECurrentBest...
                (low_bou,up_bou,X_best_pop,self.scaling_factor,1);

            % differ evolution crossover
            X_new_R1 = DECrossover...
                (low_bou,up_bou,X_best_pop,X_new_R1,self.cross_rate);
            X_new_R2 = DECrossover...
                (low_bou,up_bou,X_best_pop,X_new_R2,self.cross_rate);
            X_new_CR = DECrossover...
                (low_bou,up_bou,X_best_pop,X_new_CR,self.cross_rate);
            X_new_CB = DECrossover...
                (low_bou,up_bou,X_best_pop,X_new_CB,self.cross_rate);

            % find global infill point base kriging model from offspring X
            X_DE = [X_new_R1;X_new_R2;X_new_CR;X_new_CB];

            % step 4
            % construct RBF model
            % base on distance to x_list to repace predict variance
            [self.obj_fcn_surr,self.con_fcn_surr,output_model]=self.getSurrogateFcn...
                (X_surr,Obj_surr,Con_surr,Coneq_surr);
            self.obj_surr=output_model.RBF_obj;
            self.con_surr_list=output_model.RBF_con_list;
            self.coneq_surr_list=output_model.RBF_coneq_list;

            % evaluate each x_offspring obj and constraints
            [obj_pred_DE_list] = self.obj_fcn_surr(X_DE);

            % variance were replaced by x distance to exist x
            % distance scale
            dis = zeros(size(X_DE,1),size(X_surr,1));
            for vari_idx = 1:vari_num
                dis = dis+(X_DE(:,vari_idx)-X_surr(:,vari_idx)').^2;
            end
            D = min(sqrt(dis),[],2);
            D_min = min(D); D_max = max(D);
            obj_var_DE_list = (D_max-D)./(D_max-D_min);
            if self.expensive_con_flag

                if ~isempty(self.con_fcn_surr)
                    [con_pred_DE_list,coneq_pred_DE_list] = self.con_fcn_surr(X_DE);
                    con_var_DE_list = obj_var_DE_list;
                    coneq_var_DE_list = obj_var_DE_list;
                end

                vio_DE_list = zeros(4*pop_num,1);
                if ~isempty(Con_surr)
                    vio_DE_list = vio_DE_list+sum(max(con_pred_DE_list-self.con_torl,0),2);
                end
                if ~isempty(Coneq_surr)
                    vio_DE_list = vio_DE_list+sum((abs(coneq_pred_DE_list)-self.con_torl),2);
                end

                feasi_boolean_DE_list = vio_DE_list == 0;
            else
                feasi_boolean_DE_list = true(ones(1,4*pop_num));
            end

            % if have feasiable_idx_list,only use feasiable to choose
            if all(~feasi_boolean_DE_list)
                % base on constaints improve select global infill
                % lack process of equal constraints
                con_nomlz_base = max(min(Con_surr,[],1),0);
                Con_impove_prob = sum(...
                    normcdf((con_nomlz_base-con_pred_DE_list)./sqrt(con_var_DE_list)),2);
                [~,con_best_idx] = max(Con_impove_prob);
                con_best_idx = con_best_idx(1);
                x_infill = X_DE(con_best_idx,:);
            else
                % base on fitness DE point to select global infill
                if self.expensive_con_flag
                    X_DE = X_DE(feasi_boolean_DE_list,:);
                    obj_pred_DE_list = obj_pred_DE_list(feasi_boolean_DE_list);
                    obj_var_DE_list = obj_var_DE_list(feasi_boolean_DE_list);
                end

                obj_DE_min = min(obj_pred_DE_list,[],1);
                obj_DE_max = max(obj_pred_DE_list,[],1);

                DE_fitness_list = -(obj_pred_DE_list-obj_DE_min)/(obj_DE_max-obj_DE_min)+...
                    obj_var_DE_list;
                [~,fitness_best_idx] = max(DE_fitness_list);
                fitness_best_idx = fitness_best_idx(1);
                x_infill = X_DE(fitness_best_idx,:);
            end

            function X_new = DERand(low_bou,up_bou,X,F,x_number,rand_number)
                [x_number__,variable_number__] = size(X);
                X_new = zeros(x_number,variable_number__);
                for x_idx__ = 1:x_number
                    idx__ = randi(x_number__,2*rand_number+1,1);
                    X_new(x_idx__,:) = X(idx__(1),:);
                    for rand_idx__ = 1:rand_number
                        X_new(x_idx__,:) = X_new(x_idx__,:)+...
                            F*(X(idx__(2*rand_idx__),:)-X(idx__(2*rand_idx__+1),:));
                        X_new(x_idx__,:) = max(X_new(x_idx__,:),low_bou);
                        X_new(x_idx__,:) = min(X_new(x_idx__,:),up_bou);
                    end
                end
            end

            function X_new = DECurrentRand(low_bou,up_bou,X,F)
                [x_number__,variable_number__] = size(X);
                X_new = zeros(x_number__,variable_number__);
                for x_idx__ = 1:x_number__
                    idx__ = randi(x_number__,3,1);
                    X_new(x_idx__,:) = X(x_idx__,:)+...
                        F*(X(idx__(1),:)-X(x_idx__,:)+...
                        X(idx__(2),:)-X(idx__(3),:));
                    X_new(x_idx__,:) = max(X_new(x_idx__,:),low_bou);
                    X_new(x_idx__,:) = min(X_new(x_idx__,:),up_bou);
                end
            end

            function X_new = DECurrentBest(low_bou,up_bou,X,F,x_best_idx)
                [x_number__,variable_number__] = size(X);
                X_new = zeros(x_number__,variable_number__);
                for x_idx__ = 1:x_number__
                    idx__ = randi(x_number__,2,1);
                    X_new(x_idx__,:) = X(x_idx__,:)+...
                        F*(X(x_best_idx,:)-X(x_idx__,:)+...
                        X(idx__(1),:)-X(idx__(2),:));
                    X_new(x_idx__,:) = max(X_new(x_idx__,:),low_bou);
                    X_new(x_idx__,:) = min(X_new(x_idx__,:),up_bou);
                end
            end

            function X_new = DECrossover(low_bou,up_bou,X,V,C_R)
                [x_number__,variable_number__] = size(X);
                X_new = X;
                rand_number = rand(x_number__,variable_number__);
                idx__ = find(rand_number < C_R);
                X_new(idx__) = V(idx__);
                for x_idx__ = 1:x_number__
                    X_new(x_idx__,:) = max(X_new(x_idx__,:),low_bou);
                    X_new(x_idx__,:) = min(X_new(x_idx__,:),up_bou);
                end
            end

        end

        function x_infill = searchLocal...
                (self,X_surr,Obj_surr,Con_surr,Coneq_surr,...
                vari_num,low_bou,up_bou,con_fcn_cheap,...
                pop_num,RBF_num,fmincon_options)
            % find local infill point function
            %

            % sort X potential by Obj
            [X_surr,Obj_surr,Con_surr,Coneq_surr]=self.rankData...
                (X_surr,Obj_surr,Con_surr,Coneq_surr,con_fcn_cheap,self.con_torl);

            % step 8
            % rand select initial local point from x_list
            x_idx = randi(pop_num);
            x_initial = X_surr(x_idx,:);

            % select nearest point to construct RBF
            RBF_num = min(RBF_num,size(X_surr,1));
            dist = sum(((x_initial-X_surr)./(up_bou-low_bou)).^2,2);
            [~,idx_list] = sort(dist);
            idx_list = idx_list(1:RBF_num);
            X_RBF = X_surr(idx_list,:);
            Obj_RBF = Obj_surr(idx_list,:);
            if ~isempty(Con_surr)
                Con_RBF = Con_surr(idx_list,:);
            else
                Con_RBF = [];
            end
            if ~isempty(Coneq_surr)
                Coneq_RBF = Coneq_surr(idx_list,:);
            else
                Coneq_RBF = [];
            end

            % modify
            % get RBF model and function
            [self.obj_fcn_surr,self.con_fcn_surr,output_model]=self.getSurrogateFcn...
                (X_RBF,Obj_RBF,Con_RBF,Coneq_RBF);
            self.obj_surr=output_model.RBF_obj;
            self.con_surr_list=output_model.RBF_con_list;
            self.coneq_surr_list=output_model.RBF_coneq_list;

            % get local infill point
            % obtian total constraint function
            if ~isempty(self.con_fcn_surr) || ~isempty(con_fcn_cheap)
                con_fcn = @(x) self.conFcnTotal...
                    (x,self.con_fcn_surr,con_fcn_cheap);
            else
                con_fcn = [];
            end
            low_bou_local = min(X_RBF,[],1);
            up_bou_local = max(X_RBF,[],1);

            [x_infill,obj_infill,exit_flag,output_fmincon] = fmincon(self.obj_fcn_surr,x_initial,[],[],[],[],...
                low_bou_local,up_bou_local,con_fcn,fmincon_options);
        end

        function [obj_fcn_surr,con_fcn_surr,output]=getSurrogateFcn...
                (self,x_list,obj_list,con_list,coneq_list)
            % base on lib_data to create radialbasis objcon_fcn and function
            % if input objcon_fcn,function will updata objcon_fcn
            % object_function is single obj output
            % nonlcon_function is normal nonlcon_function which include con,coneq
            % con is colume vector,coneq is colume vector
            % var_function is same
            %
            [pred_fcn_obj,RBF_obj]=self.interpRadialBasisPreModel...
                (x_list,obj_list);

            if ~isempty(con_list)
                pred_fcn_con_list=cell(size(con_list,2),1);
                RBF_con_list=cell(size(con_list,2),1);
                for con_idx=1:size(con_list,2)
                    [pred_fcn_con_list{con_idx},RBF_con_list{con_idx}]=self.interpRadialBasisPreModel...
                        (x_list,con_list(:,con_idx));
                end
            else
                pred_fcn_con_list=[];
                RBF_con_list=[];
            end

            if ~isempty(coneq_list)
                pred_fcn_coneq_list=cell(size(coneq_list,2),1);
                RBF_coneq_list=cell(size(coneq_list,2),1);
                for coneq_idx=1:size(coneq_list,2)
                    [pred_fcn_coneq_list{coneq_idx},RBF_coneq_list{coneq_idx}]=self.interpRadialBasisPreModel...
                        (x_list,coneq_list(:,coneq_idx));
                end
            else
                pred_fcn_coneq_list=[];
                RBF_coneq_list=[];
            end

            obj_fcn_surr=@(X_pred) objFcnSurr(X_pred,pred_fcn_obj);
            if isempty(pred_fcn_con_list) && isempty(pred_fcn_coneq_list)
                con_fcn_surr=[];
            else
                con_fcn_surr=@(X_pred) conFcnSurr(X_pred,pred_fcn_con_list,pred_fcn_coneq_list);
            end

            output.RBF_obj=RBF_obj;
            output.RBF_con_list=RBF_con_list;
            output.RBF_coneq_list=RBF_coneq_list;
            function obj=objFcnSurr...
                    (X_pred,pred_fcn_obj)
                % connect all predict favl
                %
                obj=pred_fcn_obj(X_pred);
            end

            function [con,coneq]=conFcnSurr...
                    (X_pred,pred_fcn_con,pred_fcn_coneq)
                % connect all predict con and coneq
                %
                if isempty(pred_fcn_con)
                    con=[];
                else
                    con=zeros(size(X_pred,1),length(pred_fcn_con));
                    for con_idx__=1:length(pred_fcn_con)
                        con(:,con_idx__)=....
                            pred_fcn_con{con_idx__}(X_pred);
                    end
                end
                if isempty(pred_fcn_coneq)
                    coneq=[];
                else
                    coneq=zeros(size(X_pred,1),length(pred_fcn_coneq));
                    for coneq_idx__=1:length(pred_fcn_coneq)
                        coneq(:,coneq_idx__)=...
                            pred_fcn_coneq{coneq_idx__}(X_pred);
                    end
                end
            end
        end

    end

    methods(Static) % auxiliary function
        function [X_surr,Obj_surr,Con_surr,Coneq_surr,Vio_model,Ks_model,...
                obj_max,con_max_list,coneq_max_list,vio_max_list,ks_max_list]=getModelData...
                (X,Obj,Con,Coneq,Vio,Ks,nomlz_obj)

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

        function [con,coneq] = conFcnTotal...
                (x,nonlcon_function,con_fcn_cheap)
            con = [];
            coneq = [];
            if ~isempty(nonlcon_function)
                [expencon,expenconeq] = nonlcon_function(x);
                con = [con;expencon];
                coneq = [coneq;expenconeq];
            end
            if ~isempty(con_fcn_cheap)
                [expencon,expenconeq] = con_fcn_cheap(x);
                con = [con;expencon];
                coneq = [coneq;expenconeq];
            end
        end

        function [X,Obj,Con,Coneq,vio_list] = rankData...
                (X,Obj,Con,Coneq,con_fcn_cheap,con_torl)
            % rank data base on feasibility rule
            % infeasible is rank by sum of constraint
            % torlance to con_fcn_cheap is 0
            %
            if nargin < 6 || isempty(con_torl)
                con_torl = 0;
            end
            if nargin < 5
                con_fcn_cheap = [];
            end

            [x_number,~] = size(X);
            vio_list = zeros(x_number,1);
            if ~isempty(Con)
                vio_list = vio_list+sum(max(Con-con_torl,0),2);
            end
            if ~isempty(Coneq)
                vio_list = vio_list+sum((abs(Coneq)-con_torl),2);
            end

            % add cheap con
            for x_index = 1:size(X,1)
                if ~isempty(con_fcn_cheap)
                    [con,coneq] = con_fcn_cheap(X(x_index,:));
                    vio_list(x_index) = vio_list(x_index)+...
                        sum(max(con,0))+sum(max(abs(coneq),0));
                end
            end

            % rank data
            % infeasible data rank by violation, feasible data rank by obj
            Bool_feas = vio_list <= 0;
            all = 1:x_number;
            feasi_index_list = all(Bool_feas);
            infeasi_index_list = all(~Bool_feas);
            [~,index_list] = sort(Obj(feasi_index_list));
            feasi_index_list = feasi_index_list(index_list);
            [~,index_list] = sort(vio_list(infeasi_index_list));
            infeasi_index_list = infeasi_index_list(index_list);
            index_list = [feasi_index_list,infeasi_index_list];

            % rank by index_list
            X = X(index_list,:);
            Obj = Obj(index_list);
            if ~isempty(Con)
                Con = Con(index_list,:);
            end
            if ~isempty(Coneq)
                Coneq = Coneq(index_list,:);
            end
            vio_list = vio_list(index_list);

        end

    end

    methods(Static) % surrogate model
        function [predict_function,radialbasis_model] = interpRadialBasisPreModel...
                (X,Y,basis_function)
            % radial basis function interp pre model function
            % input initial data X,Y,which are real data
            % X,Y are x_number x variable_number matrix
            % aver_X,stdD_X is 1 x x_number matrix
            % output is a radial basis model,include X,Y,base_function
            % and predict_function
            %
            % Copyright 2023 Adel
            %
            if nargin < 3
                basis_function = [];
            end

            [x_number,variable_number] = size(X);

            % normalize data
            aver_X = mean(X);
            stdD_X = std(X);
            aver_Y = mean(Y);
            stdD_Y = std(Y);
            idx__ = find(stdD_X == 0);
            if ~isempty(idx__),stdD_X(idx__) = 1;end
            idx__ = find(stdD_Y == 0);
            if ~isempty(idx__),stdD_Y(idx__) = 1;end
            X_nomlz = (X-aver_X)./stdD_X;
            Y_nomlz = (Y-aver_Y)./stdD_Y;

            if isempty(basis_function)
                % c = (prod(max(X_nomlz)-min(X_nomlz))/x_number)^(1/variable_number);
                % basis_function = @(r) exp(-(r.^2)/c);
                basis_function = @(r) r.^3;
            end

            % initialization distance of all X
            X_dis = zeros(x_number,x_number);
            for variable_idx = 1:variable_number
                X_dis = X_dis+(X_nomlz(:,variable_idx)-X_nomlz(:,variable_idx)').^2;
            end
            X_dis = sqrt(X_dis);

            [beta,rdibas_matrix] = interpRadialBasis...
                (X_dis,Y_nomlz,basis_function,x_number);

            % initialization predict function
            predict_function = @(X_predict) interpRadialBasisPredictor...
                (X_predict,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
                x_number,variable_number,beta,basis_function);

            radialbasis_model.X = X;
            radialbasis_model.Y = Y;
            radialbasis_model.radialbasis_matrix = rdibas_matrix;
            radialbasis_model.beta = beta;

            radialbasis_model.aver_X = aver_X;
            radialbasis_model.stdD_X = stdD_X;
            radialbasis_model.aver_Y = aver_Y;
            radialbasis_model.stdD_Y = stdD_Y;
            radialbasis_model.basis_function = basis_function;

            radialbasis_model.predict_function = predict_function;

            % abbreviation:
            % num: number,pred: predict,vari: variable
            function [beta,rdibas_matrix] = interpRadialBasis...
                    (X_dis,Y,basis_function,x_number)
                % interp polynomial responed surface core function
                % calculation beta
                %
                % Copyright 2022 Adel
                %
                rdibas_matrix = basis_function(X_dis);

                % stabilize matrix
                rdibas_matrix = rdibas_matrix+eye(x_number)*1e-6;

                % solve beta
                beta = rdibas_matrix\Y;
            end

            function [Y_pred] = interpRadialBasisPredictor...
                    (X_pred,X_nomlz,aver_X,stdD_X,aver_Y,stdD_Y,...
                    x_num,vari_num,beta,basis_function)
                % radial basis function interpolation predict function
                %
                [x_pred_num,~] = size(X_pred);

                % normalize data
                X_pred_nomlz = (X_pred-aver_X)./stdD_X;

                % calculate distance
                X_dis_pred = zeros(x_pred_num,x_num);
                for vari_idx = 1:vari_num
                    X_dis_pred = X_dis_pred+...
                        (X_pred_nomlz(:,vari_idx)-X_nomlz(:,vari_idx)').^2;
                end
                X_dis_pred = sqrt(X_dis_pred);

                % predict variance
                Y_pred = basis_function(X_dis_pred)*beta;

                % normalize data
                Y_pred = Y_pred*stdD_Y+aver_Y;
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
