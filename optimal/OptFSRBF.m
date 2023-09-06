classdef OptFSRBF < handle
    % FSRBF and SRBF-SVM optimization algorithm
    %
    % referance: [1] SHI R,LIU L,LONG T,et al. Sequential Radial Basis Function
    % Using Support Vector Machine for Expensive Design Optimization [J]. AIAA
    % Journal,2017,55(1): 214-27.
    % [2] Shi R, Liu L, Long T, et al. Filter-Based Sequential Radial Basis
    % Function Method for Spacecraft Multidisciplinary Design Optimization[J].
    % AIAA Journal, 2018, 57(3): 1019-31.
    %
    % Copyright 2022 Adel
    %
    % abbreviation:
    % obj: objective, con: constraint, iter: iteration, tor: torlance
    % fcn: function, surr: surrogate
    % lib: library, init: initial, rst: restart, potl: potential
    % pred: predict, var:variance, vari: variable, num: number
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

        penalty_SVM=100; % SVM parameter
        CF_m=2; % clustering parameter
        FROI_min=1e-3; % min boundary of interest sample area

        str_data_file='result_total.txt'
    end

    properties
        % problem parameter
        expensive_con_flag;
        KRG_simplify_flag=true(1);

        obj_fcn_surr=[];
        con_fcn_surr=[];

        obj_surr;
        con_surr_list;
        coneq_surr_list;

        model_SVM=[];
    end

    % main function
    methods
        function self=OptFSRBF(NFE_max,iter_max,obj_torl,con_torl,data_lib,save_description)
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
            sample_num_init=min((vari_num+1)*(vari_num+2)/2,5*vari_num);
            sample_num_add=vari_num;
            sample_num_data=100*sample_num_init;
            eta=1/vari_num; % space decrease coefficient

            % NFE and iteration setting
            if isempty(self.NFE_max)
                self.NFE_max=10+10*vari_num;
            end

            if isempty(self.iter_max)
                self.iter_max=20+20*vari_num;
            end

            done=0;NFE=0;iter=0;

            x_data_list=lhsdesign(sample_num_data,vari_num).*(up_bou-low_bou)+low_bou;

            % step 2
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

            ga_option=optimoptions('ga','Display','none','ConstraintTolerance',0,'MaxGenerations',10,'HybridFcn','fmincon');

            result_x_best=zeros(self.iter_max,vari_num);
            result_obj_best=zeros(self.iter_max,1);

            iter=iter+1;

            while ~done
                % step 3
                [X,Obj,Con,Coneq,Vio,Ks]=dataLoad(self.data_lib);

                [X_surr,Obj_surr,Con_surr,Coneq_surr,Vio_surr,Ks_surr,...
                    obj_max,con_max_list,coneq_max_list,vio_max_list,ks_max_list]=self.getModelData...
                    (X,Obj,Con,Coneq,Vio,Ks,self.nomlz_value);

                % get local infill point, construct surrogate objcon_fcn
                [self.obj_fcn_surr,self.con_fcn_surr,output_model]=self.getSurrFcnRBF...
                    (X_surr,Obj_surr,Con_surr,Coneq_surr);
                self.obj_surr=output_model.obj_surr;
                self.con_surr_list=output_model.con_surr_list;
                self.coneq_surr_list=output_model.coneq_surr_list;

                % construct surrogate objcon_fcn
                [self.obj_fcn_surr,self.con_fcn_surr,output_model]=self.getSurrFcnRBF...
                    (X_surr,Obj_surr,Con_surr,Coneq_surr);
                self.obj_surr=output_model.obj_surr;
                self.con_surr_list=output_model.con_surr_list;
                self.coneq_surr_list=output_model.coneq_surr_list;

                % MSP guideline to obtain x_adapt
                con_fcn_total=@(x) self.conFcnTotal(x,self.con_fcn_surr,con_fcn_cheap);
                [x_infill,obj_infill_pred,exit_flag,output_ga]=ga...
                    (self.obj_fcn_surr,vari_num,[],[],[],[],low_bou,up_bou,con_fcn_total,ga_option);

                if exit_flag == -2
                    % step 4
                    % optimal feasiblilty if do not exist feasible point
                    vio_fcn_surr=@(x) self.vioFcnSurr(x,self.con_fcn_surr);
                    [x_infill,obj_infill_pred,exit_flag,output_ga]=ga...
                        (vio_fcn_surr,vari_num,[],[],[],[],low_bou,up_bou,con_fcn_cheap,ga_option);
                end

                % updata infill point
                [self.data_lib,x_infill,obj_infill,con_infill,coneq_infill,vio_infill,ks_infill,repeat_idx,NFE_updata]=...
                    dataUpdata(self.data_lib,x_infill,self.protect_range,self.WRIRE_FILE_FLAG,self.str_data_file);
                NFE=NFE+NFE_updata;

                % process error
                if isempty(x_infill)
                    % continue;
                    x_infill=X(repeat_idx,:);
                    obj_infill=Obj(repeat_idx,:);
                    if ~isempty(Con)
                        con_infill=Con(repeat_idx,:);
                    end
                    if ~isempty(Coneq)
                        coneq_infill=Coneq(repeat_idx,:);
                    end
                    if ~isempty(Vio)
                        vio_local_infill=Vio(repeat_idx,:);
                    end
                else
                    if self.expensive_con_flag
                        Bool_feas=[Bool_feas;vio_infill == 0];
                    end
                end

                if self.DRAW_FIGURE_FLAG && vari_num < 3
                    interpVisualize(RBF_model_fval,low_bou,up_bou);
                    line(x_infill(1),x_infill(2),fval_infill/fval_max*nomlz_fval,'Marker','o','color','r','LineStyle','none')
                end

                % step 6
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
                    fprintf('fval:    %f    violation:    %f    NFE:    %-3d\n',obj_best,vio_best,NFE);
                    %         fprintf('iteration:          %-3d    NFE:    %-3d\n',iteration,NFE);
                    %         fprintf('current x:          %s\n',num2str(x_infill));
                    %         fprintf('current value:      %f\n',fval_infill);
                    %         fprintf('current violation:  %s  %s\n',num2str(con_infill),num2str(coneq_infill));
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
                    if ( iter > 2 && ...
                            abs((obj_best-obj_best_old)/obj_best_old) < self.obj_torl && ...
                            ((~isempty(vio_best) && vio_best == 0) || isempty(vio_best)) )
                        done=1;
                    end
                end

                % Interest space sampling
                if ~done
                    % step 7-1.2
                    % using SVM to identify area which is interesting
                    [pred_fcn_SVM,x_pareto_center,pareto_idx_list]=self.trainFilter();

                    if self.DRAW_FIGURE_FLAG && vari_num < 3
                        classifyVisualization(self.model_SVM,low_bou,up_bou);
                    end

                    % get data to obtain clustering center
                    x_sup_list=x_data_list(pred_fcn_SVM(x_data_list) == 1,:);

                    % step 7-3
                    if isempty(x_sup_list)
                        % no sup point found use x_pareto_center
                        x_center=x_pareto_center;
                    else
                        [x_center,FC_model]=self.clusteringFuzzy...
                            (x_sup_list,1,self.CF_m);
                    end

                    % updata FROI
                    x_infill_nomlz=(x_infill-low_bou)./(up_bou-low_bou);
                    x_center_nomlz=(x_center-low_bou)./(up_bou-low_bou);
                    FROI_range_nomlz=eta*norm(x_infill_nomlz-x_center_nomlz,2);
                    if FROI_range_nomlz < self.FROI_min
                        FROI_range_nomlz=bou_min;
                    end
                    FROI_range=FROI_range_nomlz.*(up_bou-low_bou);
                    low_bou_FROI=x_infill-FROI_range;
                    low_bou_FROI=max(low_bou_FROI,low_bou);
                    up_bou_FROI=x_infill+FROI_range;
                    up_bou_FROI=min(up_bou_FROI,up_bou);

                    if self.DRAW_FIGURE_FLAG && vari_num < 3
                        bou_line=[low_bou_FROI;[low_bou_FROI(1),up_bou_FROI(2)];up_bou_FROI;[up_bou_FROI(1),low_bou_FROI(2)];low_bou_FROI];
                        line(bou_line(:,1),bou_line(:,2));
                        line(x_infill(1),x_infill(2),'Marker','x')
                    end

                    % sampling in ISR
                    X_add=lhsdesign(sample_num_add,vari_num).*(up_bou_FROI-low_bou_FROI)+low_bou_FROI;

                    % updata data lib
                    [self.data_lib,~,~,~,~,Vio_add,~,~,NFE_updata]=...
                        dataUpdata(self.data_lib,X_add,self.protect_range,self.WRIRE_FILE_FLAG,self.str_data_file);
                    NFE=NFE+NFE_updata;
                    Bool_feas=[Bool_feas;Vio_add == 0];
                end

                x_infill_old=x_infill;
                obj_infill_old=obj_infill;
                obj_best_old=obj_best;

                % save iteration
                if ~isempty(self.save_description)
                    data_lib=self.data_lib;
                    save(self.save_description,'data_lib');
                end
            end

            result_x_best=result_x_best(1:iter-1,:);
            result_obj_best=result_obj_best(1:iter-1);

            output.result_x_best=result_x_best;
            output.result_obj_best=result_obj_best;
            output.data_lib=self.data_lib;
        end

        function [pred_func_SVM,x_pareto_center,pareto_idx_list]=trainFilter(self)
            % train filter of gaussian process classifer
            %
            [X,Obj,~,~,~,Ks]=dataLoad(self.data_lib);

            if self.expensive_con_flag
                % base on filter to decide which x should be choose
                %     pareto_idx_list=self.getParetoFront([Obj(~Bool_feas),Ks(~Bool_feas)]);
                pareto_idx_list=self.getParetoFront(Obj,Ks);

                Class=ones(size(X,1),1);
                Class(pareto_idx_list)=-1;
                %     Class(Bool_feas)=-1; % can go into feasiable area

                x_pareto_center=sum(X(pareto_idx_list,:),1)/length(pareto_idx_list);

                [pred_func_SVM,self.model_SVM]=self.classifySupportVectorMachine(X,Class,self.penalty_SVM);
            else
                obj_threshold=prctile(Obj,50-40*sqrt(NFE/self.NFE_max));

                Class=ones(size(X,1),1);
                Class(Obj < obj_threshold)=-1;

                x_pareto_center=sum(X(Obj < obj_threshold,:),1)/sum(Obj < obj_threshold);

                [pred_func_SVM,self.model_SVM]=self.classifySupportVectorMachine(X,Class,self.penalty_SVM);
            end
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
        function vio=vioFcnSurr(x,con_fcn_surr)
            % calculate violation by con_fcn_surr(x)
            %
            [con,coneq]=con_fcn_surr(x);
            vio=0;
            if ~isempty(con)
                vio=vio+sum(max(con,0));
            end
            if ~isempty(coneq)
                vio=vio+sum(abs(coneq));
            end
        end

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
        function [center_list,FC_model]=clusteringFuzzy(X,classify_number,m)
            % get fuzzy cluster model
            % X is x_number x variable_number matrix
            % center_list is classify_number x variable_number matrix
            %
            iteration_max=100;
            torlance=1e-6;

            [x_number,variable_number]=size(X);

            % normalize data
            aver_X=mean(X);
            stdD_X=std(X);
            index__=find(stdD_X == 0);
            if  ~isempty(index__),stdD_X(index__)=1; end
            X_nomlz=(X-aver_X)./stdD_X;

            % if x_number equal 1,clustering cannot done
            if x_number == 1
                center_list=X;
                FC_model.X=X;
                FC_model.X_normalize=X_nomlz;
                FC_model.center_list=X;
                FC_model.fval_loss_list=[];
                return;
            end

            U=zeros(classify_number,x_number);
            center_list=rand(classify_number,variable_number)*4-2;
            iteration=0;
            done=0;
            fval_loss_list=zeros(iteration_max,1);

            % get X_center_dis_sq
            X_center_dis_sq=zeros(classify_number,x_number);
            for classify_index=1:classify_number
                for x_index=1:x_number
                    X_center_dis_sq(classify_index,x_index)=...
                        getSq((X_nomlz(x_index,:)-center_list(classify_index,:)));
                end
            end

            while ~done
                % updata classify matrix U
                for classify_index=1:classify_number
                    for x_index=1:x_number
                        U(classify_index,x_index)=...
                            1/sum((X_center_dis_sq(classify_index,x_index)./X_center_dis_sq(:,x_index)).^(1/(m-1)));
                    end
                end

                % updata center_list
                center_list_old=center_list;
                for classify_index=1:classify_number
                    center_list(classify_index,:)=...
                        sum((U(classify_index,:)').^m.*X_nomlz,1)./...
                        sum((U(classify_index,:)').^m,1);
                end

                % updata X_center_dis_sq
                X_center_dis_sq=zeros(classify_number,x_number);
                for classify_index=1:classify_number
                    for x_index=1:x_number
                        X_center_dis_sq(classify_index,x_index)=...
                            getSq((X_nomlz(x_index,:)-center_list(classify_index,:)));
                    end
                end

                %     scatter(X_nomlz(:,1),X_nomlz(:,2));
                %     line(center_list(:,1),center_list(:,2),'Marker','o','LineStyle','None','Color','r');

                % forced interrupt
                if iteration > iteration_max
                    done=1;
                end

                % convergence judgment
                if sum(sum(center_list_old-center_list).^2)<torlance
                    done=1;
                end

                iteration=iteration+1;
                fval_loss_list(iteration)=sum(sum(U.^m.*X_center_dis_sq));
            end
            fval_loss_list(iteration+1:end)=[];

            % normalize
            center_list=center_list.*stdD_X+aver_X;

            FC_model.X=X;
            FC_model.X_normalize=X_nomlz;
            FC_model.center_list=center_list;
            FC_model.fval_loss_list=fval_loss_list;

            function sq=getSq(dx)
                % dx is 1 x variable_number matrix
                %
                sq=dx*dx';
            end

        end

        function [pred_fcn,model_SVM]=classifySupportVectorMachine...
                (X,Class,C,kernel_fcn)
            % generate support vector machine model
            % use fmincon to get alpha
            %
            % input:
            % X(x_num x vari_num matrix), Class(x_num x 1 matrix),...
            % C(penalty factor, default is empty), kernel_fcn(default is gauss kernal function)
            %
            % output:
            % pred_fcn, model_GPC
            %
            % abbreviation:
            % pred: predicted, nomlz: normalization,num: number
            % var: variance, fcn: function
            %
            if nargin < 4
                kernel_fcn=[];
                if nargin < 3
                    C=[];
                end
            end

            [x_number,variable_number]=size(X);

            % normalization data
            aver_X=mean(X);
            stdD_X=std(X);
            index__=find(stdD_X == 0);
            if  ~isempty(index__),  stdD_X(index__)=1; end
            X_nomlz=(X-aver_X)./stdD_X;

            Y=Class;

            % default kernal function
            if isempty(kernel_fcn)
                % notice after standard normal distribution normalize
                % X usually distribution in -2 to 2, so divide by 16
                %     sigma=100*log(sqrt(x_number))/variable_number^2/16;
                sigma=1/variable_number;
                kernel_fcn=@(U,V) kernelGaussian(U,V,sigma);
            end

            % initialization kernal function process X_cov
            K=kernel_fcn(X_nomlz,X_nomlz);

            % min SVM object function to get alpha
            object_function=@(alpha) -objectFunction(alpha,K,Y);
            alpha=ones(x_number,1)*0.5;
            low_bou_fmincon=0*ones(x_number,1);
            if isempty(C) || C==0
                up_bou_fmincon=[];
            else
                up_bou_fmincon=C*ones(x_number,1);
            end
            Aeq=Y';
            fmincon_options=optimoptions('fmincon','Display','none','Algorithm','sqp');
            alpha=fmincon(object_function,alpha,...
                [],[],Aeq,0,low_bou_fmincon,up_bou_fmincon,[],fmincon_options);

            % obtain other paramter
            alpha_Y=alpha.*Y;

            w=sum(alpha_Y.*X_nomlz);
            index_list=find(alpha > 1e-6); % support vector
            alpha_Y_cov=K*alpha_Y;
            b=sum(Y(index_list)-alpha_Y_cov(index_list))/length(index_list);

            % generate predict function
            pred_fcn=@(x) classifySupportVectorMachinePredictor...
                (x,X_nomlz,alpha_Y,b,aver_X,stdD_X,kernel_fcn);

            % output model
            model_SVM.X=X;
            model_SVM.Class=Class;
            model_SVM.Y=Y;
            model_SVM.X_nomlz=X_nomlz;
            model_SVM.aver_X=aver_X;
            model_SVM.stdD_X=stdD_X;
            model_SVM.alpha=alpha;
            model_SVM.w=w;
            model_SVM.b=b;
            model_SVM.kernel_function=kernel_fcn;
            model_SVM.predict_function=pred_fcn;

            function fval=objectFunction(alpha,K,Y)
                % support vector machine maximum object function
                %
                alpha=alpha(:);
                alpha_Y__=alpha.*Y;
                fval=sum(alpha)-alpha_Y__'*K*alpha_Y__/2;
            end

            function [Class_pred,Probability]=classifySupportVectorMachinePredictor...
                    (X_pred,X_nomlz,alpha_Y,b,aver_X,stdD_X,kernel_function)
                % predict_fval is 1 or -1, predict_class is 1 or 0
                %
                % x input is colume vector
                %
                X_pred_nomlz=(X_pred-aver_X)./stdD_X;
                K_pred=kernel_function(X_pred_nomlz,X_nomlz);
                Probability=K_pred*alpha_Y+b;
                Probability=1./(1+exp(-Probability));
                Class_pred=Probability > 0.5;
            end

            function K=kernelGaussian(U,V,sigma)
                % gaussian kernal function
                %
                K=zeros(size(U,1),size(V,1));
                vari_num=size(U,2);
                for vari_index=1:vari_num
                    K=K+(U(:,vari_index)-V(:,vari_index)').^2;
                end
                K=exp(-K*sigma);
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
