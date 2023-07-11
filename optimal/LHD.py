import numpy as np
import matplotlib.pyplot as plt
import copy
import random

def getLatinHypercube(sample_num, vari_num,low_bou = [],up_bou = [],X_exist = [],cheapcon_fcntion = None): 
    # generate latin hypercube desgin
    # more uniform point distribution by simulating particle motion

    # input:
    # sample num(new point to sample),vari_num
    # low_bou,up_bou,X_exist,cheapcon_fcntion

    # output:
    # X_sample,dist_min_nomlz(min distance of normalize data)
    # X_total,include all data in area

    # Copyright 2023 3 Adel

    if len(up_bou)==0:
        up_bou = np.ones((1,vari_num))
    if len(low_bou)==0:
        low_bou = np.zeros((1,vari_num))

    iteration_max = 100
    X_exist=np.array(X_exist).reshape(len(X_exist),vari_num)
    if not len(X_exist)==0:
        X_exist_nomlz = (X_exist - low_bou) / (up_bou - low_bou)
    else:
        X_exist_nomlz = np.array([]).reshape(0,vari_num)

    exist_num = X_exist.shape[0]
    total_num = sample_num + exist_num
    if sample_num < 0:
        X_total = X_exist
        X_sample = []
        dist_min_nomlz = calMinDistance(X_exist_nomlz)
        return X_sample,dist_min_nomlz,X_total

    low_bou_nomlz = np.zeros((1,vari_num))
    up_bou_nomlz = np.ones((1,vari_num))
    # obtain initial point
    if cheapcon_fcntion is not None :
        # obtian feasible point
        X_quasi_nomlz = np.array([[]]).reshape(0,vari_num)
        # check if have enough X_supply_nomlz
        iteration = 0
        while X_quasi_nomlz.shape[0] < 10 * sample_num and iteration < 500:

            X_quasi_nomlz_initial = np.random.rand(10 * sample_num,vari_num)
            qusai_idx = []
            for x_idx in range(len(X_quasi_nomlz_initial)):
                con=cheapcon_fcntion((X_quasi_nomlz_initial[x_idx,:]*(up_bou - low_bou) + low_bou).reshape(vari_num))
                if (con <= 0).all():
                    qusai_idx.append(x_idx)
            X_quasi_nomlz = np.concatenate((X_quasi_nomlz,X_quasi_nomlz_initial[qusai_idx]),axis=0)
            iteration = iteration + 1

        if iteration == 500 and X_quasi_nomlz.shape[0] < sample_num:
            raise Exception('getLatinHypercube: feasible quasi point cannot be found')
        # use fuzzy clustering to get feasible point center
        X_sample_nomlz,_ = clusteringFuzzy(X_quasi_nomlz,sample_num,2)
        X_feas_center_nomlz = X_sample_nomlz
        #     scatter(X_quasi_nomlz(:,1),X_quasi_nomlz(:,2));
        #     hold on;
        #     scatter(X_sample_nomlz(:,1),X_sample_nomlz(:,2),'red');
        #     hold off;
    else:
        X_sample_nomlz = np.random.rand(sample_num,vari_num)

    # pic_num = 1;

    iteration = 0
    fval_list = np.zeros((sample_num,1))
    gradient_list = np.zeros((sample_num,vari_num))
    while iteration < iteration_max:

        # change each x place by newton methods
        for x_idx in range(sample_num):
            # get gradient
            X_sample_exist=np.concatenate((X_sample_nomlz[:x_idx-1],X_sample_nomlz[x_idx+1:]),axis=0)
            X_sample_exist=np.concatenate((X_sample_exist,X_exist_nomlz),axis=0)
            fval_list[x_idx,:],gradient_list[x_idx,:] = calParticleEnergy(X_sample_nomlz[x_idx,:],X_sample_exist,sample_num,vari_num,low_bou_nomlz - 1 / 2 / sample_num,up_bou_nomlz + 1 / 2 / sample_num)
            
        C = (1 - iteration / iteration_max) * 0.5
        # updata partical location
        for x_idx in range(sample_num):
            x = X_sample_nomlz[x_idx,:].reshape(1,vari_num)
            gradient = gradient_list[x_idx,:].reshape(1,vari_num)
            # check if feasible
            if cheapcon_fcntion is not None:
                con = cheapcon_fcntion((x*(up_bou - low_bou) + low_bou).reshape(vari_num))
                # if no feasible,move point to close point
                if con > 0:
                    gradient = x - X_feas_center_nomlz[x_idx,:]
            gradient[gradient>0.5]=0.5
            gradient[gradient<-0.5]=-0.5
            x = x - gradient * C
            boolean = (x < low_bou_nomlz).reshape(1,vari_num)
            x[boolean] = - x[boolean]
            boolean = (x > up_bou_nomlz).reshape(1,vari_num)
            x[boolean] = 2 - x[boolean]
            X_sample_nomlz[x_idx,:] = x
        iteration = iteration + 1

    # process out of boundary point
    for x_idx in range(sample_num):
        x = X_sample_nomlz[x_idx,:].reshape(1,vari_num)
        # check if feasible
        if cheapcon_fcntion is not None :
            con = cheapcon_fcntion((x*(up_bou - low_bou) + low_bou).reshape(vari_num))
            # if no feasible,move point to close point
            if (con > 0).all():
                # search closest point
                dx_center = x - X_feas_center_nomlz
                idx = np.argmax(sum(abs(dx_center),1))
                gradient = dx_center[idx]
        x = x - gradient * C
        boolean = (x < low_bou_nomlz).reshape(1,vari_num)
        x[boolean] = - x[boolean]
        boolean = (x > up_bou_nomlz).reshape(1,vari_num)
        x[boolean] = 2 - x[boolean]
        X_sample_nomlz[x_idx,:] = x

    X_sample_nomlz = np.maximum(X_sample_nomlz,low_bou_nomlz)
    X_sample_nomlz = np.minimum(X_sample_nomlz,up_bou_nomlz)
    dist_min_nomlz = calMinDistance(np.concatenate((X_sample_nomlz,X_exist_nomlz),axis=0))
    X_sample = X_sample_nomlz*(up_bou - low_bou) + low_bou
    X_total = np.concatenate((X_sample,X_exist),axis=0)
    
    return X_sample,dist_min_nomlz,X_total


def getNestedHypercube(X_base,sample_num, vari_num,low_bou = [],up_bou = [],X_exist = []): 
    # generate nested latin hypercube design
    # SLE method is used(sample and iteration, select max min distance group)
    # election combination mode of point and find best combination

    # input:
    # X_base(which will be sample), sample num(new point to sample), ...
    # vari_num, low_bou, up_bou, X_exist(exist point)

    # output:
    # X_sample, dist_min_nomlz(min distance of normalize data)
    # X_total(include all data in area)

    if len(up_bou)==0:
        up_bou = np.ones((1,vari_num))
    if len(low_bou)==0:
        low_bou = np.zeros((1,vari_num))

    iteration_max = 100 * sample_num
    # check X_exist if meet boundary
    X_exist=np.array(X_exist).reshape(len(X_exist),vari_num)
    if not len(X_exist)==0 :
        X_exist_nomlz = (X_exist - low_bou) / (up_bou - low_bou)
    else:
        X_exist_nomlz = np.array([]).reshape(0,vari_num)

    exist_num = X_exist_nomlz.shape[0]
    total_num = sample_num + exist_num
    if sample_num <= 0:
        X_total = X_exist
        X_sample = []
        dist_min_nomlz = calMinDistance(X_exist_nomlz)
        return X_sample,dist_min_nomlz,X_total

    # get quasi-feasible point
    X_base_nomlz = (X_base - low_bou) / (up_bou - low_bou)
    # iterate and get final x_supply_list
    iteration = 0
    x_supply_quasi_num = X_base_nomlz.shape[1-1]
    dist_min_nomlz = 0
    X_sample_nomlz = []
    # dist_min_nomlz_result = zeros(1,iteration);
    while iteration <= iteration_max:

        # random select x_new_num X to X_trial_nomlz
        x_select_idx = random.sample(range(x_supply_quasi_num),sample_num)
        # get distance min itertion X_
        distance_min_iteration = calMinDistanceIter(X_base_nomlz[x_select_idx,:],X_exist_nomlz)
        # if distance_min_iteration is large than last time
        if distance_min_iteration > dist_min_nomlz:
            dist_min_nomlz = distance_min_iteration
            X_sample_nomlz = X_base_nomlz[x_select_idx,:]
        iteration = iteration + 1
        #     dist_min_nomlz_result(iteration) = dist_min_nomlz;


    dist_min_nomlz = calMinDistance(np.concatenate((X_sample_nomlz,X_exist_nomlz),axis=0))
    X_sample = X_sample_nomlz*(up_bou - low_bou) + low_bou
    X_total = np.concatenate((X_sample,X_exist),axis=0)

    return X_sample,dist_min_nomlz,X_total


def calMinDistance(X): 
    # get distance min from X
    X=np.array(X)
    
    if len(X)==0:
        distance_min = []
        return distance_min

    # sort x_supply_list_initial to decrese distance calculate times
    X = np.array(sorted(X, key = lambda x:x[0]))
    
    sample_num,vari_num = X.shape
    distance_min = vari_num
    for x_idx in range(sample_num):
        x_curr = X[x_idx,:]
        x_next_idx = x_idx + 1
        # first dimension only search in min_distance
        search_range = vari_num
        while x_next_idx < sample_num and (X[x_next_idx,0] - X[x_idx,0]) ** 2 < search_range:

            x_next = X[x_next_idx]
            distance_temp = sum((x_next - x_curr) ** 2)
            if distance_temp < distance_min:
                distance_min = distance_temp
            if distance_temp < search_range:
                search_range = distance_temp
            x_next_idx = x_next_idx + 1

    distance_min = np.sqrt(distance_min)
    
    return distance_min


def calMinDistanceIter(X,X_exist): 
    # get distance min from X
    # this version do not consider distance between x exist
    X=np.array(X)
    X_exist=np.array(X_exist)
    
    # sort x_supply_list_initial to decrese distance calculate times
    X = np.array(sorted(X, key = lambda x:x[0]))
    
    sample_num,vari_num = X.shape
    exist_num=X_exist.shape[0]
    
    distance_min = vari_num
    for x_idx in range(sample_num):
        x_curr = X[x_idx,:]
        x_next_idx = x_idx + 1
        # only search in min_distance(X had been sort)
        search_range__ = vari_num
        while x_next_idx < sample_num and (X[x_next_idx,0] - X[x_idx,0]) ** 2 < search_range__:

            x_next__ = X[x_next_idx,:]
            distance_temp__ = sum((x_next__ - x_curr) ** 2)
            if distance_temp__ < distance_min:
                distance_min = distance_temp__
            if distance_temp__ < search_range__:
                search_range__ = distance_temp__
            x_next_idx = x_next_idx + 1

        for x_exist_idx in range(exist_num):
            x_next__ = X_exist[x_exist_idx,:]
            distance_temp__ = sum((x_next__ - x_curr) ** 2)
            if distance_temp__ < distance_min:
                distance_min = distance_temp__

    distance_min = np.sqrt(distance_min)
    
    return distance_min


def calParticleEnergy(x = None,X_surplus = None,sample_num = None,vari_num = None,low_bou = None,up_bou = None): 
    # fcntion describe distance between X and X_supply
    # x is colume vector and X_surplus is matrix which is num-1 x var
    # low_bou_limit__ and up_bou_limit__ is colume vector
    # variable in colume

    a__ = 10
    a_bou__ = 10
    sign__ = ((x > X_surplus) - 0.5) * 2
    xi__ = np.multiply(- a__ * (x - X_surplus),sign__)
    psi__ = a_bou__ * (low_bou - x)
    zeta__ = a_bou__ * (x - up_bou)
    exp_psi__ = np.exp(psi__)
    exp_zeta__ = np.exp(zeta__)

    sum_xi__ = (np.sum(xi__, 1) / vari_num).reshape((X_surplus.shape[0],1))
    exp_sum_xi__ = np.exp(sum_xi__)
    exp_xi__ = np.exp(xi__)
    sum_exp_xi__ = (np.sum(exp_xi__, 1) / vari_num).reshape((X_surplus.shape[0],1))
    # get fval
    fval = (np.sum(sum_exp_xi__, 0) + np.sum(exp_sum_xi__, 0)) / 2 / sample_num + np.sum(exp_psi__ + exp_zeta__, 1) / vari_num * 0.1

    # get gradient
    gradient = (np.sum(np.multiply(- a__ * sign__,exp_sum_xi__), 0) + np.sum(np.multiply(- a__ * sign__,exp_xi__), 0)) / 2 / vari_num / sample_num + (- a_bou__ * exp_psi__ + a_bou__ * exp_zeta__) / vari_num * 0.1
    return fval,gradient


def clusteringFuzzy(X,classify_num,m = 2): 
    # get fuzzy cluster model
    # X is x_num x vari_num matrix
    # center_list is classify_num x vari_num matrix

    iteration_max = 100
    torlance = 1e-06
    x_num,vari_num = X.shape
    # normalization data
    aver_X = np.mean(X,axis=0).reshape((1,vari_num))
    stdD_X = np.std(X,ddof =1,axis=0).reshape((1,vari_num))
    stdD_X[stdD_X==0]=1

    X_nomlz = (copy.deepcopy(X) - aver_X) / stdD_X
    # X_nomlz=X
    # if x_num equal 1,clustering cannot done
    if x_num == 1:
        FC_model={'X':X,'X_normalize':X_nomlz,'center_list' : X,'fval_loss_list':[]}
        return center_list,FC_model

    U = np.zeros((classify_num,x_num))
    center_list = np.random.rand(classify_num,vari_num) * 4-2
    iteration = 0
    done = 0
    fval_loss_list = np.zeros((iteration_max+1,1))
    # get X_center_dis_sq
    X_center_dis_sq = np.zeros((classify_num,x_num))
    for classify_idx in range(classify_num):
        for x_idx in range(x_num):
            dx_center=X_nomlz[x_idx,:] - center_list[classify_idx]
            X_center_dis_sq[classify_idx,x_idx] = np.dot(dx_center,dx_center)

    while not done :
        # updata classify matrix U
        for classify_idx in range(classify_num):
            for x_idx in range(x_num):
                U[classify_idx,x_idx] = 1 / sum((X_center_dis_sq[classify_idx,x_idx] / X_center_dis_sq[:,x_idx]) ** (1 / (m - 1)))
        # updata center_list
        center_list_old = copy.deepcopy(center_list)
        for classify_idx in range(classify_num):
            U_tran=np.transpose(U[classify_idx]).reshape(x_num,1)
            up=np.sum(((U_tran) ** m*X_nomlz), 0)
            down=np.sum((U_tran) ** m, 0)
            center_list[classify_idx] =  up/ down
        # updata X_center_dis_sq
        X_center_dis_sq = np.zeros((classify_num,x_num))
        for classify_idx in range(classify_num):
            for x_idx in range(x_num):
                dx_center=X_nomlz[x_idx,:] - center_list[classify_idx]
                X_center_dis_sq[classify_idx,x_idx] = np.dot(dx_center,dx_center)
        #     plot(center_list(:,1),center_list(:,2));
        # forced interrupt
        if iteration >= iteration_max:
            done = 1
        # convergence judgment
        if sum(sum(center_list_old - center_list) ** 2) < torlance:
            done = 1
        fval_loss_list[iteration] = sum(sum(np.multiply(U ** m,X_center_dis_sq)))
        iteration = iteration + 1

    fval_loss_list=fval_loss_list[:iteration]
    center_list = np.multiply(center_list,stdD_X) + aver_X
    FC_model={'X':X,'X_normalize':X_nomlz,'center_list' : X,'fval_loss_list':[]}
    
    return center_list,FC_model


if __name__ == '__main__':
    X_LF,_,_=getLatinHypercube(20, 2,low_bou = [],up_bou = [],X_exist = [],cheapcon_fcntion = lambda x:sum(x**2)-1)
    X_HF,_,_=getNestedHypercube(X_LF,5,2)
    
    plt.scatter(X_LF[:,0],X_LF[:,1])
    plt.scatter(X_HF[:,0],X_HF[:,1],marker='x')
    plt.show()
    