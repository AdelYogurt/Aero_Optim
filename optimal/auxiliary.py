import numpy as np
import random
from scipy import io
import time
import copy

from GP import GPC, GPCMF
from LHD import getLatinHypercube, getNestedHypercube
from SM import RBF, RBFMF
import matplotlib.pyplot as plt

# GPC inf

def conGPCFcn(x, model_GPC):
    # fcn to obtain probability pred fcn
    Y_pred, possibility, miu_pred, var_pred = model_GPC.predict(x)
    con = possibility
    return con

# other function


def calMinDistanceIter(XMat, X_exist):
    # get distance min from XMat
    # this version do not consider distance between x exist
    XMat = np.array(XMat)
    X_exist = np.array(X_exist)

    # sort x_supply_list_initial to decrese distance calculate times
    XMat = np.array(sorted(XMat, key=lambda x: x[0]))

    sample_num, vari_num = XMat.shape
    exist_num = X_exist.shape[0]

    distance_min = vari_num
    for x_idx in range(sample_num):
        x_curr = XMat[x_idx, :]
        x_next_idx = x_idx + 1
        # only search in min_distance(XMat had been sort)
        search_range__ = vari_num
        while x_next_idx < sample_num and (XMat[x_next_idx, 0] - XMat[x_idx, 0]) ** 2 < search_range__:

            x_next__ = XMat[x_next_idx, :]
            distance_temp__ = sum((x_next__ - x_curr) ** 2)
            if distance_temp__ < distance_min:
                distance_min = distance_temp__
            if distance_temp__ < search_range__:
                search_range__ = distance_temp__
            x_next_idx = x_next_idx + 1

        for x_exist_idx in range(exist_num):
            x_next__ = X_exist[x_exist_idx, :]
            distance_temp__ = sum((x_next__ - x_curr) ** 2)
            if distance_temp__ < distance_min:
                distance_min = distance_temp__

    distance_min = np.sqrt(distance_min)

    return distance_min


def calViolation(con_list, coneq_list, con_torl):
    # calculate violation of data

    if con_list.shape[1] == 0 and coneq_list.shape[1] == 0:
        vio_list = np.array([[]])
    else:
        vio_list = np.zeros(((max(con_list.shape[0], coneq_list.shape[0])), 1))
        if not con_list.shape[1] == 0:
            con_list_bias = np.array(con_list) - con_torl
            con_list_bias[con_list_bias < 0] = 0
            vio_list = vio_list + \
                np.sum(con_list_bias, axis=1)
        if not coneq_list.shape[1] == 0:
            coneq_list_bias = np.array(coneq_list) - con_torl
            vio_list = vio_list + \
                np.sum(np.abs(coneq_list_bias), axis=1)

    return vio_list.reshape((-1, 1))


def getParetoFront(DataMat):
    # distinguish pareto front of data list
    # DataMat is x_num x data_number matrix
    # notice if all data of x1 is less than x2,x1 domain x2
    DataMat = np.array(DataMat)
    x_num = len(DataMat)
    pareto_idx_list = []

    # select no domain filter
    for x_idx in range(x_num):
        data = DataMat[x_idx, :]
        pareto_idx = 0
        add_filter_flag = True
        while pareto_idx < len(pareto_idx_list):

            # compare x with exit pareto front point
            x_pareto_idx = pareto_idx_list[pareto_idx]
            # contain constraint of x_filter
            data_pareto = DataMat[x_pareto_idx, :]
            # compare x with x_pareto
            judge = data >= data_pareto
            if np.all(judge):
                add_filter_flag = False
                break
            # if better or equal than exit pareto point,reject pareto point
            judge = data <= data_pareto
            if np.all(judge):
                pareto_idx_list.pop(pareto_idx)
                pareto_idx = pareto_idx - 1

            pareto_idx = pareto_idx + 1

        # add into pareto list if possible
        if add_filter_flag:
            pareto_idx_list.append(x_idx)

    return pareto_idx_list


def totalConFcn(x, con_fcn, cheapcon_fcn, con_GPC_fcn):
    con = np.array([[]])
    coneq = np.array([[]])
    if con_fcn is not None:
        expencon, expenconeq = con_fcn(x)
        con = np.append(con, expencon)
        coneq = np.append(coneq, expenconeq)

    if cheapcon_fcn is not None:
        expencon, expenconeq = cheapcon_fcn(x)
        con = np.append(con, expencon)
        coneq = np.append(coneq, expenconeq)

    if con_GPC_fcn is not None:
        expencon, expenconeq = con_GPC_fcn(x)
        con = np.append(con, expencon)
        coneq = np.append(coneq, expenconeq)
    con = np.reshape(con, (1, -1))
    coneq = np.reshape(coneq, (1, -1))
    return con, coneq


# inf of minimize
def objFcnInf(obj_fcn, x):
    obj = obj_fcn(x)
    obj = obj[0, 0]
    return obj


def conFcnInf(total_con_Fcn, x):
    if total_con_Fcn is None:
        con_list = [0]
    else:
        con, coneq = total_con_Fcn(x)
        con_list = list(con[0, :])
    return con_list


def coneqFcnInf(total_con_Fcn, x):
    if total_con_Fcn is None:
        coneq_list = [0]
    else:
        con, coneq = total_con_Fcn(x)
        coneq_list = list(coneq[0, :])
    return coneq_list

# filter


def trainFilter(data_lib, xMat_infill, hyp, train_num, expensive_con_flag, Bool_conv):
    # train filter of gaussian process classifer
    Bool_conv = np.array(copy.deepcopy(Bool_conv))
    XMat, ObjMat, __, __, VioMat, KsMat = data_lib.dataLoad()
    # base on distance select usable point
    x_dist = np.sum(np.abs(XMat - xMat_infill), axis=1)
    idx = np.argsort(x_dist)
    ObjMat = ObjMat[idx[range(train_num)], 0].reshape((-1, 1))
    KsMat = KsMat[idx[range(train_num)], 0].reshape((-1, 1))
    XMat = XMat[idx[range(train_num)], :]
    VioMat = VioMat[idx[range(train_num)], 0].reshape((-1, 1))
    Bool_conv = Bool_conv[idx[range(train_num)]]
    Bool_feas = VioMat == 0

    if expensive_con_flag:
        # base on filter to decide which x should be choose
        # pareto_idx_list = getParetoFront([ObjMat(~Bool_feas),KsMat(~Bool_feas)]);
        pareto_idx_list = getParetoFront(
            np.concatenate((ObjMat, KsMat), axis=1))
        Y = np.ones((len(XMat), 1))
        Y[pareto_idx_list] = - 1
        #     Y(Bool_feas) = -1; # can go into feasiable area
        Y[Bool_conv] = 1

        x_pareto_center = (np.sum(
            XMat[pareto_idx_list, :], axis=0) / len(pareto_idx_list)).reshape((-1, xMat_infill.shape[1]))

    else:
        obj_threshold = np.quantile(ObjMat(not Bool_feas), 0.25)
        Y = np.ones((XMat.shape[1-1], 1))
        Y[ObjMat < obj_threshold] = - 1
        Y[Bool_feas] = - 1
        Y[Bool_conv] = 1
        x_pareto_center = (np.sum(
            XMat[pareto_idx_list, :], axis=0) / len(pareto_idx_list)).reshape((-1, xMat_infill.shape[1]))

    model_GPC = GPC(XMat, Y, hyp)
    model_GPC.train()

    # draw_X, draw_Y = np.meshgrid(np.linspace(0, 100, 21), np.linspace(0, 100, 21))
    # draw_Point = np.concatenate(
    #     (draw_X.reshape((441, 1)), draw_Y.reshape((441, 1))), axis=1)
    # draw_Z,possibility,miu_pred,var_pred = model_GPC.pred(draw_Point)
    # possibility=possibility.reshape((21, 21))

    # plt.contour(draw_X, draw_Y, possibility)
    # plt.show()

    hyp = model_GPC.hyp
    def pred_fcn_GPC(X_pred): return model_GPC.predict(X_pred)

    return pred_fcn_GPC, model_GPC, x_pareto_center, hyp


def trainFilterMF(data_lib_HF, data_lib_LF, Ratio, xMat_infill, hyp_MF, train_num, expensive_con_flag, Bool_conv_HF, Bool_conv_LF):
    # train filter of gaussian process classifer
    Bool_conv_HF = np.array(copy.deepcopy(Bool_conv_HF))
    Bool_conv_LF = np.array(copy.deepcopy(Bool_conv_LF))

    ratio_HF = Ratio[0]
    ratio_LF = Ratio[1]
    XMat_HF, ObjMat_HF, __, __, VioMat_HF, KsMat_HF = data_lib_HF.dataLoad()
    train_num_HF = min(train_num * ratio_HF, XMat_HF.shape[0])
    # base on dist select usable point
    x_dist = np.sum(np.abs(XMat_HF - xMat_infill), 1)
    idx = np.argsort(x_dist)
    ObjMat_HF = ObjMat_HF[idx[range(train_num_HF)], 0].reshape((-1, 1))
    KsMat_HF = KsMat_HF[idx[range(train_num_HF)], 0].reshape((-1, 1))
    XMat_HF = XMat_HF[idx[range(train_num_HF)], :]
    VioMat_HF = VioMat_HF[idx[range(train_num_HF)], 0].reshape((-1, 1))
    Bool_conv_HF = Bool_conv_HF[idx[range(train_num_HF)]]
    Bool_feas_HF = VioMat_HF == 0

    XMat_LF, ObjMat_LF, __, __, VioMat_LF, KsMat_LF = data_lib_LF.dataLoad()
    train_num_LF = min(train_num * ratio_LF, XMat_LF.shape[1-1])
    # base on dist select usable point
    x_dist = np.sum(np.abs(XMat_LF - xMat_infill), 1)
    idx = np.argsort(x_dist)
    ObjMat_LF = ObjMat_LF[idx[range(train_num_LF)], 0].reshape((-1, 1))
    KsMat_LF = KsMat_LF[idx[range(train_num_LF)], 0].reshape((-1, 1))
    XMat_LF = XMat_LF[idx[range(train_num_LF)], :]
    VioMat_LF = VioMat_LF[idx[range(train_num_LF)], 0].reshape((-1, 1))
    Bool_conv_LF = Bool_conv_LF[idx[range(train_num_LF)]]
    Bool_feas_LF = VioMat_LF == 0

    if expensive_con_flag:
        # base on filter to decide which x should be choose
        pareto_idx_list = getParetoFront(
            np.concatenate((ObjMat_HF, KsMat_HF), axis=1))
        Y_HF = np.ones((len(XMat_HF), 1))
        Y_HF[pareto_idx_list] = - 1
        #     Class(Bool_feas) = -1; # can go into feasiable area
        Y_HF[Bool_conv_HF] = 1

        x_pareto_center = np.sum(
            XMat_HF[pareto_idx_list, :], 1-1) / len(pareto_idx_list)

        pareto_idx_list = getParetoFront(
            np.concatenate((ObjMat_LF, KsMat_LF), axis=1))

        Y_LF = np.ones((len(XMat_LF), 1))
        Y_LF[pareto_idx_list] = - 1
        #     Class(Bool_feas) = -1; # can go into feasiable area
        Y_LF[Bool_conv_LF] = 1
    else:
        fval_threshold = np.quantile(ObjMat_HF(not Bool_feas_HF), 0.25)

        Y_HF = np.ones((XMat_HF.shape[1-1], 1))
        Y_HF[ObjMat_HF < fval_threshold] = - 1
        Y_HF[Bool_conv_HF] = 1

        Y_LF = np.ones((XMat_LF.shape[1-1], 1))
        Y_LF[ObjMat_LF < fval_threshold] = - 1
        Y_LF[Bool_conv_LF] = 1

        x_pareto_center = np.sum(
            XMat_HF[ObjMat_LF < fval_threshold, :], axis=0) / np.sum(ObjMat_HF < fval_threshold)

    model_GPCMF = GPCMF(XMat_HF, Y_HF, XMat_LF, Y_LF, hyp_MF)
    model_GPCMF.train()

    def pred_fcn_GPCMF(X_pred): return model_GPCMF.predict(X_pred)

    hyp_MF = model_GPCMF.hyp
    return pred_fcn_GPCMF, model_GPCMF, x_pareto_center, hyp_MF

# model function


def getModelData(XMat, ObjMat, Con, Coneq, VioMat, KsMat, nomlz_obj):
    X_model = XMat
    obj_max = np.amax(np.abs(ObjMat), axis=0)
    Obj_model = ObjMat / obj_max * nomlz_obj
    if not len(Con) == 0:
        con_max_list = np.amax(np.abs(Con), axis=0)
        Con_model = Con / con_max_list * nomlz_obj
    else:
        con_max_list = np.array([[]])
        Con_model = np.array([[]])

    if not len(Coneq) == 0:
        coneq_max_list = np.amax(np.abs(Coneq), axis=0)
        Coneq_model = Coneq / coneq_max_list * nomlz_obj
    else:
        coneq_max_list = np.array([[]])
        Coneq_model = np.array([[]])

    if not len(VioMat) == 0:
        vio_max_list = np.amax(np.abs(VioMat), axis=0)
        Vio_model = VioMat / vio_max_list * nomlz_obj
    else:
        vio_max_list = np.array([[]])
        Vio_model = np.array([[]])

    if not len(KsMat) == 0:
        ks_max_list = np.amax(np.abs(KsMat), axis=0)
        Ks_model = KsMat / ks_max_list * nomlz_obj
    else:
        ks_max_list = np.array([[]])
        Ks_model = np.array([[]])

    return X_model, Obj_model, Con_model, Coneq_model, Vio_model, Ks_model, obj_max, con_max_list, coneq_max_list, vio_max_list, ks_max_list


def getModelDataMF(fidelity_num, X_MF, Obj_MF, Con_MF, Coneq_MF, Vio_MF, Ks_MF, nomlz_value):
    # normalize data to construct surrogate model

    X_model_MF = X_MF
    Obj_model_MF, fval_max = calNormalizeDataMF(
        fidelity_num, Obj_MF, nomlz_value)
    if not len(Con_MF) == 0:
        Con_model_MF, con_max_list = calNormalizeDataMF(
            fidelity_num, Con_MF, nomlz_value)
    else:
        Con_model_MF = []
        con_max_list = np.array([[]])

    if not len(Coneq_MF) == 0:
        Coneq_model_MF, coneq_max_list = calNormalizeDataMF(
            fidelity_num, Coneq_MF, nomlz_value)
    else:
        Coneq_model_MF = []
        coneq_max_list = np.array([[]])

    if not len(Vio_MF) == 0:
        Vio_model_MF, vio_max = calNormalizeDataMF(
            fidelity_num, Vio_MF, nomlz_value)
    else:
        Vio_model_MF = []
        vio_max = np.array([[]])

    if not len(Ks_MF) == 0:
        Ks_model_MF, ks_max = calNormalizeDataMF(
            fidelity_num, Ks_MF, nomlz_value)
    else:
        Ks_model_MF = []
        ks_max = np.array([[]])

    return X_model_MF, Obj_model_MF, Con_model_MF, Coneq_model_MF, Vio_model_MF, Ks_model_MF, fval_max, con_max_list, coneq_max_list, vio_max, ks_max


def calNormalizeDataMF(fidelity_num, Data_MF, nomlz_value):
    Data_total = Data_MF[0]
    for fidelity_idx in range(1, fidelity_num):
        Data_total = np.concatenate(
            (Data_total, Data_MF[fidelity_idx]), axis=0)
    data_max = np.amax(np.abs(Data_total), axis=0)
    Data_model_MF = []
    for fidelity_idx in range(0, fidelity_num):
        Data_model_MF.append(Data_MF[fidelity_idx] / data_max * nomlz_value)

    return Data_model_MF, data_max


def getSurrogatefcn(XMat, ObjMat, con_list, coneq_list):
    # base on lib_data to create radialbasis model and fcn
    # if input model,fcn will updata model
    # obj_fcn is single obj output
    # con_fcn is normal con_fcn which include con,coneq
    # con is colume vector,coneq is colume vector
    # var_fcn is same
    con_list = np.array(con_list)
    coneq_list = np.array(coneq_list)
    RBF_obj = RBF(XMat, ObjMat)
    RBF_obj.train()

    RBF_con_list = []
    if not len(con_list) == 0:
        for con_idx in range(con_list.shape[1]):
            RBF_con = RBF(XMat, con_list[:, con_idx])
            RBF_con.train()
            RBF_con_list.append(RBF_con)

    RBF_coneq_list = []
    if not len(coneq_list) == 0:
        for coneq_idx in range(coneq_list.shape[1]):
            RBF_coneq = RBF(XMat, coneq_list[:, coneq_idx])
            RBF_coneq.train()
            RBF_coneq_list.append(RBF_coneq)

    def obj_fcn_surrogate(X_pred): return objfcnSurrogate(X_pred, RBF_obj)
    if len(RBF_con_list) == 0 and len(RBF_coneq_list) == 0:
        con_fcn_surrogate = None
    else:
        def con_fcn_surrogate(X_pred): return confcnSurrogate(
            X_pred, RBF_con_list, RBF_coneq_list)

    output = {}
    output['RBF_obj'] = RBF_obj
    output['RBF_con_list'] = RBF_con_list
    output['RBF_coneq_list'] = RBF_coneq_list

    return obj_fcn_surrogate, con_fcn_surrogate, output


def getSurrogateFcnMF(X_MF, Obj_MF, Con_MF, Coneq_MF):
    # base on input data to generate surrogate pred fcn
    # nonlcon_fcn_surrogate if format of nonlcon fcn in fmincon
    # judge MF and SF quality and select best one

    model_obj, type_obj = getBestModel(X_MF, Obj_MF)

    model_con_list = []
    type_con_list = []
    if not len(Con_MF) == 0:
        con_num = Con_MF[0].shape[1]

        for con_idx in range(con_num):
            Con_MF_unit = {
                0: Con_MF[0][:, con_idx], 1: Con_MF[1][:, con_idx]}
            model_con, type_con = getBestModel(X_MF, Con_MF_unit)
            model_con_list.append(model_con)
            type_con_list.append(type_con)

    model_coneq_list = []
    type_coneq_list = []
    if not len(Coneq_MF) == 0:
        coneq_num = Coneq_MF[0].shape[1]

        for coneq_idx in range(coneq_num):
            Coneq_MF_unit = {
                0: Coneq_MF[0][:, coneq_idx], 1: Coneq_MF[1][:, coneq_idx]}
            model_coneq, type_coneq = getBestModel(X_MF, Coneq_MF_unit)
            model_coneq_list.append(model_coneq)
            type_coneq_list.append(type_coneq)

    def obj_fcn_surrogate(X_pred): return model_obj.predict(X_pred)
    if len(model_con_list) == 0 and len(model_coneq_list) == 0:
        con_fcn_surrogate = None
    else:
        def con_fcn_surrogate(X_pred): return confcnSurrogate(
            X_pred, model_con_list, model_coneq_list)

    output = {}
    output['model_obj'] = model_obj
    output['model_con_list'] = model_con_list
    output['model_coneq_list'] = model_coneq_list
    output['type_obj']=type_obj
    output['type_con_list']=type_con_list
    output['type_coneq_list']=type_coneq_list

    return obj_fcn_surrogate, con_fcn_surrogate, output


def getBestModel(X_MF, Obj_MF):
    # judge use single fidelity of mulit fidelity by R^2

    x_HF_num = X_MF[0].shape[0]

    RBFMF_obj = RBFMF(X_HF=X_MF[0], Y_HF=Obj_MF[0],
                      X_LF=X_MF[1], Y_LF=Obj_MF[1])
    RBFMF_obj.train()
    error_MF = (RBFMF_obj.beta / np.diag(np.linalg.solve(RBFMF_obj.H[:, 0:x_HF_num], np.eye(x_HF_num))).reshape(RBFMF_obj.beta.shape) + RBFMF_obj.alpha / np.diag(
        np.linalg.solve(RBFMF_obj.H[:, x_HF_num:], np.eye(x_HF_num))).reshape(RBFMF_obj.beta.shape)) * RBFMF_obj.stdD_Y
    Rsq_MF = 1 - sum(error_MF ** 2) / \
        sum((np.mean(Obj_MF[0]) - Obj_MF[0]) ** 2)

    RBF_obj = RBF(X_MF[0], Obj_MF[0])
    RBF_obj.train()
    error_SF = (RBF_obj.beta /
                np.diag(RBF_obj.inv_rdibas_matrix).reshape(RBF_obj.beta.shape)) * RBF_obj.stdD_Y
    Rsq_SF = 1 - sum(error_SF ** 2) / \
        sum((np.mean(Obj_MF[0]) - Obj_MF[0]) ** 2)

    if Rsq_MF > Rsq_SF:
        model = RBFMF_obj
        model_type = 'MF'
    else:
        model = RBF_obj
        model_type = 'SF'

    return model, model_type


def objfcnSurrogate(X_pred, RBF_obj):
    # connect all pred favl
    obj = RBF_obj.predict(X_pred)
    return obj


def confcnSurrogate(X_pred, RBF_con_list, RBF_coneq_list):
    # connect all pred con and coneq
    con = np.array([[]])
    coneq = np.array([[]])
    con_num = len(RBF_con_list)
    coneq_num = len(RBF_coneq_list)
    if con_num != 0:
        con = RBF_con_list[0].predict(X_pred)
        for con_idx in range(1, con_num):
            con = np.concatenate(
                (con, RBF_con_list[con_idx].predict(X_pred)), axis=1)

    if coneq_num != 0:
        coneq = RBF_coneq_list[0].predict(X_pred)
        for coneq_idx in range(1, coneq_num):
            coneq = np.coneqcatenate(
                (coneq, RBF_coneq_list[coneq_idx].predict(X_pred)), axis=1)

    return con, coneq


class DataLibrary():
    # data lib fcn
    def __init__(self, model, vari_num, low_bou, up_bou, con_torl=0, str_data_file: str = 'result_total.txt', write_file_flag: bool = True):
        self.model = model
        self.vari_num = vari_num
        self.low_bou = low_bou
        self.up_bou = up_bou
        self.con_torl = con_torl

        self.result_best_idx = list()
        self.write_file_flag = write_file_flag
        self.str_data_file = str_data_file

        self.X = np.array([[]])
        self.Obj = np.array([[]])
        self.Con = np.array([[]])
        self.Coneq = np.array([[]])
        self.Ks = np.array([[]])
        self.Vio = np.array([[]])

        if self.write_file_flag:
            with open(str_data_file, 'a') as file_data:
                file_data.write('%s\n' % time.asctime())

        pass

    def dataUpdata(self, xMat_origin_new, protect_range=0):
        # updata data lib
        # updata format:
        # vari_num,obj_number,con_number,coneq_number
        # xMat,objMat,conMat,coneqMat
        if isinstance(xMat_origin_new, list):
            xMat_origin_new = np.array(xMat_origin_new).reshape(
                (len(xMat_origin_new), -1))

        x_new_num = len(xMat_origin_new)
        xMat_new = np.array([[]])
        objMat_new = np.array([[]])
        conMat_new = np.array([[]])
        coneqMat_new = np.array([[]])
        vioMat_new = np.array([[]])
        ksMat_new = np.array([[]])
        repeat_idx = list()
        NFE = 0
        if self.write_file_flag:
            file_data = open('result_total.txt', 'a')

        # updata format:
        # vari_num,obj_number,con_number,coneq_number
        # xMat,objMat,conMat,coneqMat
        for x_idx in range(x_new_num):
            xMat = xMat_origin_new[x_idx, :].reshape((-1, self.vari_num))
            if protect_range != 0:
                # updata data with same_point_avoid protect
                # check xMat_potc if exist in data lib
                # if exist, jump updata
                distance = np.sum(
                    (np.abs(xMat - self.X) / (self.up_bou - self.low_bou)), axis=1)
                distance_min = np.amin(distance)
                min_idx = np.argmin(distance)
                if distance_min < self.vari_num * protect_range:
                    # distance to exist point of point to add is small than protect_range
                    repeat_idx.append(min_idx)
                    continue
            objMat, conMat, coneqMat = self.model(xMat)
            objMat = np.array(objMat).reshape((1, 1))
            conMat = np.array(conMat).reshape((1, -1))
            coneqMat = np.array(coneqMat).reshape((1, -1))
            NFE = NFE + 1
            # calculate vioMat
            if conMat.shape[1] == 0 and coneqMat.shape[1] == 0:
                vioMat = np.array([]).reshape((1, -1))
                ksMat = np.array([]).reshape((1, -1))
            else:
                vioMat = self.calViolation(conMat, coneqMat)
                ksMat = np.amax(np.concatenate(
                    (conMat, coneqMat), axis=1)).reshape((1, 1))
            xMat_new = np.append(xMat_new, xMat).reshape((-1, xMat.shape[1]))
            objMat_new = np.append(objMat_new, objMat).reshape(
                (-1, objMat.shape[1]))
            if not conMat.shape[1] == 0:
                conMat_new = np.append(conMat_new, conMat).reshape(
                    (-1, conMat.shape[1]))
            if not coneqMat.shape[1] == 0:
                coneqMat_new = np.append(coneqMat_new, coneqMat).reshape(
                    (-1, coneqMat.shape[1]))
            if not vioMat.shape[1] == 0:
                vioMat_new = np.append(vioMat_new, vioMat).reshape(
                    (-1, vioMat.shape[1]))
            if not ksMat.shape[1] == 0:
                ksMat_new = np.append(ksMat_new, ksMat).reshape(
                    (-1, ksMat.shape[1]))

            if self.write_file_flag:
                # write data to txt_result
                file_data.write('%d ' % (self.vari_num))
                file_data.write('%d ' % (len(objMat)))
                file_data.write('%d ' % (conMat.shape[1]))
                file_data.write('%d ' % (coneqMat.shape[1]))

                file_data.write(list(xMat))
                file_data.write(objMat)
                file_data.write(list(conMat))
                file_data.write(list(coneqMat))

            self.dataJoin(xMat, objMat, conMat, coneqMat, vioMat, ksMat)
            # record best
            if len(self.result_best_idx) == 0:
                self.result_best_idx.append(1)
            else:
                if vioMat is None or vioMat == 0:
                    if objMat <= self.Obj[self.result_best_idx[-1]]:
                        self.result_best_idx.append(len(self.X)-1)
                    else:
                        self.result_best_idx.append(
                            self.result_best_idx[-1])
                else:
                    if vioMat <= self.Vio[self.result_best_idx[-1]]:
                        self.result_best_idx.append(len(self.X)-1)
                    else:
                        self.result_best_idx.append(
                            self.result_best_idx[-1])

        if self.write_file_flag:
            file_data.close()

        return xMat_new, objMat_new, conMat_new, coneqMat_new, vioMat_new, ksMat_new, repeat_idx, NFE

    def dataLoad(self, low_bou=float("-inf"), up_bou=float("inf")):
        # updata data to exist data lib

        idx = []
        for x_idx in range(self.X.shape[0]):
            xMat = self.X[x_idx, :]
            if np.all(xMat > low_bou) and np.all(xMat < up_bou):
                idx.append(x_idx)

        XMat = self.X[idx, :].reshape((-1, self.X.shape[1]))
        ObjMat = self.Obj[idx, :].reshape((-1, self.Obj.shape[1]))
        if not self.Con.shape[1] == 0:
            con_list = self.Con[idx, :].reshape((-1, self.Con.shape[1]))
        else:
            con_list = np.array([[]])

        if not self.Coneq.shape[1] == 0:
            coneq_list = self.Coneq[idx, :].reshape(
                (-1, self.Coneq.shape[1]))
        else:
            coneq_list = np.array([[]])

        if not self.Vio.shape[1] == 0:
            vio_list = self.Vio[idx, :].reshape((-1, self.Vio.shape[1]))
        else:
            vio_list = np.array([[]])

        if not self.Ks.shape[1] == 0:
            ks_list = self.Ks[idx, :].reshape((-1, self.Ks.shape[1]))
        else:
            ks_list = np.array([[]])

        return XMat, ObjMat, con_list, coneq_list, vio_list, ks_list

    def dataJoin(self, xMat, objMat, conMat, coneqMat, vioMat, ksMat):
        # updata data to exist data lib

        self.X = np.append(self.X, xMat).reshape((-1, xMat.shape[1]))
        self.Obj = np.append(self.Obj, objMat).reshape(
            (-1, objMat.shape[1]))
        if (not self.Con.shape[1] == 0) or (not conMat.shape[1] == 0):
            self.Con = np.append(self.Con, conMat).reshape(
                (-1, conMat.shape[1]))
        if (not self.Coneq.shape[1] == 0) or (not coneqMat.shape[1] == 0):
            self.Coneq = np.append(self.Coneq, coneqMat).reshape(
                (-1, coneqMat.shape[1]))
        if (not self.Vio.shape[1] == 0) or (not vioMat.shape[1] == 0):
            self.Vio = np.append(self.Vio, vioMat).reshape(
                (-1, vioMat.shape[1]))
        if (not self.Ks.shape[1] == 0) or (not ksMat.shape[1] == 0):
            self.Ks = np.append(self.Ks, ksMat).reshape(
                (-1, ksMat.shape[1]))

        return self

    def dataSave(self, data_librart: str = 'data_lib'):
        data_lib = {'model': str(self.model), 'vari_num': self.vari_num, 'low_bou': self.low_bou, 'up_bou': self.up_bou, 'con_torl': self.con_torl, 'str_data_file': self.str_data_file,
                    'write_file_flag': self.write_file_flag, 'X': list(self.X), 'Obj': list(self.Obj), 'Con': list(self.Con), 'Coneq': list(self.Coneq), 'Vio': list(self.Vio), 'Ks': list(self.Ks), 'result_best_idx': self.result_best_idx}
        io.savemat(data_librart+'.mat', {'data_lib': data_lib})

    def calViolation(self, con_list, coneq_list):
        # calculate violation of data

        if con_list.shape[1] == 0 and coneq_list.shape[1] == 0:
            vio_list = np.array([[]])
        else:
            vio_list = np.zeros(
                ((max(con_list.shape[0], coneq_list.shape[0])), 1))
            if not con_list.shape[1] == 0:
                con_list_bias = np.array(con_list) - self.con_torl
                con_list_bias[con_list_bias < 0] = 0
                vio_list = vio_list + \
                    np.sum(con_list_bias, axis=1)
            if not coneq_list.shape[1] == 0:
                coneq_list_bias = np.array(coneq_list) - self.con_torl
                vio_list = vio_list + \
                    np.sum(np.abs(coneq_list_bias), axis=1)

        return vio_list.reshape((-1, 1))


if __name__ == "__main__":

    point_list = np.random.random((10, 2))

    # plt.scatter(point_list[:,0],point_list[:,1])
    # plt.show()
    # pareto_idx_list = getParetoFront(point_list)
    # point_select=point_list[pareto_idx_list,:]

    # plt.scatter(point_list[:,0],point_list[:,1])
    # plt.scatter(point_select[:,0],point_select[:,1],marker='x')
    # plt.show()
