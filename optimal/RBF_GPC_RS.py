import sys
import os
import numpy as np
from scipy import io
from pyDOE2 import lhs
from scipy.optimize import minimize, Bounds, NonlinearConstraint
import random
import time
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from GP import GPC
from LHD import getLatinHypercube, getNestedHypercube
from SM import RBF, RBFMF, interpVisualize
from auxiliary import *


class RBF_GPC_RS():
    def __init__(self, NFE_max=None, iter_max=None, torl=1e-6, con_torl=1e-3, data_lib:dict=None):
        self.NFE_max = NFE_max
        self.iter_max = iter_max
        self.torl = torl
        self.con_torl = con_torl

        self.data_lib = data_lib

        pass

    def optimize(self, model, vari_num, low_bou, up_bou, cheapcon_fcn=None):
        # RBF-AL-KTR optimization algorithm
        # potc: potential

        # Copyright 2023 4 Adel
        low_bou = np.array(low_bou).reshape((1, vari_num))
        up_bou = np.array(up_bou).reshape((1, vari_num))

        DRAW_FIGURE_FLAG = 0
        INFORMATION_FLAG = 1
        CONVERGENCE_JUDGMENT_FLAG = 0
        WRIRE_FILE_FLAG = 0

        # hyper parameter
        sample_num_initial = 6 + 3 * vari_num
        sample_num_restart = sample_num_initial
        sample_num_add = int(np.ceil(np.log(sample_num_initial)/2))

        min_bou_interest = 0.001
        trial_num = np.minimum(100 * vari_num, 100)
        nomlz_value = 10

        protect_range = 1e-6

        identiy_torl = 0.001

        # NFE and iter setting
        if self.NFE_max is None:
            self.NFE_max = 10 + 10 * vari_num

        if self.iter_max is None:
            self.iter_max = 20 + 20 * vari_num

        done = 0
        NFE = 0
        iter = 0
        # step 1
        # generate initial data library
        if self.data_lib is not None:
            data_lib = DataLibrary(model=model, vari_num=vari_num, low_bou=low_bou, up_bou=up_bou, con_torl=self.con_torl,
                                   write_file_flag=WRIRE_FILE_FLAG)
            data_lib.X = self.data_lib['X']
            data_lib.Obj = self.data_lib['Obj']
            data_lib.Con = self.data_lib['Con']
            data_lib.Coneq = self.data_lib['Coneq']
            data_lib.Ks = self.data_lib['Ks']
            data_lib.Vio = self.data_lib['Vio']
            if (data_lib.Con.shape[1] != 0) or (data_lib.Coneq.shape[1] != 0):
                expensive_con_flag = True
            else:
                expensive_con_flag = False

        else:
            XMat_updata = lhs(vari_num, sample_num_initial) * \
                (up_bou - low_bou) + low_bou
            data_lib = DataLibrary(model=model, vari_num=vari_num, low_bou=low_bou, up_bou=up_bou, con_torl=self.con_torl,
                                   write_file_flag=WRIRE_FILE_FLAG)
            # detech expensive constraints and initialize data library
            xMat_new, objMat_new, conMat_new, coneqMat_new, vioMat_new, ksMat_new, repeat_idx, NFE_updata = data_lib.dataUpdata(
                [XMat_updata[0]], 0)
            NFE = NFE + NFE_updata
            if not len(data_lib.Vio) == 0:
                expensive_con_flag = True
            else:
                expensive_con_flag = False

            # updata data library by XMat
            xMat_new, objMat_new, conMat_new, coneqMat_new, vioMat_new, ksMat_new, repeat_idx, NFE_updata = data_lib.dataUpdata(
                XMat_updata[1:], 0)
            NFE = NFE + NFE_updata
            # find fesiable data in current data library

        if expensive_con_flag:
            Bool_feas = [(data_lib.Vio[idx, 0] == 0)
                         for idx in range(data_lib.X.shape[0])]
        Bool_conv = [False]*data_lib.X.shape[0]

        hyp = {'mean': [0], 'cov': np.array([0, 0])}

        XMat_local_best = np.array([[]])
        ObjMat_local_best = np.array([[]])
        XMat_potc = np.array([[]])
        ObjMat_potc = np.array([[]])
        VioMat_potc = np.array([[]])
        detect_local_flag = True

        SLSQP_options = {'ftol': 1e-6, 'disp': False, 'maxiter': 100}
        result_x_best = list()
        result_obj_best = list()
        iter = iter + 1
        while not done:
            # step 2
            # nomalization all data by max objMat and to create surrogate model
            XMat, ObjMat, ConMat, ConeqMat, VioMat, KsMat = data_lib.dataLoad()
            XMat_model, ObjMat_model, ConMat_model, ConeqMat_model, VioMat_model, KsMat_model, obj_max, con_max_list, coneq_max_list, vio_max_list, ks_max_list = getModelData(
                XMat, ObjMat, ConMat, ConeqMat, VioMat, KsMat, nomlz_value)
            # get local infill point, construct surrogate model
            obj_fcn_surrogate, con_fcn_surrogate, output_model = getSurrogatefcn(
                XMat_model, ObjMat_model, ConMat_model, ConeqMat_model)
            RBF_obj = output_model['RBF_obj']
            RBF_con_list = output_model['RBF_con_list']
            RBF_coneq_list = output_model['RBF_coneq_list']
            if (con_fcn_surrogate is not None) or (cheapcon_fcn is not None):
                def con_fcn_total(xMat): return totalConFcn(
                    xMat, con_fcn_surrogate, cheapcon_fcn, None)
            else:
                con_fcn_total = None

            # step 3
            if detect_local_flag:
                # detech potc local best point
                # for x_idx in range(1):
                for x_idx in range(len(XMat_model)):
                    x_initial = XMat_model[x_idx]

                    def obj_fcn_inf(xMat): return objFcnInf(
                        obj_fcn_surrogate, xMat)

                    def con_fcn_inf(xMat): return conFcnInf(
                        con_fcn_total, xMat)
                    opt_res = minimize(obj_fcn_inf, x_initial, method='SLSQP', constraints=NonlinearConstraint(con_fcn_inf, -np.inf, 0), bounds=Bounds(low_bou[0, :].reshape((vari_num,)), up_bou[0, :].reshape((vari_num,))),
                                       options=SLSQP_options)
                    xMat_potc = np.array(opt_res.x).reshape((1, -1))
                    objMat_potc_pred = np.array(opt_res.fun).reshape((1, -1))

                    conMat_potc, coneqMat_potc = con_fcn_surrogate(xMat_potc)
                    if not conMat_potc.shape[1]:
                        conMat_potc = conMat_potc / nomlz_value*con_max_list
                    if not coneqMat_potc.shape[1]:
                        coneqMat_potc = coneqMat_potc / nomlz_value*coneq_max_list
                    vioMat_potc = calViolation(
                        conMat_potc, coneqMat_potc, self.con_torl)

                    if opt_res.success:
                        # check if xMat_potc have existed
                        add_flag = True
                        if XMat_potc.shape[1] != 0:
                            for x_check_idx in range(len(XMat_potc)):
                                if np.sum(np.abs(XMat_potc[x_check_idx, :] - xMat_potc)) / vari_num < identiy_torl:
                                    add_flag = False
                                    break
                        if XMat_local_best.shape[1] != 0:
                            for x_check_idx in range(len(XMat_local_best)):
                                if np.sum(np.abs(XMat_local_best[x_check_idx, :] - xMat_potc)) / vari_num < identiy_torl:
                                    add_flag = False
                                    break

                        # updata into XMat_potc
                        if add_flag:
                            if XMat_potc.shape[1] == 0:
                                XMat_potc = xMat_potc
                                ObjMat_potc = objMat_potc_pred
                                VioMat_potc = vioMat_potc
                            else:
                                XMat_potc = np.concatenate(
                                    (XMat_potc, xMat_potc), axis=0)
                                ObjMat_potc = np.concatenate(
                                    (ObjMat_potc, objMat_potc_pred/nomlz_value*obj_max), axis=0)
                                VioMat_potc = np.concatenate(
                                    (VioMat_potc, vioMat_potc), axis=0)

                # if XMat_potc is empty, try to use KS surrogate as xMat potc
                if XMat_potc.shape[1] == 0:
                    model_ks = RBF(XMat_model, KsMat_model)
                    model_ks.train()

                    def ks_fcn_surrogate(
                        X_pred): return model_ks.predict(X_pred)

                    x_idx = np.argmin(VioMat_model)
                    x_initial = XMat_model[x_idx]
                    def obj_fcn_inf(xMat): return objFcnInf(
                        ks_fcn_surrogate, xMat)

                    def con_fcn_inf(xMat): return conFcnInf(cheapcon_fcn, xMat)
                    opt_res = minimize(con_fcn_inf, x_initial, method='SLSQP', constraints=NonlinearConstraint(obj_fcn_inf, -np.inf, 0), bounds=Bounds(low_bou[0, :].reshape((vari_num,)), up_bou[0, :].reshape((vari_num,))),
                                       options=SLSQP_options)

                    xMat_potc = np.array(opt_res.x).reshape((1, -1))
                    objMat_potc_pred = obj_fcn_surrogate(xMat_potc)
                    conMat_potc, coneqMat_potc = con_fcn_surrogate(xMat_potc)
                    if not conMat_potc.shape[1]:
                        conMat_potc = conMat_potc / nomlz_value*con_max_list
                    if not coneqMat_potc.shape[1]:
                        coneqMat_potc = coneqMat_potc / nomlz_value*coneq_max_list

                    vioMat_potc = calViolation(
                        conMat_potc, coneqMat_potc, self.con_torl)

                    # updata into XMat_potc
                    XMat_potc = xMat_potc
                    ObjMat_potc = objMat_potc_pred
                    VioMat_potc = vioMat_potc

                detect_local_flag = False
            else:
                # updata XMat potc
                # for x_idx in range(1):
                for x_idx in range(len(XMat_potc)):
                    xMat_potc = XMat_potc[x_idx]

                    def obj_fcn_inf(xMat): return objFcnInf(
                        obj_fcn_surrogate, xMat)

                    def con_fcn_inf(xMat): return conFcnInf(
                        con_fcn_total, xMat)
                    opt_res = minimize(obj_fcn_inf, x_initial, method='SLSQP', constraints=NonlinearConstraint(con_fcn_inf, -np.inf, 0), bounds=Bounds(low_bou[0, :].reshape((vari_num,)), up_bou[0, :].reshape((vari_num,))),
                                       options=SLSQP_options)

                    xMat_potc = np.array(opt_res.x).reshape((1, -1))
                    objMat_potc_pred = np.array(opt_res.fun).reshape((1, -1))
                    conMat_potc, coneqMat_potc = con_fcn_surrogate(xMat_potc)
                    if not conMat_potc.shape[1]:
                        conMat_potc = conMat_potc / nomlz_value*con_max_list
                    if not coneqMat_potc.shape[1]:
                        coneqMat_potc = coneqMat_potc / nomlz_value*coneq_max_list
                    vioMat_potc = calViolation(
                        conMat_potc, coneqMat_potc, self.con_torl)

                    XMat_potc[x_idx, :] = xMat_potc
                    ObjMat_potc[x_idx, :] = objMat_potc_pred / \
                        nomlz_value * obj_max
                    VioMat_potc[x_idx, :] = vioMat_potc

                # merge XMat potc
                # Upward merge
                order_list = list(range(1, len(XMat_potc)-1))
                order_list = order_list[::-1]
                for x_idx in order_list:
                    xMat_potc = XMat_potc[x_idx]
                    # check if xMat_potc have existed
                    merge_flag = False
                    for x_check_idx in range(x_idx-1):
                        if np.sum(np.abs(XMat_potc[x_check_idx] - xMat_potc)) / vari_num < identiy_torl:
                            merge_flag = True
                            break
                    # updata into XMat_potc
                    if merge_flag:
                        XMat_potc = np.delete(XMat_potc, obj=x_idx, axis=0)
                        ObjMat_potc = np.delete(ObjMat_potc, obj=x_idx, axis=0)
                        VioMat_potc = np.delete(VioMat_potc, obj=x_idx, axis=0)

            # sort XMat potc by VioMat
            idx = np.argsort(VioMat_potc[:, 0])
            VioMat_potc = VioMat_potc[idx, :]
            ObjMat_potc = ObjMat_potc[idx, :]
            XMat_potc = XMat_potc[idx, :]

            # sort XMat potc by ObjMat
            if all(VioMat_potc == 0):
                idx = np.argsort(ObjMat_potc[:, 0])
                ObjMat_potc = ObjMat_potc[idx, :]
                VioMat_potc = VioMat_potc[idx, :]
                XMat_potc = XMat_potc[idx, :]
            else:
                index = [idx for idx in range(
                    len(XMat_potc)) if VioMat_potc[idx] == 0]
                if len(index) == 0:
                    idx = np.argsort(VioMat_potc[:, 0])
                    ObjMat_potc = ObjMat_potc[idx, :]
                    VioMat_potc = VioMat_potc[idx, :]
                    XMat_potc = XMat_potc[idx, :]
                else:
                    flag = max(index)

                    idx_feas = np.argsort(ObjMat_potc[:flag, 0])
                    idx_infeas = np.argsort(ObjMat_potc[flag:, 0])
                    idx = np.concatenate((idx_feas, idx_infeas + flag), axis=0)
                    ObjMat_potc = ObjMat_potc[idx, :]
                    VioMat_potc = VioMat_potc[idx, :]
                    XMat_potc = XMat_potc[idx, :]

            # step 4
            # select best potc point as xMat_infill
            xMat_infill = XMat_potc[0, :].reshape((1, -1))
            obj_infill_pred = ObjMat_potc[0, :].reshape((1, -1))
            # updata infill point

            xMat_infill, objMat_infill, conMat_infill, coneqMat_infill, vioMat_infill, ksMat_infill, repeat_idx, NFE_updata = data_lib.dataUpdata(
                xMat_infill, protect_range)
            NFE = NFE + NFE_updata
            if xMat_infill.shape[1] == 0:
                # process error
                xMat_infill = data_lib.X[repeat_idx]
                objMat_infill = data_lib.Obj[repeat_idx]
                if not ConMat.shape[1] == 0:
                    conMat_infill = data_lib.Con[repeat_idx]
                if not ConeqMat.shape[1] == 0:
                    coneqMat_infill = data_lib.Coneq[repeat_idx]
                if not VioMat.shape[1] == 0:
                    vioMat_infill = data_lib.Vio[repeat_idx]
            else:
                if not vioMat_infill.shape[1] == 0 and vioMat_infill > 0:
                    Bool_feas.append(False)
                else:
                    Bool_feas.append(True)
                Bool_conv.append(False)

            ObjMat_potc[0, :] = objMat_infill
            VioMat_potc[0, :] = vioMat_infill

            # find best result to record
            XMat, ObjMat, __, __, VioMat, __ = data_lib.dataLoad()
            XMat_unconv = np.array(
                [XMat[idx] for idx in range(len(XMat)) if not Bool_conv[idx]])
            ObjMat_unconv = np.array([ObjMat[idx]
                                      for idx in range(len(ObjMat)) if not Bool_conv[idx]])
            if not VioMat.shape[0] == 0:
                VioMat_unconv = np.array([VioMat[idx]
                                          for idx in range(len(VioMat)) if not Bool_conv[idx]])
            else:
                VioMat_unconv = np.array([[]])

            idx = [idx_unit for idx_unit in range(
                len(XMat_unconv)) if VioMat_unconv[idx_unit] == 0]
            if len(idx) == 0:
                min_idx = np.argmin(VioMat_unconv[:, 0])
                vio_best = VioMat_unconv[min_idx][0]
                obj_best = ObjMat_unconv[min_idx][0]
                x_best = list(XMat_unconv[min_idx])
            else:
                min_idx = np.argmin(ObjMat_unconv[idx, 0])
                obj_best = ObjMat_unconv[idx[min_idx]]
                vio_best = 0
                x_best = list(XMat_unconv[idx[min_idx]])

            if INFORMATION_FLAG:
                print('________________current best data________________')
                print('objMat:    %f    violation:    %f    NFE:    %-3d' %
                      (obj_best, vio_best, NFE))
                print(x_best)
                print('\n')

            result_obj_best.append(obj_best)
            iter = iter + 1
            # forced interrupt
            if iter > self.iter_max or NFE >= self.NFE_max:
                done = 1
            # convergence judgment
            if CONVERGENCE_JUDGMENT_FLAG:
                if (((iter > 2) and (np.abs((objMat_infill - objMat_infill_old) / objMat_infill_old) < self.torl)) and ((not vioMat_infill.shape[1] == 0 and vioMat_infill[0, 0] == 0) or vioMat_infill.shape[1] == 0)):
                    done = 1

            if not done:
                XMat, ObjMat, ConMat, ConeqMat, VioMat, KsMat = data_lib.dataLoad()
                # check if converage
                if (((iter > 2) and (np.abs((objMat_infill - objMat_infill_old) / objMat_infill_old) < self.torl)) and ((not len(vioMat_infill) == 0 and vioMat_infill == 0) or len(vioMat_infill) == 0)):
                    conv_con_GPC_fcn = None
                    # resample LHD

                    # step 6.1
                    # updata local best data
                    XMat_local_best = np.append(
                        XMat_local_best, xMat_infill).reshape((-1, vari_num))
                    ObjMat_local_best = np.append(
                        ObjMat_local_best, objMat_infill).reshape((-1, 1))

                    XMat_potc = XMat_potc[1:, :]
                    ObjMat_potc = ObjMat_potc[1:, :]
                    VioMat_potc = VioMat_potc[1:, :]
                    if len(XMat_potc) == 0:
                        detect_local_flag = True

                    # step 6.2
                    # detech converage
                    for x_idx in range(len(XMat)):
                        if not Bool_conv[x_idx]:
                            def obj_fcn_inf(xMat): return objFcnInf(
                                obj_fcn_surrogate, xMat)
                            def con_fcn_inf(xMat): return conFcnInf(
                                con_fcn_total, xMat)
                            opt_res = minimize(obj_fcn_inf, XMat[x_idx], method='SLSQP', constraints=NonlinearConstraint(con_fcn_inf, -np.inf, 0), bounds=Bounds(low_bou[0, :].reshape((vari_num,)), up_bou[0, :].reshape((vari_num,))),
                                               options=SLSQP_options)
                            x_single_pred = opt_res.x
                            converage_flag = False
                            for x_check_idx in range(len(XMat_local_best)):
                                if np.sum(np.abs(XMat_local_best[x_check_idx] - x_single_pred)) / vari_num < identiy_torl:
                                    converage_flag = True
                                    break
                            if converage_flag:
                                # if converage to local minimum, set to infeasible
                                Bool_conv[x_idx] = True

                    # step 6.3
                    # use GPC to limit do not converage to exist local best
                    if not np.all(Bool_conv):
                        Class = - 1 * np.ones((XMat.shape[0], 1))
                        Class[Bool_conv] = 1
                        model_GPC = GPC(XMat, Class, hyp)
                        model_GPC.train()

                        def conv_con_GPC_fcn(
                            xMat): return conGPCFcn(xMat, model_GPC)
                    else:
                        conv_con_GPC_fcn = None

                    # step 6.4
                    # resample latin hypercubic and updata into data library
                    try:
                        XMat_add, _, _ = getLatinHypercube(min(
                            sample_num_restart, self.NFE_max - NFE), vari_num, low_bou, up_bou, XMat, conv_con_GPC_fcn)
                    except:
                        XMat_add = lhs(vari_num, sample_num_initial) * \
                            (up_bou - low_bou) + low_bou

                    conv_con_GPC_fcn = None

                    # end select point RS
                else:
                    # step 5.1
                    # check if improve
                    improve = False
                    if len(repeat_idx) == 0:
                        VioMat_comp = np.array([VioMat[idx] for idx in range(
                            len(Bool_conv)-1) if (not Bool_conv[idx]) and Bool_feas[idx]])
                        ObjMat_comp = np.array([ObjMat[idx] for idx in range(
                            len(Bool_conv)-1) if (not Bool_conv[idx]) and Bool_feas[idx]])

                        
                        Bool_comp = np.logical_and(Bool_conv, Bool_feas)
                        if expensive_con_flag:
                            min_vio = min(VioMat_unconv[:, 0])
                            if not np.any(Bool_comp):
                                min_obj = None
                            else:
                                min_obj = min(ObjMat_comp[:, 0])

                            # if all point is infeasible,violation of point infilled is
                            # less than min violation of all point means improve.if
                            # feasible point exist,objMat of point infilled is less than min
                            # objMat means improve
                            if vioMat_infill[0, 0] == 0 or vioMat_infill < min_vio:
                                if min_obj is not None:
                                    if objMat_infill < min_obj:
                                        # improve, continue local search
                                        improve = True
                                else:
                                    # improve, continue local search
                                    improve = True
                        else:
                            ObjMat_comp = [ObjMat[idx] for idx in range(
                                len(Bool_conv-1)) if not Bool_conv[idx]]
                            min_obj = min(ObjMat_comp)
                            # objMat of point infilled is less than min objMat means improve
                            if objMat_infill < min_obj:
                                # imporve, continue local search
                                improve = True

                        # end check

                    # step 5.2
                    # if objMat no improve, use GPC to identify interest area
                    # than, imporve interest area surrogate quality
                    if not improve:
                        # construct GPC
                        train_num = min(
                            data_lib.X.shape[0], 11 * vari_num - 1 + 25)
                        pred_fcn_GPC, model_GPC, x_pareto_center, hyp = trainFilter(
                            data_lib, xMat_infill, hyp, train_num, expensive_con_flag, Bool_conv)

                        # step 5.3
                        # identify interest area
                        def con_GPC_fcn(xMat): return conGPCFcn(
                            xMat, model_GPC)

                        def obj_fcn_inf(xMat): return objFcnInf(
                            con_GPC_fcn, xMat)
                        def con_fcn_inf(xMat): return conFcnInf(
                            cheapcon_fcn, xMat)
                        opti_res = minimize(obj_fcn_inf, x_pareto_center.reshape((vari_num,)), method='SLSQP', constraints=NonlinearConstraint(con_fcn_inf, -np.inf, 0), bounds=Bounds(low_bou[0, :].reshape((vari_num,)), up_bou[0, :].reshape((vari_num,))),
                                            options=SLSQP_options)
                        center_point = opti_res.x

                        bou_interest = np.abs(center_point - xMat_infill)
                        bou_interest = np.maximum(np.multiply(
                            min_bou_interest, (up_bou - low_bou)), bou_interest)
                        low_bou_interest = xMat_infill - bou_interest
                        up_bou_interest = xMat_infill + bou_interest
                        low_bou_interest = np.maximum(
                            low_bou_interest, low_bou)
                        up_bou_interest = np.minimum(up_bou_interest, up_bou)

                        # generate trial point
                        trial_point = np.tile(xMat_infill, (trial_num, 1))
                        for vari_idx in range(vari_num):
                            trial_point[:, vari_idx] = trial_point[:, vari_idx] + np.random.normal(
                                loc=0, scale=bou_interest[0][vari_idx], size=(trial_num,))

                        trial_point = np.maximum(trial_point, low_bou)
                        trial_point = np.minimum(trial_point, up_bou)

                        Bool_negetive = pred_fcn_GPC(trial_point) == - 1
                        if np.sum(Bool_negetive) < sample_num_add:
                            value = con_GPC_fcn(trial_point)
                            thres = np.quantile(value, 0.25)
                            Bool_negetive = value < thres

                        trial_point = trial_point[Bool_negetive[:, 0], :]

                        # step 5.4
                        # select point
                        select_num = min(sample_num_add, self.NFE_max - NFE)
                        if trial_point.shape[0] <= select_num:
                            XMat_add = trial_point
                        else:
                            XMat, __, __, __, __, __ = data_lib.dataLoad()
                            max_dist = 0.0
                            iter_select = 1
                            while iter_select < 100:

                                select_idx = random.sample(
                                    range(trial_point.shape[0]), select_num)
                                dist = calMinDistanceIter(
                                    trial_point[select_idx, :], XMat)
                                if max_dist < dist:
                                    XMat_add = trial_point[select_idx, :]
                                    max_dist = dist
                                iter_select = iter_select + 1
                    else:
                        XMat_add = np.array([[]])
                    # end select point GPC

                # end select point
                # updata tp data library
                if XMat_add.shape[1] != 0:
                    XMat_add, __, __, __, VioMat_add, __, __, NFE_updata = data_lib.dataUpdata(
                        XMat_add, protect_range)
                    NFE = NFE + NFE_updata
                    if NFE_updata > 0:
                        Bool_feas_new = [VioMat_add[idx][0] ==
                                        0 for idx in range(len(VioMat_add))]
                        Bool_feas = Bool_feas+Bool_feas_new
                        Bool_conv = Bool_conv+[False]*len(XMat_add)

                # forced interrupt
                if iter > self.iter_max or NFE >= self.NFE_max:
                    done = 1

            # end sample
            obj_best_old = obj_best
            objMat_infill_old = objMat_infill
            con_infill_old = conMat_infill
            coneq_infill_old = coneqMat_infill
            vio_infill_old = vioMat_infill

            XMat_local_best = np.array([[]])
            ObjMat_local_best = np.array([[]])
            XMat_potc = np.array([[]])
            ObjMat_potc = np.array([[]])
            VioMat_potc = np.array([[]])
            detect_local_flag = True

            io.savemat('iter.mat', {
                'data_lib': data_lib, 'X_local_best': XMat_local_best, 'Obj_local_best': ObjMat_local_best, 'X_potential': XMat_potc, 'Obj_potential': ObjMat_potc, 'Vio_potential': VioMat_potc, 'detect_local_flag': detect_local_flag})

        # end algoritm
        # find best result to record
        x_best = data_lib.X[data_lib.result_best_idx[-1]]
        obj_best = data_lib.Obj[data_lib.result_best_idx[-1]]

        output = {}
        output['result_x_best'] = result_x_best
        output['result_obj_best'] = result_obj_best
        output['x_local_best'] = XMat_local_best
        output['obj_local_best'] = ObjMat_local_best
        output['data_lib'] = data_lib

        io.savemat('result.mat', {
            'data_lib': data_lib, 'X_local_best': XMat_local_best, 'Obj_local_best': ObjMat_local_best})

        return x_best, obj_best, NFE, output


class RBF_GPC_RS_MF():
    def __init__(self, NFE_max=None, iter_max=None, torl=1e-6, con_torl=1e-3, data_lib_HF:dict=None, data_lib_LF:dict=None):
        self.NFE_max = NFE_max
        self.iter_max = iter_max
        self.torl = torl
        self.con_torl = con_torl

        self.data_lib_HF = data_lib_HF
        self.data_lib_LF = data_lib_LF

        pass

    def optimize(self, MF_model, Cost, Ratio, vari_num, low_bou, up_bou, cheapcon_fcn=None):
        # RBF-AL-KTR optimization algorithm
        # potc: potential

        # Copyright 2023 4 Adel
        low_bou = np.array(low_bou).reshape((1, vari_num))
        up_bou = np.array(up_bou).reshape((1, vari_num))

        DRAW_FIGURE_FLAG = 0
        INFORMATION_FLAG = 1
        CONVERGENCE_JUDGMENT_FLAG = 0
        WRIRE_FILE_FLAG = 0

        fidelity_num = len(MF_model)
        HF_model = MF_model[0]
        LF_model = MF_model[1]
        cost_HF = Cost[0]
        cost_LF = Cost[1]
        ratio_HF = Ratio[0]
        ratio_LF = Ratio[1]

        # hyper parameter
        sample_num_initial = 4 + 2 * vari_num
        sample_num_restart = sample_num_initial
        sample_num_add = int(np.ceil(np.log(sample_num_initial)/2))

        min_bou_interest = 0.001
        trial_num = np.minimum(100 * vari_num, 100)
        nomlz_value = 10
        protect_range = 1e-6
        identiy_torl = 0.001

        # NFE and iter setting
        if self.NFE_max is None:
            self.NFE_max = 10 + 10 * vari_num

        if self.iter_max is None:
            self.iter_max = 20 + 20 * vari_num

        done = 0
        NFE = 0
        iter = 0
        NFE_list = [0, 0]
        # step 1
        # generate initial data library
        
        if self.data_lib_LF is not None:
            data_lib_LF = DataLibrary(model=LF_model, vari_num=vari_num, low_bou=low_bou, up_bou=up_bou, con_torl=self.con_torl,
                                      write_file_flag=WRIRE_FILE_FLAG)
            data_lib_LF.X = self.data_lib_LF['X']
            data_lib_LF.Obj = self.data_lib_LF['Obj']
            data_lib_LF.Con = self.data_lib_LF['Con']
            data_lib_LF.Coneq = self.data_lib_LF['Coneq']
            data_lib_LF.Ks = self.data_lib_LF['Ks']
            data_lib_LF.Vio = self.data_lib_LF['Vio']
            if (data_lib_LF.Con.shape[1] != 0) or (data_lib_LF.Coneq.shape[1] != 0):
                expensive_con_flag = True
            else:
                expensive_con_flag = False
        else:
            sample_num_initial_LF = int(np.ceil(sample_num_initial*ratio_LF))
            XMat_LF_updata = lhs(
                vari_num, sample_num_initial_LF)*(up_bou - low_bou) + low_bou
            
            data_lib_LF = DataLibrary(model=LF_model, vari_num=vari_num, low_bou=low_bou, up_bou=up_bou, con_torl=self.con_torl,
                                      write_file_flag=WRIRE_FILE_FLAG)
            # detech expensive constraints and initialize data library
            xMat_new, objMat_new, conMat_new, coneqMat_new, vioMat_new, ksMat_new, repeat_idx, NFE_updata = data_lib_LF.dataUpdata(
                [XMat_LF_updata[0]], 0)
            NFE = NFE + NFE_updata * cost_LF
            NFE_list[1] = NFE_list[1] + NFE_updata
            if not data_lib_LF.Vio.shape[1] == 0:
                expensive_con_flag = True
            else:
                expensive_con_flag = False

            # updata data library by XMat
            xMat_new, objMat_new, conMat_new, coneqMat_new, vioMat_new, ksMat_new, repeat_idx, NFE_updata = data_lib_LF.dataUpdata(
                XMat_LF_updata[1:], 0)
            NFE = NFE + NFE_updata * cost_LF
            NFE_list[1] = NFE_list[1] + NFE_updata
     
            
        if self.data_lib_HF is not None:
            data_lib_HF = DataLibrary(model=HF_model, vari_num=vari_num, low_bou=low_bou, up_bou=up_bou, con_torl=self.con_torl,
                                      write_file_flag=WRIRE_FILE_FLAG)
            data_lib_HF.X = self.data_lib_HF['X']
            data_lib_HF.Obj = self.data_lib_HF['Obj']
            data_lib_HF.Con = self.data_lib_HF['Con']
            data_lib_HF.Coneq = self.data_lib_HF['Coneq']
            data_lib_HF.Ks = self.data_lib_HF['Ks']
            data_lib_HF.Vio = self.data_lib_HF['Vio']
        else:
            sample_num_initial_HF = int(np.ceil(sample_num_initial*ratio_HF))

            XMat_HF_updata,_,_ = getNestedHypercube(
                data_lib_LF.X, sample_num_initial_HF * ratio_HF, vari_num, low_bou, up_bou)

            data_lib_HF = DataLibrary(model=HF_model, vari_num=vari_num, low_bou=low_bou, up_bou=up_bou, con_torl=self.con_torl,
                            write_file_flag=WRIRE_FILE_FLAG)       
            xMat_new, objMat_new, conMat_new, coneqMat_new, vioMat_new, ksMat_new, repeat_idx, NFE_updata = data_lib_HF.dataUpdata(
                XMat_HF_updata, 0)
            NFE = NFE + NFE_updata * cost_HF
            NFE_list[0] = NFE_list[0] + NFE_updata


        # find fesiable data in current data library
        if expensive_con_flag:
            Bool_feas_HF = [(data_lib_HF.Vio[idx, 0] == 0)
                            for idx in range(data_lib_HF.X.shape[0])]
            Bool_feas_LF = [(data_lib_LF.Vio[idx, 0] == 0)
                            for idx in range(data_lib_LF.X.shape[0])]

        Bool_conv_HF = [False]*data_lib_HF.X.shape[0]
        Bool_conv_LF = [False]*data_lib_LF.X.shape[0]

        hyp_MF = {'mean': [0], 'cov': np.array([0, 0, 0, 0, 0])}
        hyp_SF = {'mean': [0], 'cov': np.array([0, 0])}

        XMat_local_best = np.array([[]])
        ObjMat_local_best = np.array([[]])
        XMat_potc = np.array([[]])
        ObjMat_potc = np.array([[]])
        VioMat_potc = np.array([[]])
        detect_local_flag = True

        SLSQP_options = {'ftol': 1e-6, 'disp': False, 'maxiter': 100}
        result_x_best = list()
        result_obj_best = list()
        add_LF_flag = True
        iter = iter + 1
        while not done:
            # step 2
            # nomalization all data by max objMat and to create surrogate model
            XMat_HF, ObjMat_HF, Con_HF, Coneq_HF, VioMat_HF, KsMat_HF = data_lib_HF.dataLoad()
            XMat_LF, ObjMat_LF, Con_LF, Coneq_LF, VioMat_LF, KsMat_LF = data_lib_LF.dataLoad()
            # get local infill point, construct surrogate model
            XMat_MF = [XMat_HF, XMat_LF]
            ObjMat_MF = [ObjMat_HF, ObjMat_LF]
            if Con_HF.shape[1] != 0:
                ConMat_MF = [Con_HF, Con_LF]
            else:
                ConMat_MF = []
            if Coneq_HF.shape[1] != 0:
                ConeqMat_MF = [Coneq_HF, Coneq_LF]
            else:
                ConeqMat_MF = []
            if VioMat_HF.shape[1] != 0:
                VioMat_MF = [VioMat_HF, VioMat_LF]
            else:
                VioMat_MF = []
            if KsMat_HF.shape[1] != 0:
                KsMat_MF = [KsMat_HF, KsMat_LF]
            else:
                KsMat_MF = []

            XMat_model_MF, ObjMat_model_MF, ConMat_model_MF, ConeqMat_model_MF, VioMat_model_MF, KsMat_model_MF, obj_max, con_max_list, coneq_max_list, vio_max, ks_max = getModelDataMF(
                fidelity_num, XMat_MF, ObjMat_MF, ConMat_MF, ConeqMat_MF, VioMat_MF, KsMat_MF, nomlz_value)
            obj_fcn_surrogate, con_fcn_surrogate, output_model = getSurrogateFcnMF(
                XMat_model_MF, ObjMat_model_MF, ConMat_model_MF, ConeqMat_model_MF)

            model_obj = output_model['model_obj']
            model_con_list = output_model['model_con_list']
            model_coneq_list = output_model['model_coneq_list']
            type_obj = output_model['type_obj']
            type_con_list = output_model['type_con_list']
            type_coneq_list = output_model['type_coneq_list']

            if (con_fcn_surrogate is not None) or (cheapcon_fcn is not None):
                def con_fcn_total(xMat): return totalConFcn(
                    xMat, con_fcn_surrogate, cheapcon_fcn, None)
            else:
                con_fcn_total = None

            # check if all model is SF
            all_SF_flag = True
            if type_obj == 'MF':
                all_SF_flag = False
            else:
                if not len(model_con_list) == 0:
                    for con_idx in range(len(model_con_list)):
                        if type_con_list[con_idx] == 'MF':
                            all_SF_flag = False
                            break
                if not len(model_coneq_list) == 0:
                    for coneq_idx in range(len(model_coneq_list)):
                        if type_coneq_list[coneq_idx] == 'MF':
                            all_SF_flag = False
                            break
            if all_SF_flag:
                add_LF_flag = False

            # step 3
            if detect_local_flag:
                # detech potc local best point
                # for x_idx in range(1):
                for x_idx in range(len(XMat_HF)):
                    x_initial = XMat_HF[x_idx]

                    def obj_fcn_inf(xMat): return objFcnInf(
                        obj_fcn_surrogate, xMat)

                    def con_fcn_inf(xMat): return conFcnInf(
                        con_fcn_total, xMat)
                    opt_res = minimize(obj_fcn_inf, x_initial, method='SLSQP', constraints=NonlinearConstraint(con_fcn_inf, -np.inf, 0), bounds=Bounds(low_bou[0, :].reshape((vari_num,)), up_bou[0, :].reshape((vari_num,))),
                                       options=SLSQP_options)
                    xMat_potc = np.array(opt_res.x).reshape((1, -1))
                    objMat_potc_pred = np.array(opt_res.fun).reshape((1, -1))

                    conMat_potc, coneqMat_potc = con_fcn_surrogate(xMat_potc)
                    if not conMat_potc.shape[1]:
                        conMat_potc = conMat_potc / nomlz_value*con_max_list
                    if not coneqMat_potc.shape[1]:
                        coneqMat_potc = coneqMat_potc / nomlz_value*coneq_max_list
                    vioMat_potc = calViolation(
                        conMat_potc, coneqMat_potc, self.con_torl)

                    if opt_res.success:
                        # check if xMat_potc have existed
                        add_flag = True
                        if XMat_potc.shape[1] != 0:
                            for x_check_idx in range(len(XMat_potc)):
                                if np.sum(np.abs(XMat_potc[x_check_idx, :] - xMat_potc)) / vari_num < identiy_torl:
                                    add_flag = False
                                    break
                        if XMat_local_best.shape[1] != 0:
                            for x_check_idx in range(len(XMat_local_best)):
                                if np.sum(np.abs(XMat_local_best[x_check_idx, :] - xMat_potc)) / vari_num < identiy_torl:
                                    add_flag = False
                                    break

                        # updata into XMat_potc
                        if add_flag:
                            if XMat_potc.shape[1] == 0:
                                XMat_potc = xMat_potc
                                ObjMat_potc = objMat_potc_pred
                                VioMat_potc = vioMat_potc
                            else:
                                XMat_potc = np.concatenate(
                                    (XMat_potc, xMat_potc), axis=0)
                                ObjMat_potc = np.concatenate(
                                    (ObjMat_potc, objMat_potc_pred/nomlz_value*obj_max), axis=0)
                                VioMat_potc = np.concatenate(
                                    (VioMat_potc, vioMat_potc), axis=0)

                # if XMat_potc is empty, try to use KS surrogate as xMat potc
                if XMat_potc.shape[1] == 0:
                    model_ks, model_type = getBestModel(
                        XMat_model_MF, KsMat_model_MF)

                    def ks_fcn_surrogate(
                        X_pred): return model_ks.predict(X_pred)

                    x_idx = np.argmin(VioMat_HF)
                    x_initial = XMat_HF[x_idx]
                    def obj_fcn_inf(xMat): return objFcnInf(
                        ks_fcn_surrogate, xMat)

                    def con_fcn_inf(xMat): return conFcnInf(cheapcon_fcn, xMat)
                    opt_res = minimize(con_fcn_inf, x_initial, method='SLSQP', constraints=NonlinearConstraint(obj_fcn_inf, -np.inf, 0), bounds=Bounds(low_bou[0, :].reshape((vari_num,)), up_bou[0, :].reshape((vari_num,))),
                                       options=SLSQP_options)

                    xMat_potc = np.array(opt_res.x).reshape((1, -1))
                    objMat_potc_pred = obj_fcn_surrogate(xMat_potc)
                    conMat_potc, coneqMat_potc = con_fcn_surrogate(xMat_potc)
                    if not conMat_potc.shape[1]:
                        conMat_potc = conMat_potc / nomlz_value*con_max_list
                    if not coneqMat_potc.shape[1]:
                        coneqMat_potc = coneqMat_potc / nomlz_value*coneq_max_list

                    vioMat_potc = calViolation(
                        conMat_potc, coneqMat_potc, self.con_torl)

                    # updata into XMat_potc
                    XMat_potc = xMat_potc
                    ObjMat_potc = objMat_potc_pred
                    VioMat_potc = vioMat_potc

                detect_local_flag = False
            else:
                # updata XMat potc
                # for x_idx in range(1):
                for x_idx in range(len(XMat_potc)):
                    xMat_potc = XMat_potc[x_idx]

                    def obj_fcn_inf(xMat): return objFcnInf(
                        obj_fcn_surrogate, xMat)

                    def con_fcn_inf(xMat): return conFcnInf(
                        con_fcn_total, xMat)
                    opt_res = minimize(obj_fcn_inf, x_initial, method='SLSQP', constraints=NonlinearConstraint(con_fcn_inf, -np.inf, 0), bounds=Bounds(low_bou[0, :].reshape((vari_num,)), up_bou[0, :].reshape((vari_num,))),
                                       options=SLSQP_options)

                    xMat_potc = np.array(opt_res.x).reshape((1, -1))
                    objMat_potc_pred = np.array(opt_res.fun).reshape((1, -1))
                    conMat_potc, coneqMat_potc = con_fcn_surrogate(xMat_potc)
                    if not conMat_potc.shape[1]:
                        conMat_potc = conMat_potc / nomlz_value*con_max_list
                    if not coneqMat_potc.shape[1]:
                        coneqMat_potc = coneqMat_potc / nomlz_value*coneq_max_list
                    vioMat_potc = calViolation(
                        conMat_potc, coneqMat_potc, self.con_torl)

                    XMat_potc[x_idx, :] = xMat_potc
                    ObjMat_potc[x_idx, :] = objMat_potc_pred / \
                        nomlz_value * obj_max
                    VioMat_potc[x_idx, :] = vioMat_potc

                # merge XMat potc
                # Upward merge
                order_list = list(range(1, len(XMat_potc)-1))
                order_list = order_list[::-1]
                for x_idx in order_list:
                    xMat_potc = XMat_potc[x_idx]
                    # check if xMat_potc have existed
                    merge_flag = False
                    for x_check_idx in range(x_idx-1):
                        if np.sum(np.abs(XMat_potc[x_check_idx] - xMat_potc)) / vari_num < identiy_torl:
                            merge_flag = True
                            break
                    # updata into XMat_potc
                    if merge_flag:
                        XMat_potc = np.delete(XMat_potc, obj=x_idx, axis=0)
                        ObjMat_potc = np.delete(ObjMat_potc, obj=x_idx, axis=0)
                        VioMat_potc = np.delete(VioMat_potc, obj=x_idx, axis=0)

            # sort XMat potc by VioMat
            idx = np.argsort(VioMat_potc[:, 0])
            VioMat_potc = VioMat_potc[idx, :]
            ObjMat_potc = ObjMat_potc[idx, :]
            XMat_potc = XMat_potc[idx, :]

            # sort XMat potc by ObjMat
            if all(VioMat_potc == 0):
                idx = np.argsort(ObjMat_potc[:, 0])
                ObjMat_potc = ObjMat_potc[idx, :]
                VioMat_potc = VioMat_potc[idx, :]
                XMat_potc = XMat_potc[idx, :]
            else:
                index = [idx for idx in range(
                    len(XMat_potc)) if VioMat_potc[idx] == 0]
                if len(index) == 0:
                    idx = np.argsort(VioMat_potc[:, 0])
                    ObjMat_potc = ObjMat_potc[idx, :]
                    VioMat_potc = VioMat_potc[idx, :]
                    XMat_potc = XMat_potc[idx, :]
                else:
                    flag = max(index)

                    idx_feas = np.argsort(ObjMat_potc[:flag, 0])
                    idx_infeas = np.argsort(ObjMat_potc[flag:, 0])
                    idx = np.concatenate((idx_feas, idx_infeas + flag), axis=0)
                    ObjMat_potc = ObjMat_potc[idx, :]
                    VioMat_potc = VioMat_potc[idx, :]
                    XMat_potc = XMat_potc[idx, :]

            # step 4
            # select best potc point as xMat_infill
            xMat_infill = XMat_potc[0, :].reshape((1, -1))
            obj_infill_pred = ObjMat_potc[0, :].reshape((1, -1))

            # updata infill point
            xMat_infill, objMat_infill, conMat_infill, coneqMat_infill, vioMat_infill, ksMat_infill, repeat_idx, NFE_updata = data_lib_HF.dataUpdata(
                xMat_infill, protect_range)
            NFE = NFE + NFE_updata * cost_HF
            NFE_list[0] = NFE_list[0] + NFE_updata

            if xMat_infill.shape[1] == 0:
                # process error
                xMat_infill = data_lib_HF.X[repeat_idx]
                objMat_infill = data_lib_HF.Obj[repeat_idx]
                if not ConMat_HF.shape[1] == 0:
                    conMat_infill = data_lib_HF.Con[repeat_idx]
                if not ConeqMat_HF.shape[1] == 0:
                    coneqMat_infill = data_lib_HF.Coneq[repeat_idx]
                if not VioMat_HF.shape[1] == 0:
                    vioMat_infill = data_lib_HF.Vio[repeat_idx]
            else:
                if not vioMat_infill.shape[1] == 0 and vioMat_infill > 0:
                    Bool_feas_HF.append(False)
                else:
                    Bool_feas_HF.append(True)
                Bool_conv_HF.append(False)

            ObjMat_potc[0, :] = objMat_infill
            VioMat_potc[0, :] = vioMat_infill

            # find best result to record
            XMat_HF, ObjMat_HF, __, __, VioMat_HF, __ = data_lib_HF.dataLoad()
            XMat_HF_unconv = np.array(
                [XMat_HF[idx] for idx in range(len(XMat_HF)) if not Bool_conv_HF[idx]])
            ObjMat_HF_unconv = np.array([ObjMat_HF[idx]
                                         for idx in range(len(ObjMat_HF)) if not Bool_conv_HF[idx]])
            if not VioMat_HF.shape[0] == 0:
                VioMat_HF_unconv = np.array([VioMat_HF[idx]
                                             for idx in range(len(VioMat_HF)) if not Bool_conv_HF[idx]])
            else:
                VioMat_HF_unconv = np.array([[]])

            idx = [idx_unit for idx_unit in range(
                len(XMat_HF_unconv)) if VioMat_HF_unconv[idx_unit] == 0]
            if len(idx) == 0:
                min_idx = np.argmin(VioMat_HF_unconv[:, 0])
                vio_best = VioMat_HF_unconv[min_idx][0]
                obj_best = ObjMat_HF_unconv[min_idx][0]
                x_best = list(XMat_HF_unconv[min_idx])
            else:
                min_idx = np.argmin(ObjMat_HF_unconv[idx, 0])
                obj_best = ObjMat_HF_unconv[idx[min_idx]]
                vio_best = 0
                x_best = list(XMat_HF_unconv[idx[min_idx]])

            if INFORMATION_FLAG:
                print('________________current best data________________')
                print('objMat:    %f    violation:    %f    NFE:    %-3d' %
                      (obj_best, vio_best, NFE))
                print(x_best)
                print('\n')

            result_obj_best.append(obj_best)
            iter = iter + 1
            # forced interrupt
            if iter > self.iter_max or NFE >= self.NFE_max:
                done = 1
            # convergence judgment
            if CONVERGENCE_JUDGMENT_FLAG:
                if (((iter > 2) and (np.abs((objMat_infill - objMat_infill_old) / objMat_infill_old) < self.torl)) and ((not vioMat_infill.shape[1] == 0 and vioMat_infill[0, 0] == 0) or vioMat_infill.shape[1] == 0)):
                    done = 1

            if not done:
                XMat_HF, ObjMat_HF, ConMat_HF, ConeqMat_HF, VioMat_HF, KsMat_HF = data_lib_HF.dataLoad()
                XMat_LF, ObjMat_LF, ConMat_LF, ConeqMat_LF, VioMat_LF, KsMat_LF = data_lib_LF.dataLoad()
                # check if converage
                if (((iter > 2) and (np.abs((objMat_infill - objMat_infill_old) / objMat_infill_old) < self.torl)) and ((not len(vioMat_infill) == 0 and vioMat_infill == 0) or len(vioMat_infill) == 0)):
                    conv_con_GPC_fcn = None

                    # step 6
                    # resample LHD
                    XMat_local_best = np.append(
                        XMat_local_best, xMat_infill).reshape((-1, vari_num))
                    ObjMat_local_best = np.append(
                        ObjMat_local_best, objMat_infill).reshape((-1, 1))

                    XMat_potc = XMat_potc[1:, :]
                    ObjMat_potc = ObjMat_potc[1:, :]
                    VioMat_potc = VioMat_potc[1:, :]
                    if len(XMat_potc) == 0:
                        detect_local_flag = True

                    # step 6.1
                    # detech converage
                    # HF
                    for x_idx in range(len(XMat_HF)):
                        if not Bool_conv_HF[x_idx]:
                            def obj_fcn_inf(xMat): return objFcnInf(
                                obj_fcn_surrogate, xMat)

                            def con_fcn_inf(xMat): return conFcnInf(
                                con_fcn_total, xMat)
                            opt_res = minimize(obj_fcn_inf, XMat_HF[x_idx], method='SLSQP', constraints=NonlinearConstraint(con_fcn_inf, -np.inf, 0), bounds=Bounds(low_bou[0, :].reshape((vari_num,)), up_bou[0, :].reshape((vari_num,))),
                                               options=SLSQP_options)
                            x_single_pred = opt_res.x
                            converage_flag = False
                            for x_check_idx in range(len(XMat_local_best)):
                                if np.sum(np.abs(XMat_local_best[x_check_idx] - x_single_pred)) / vari_num < identiy_torl:
                                    converage_flag = True
                                    break
                            if converage_flag:
                                # if converage to local minimum, set to infeasible
                                Bool_conv_HF[x_idx] = True

                    # LF
                    for x_idx in range(len(XMat_LF)):
                        if not Bool_conv_LF[x_idx]:
                            def obj_fcn_inf(xMat): return objFcnInf(
                                obj_fcn_surrogate, xMat)

                            def con_fcn_inf(xMat): return conFcnInf(
                                con_fcn_total, xMat)
                            opt_res = minimize(obj_fcn_inf, XMat_LF[x_idx], method='SLSQP', constraints=NonlinearConstraint(con_fcn_inf, -np.inf, 0), bounds=Bounds(low_bou[0, :].reshape((vari_num,)), up_bou[0, :].reshape((vari_num,))),
                                               options=SLSQP_options)
                            x_single_pred = opt_res.x
                            converage_flag = False
                            for x_check_idx in range(len(XMat_local_best)):
                                if np.sum(np.abs(XMat_local_best[x_check_idx] - x_single_pred)) / vari_num < identiy_torl:
                                    converage_flag = True
                                    break
                            if converage_flag:
                                # if converage to local minimum, set to infeasible
                                Bool_conv_LF[x_idx] = True

                    # step 6.2
                    # use GPC to limit do not converage to exist local best
                    if not np.all(Bool_conv_HF):
                        Y_HF = - 1 * np.ones((XMat_HF.shape[0], 1))
                        Y_HF[Bool_conv_HF] = 1
                        Y_LF = - 1 * np.ones((XMat_LF.shape[0], 1))
                        Y_LF[Bool_conv_LF] = 1

                        if add_LF_flag:
                            model_GPCMF = GPCMF(
                                XMat_HF, Y_HF, XMat_LF, Y_LF, hyp_MF)
                            model_GPCMF.train()

                            def conv_con_GPC_fcn(
                                xMat): return conGPCFcn(xMat, model_GPCMF)
                        else:
                            model_GPCSF = GPC(XMat_HF, Y_HF, hyp_SF)
                            model_GPCSF.train()

                            def conv_con_GPC_fcn(
                                xMat): return conGPCFcn(xMat, model_GPCSF)
                    else:
                        conv_con_GPC_fcn = None

                    # step 6.3
                    # resample latin hypercubic and updata into data library
                    usable_NFE = self.NFE_max - NFE
                    sample_num = min(
                        sample_num_restart, int(np.ceil(usable_NFE / sum(Ratio))))
                    sample_num_HF = int(np.ceil(sample_num*ratio_HF))
                    sample_num_LF = int(np.ceil(sample_num*ratio_LF))

                    if add_LF_flag:
                        try:
                            XMat_add_LF, _, _ = getLatinHypercube(min(
                                sample_num_LF, self.NFE_max - NFE), vari_num, low_bou, up_bou, XMat_LF, conv_con_GPC_fcn)
                        except:
                            XMat_add_LF = lhs(
                                vari_num, sample_num_LF) * (up_bou - low_bou) + low_bou

                        XMat_add_HF,_,_ = getNestedHypercube(
                            XMat_add_LF, sample_num_HF, vari_num, low_bou, up_bou, XMat_HF)
                    else:
                        XMat_add_LF = np.array([[]])

                        try:
                            XMat_add_HF, _, _ = getLatinHypercube(min(
                                sample_num_HF, self.NFE_max - NFE), vari_num, low_bou, up_bou, XMat_HF, conv_con_GPC_fcn)
                        except:
                            XMat_add_HF = lhs(
                                vari_num, sample_num_HF) * (up_bou - low_bou) + low_bou

                    conv_con_GPC_fcn = None

                    # end select add point RS
                else:
                    # step 5.1
                    # check if improve
                    XMat_add_HF=np.array([[]])
                    XMat_add_LF=np.array([[]])
                    improve = False
                    if len(repeat_idx) == 0:
                        VioMat_HF_comp = np.array([VioMat_HF[idx] for idx in range(
                            len(Bool_conv_HF)-1) if (not Bool_conv_HF[idx]) and Bool_feas_HF[idx]])
                        ObjMat_HF_comp = np.array([ObjMat_HF[idx] for idx in range(
                            len(Bool_conv_HF)-1) if (not Bool_conv_HF[idx]) and Bool_feas_HF[idx]])

                        Bool_HF_comp = np.logical_and(Bool_conv_HF, Bool_feas_HF)
                        if expensive_con_flag:
                            min_vio = min(VioMat_HF_unconv[:, 0])
                            if not np.any(Bool_HF_comp):
                                min_obj = None
                            else:
                                min_obj = min(ObjMat_HF_comp[:, 0])

                            # if all point is infeasible,violation of point infilled is
                            # less than min violation of all point means improve.if
                            # feasible point exist,objMat of point infilled is less than min
                            # objMat means improve
                            if vioMat_infill[0, 0] == 0 or vioMat_infill < min_vio:
                                if min_obj is not None:
                                    if objMat_infill < min_obj:
                                        # improve, continue local search
                                        improve = True
                                else:
                                    # improve, continue local search
                                    improve = True
                        else:
                            ObjMat_HF_comp = [ObjMat_HF[idx] for idx in range(
                                len(Bool_conv_HF-1)) if not Bool_conv_HF[idx]]
                            min_obj = min(ObjMat_HF_comp)
                            # objMat of point infilled is less than min objMat means improve
                            if objMat_infill < min_obj:
                                # imporve, continue local search
                                improve = True

                    # step 5.2
                    # if objMat no improve, use GPC to identify interest area
                    # than, imporve interest area surrogate quality
                    if not improve:
                        # construct GPC
                        train_num = min(
                            data_lib_HF.X.shape[0], 11 * vari_num - 1 + 25)
                        if add_LF_flag:
                            pred_fcn_GPC, model_GPCMF, x_pareto_center, hyp_MF = trainFilterMF(
                                data_lib_HF, data_lib_LF, Ratio, xMat_infill, hyp_MF, train_num, expensive_con_flag, Bool_conv_HF, Bool_conv_LF)

                            def con_GPC_fcn(xMat): return conGPCFcn(
                                xMat, model_GPCMF)
                        else:
                            pred_fcn_GPC, model_GPC, x_pareto_center, hyp_SF = trainFilter(
                                data_lib_HF, xMat_infill, hyp_SF, train_num, expensive_con_flag, Bool_conv_HF)
                            def con_GPC_fcn(xMat): return conGPCFcn(
                                xMat, model_GPC)

                        # step 5.3
                        # identify interest area
                        def obj_fcn_inf(xMat): return objFcnInf(
                            con_GPC_fcn, xMat)
                        def con_fcn_inf(xMat): return conFcnInf(
                            cheapcon_fcn, xMat)
                        opti_res = minimize(obj_fcn_inf, x_pareto_center.reshape((vari_num,)), method='SLSQP', constraints=NonlinearConstraint(con_fcn_inf, -np.inf, 0), bounds=Bounds(low_bou[0, :].reshape((vari_num,)), up_bou[0, :].reshape((vari_num,))),
                                            options=SLSQP_options)
                        center_point = opti_res.x

                        bou_interest = np.abs(center_point - xMat_infill)
                        bou_interest = np.maximum(np.multiply(
                            min_bou_interest, (up_bou - low_bou)), bou_interest)
                        low_bou_interest = xMat_infill - bou_interest
                        up_bou_interest = xMat_infill + bou_interest
                        low_bou_interest = np.maximum(
                            low_bou_interest, low_bou)
                        up_bou_interest = np.minimum(up_bou_interest, up_bou)

                        # generate trial point
                        usable_NFE = self.NFE_max - NFE
                        sample_num = min(
                            sample_num_add, int(np.ceil(usable_NFE / sum(Ratio))))
                        sample_num_HF = int(np.ceil(sample_num*ratio_HF))
                        sample_num_LF = int(np.ceil(sample_num*ratio_LF))
                        
                        trial_point = np.tile(xMat_infill, (trial_num, 1))
                        for vari_idx in range(vari_num):
                            trial_point[:, vari_idx] = trial_point[:, vari_idx] + np.random.normal(
                                loc=0, scale=bou_interest[0][vari_idx], size=(trial_num,))

                        trial_point = np.maximum(trial_point, low_bou)
                        trial_point = np.minimum(trial_point, up_bou)

                        Bool_negetive = pred_fcn_GPC(trial_point) == - 1
                        if np.sum(Bool_negetive) < sample_num_LF:
                            value = con_GPC_fcn(trial_point)
                            thres = np.quantile(value, 0.25)
                            Bool_negetive = value < thres

                        trial_point = trial_point[Bool_negetive[:, 0], :]
                        
                        # step 5.4
                        # select point
                        if add_LF_flag:
                            if trial_point.shape[0] <= sample_num_LF:
                                XMat_add_LF = trial_point
                            else:
                                # sample LF first
                                max_dist = 0.0
                                iter_select = 1
                                while iter_select < 100:
                                    select_idx = random.sample(
                                        range(trial_point.shape[0]), sample_num_LF)
                                    dist = calMinDistanceIter(
                                        trial_point[select_idx, :], XMat_LF)
                                    if max_dist < dist:
                                        XMat_add_LF = trial_point[select_idx, :]
                                        max_dist = dist
                                    iter_select = iter_select + 1

                                # sample HF from LF
                                max_dist = 0.0
                                iter_select = 1
                                while iter_select < 100:
                                    select_idx = random.sample(
                                        range(XMat_add_LF.shape[0]), sample_num_HF)
                                    dist = calMinDistanceIter(
                                        XMat_add_LF[select_idx, :], XMat_HF)
                                    if max_dist < dist:
                                        XMat_add_HF = XMat_add_LF[select_idx, :]
                                        max_dist = dist
                                    iter_select = iter_select + 1

                        else:
                            XMat_add_LF=np.array([[]])
                            if trial_point.shape[0] <= sample_num_HF:
                                XMat_add_HF = trial_point
                            else:
                                # sample HF first
                                max_dist = 0.0
                                iter_select = 1
                                while iter_select < 100:
                                    select_idx = random.sample(
                                        range(trial_point.shape[0]), sample_num_HF)
                                    dist = calMinDistanceIter(
                                        trial_point[select_idx, :], XMat_HF)
                                    if max_dist < dist:
                                        XMat_add_HF = trial_point[select_idx, :]
                                        max_dist = dist
                                    iter_select = iter_select + 1
                        # end select add point GPC

                # end select add point
                # updata to HF library
                if XMat_add_HF.shape[1] != 0:
                    XMat_add_HF, __, __, __, VioMat_add_HF, __, __, NFE_updata = data_lib_HF.dataUpdata(
                        XMat_add_HF, protect_range)
                    NFE = NFE + NFE_updata * cost_HF
                    NFE_list[0] = NFE_list[0] + NFE_updata
                    if NFE_updata > 0:
                        Bool_feas_new_HF = [XMat_add_HF[idx][0] ==
                                            0 for idx in range(len(XMat_add_HF))]
                        Bool_feas_HF = Bool_feas_HF+Bool_feas_new_HF
                        Bool_conv_HF = Bool_conv_HF + \
                            [False]*XMat_add_HF.shape[0]

                # updata to LF library
                if XMat_add_LF.shape[1] != 0:
                    XMat_add_LF, __, __, __, VioMat_add_LF, __, __, NFE_updata = data_lib_LF.dataUpdata(
                        XMat_add_LF, protect_range)
                    NFE = NFE + NFE_updata * cost_LF
                    NFE_list[1] = NFE_list[1] + NFE_updata
                    if NFE_updata > 0:
                        Bool_feas_new_LF = [XMat_add_LF[idx][0] ==
                                            0 for idx in range(len(XMat_add_LF))]
                        Bool_feas_LF = Bool_feas_LF+Bool_feas_new_LF
                        Bool_conv_LF = Bool_conv_LF + \
                            [False]*XMat_add_LF.shape[0]

                # end add point
                # forced interrupt
                if iter > self.iter_max or NFE >= self.NFE_max:
                    done = 1

            # end sample
            obj_best_old = obj_best
            objMat_infill_old = objMat_infill
            con_infill_old = conMat_infill
            coneq_infill_old = coneqMat_infill
            vio_infill_old = vioMat_infill

            XMat_local_best = np.array([[]])
            ObjMat_local_best = np.array([[]])
            XMat_potc = np.array([[]])
            ObjMat_potc = np.array([[]])
            VioMat_potc = np.array([[]])
            detect_local_flag = True

            io.savemat('iter.mat', {
                'data_lib_HF': data_lib_HF, 'data_lib_LF': data_lib_LF, 'X_local_best': XMat_local_best, 'Obj_local_best': ObjMat_local_best, 'X_potential': XMat_potc, 'Obj_potential': ObjMat_potc, 'Vio_potential': VioMat_potc, 'detect_local_flag': detect_local_flag})

        # end algoritm
        # find best result to record
        x_best = data_lib_HF.X[data_lib_HF.result_best_idx[-1]]
        obj_best = data_lib_HF.Obj[data_lib_HF.result_best_idx[-1]]

        output = {}
        output['result_x_best'] = result_x_best
        output['result_obj_best'] = result_obj_best
        output['x_local_best'] = XMat_local_best
        output['obj_local_best'] = ObjMat_local_best
        output['data_lib_HF'] = data_lib_HF
        output['data_lib_LF'] = data_lib_LF

        io.savemat('result.mat', {
            'data_lib_HF': data_lib_HF, 'data_lib_LF': data_lib_LF, 'X_local_best': XMat_local_best, 'Obj_local_best': ObjMat_local_best})

        return x_best, obj_best, NFE, output


if __name__ == '__main__':
    import problem

    # ________________SF test________________
    # prob = problem.G6()

    # def model(xMat):
    #     objMat = prob.objfcn(xMat)
    #     conMat = prob.confcn(xMat)
    #     coneqMat = []
    #     return objMat, conMat, coneqMat
    # matBnd = prob.matBnd
    # low_bou = np.transpose(matBnd[:, 0])
    # up_bou = np.transpose(matBnd[:, 1])
    # vari_num = prob.options['numVar']

    # def model(x):
    #     obj=x+1
    #     con=x-5
    #     coneq=[]
    #     return x,con,coneq
    # vari_num=1
    # low_bou=[1]
    # up_bou=[2]

    # max_NFE = 50
    # optimizer = RBF_GPC_RS(max_NFE, 300)
    
    # initia_datal=io.loadmat('iter.mat')['data_lib']
    # inital_data_lib={'X':initia_datal['X'][0][0],'Obj':initia_datal['Obj'][0][0],'Con':initia_datal['Con'][0][0],'Coneq':initia_datal['Coneq'][0][0],'Ks':initia_datal['Ks'][0][0],'Vio':initia_datal['Vio'][0][0]}
    # optimizer.data_lib=inital_data_lib
    
    # x_best, obj_best, NFE, output = optimizer.optimize(
    #     model=model, vari_num=vari_num, low_bou=low_bou, up_bou=up_bou)

    # print(x_best)
    # print(obj_best)

    # ________________MF test________________
    # prob = problem.G7MODMF()

    # def HF_model(xMat):
    #     objMat = prob.objfcn(xMat)
    #     conMat = prob.confcn(xMat)
    #     coneqMat = []
    #     return objMat, conMat, coneqMat

    # def LF_model(xMat):
    #     objMat = prob.objfcn_lf(xMat)
    #     conMat = prob.confcn_lf(xMat)
    #     coneqMat = []
    #     return objMat[0], conMat[0], coneqMat
    # MF_model = [HF_model, LF_model]
    # Cost = [1, 0.25]
    # Ratio = [1, 4]

    # matBnd = prob.matBnd
    # low_bou = np.transpose(matBnd[:, 0])
    # up_bou = np.transpose(matBnd[:, 1])
    # vari_num = prob.options['numVar']

    def HF_model(xMat):
        objMat=sum(sum(xMat+5))**2
        conMat=sum(sum(xMat-2))
        coneqMat=[]
        return objMat, conMat, coneqMat

    def LF_model(xMat):
        objMat=sum(sum(xMat+5))**2*0.8
        conMat=sum(sum(xMat-2))*0.9
        coneqMat=[]
        return objMat, conMat, coneqMat
    
    MF_model = [HF_model, LF_model]
    Cost = [1, 0.25]
    Ratio = [1, 4]
    
    low_bou = [-4]*5
    up_bou = [0]*5
    vari_num = 5
    
    max_NFE = 100
    inital_data_lib_HF=None
    inital_data_lib_LF=None
    
    optimizer = RBF_GPC_RS_MF(max_NFE, 300,data_lib_HF=inital_data_lib_HF,data_lib_LF=inital_data_lib_LF)

    x_best, obj_best, NFE, output = optimizer.optimize(
        MF_model=MF_model, Cost=Cost, Ratio=Ratio, vari_num=vari_num, low_bou=low_bou, up_bou=up_bou)

    print(x_best)
    print(obj_best)
    
    