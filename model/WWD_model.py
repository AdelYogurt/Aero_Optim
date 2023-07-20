import os
import sys
import numpy as np
from multiprocessing import current_process
from scipy import io
from stl import mesh

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, proj3d

from model.util import MyLogger, safeCopyDirs, safeMakeDirs, safeSplitStrFile
from geometry.Waverider import WaveriderWingDia
from model.EvalSU2 import runSU2CFDRestart, runSU2CFD, runSU2DEF, readHistory
from model.EvalSHABP import runSHABP
from model.EvalThermo import infStagHeat
from model.EvalAtmos import getAtmosphereEnvironment

WDB_logger = MyLogger("WDB_prblem.log")

class WWDProblem():
    def __init__(self, SU2_partitions, AOA=8, SIDESLIP_ANGLE=0, Ma=13.8, P=143.9885, T=264.9206, T_w=300):
        self.partitions = SU2_partitions
        self.AOA = AOA
        self.SIDESLIP_ANGLE = SIDESLIP_ANGLE
        self.Ma = Ma
        self.T = T
        self.P = P
        self.T_w = T_w

        # problem basic problem define
        self.vari_num = 15
        self.con_num = 3
        self.coneq_num = 0

        # problem constraint define
        self.HF_max = 8e6
        self.Cl_min = 0.04
        self.V_eff_min = 0.23

        self.low_bou=[0.7,0.7,2.4,0.6,0.5,0.5,0.1,0.1,0.005,0.25,0.6,0.6,0.5,0.01,0.01]
        self.up_bou=[1.0,1.0,3.2,0.8,2.0,2.0,0.5,0.5,0.02,0.35,0.8,0.8,0.7,0.05,0.05]

        # temp file
        safeMakeDirs('model/SU2_temp', WDB_logger)
        safeMakeDirs('model/SHABP_temp', WDB_logger)

    def proPoint(self, xMat, str_data_file=None, dir_point_data=None):
        total_length = 4

        xMat=np.array(xMat).reshape((self.vari_num,))
        # length parameter
        par_M_up = xMat[0]
        par_M_low = xMat[1]
        # width parameter
        par_width = xMat[2]
        par_T = xMat[3]
        par_N_up = xMat[4]
        par_N_low = xMat[5]
        # height parameter
        par_hight_up = xMat[6]
        par_hight_low = xMat[7]
        par_R = xMat[8]
        # wing parameter
        par_rho1 = xMat[9]
        par_rho12 = xMat[10]
        par_WS1 = xMat[11]
        par_WS2 = xMat[12]
        par_TWU = xMat[13]
        par_TWL = xMat[14]

        # gen object
        waverider_wing = WaveriderWingDia(total_length, par_width, par_hight_up, par_hight_low, par_T, par_M_up,
                                          par_M_low, par_N_up, par_N_low, par_R, par_rho1, par_rho12, par_WS1, par_WS2, par_TWU, par_TWL)
        surface_name_list = list(waverider_wing.surface_list.keys())
        total_point_list = list()

        def stag_deform_Y_up(Y, Z): return (np.minimum(
            (Y / (par_R*3)) ** 2 + ((Z-par_R)/(par_R))**2, 1)) ** (1.0 - par_T / 0.7)

        def stag_deform_Y_low(Y, Z): return (np.minimum(
            (Y / (par_R*3)) ** 2 + ((Z+par_R)/(par_R))**2, 1)) ** (1.0 - par_T / 0.7)
        
        def stag_deform_Y_mid(Y, Z): return (np.minimum(
            (Y / (par_R*3)) ** 2 + ((Z)/(par_R))**2, 1)) ** (1.0 - par_T / 0.7)

        for surface_index in range(len(surface_name_list)):
            surface_name = surface_name_list[surface_index]
            data_file_name = surface_name+'_local_coord.txt'
            with open(os.path.join(dir_point_data, data_file_name)) as file_data:
                str_data = file_data.readlines()

            data = [data_unit.split() for data_unit in str_data]
            index = [int(data_rank[0]) for data_rank in data]
            XI = np.array([float(data_rank[1]) for data_rank in data])
            PSI = np.array([float(data_rank[2]) for data_rank in data])

            X, Y, Z = waverider_wing.surface_list[surface_name].calPoint(
                XI, PSI)

            if str(surface_name) == str('head_side'):

                data_file_name = 'symmetry.txt'
                with open(os.path.join(dir_point_data, data_file_name)) as file_data:
                    str_data = file_data.readlines()

                data = [data_unit.split() for data_unit in str_data]
                symmetry_index = [int(data_rank[0])
                                  for data_rank in data if data_rank]

                sam_index = np.intersect1d(index, symmetry_index)
                edge_index = [equal for equal in range(
                    len(index)) if index[equal] in sam_index]
                y_edge = Y[edge_index]
                z_edge = Z[edge_index]

                def y_edge_func(Z=None): return self.interpLinear(
                    Z, z_edge, y_edge)

                PSI = Y / waverider_wing.y_cut
                Y = Y - (1 - PSI)*y_edge_func(Z)

                # deform of stag
                Y = Y*stag_deform_Y_up(Y, Z)
                Y = Y*stag_deform_Y_low(Y, Z)
                Y = Y*stag_deform_Y_mid(Y,Z)
            if str(surface_name) == str('head_up') or str(surface_name) == str('head_low'):
                # deform of stag
                Y = Y*stag_deform_Y_up(Y, Z)
                Y = Y*stag_deform_Y_low(Y, Z)
                Y = Y*stag_deform_Y_mid(Y,Z)

            total_point_list.append([index, X, Y, Z])

        # process repeat
        for base_index in range(len(surface_name_list)):
            point_base = total_point_list[base_index]
            index_base = point_base[0]

            for check_index in range(base_index + 1, len(surface_name_list)):
                point_check = total_point_list[check_index]
                index_check = point_check[0]

                index_equal = np.intersect1d(index_base, index_check)
                posi_base = [equal for equal in range(
                    len(index_base)) if index_base[equal] in index_equal]
                posi_check = [equal for equal in range(
                    len(index_check)) if index_check[equal] in index_equal]
                if not len(index_equal) == 0:
                    X_equal = (point_base[1][posi_base] +
                               point_check[1][posi_check]) / 2
                    Y_equal = (point_base[2][posi_base] +
                               point_check[2][posi_check]) / 2
                    Z_equal = (point_base[3][posi_base] +
                               point_check[3][posi_check]) / 2
                    point_base[1][posi_base] = X_equal
                    point_check[1][posi_check] = X_equal
                    point_base[2][posi_base] = Y_equal
                    point_check[2][posi_check] = Y_equal
                    point_base[3][posi_base] = Z_equal
                    point_check[3][posi_check] = Z_equal

                total_point_list[check_index] = point_check

            total_point_list[base_index] = point_base

        # io.savemat('data.mat',{'point':total_point_list})

        return total_point_list

    def interpLinear(self, X_pred, X_origin, Y_origin):
        X = sorted(X_origin)
        index = sorted(range(len(X_origin)), key=lambda k: X_origin[k])
        Y = Y_origin[index]
        Y_pred = list()
        for x_index in range(len(X_pred)):
            x_pred = X_pred[x_index]
            num = len(X)
            index = num-1
            while ((index > 0) and (X[index] > x_pred)):
                index = index - 1

            if (index == num-1):
                # out interp
                Y_pred.append((Y[-1] - Y[-2]) / (X[-1] - X[-2])
                              * (x_pred - X[-1]) + Y[-1])
            else:
                if (index == 0):
                    Y_pred.append((Y[2] - Y[1]) / (X[2] - X[1])
                                  * (x_pred - X[1]) + Y[1])
                else:
                    # linear interpolation
                    Y_pred.append(Y[index] + (Y[index + 1] - Y[index])
                                  * (x_pred - X[index]) / (X[index + 1] - X[index]))

        return Y_pred

    def writePoint(self, xMat, str_data_file, dir_point_data):
        total_point_list = self.proPoint(xMat, str_data_file, dir_point_data)

        point_out = open(str_data_file, 'w')
        for surface_index in range(len(total_point_list)):
            point_list = total_point_list[surface_index]
            for point_index in range(len(point_list[0])):
                point_out.write('%d %f %f %f\n' % (
                    point_list[0][point_index]-1, point_list[1][point_index], point_list[2][point_index], point_list[3][point_index]))

        point_out.close()

        return total_point_list

    def _HF_model(self, xMat, description: str = None, SA_retry_flag: bool = True):

        str_data_file = 'model/WDB_deform.dat'
        dir_point_data = 'coord_data'

        # rebulid surface point
        self.writePoint(xMat, str_data_file, dir_point_data)

        # generate new grid
        str_mesh_file = 'initial_data/WDB.cgns'
        str_dat_file = 'model/WDB_deform.dat'
        str_cfg_file = 'initial_data/WDB_deform.cfg'
        str_mesh_out = 'WDB_deform.su2'
        runSU2DEF(str_mesh_file, str_cfg_file, str_dat_file, str_mesh_out, self.partitions,
                  description=description)

        procid = int(current_process().pid)
        dir_temp = os.path.join(os.path.dirname(
            os.path.abspath(__file__)), "model/SU2_temp")
        dir_work = os.path.join(dir_temp, description+"_{}".format(procid))

        try:
            str_cfg_file = "initial_data/WDB.cfg"
            str_mesh_file = 'WDB_deform.su2'
            if os.path.exists('initial_data/restart_flow.dat'):
                str_restart_file = 'initial_data/restart_flow.dat'
                SU2_out, SU2_history = runSU2CFDRestart(str_mesh_file, str_cfg_file, str_restart_file, self.partitions,
                                                        self.AOA, self.SIDESLIP_ANGLE, self.Ma, self.P, self.T,
                                                        description=description)
            else:
                SU2_out, SU2_history = runSU2CFD(str_mesh_file, str_cfg_file, self.partitions,
                                                    self.AOA, self.SIDESLIP_ANGLE, self.Ma, self.P, self.T,
                                                    description=description)
            SU2_RANS = 'SST'
            SU2_CONV = 'True'
            Cl = SU2_out['CL']
            CD = SU2_out['CD']
            HF = SU2_out['maxHF']
            CEff = Cl/CD
        except:
            SU2_fail = True
            WDB_logger.logger.info('SU2 SST-LM2003 analysis fail')

            # if specify retry by SA model, rerun SU2 by retry cfg
            if SA_retry_flag:
                WDB_logger.logger.info('retry SU2 SA analysis')
                try:
                    str_cfg_file = "initial_data/WDB_retry.cfg"
                    str_mesh_file = 'WDB_deform.su2'
                    SU2_out, SU2_history = runSU2CFD(str_mesh_file, str_cfg_file, self.partitions,
                                                        self.AOA, self.SIDESLIP_ANGLE, self.Ma, self.P, self.T,
                                                        description=description)
                    SU2_RANS = 'SA'
                    SU2_CONV = 'True'
                    Cl = SU2_out['CL']
                    CD = SU2_out['CD']
                    HF = SU2_out['maxHF']
                    CEff = Cl/CD
                    SU2_fail = False
                except:
                    SU2_fail = True

            if SU2_fail:
                run_SHABP = True

                str_result_file = os.path.join(dir_work, 'history.csv')
                if os.path.exists(str_result_file):
                    SU2_out, SU2_history = readHistory(str_result_file)
                    
                    if SU2_out['Inner_Iter'] > 100 and SU2_out['rms[Rho]'] > -5:
                        WDB_logger.logger.info(
                            'SU2 analysis fail beform converage, use unconverage data')
                        
                        SU2_CONV = 'False'
                        Cl = SU2_out['CL']
                        CD = SU2_out['CD']
                        HF = SU2_out['maxHF']
                        CEff = Cl/CD
                        SU2_RANS='SA'
                        
                        run_SHABP = False

                if run_SHABP:
                    WDB_logger.logger.info(
                        'all SU2 analysis fail, instead by SHABP')
                    # instead by SHABP
                    fval, con, coneq,  SHABP_out, HF = self._LF_model(xMat)
                    Cl = SHABP_out['clift']
                    Cd = SHABP_out['cdrag']
                    CEff = -2

                    SU2_out = {"data": "error"}
                    SU2_CONV = 'error'
                    SU2_history = 'error'
                    SU2_RANS = 'error'

        # evaluate geometry
        S, V, V_eff = self.calGEO(xMat)

        fval = -CEff
        con = [(HF-self.HF_max)/8e6, self.Cl_min-Cl, self.V_eff_min-V_eff]
        coneq = []

        return fval, con, coneq, SU2_out, SU2_history, SU2_RANS, SU2_CONV

    def _LF_model(self, xMat, cmethods=None, emethods=None, Q_method='Detra-Kemp-Riddell'):

        total_length = 4

        xMat=np.array(xMat).reshape((self.vari_num,))
        # length parameter
        par_M_up = xMat[0]
        par_M_low = xMat[1]
        # width parameter
        par_width = xMat[2]
        par_T = xMat[3]
        par_N_up = xMat[4]
        par_N_low = xMat[5]
        # height parameter
        par_hight_up = xMat[6]
        par_hight_low = xMat[7]
        par_R = xMat[8]
        # wing parameter
        par_rho1 = xMat[9]
        par_rho12 = xMat[10]
        par_WS1 = xMat[11]
        par_WS2 = xMat[12]
        par_TWU = xMat[13]
        par_TWL = xMat[14]

        # gen object
        waverider_wing = WaveriderWingDia(total_length, par_width, par_hight_up, par_hight_low, par_T, par_M_up,
                                          par_M_low, par_N_up, par_N_low, par_R, par_rho1, par_rho12, par_WS1, par_WS2, par_TWU, par_TWL)

        xi_grid_num_head = 30
        eta_grid_num_head = 20
        xi_grid_num_body = 20
        eta_grid_num_wing = 6

        edge_gird_num = int((par_R / 0.003))

        # get wgs
        waverider_wing.writeWGS('model/WDB_deform', xi_grid_num_head,
                                eta_grid_num_head, xi_grid_num_body, eta_grid_num_wing, edge_gird_num)

        # run SHABP to get aerodynamics
        Ma = 13.8
        AOA = 8
        str_wgs_file = 'model/WDB_deform.wgs'
        str_inp_file = None
        cmethods = 14
        emethods = 3

        SHABP_out = runSHABP(str_inp_file, str_wgs_file,
                             self.Ma, self.AOA, cmethods, emethods, dir_temp='model/SHABP_temp')
        WDB_logger.logger.info("S/HABP analysis finished")
        Cl = SHABP_out['clift']
        Cd = SHABP_out['cdrag']

        # use equation to estimate heat flux

        if par_M_up > 1:
            radius_center_up = 0
        else:
            dz_dx = (par_hight_up * (0.002) ** par_M_up) / \
                (total_length * 0.002)
            radius_center_up = par_R * dz_dx

        if par_M_low > 1:
            radius_center_low = 0
        else:
            dz_dx = (par_hight_low * (0.002) ** par_M_low) / \
                (total_length * 0.002)
            radius_center_low = par_R * dz_dx

        radius_center_head = (radius_center_up + radius_center_low) / 2
        radius_head_sq = (radius_center_head *
                          radius_center_head + par_R * par_R)
        radius_head_ZX = np.sqrt(radius_head_sq)
        # convert to y=f(xMat) can calculate xMat=0 point curvature
        radius_head_XY = 2**(1/par_T-2) * 2*par_width/(1/par_T/2)

        HF = (infStagHeat(radius_head_ZX, radius_head_ZX,
                         self.Ma, self.P, self.T, self.T_w, None, Q_method))

        S, V, V_eff = self.calGEO(xMat)

        fval = -Cl/Cd
        con = [(HF-self.HF_max)/8e6, self.Cl_min -
               Cl, self.V_eff_min-V_eff]
        coneq = []

        return fval, con, coneq, SHABP_out, HF

    def calGEO(self, xMat):
        total_length = 4

        xMat=np.array(xMat).reshape((self.vari_num,))
        # length parameter
        par_M_up = xMat[0]
        par_M_low = xMat[1]
        # width parameter
        par_width = xMat[2]
        par_T = xMat[3]
        par_N_up = xMat[4]
        par_N_low = xMat[5]
        # height parameter
        par_hight_up = xMat[6]
        par_hight_low = xMat[7]
        par_R = xMat[8]
        # wing parameter
        par_rho1 = xMat[9]
        par_rho12 = xMat[10]
        par_WS1 = xMat[11]
        par_WS2 = xMat[12]
        par_TWU = xMat[13]
        par_TWL = xMat[14]

        # gen object
        waverider_wing = WaveriderWingDia(total_length, par_width, par_hight_up, par_hight_low, par_T, par_M_up,
                                          par_M_low, par_N_up, par_N_low, par_R, par_rho1, par_rho12, par_WS1, par_WS2, par_TWU, par_TWL, False)

        xi_grid_num_head = 30
        eta_grid_num_head = 20
        xi_grid_num_body = 20
        eta_grid_num_wing = 6

        edge_gird_num = int((par_R / 0.003))

        waverider_wing.writeSTL('model/waverider_wing_dia_full', xi_grid_num_head,
                                eta_grid_num_head, xi_grid_num_body, eta_grid_num_wing, edge_gird_num)

        mymesh = mesh.Mesh.from_file('model/waverider_wing_dia_full.stl')
        V = mymesh.get_mass_properties()[0]
        S = float(sum(mymesh.areas))

        V_eff = 6*np.sqrt(np.pi) * V / S**(3/2)

        return S, V, V_eff

if __name__ == "__main__":
    # initialize
    my_logger = MyLogger("OptBaseMF.log")
    partitions = 128
    problem=WDBProblem(partitions)
    
    print(problem.calGEO((np.array(problem.low_bou)+np.array(problem.up_bou))/2))
    
    