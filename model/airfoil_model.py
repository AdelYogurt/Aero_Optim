import os
import numpy as np
from scipy import io

from model.util import safeMakeDirs
from model.EvalSU2 import runSU2CFD, runSU2DEF

from geometry.CurveBezier import CurveBezier
from geometry.CST import CST2DCurve


def importdata(filestr):
    with open(filestr) as data_file:
        str_data = data_file.readlines()
    data = [data_unit.split()[0].split(',') for data_unit in str_data]
    return data


class Airfoil():
    def __init__(self, control_point_up, control_point_low):
        class_par_Y = [0.5, 1]

        self.bezier_up = CurveBezier(control_point_up)
        self.bezier_low = CurveBezier(control_point_low)
        self.airfoil_up = CST2DCurve(1, 1, lambda u: self.bezier_up.interpPoint(u)[
            :, 1].reshape((1, -1)), class_par_Y)
        self.airfoil_low = CST2DCurve(1, -1, lambda u: self.bezier_low.interpPoint(u)[
            :, 1].reshape((1, -1)), class_par_Y)


class AirfoilProblem():
    def __init__(self, SU2_partitions, out_logger=None):
        self.partitions = SU2_partitions
        self.out_logger = out_logger

    def prePoint(self, control_point_up, control_point_low, coord_data_dir, dat_filestr=None):
        new_airfoil = Airfoil(control_point_up, control_point_low)

        airfoil_coord_up = importdata(os.path.join(
            coord_data_dir, 'airfoil_up_local_coord.txt'))
        X_up, Y_up = new_airfoil.airfoil_up.calPoint(
            [float(data[1]) for data in airfoil_coord_up])

        airfoil_coord_low = importdata(os.path.join(
            coord_data_dir, 'airfoil_low_local_coord.txt'))
        X_low, Y_low = new_airfoil.airfoil_low.calPoint(
            [float(data[1]) for data in airfoil_coord_low])

        total_point_list = [[[int(data[0]) for data in airfoil_coord_up], X_up, Y_up], [
            [int(data[0]) for data in airfoil_coord_low], X_low, Y_low]]

        # io.savemat('data.mat', {'point': total_point_list})

        if dat_filestr is not None:
            point_out = open(dat_filestr, 'w')
            for surface_index in range(len(total_point_list)):
                point_list = total_point_list[surface_index]
                for point_index in range(len(point_list[0])):
                    point_out.write('%d %f %f\n' % (
                        point_list[0][point_index]-1, point_list[1][0, point_index], point_list[2][0, point_index]))
            point_out.close()

        return total_point_list

    def solveAirfoilSU2(self, control_point_up, control_point_low, coord_data_dir,
                        AOA=0.0, SIDESLIP_ANGLE=0.0, Ma=0.0, T=273.0, P=101325, T_w=273.0,
                        description: str = None,
                        initial_data_dir='initial_data', dir_out='model', restart_filestr=None):
        dat_filestr = os.path.join(dir_out, 'airfoil_deform.dat')
        mesh_out_filestr = os.path.join(dir_out, 'airfoil_deform.su2')

        # rebulid surface point
        mesh_filestr = os.path.join(initial_data_dir, 'NACA0012.cgns')
        cfg_filestr = os.path.join(initial_data_dir, 'NACA0012_deform.cfg')
        total_point_list = self.prePoint(
            control_point_up, control_point_low, coord_data_dir, dat_filestr)
        runSU2DEF(mesh_filestr, cfg_filestr, dat_filestr, mesh_out_filestr, self.partitions,DEF_config=None,
                  description=description, out_logger=self.out_logger)

        # run SU2_CFD to get result
        cfg_filestr = os.path.join(initial_data_dir, 'airfoil_turb.cfg')
        mesh_filestr = mesh_out_filestr
        CFD_config=dict()
        CFD_config['AOA']=AOA
        CFD_config['SIDESLIP_ANGLE']=SIDESLIP_ANGLE
        CFD_config['MACH_NUMBER']=Ma
        CFD_config['FREESTREAM_TEMPERATURE']=T
        CFD_config['FREESTREAM_PRESSURE']=P
        SU2_out, SU2_history, SU2_surface = runSU2CFD(mesh_filestr, cfg_filestr, self.partitions,CFD_config,
                                                      restart_filestr,
                                                      description=description, out_logger=self.out_logger)

        return total_point_list, SU2_out, SU2_history, SU2_surface


class NACA0012AirfoilProblem(AirfoilProblem):
    def __init__(self, SU2_partitions, AOA=2.0, SIDESLIP_ANGLE=0.0, Ma=0.63, T=273.0, P=101325, T_w=273.0,
                 description: str = None,
                 initial_data_dir='initial_data', dir_out='model', restart_filestr=None):
        self.partitions = SU2_partitions
        self.AOA = AOA
        self.SIDESLIP_ANGLE = SIDESLIP_ANGLE
        self.Ma = Ma
        self.T = T
        self.P = P
        self.T_w = T_w

        # problem basic problem define
        self.vari_num = 10
        self.con_num = 2
        self.coneq_num = 0

        # problem constraint define
        self.t_mim = 0.096
        self.Cl_min = 0.268863

        self.low_bou = np.array([0.5040,2.1551,-0.4582,1.1015,0.1962,
                                0.5040,2.1551,-0.4582,1.1015,0.1962]).reshape((1, self.vari_num))*1e-2
        self.up_bou = np.array([10.5040,12.1551,9.5418,11.1015,10.1962,
                               10.5040,12.1551,9.5418,11.1015,10.1962]).reshape((1, self.vari_num))*1e-2

        self.description = description
        self.initial_data_dir = initial_data_dir
        self.coord_data_dir = 'coord_data'
        self.dir_out = dir_out
        self.restart_filestr = restart_filestr

        # temp file
        safeMakeDirs('model/SU2_temp')

    def model(self, xMat):

        xMat = np.array(xMat).reshape((1, -1))
        x_control_number = int(xMat.shape[1]/2)

        # generate Bezier controll point
        control_point_up = np.zeros((x_control_number+2, 2))
        control_point_up[0, :] = [0.000000000000, 0.067924153974]
        control_point_up[-1, :] = [1.000000000000, 0.055565107922]
        control_point_up[:, 0] = np.linspace(0, 1, x_control_number+2)
        control_point_low = np.zeros((x_control_number+2, 2))
        control_point_low[0, :] = [0.000000000000, 0.067924153974]
        control_point_low[-1, :] = [1.000000000000, 0.055565107922]
        control_point_low[:, 0] = np.linspace(0, 1, x_control_number+2)

        # add xMat
        control_point_up[1:-1, 1] = xMat[0, :x_control_number]
        control_point_low[1:-1, 1] = xMat[0, x_control_number:]

        total_point_list, SU2_out, SU2_history, SU2_surface = self.solveAirfoilSU2(control_point_up, control_point_low, self.coord_data_dir,
                                                                                   self.AOA, self.SIDESLIP_ANGLE, self.Ma, self.T, self.P, self.T_w,
                                                                                   description=self.description,
                                                                                   initial_data_dir=self.initial_data_dir, dir_out=self.dir_out, restart_filestr=self.restart_filestr)

        # evaluate geometry
        y_up = total_point_list[0][2]
        y_low = total_point_list[1][2]
        t = np.max(y_up[0, :]-y_low[0, :])

        SU2_RANS = 'SST'
        SU2_CONV = 'True'
        Cl = SU2_out['CL']
        CD = SU2_out['CD']
        CEff = Cl/CD

        fval = -CEff
        con = [self.Cl_min-Cl, self.t_mim-t]
        coneq = []

        return fval, con, coneq, SU2_out, SU2_history, SU2_surface, SU2_RANS, SU2_CONV


if __name__ == '__main__':
    porblem = NACA0012AirfoilProblem(1)
    print(porblem.model(np.array([0]*6+[0.1]*6).reshape((1, 12))))
