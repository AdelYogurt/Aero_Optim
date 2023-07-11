import os
import sys
import numpy as np
from scipy import io

from model.util import MyLogger, safeMakeDirs
from model.airfoil_model import *


def evalX(X, partitions, dir_out, AOA, Ma, my_logger, eval_range=None):
    if eval_range is None:
        eval_range = range(0, len(X))
    problem = AirfoilProblem(partitions, my_logger)

    safeMakeDirs('data', my_logger)
    safeMakeDirs('model/SU2_temp', my_logger)
    safeMakeDirs(dir_out, my_logger)

    for x_index in eval_range:
        x = X[x_index]

        coord_data_dir = 'coord_data'
        # generate Bezier controll point
        control_number = int(len(x)/2)
        control_point_up = np.zeros((control_number, 2))
        control_point_up[:, 0] = np.linspace(0, 1, control_number)
        control_point_low = np.zeros((control_number, 2))
        control_point_low[:, 0] = np.linspace(0, 1, control_number)

        # add xMat
        control_point_up[:, 1] = x[:control_number]
        control_point_low[:, 1] = x[control_number:]

        my_logger.logger.info("current data")
        my_logger.logger.info('index: '+str(x_index))
        my_logger.logger.info(str(x))

        total_point_list, SU2_out, SU2_history, SU2_surface = problem.solveAirfoilSU2(
            control_point_up, control_point_low, coord_data_dir, AOA=AOA, Ma=Ma, description=str(x_index), initial_data_dir='SU2/', dir_out=dir_out)

        my_logger.logger.info(SU2_out)
        str_mat_file = "data/"+str(x_index)+".mat"
        io.savemat(str_mat_file, {
                   'x': x, 'type': list(SU2_out.keys()), 'data': list(SU2_out.values()), 'history': list(SU2_history), 'surface': SU2_surface, 'AOA': AOA, 'Ma': Ma})


if __name__ == "__main__":

    # procid = int(current_process().pid)
    # dir_temp = os.path.join(os.path.dirname(
    #     os.path.abspath(__file__)), "model/SU2_temp")
    # if description is not None:
    #     dir_work = os.path.join(dir_temp, description+"_{}".format(procid))
    # else:
    #     dir_work = os.path.join(dir_temp, "_{}".format(procid))

    # initialize
    my_logger = MyLogger("eval_model.log")
    partitions = 4

    # problem = NACA0012AirfoilProblem(partitions, AOA=2.5, Ma=0.8)

    data = io.loadmat('LHD.mat')
    X = data['X']
    dir_out = 'model/SU2_temp/out'
    evalX(X, partitions, dir_out, 2.5, 0.8, my_logger)
    my_logger.logger.info('SU2 all done')

    # safeMakeDirs('data', my_logger)
    # safeMakeDirs('model/SU2_temp', my_logger)
    # dir_out = 'model/SU2_temp/slave1'
    # safeMakeDirs(dir_out, my_logger)

    # for x_index in range(0, 500):
    #     x = X[x_index]

    #     my_logger.logger.info("current data")
    #     my_logger.logger.info('index: '+str(x_index))
    #     my_logger.logger.info(str(x))
    #     fval, con, coneq, SU2_out, SU2_history, SU2_surface, SU2_RANS, SU2_CONV = problem.model(
    #         x, str(x_index),initial_data_dir='SU2/',dir_out=dir_out)
    #     my_logger.logger.info(str(fval))
    #     my_logger.logger.info(str(con))
    #     my_logger.logger.info(str(coneq))
    #     str_mat_file = "data/"+str(x_index)+".mat"
    #     io.savemat(str_mat_file, {
    #                'type': list(SU2_out.keys()), 'data': list(SU2_out.values()), 'history': list(SU2_history), 'surface': SU2_surface, 'conv': SU2_CONV, 'x': x, 'fval': fval, 'con': con, 'coneq': coneq})

    # my_logger.logger.info('SU2 all done')
