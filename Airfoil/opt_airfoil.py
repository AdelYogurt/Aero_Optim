import os
import sys
import numpy as np
from scipy import io

from airfoil_model import AirfoilProblem
from optimal.RBF_GPC_RS import RBF_GPC_RS

if __name__ == "__main__":
    # initialize
    partitions = 127
    problem = AirfoilProblem(partitions, AOA=2.5, Ma=0.8)

    def model(x):
        fval, con, coneq, SU2_out, SU2_history, _, _ = problem.model(x, 'Opt')
        return fval, con, coneq
    vari_num = problem.vari_num
    low_bou = problem.low_bou
    up_bou = problem.up_bou

    max_NFE = 200

    # initia_data=io.loadmat('initial_data_SF.mat')['data_lib']
    # initia_data = io.loadmat('iter.mat')['data_lib']
    # inital_data_lib = {'X': initia_data['X'][0][0], 'Obj': initia_data['Obj'][0][0], 'Con': initia_data['Con'][0]
    #                    [0], 'Coneq': initia_data['Coneq'][0][0], 'Ks': initia_data['Ks'][0][0], 'Vio': initia_data['Vio'][0][0]}
    optimizer = RBF_GPC_RS(max_NFE, 300, data_lib=None)

    x_best, obj_best, NFE, output = optimizer.optimize(
        model=model, vari_num=vari_num, low_bou=low_bou, up_bou=up_bou)

    io.savemat('opt_base_res.mat', {
        'x_best': x_best, 'data': obj_best, 'history': output})

    print('optimization all done')
