import os
import sys
import numpy as np
from scipy import io

from model.util import MyLogger, safeMakeDirs
from WDB_model import WDBProblem
from optimal.RBF_GPC_RS import RBF_GPC_RS

if __name__ == "__main__":
    # initialize
    my_logger = MyLogger("OptBaseSF.log")
    partitions = 128
    problem=WDBProblem(partitions)
    
    
    def model(x):
        fval, con, coneq, SU2_out, SU2_history, _, _ = problem._HF_model(x,'Opt1SF',SA_retry_flag=True)
        return fval, con, coneq
    vari_num=problem.vari_num
    low_bou=problem.low_bou
    up_bou=problem.up_bou
    
    max_NFE=200
    
    initia_data=io.loadmat('initial_data_SF.mat')['data_lib']
    # initia_data=io.loadmat('iter.mat')['data_lib']
    inital_data_lib={'X':initia_data['X'][0][0],'Obj':initia_data['Obj'][0][0],'Con':initia_data['Con'][0][0],'Coneq':initia_data['Coneq'][0][0],'Ks':initia_data['Ks'][0][0],'Vio':initia_data['Vio'][0][0]}
    optimizer = RBF_GPC_RS(max_NFE-inital_data_lib['X'].shape[1], 300,con_torl=1e-6,data_lib=inital_data_lib)
    
    x_best, obj_best, NFE, output = optimizer.optimize(
        model=model, vari_num=vari_num, low_bou=low_bou, up_bou=up_bou)
    
    io.savemat('opt_base_res.mat', {
                'x_best': x_best, 'data': obj_best,'history':output})
        
    my_logger.logger.info('SU2 all done')