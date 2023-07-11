import os
import sys
import numpy as np
from scipy import io

from model.util import MyLogger, safeMakeDirs
from WDB_model import WDBProblem
from optimal.RBF_GPC_RS import RBF_GPC_RS_MF

if __name__ == "__main__":
    # initialize
    my_logger = MyLogger("OptBaseMF.log")
    partitions = 128
    problem=WDBProblem(partitions)
    max_NFE=200
    
    def HF_model(x):
        fval, con, coneq, SU2_out, SU2_history, _, _ = problem._HF_model(x,'Opt1MF',SA_retry_flag=True)
        return fval, con, coneq
    def LF_model(x):
        fval, con, coneq,_,_ = problem._LF_model(x,'Opt1MF')
        return fval, con, coneq
    MF_model=[HF_model,LF_model]
    Cost=[1,0.01]
    Ratio=[1,4]
    
    vari_num=problem.vari_num
    low_bou=problem.low_bou
    up_bou=problem.up_bou
    
    initia_data_HF=io.loadmat('initial_data_MF.mat')['data_lib_HF']
    initia_data_LF=io.loadmat('initial_data_MF.mat')['data_lib_LF']
    
    # initia_data_HF=io.loadmat('iter.mat')['data_lib_HF']
    # initia_data_LF=io.loadmat('iter.mat')['data_lib_LF']
    
    inital_data_lib_HF={'X':initia_data_HF['X'][0][0],'Obj':initia_data_HF['Obj'][0][0],'Con':initia_data_HF['Con'][0][0],'Coneq':initia_data_HF['Coneq'][0][0],'Ks':initia_data_HF['Ks'][0][0],'Vio':initia_data_HF['Vio'][0][0]}
    inital_data_lib_LF={'X':initia_data_LF['X'][0][0],'Obj':initia_data_LF['Obj'][0][0],'Con':initia_data_LF['Con'][0][0],'Coneq':initia_data_LF['Coneq'][0][0],'Ks':initia_data_LF['Ks'][0][0],'Vio':initia_data_LF['Vio'][0][0]}
    
    optimizer = RBF_GPC_RS_MF(max_NFE-inital_data_lib_HF['X'].shape[0]*Cost[0]-inital_data_lib_LF['X'].shape[0]*Cost[1], 300,con_torl=1e-6,data_lib_HF=inital_data_lib_HF,data_lib_LF=inital_data_lib_LF)
    
    x_best, obj_best, NFE, output = optimizer.optimize(
        MF_model=MF_model,Cost=Cost,Ratio=Ratio, vari_num=vari_num, low_bou=low_bou, up_bou=up_bou)
    
    io.savemat('opt_base_res.mat', {
                'x_best': x_best, 'data': obj_best,'history':output})
        
    my_logger.logger.info('SU2 all done')