import os
import sys
import numpy as np
from scipy import io

from model.util import MyLogger, safeMakeDirs
from WDB_model import WDBProblem

if __name__ == "__main__":
    # initialize
    partitions = 128
    problem=WDBProblem(partitions)
    
    x=[0.7775,0.9882,3.1125,0.6000,1.5382,1.7204,0.4090,0.2252,0.0104,0.3395,0.6000,0.6613,0.6425,0.0100,0.0100]
    fval, con, coneq, SU2_out, SU2_history, SU2_RANS, SU2_CONV = problem._HF_model(x,'opt',SA_retry_flag=True)
    str_mat_file = "HF_data/HF_opt.mat"
    io.savemat(str_mat_file, {
                'type': list(SU2_out.keys()), 'data': list(SU2_out.values()),'history':list(SU2_history),'conv':SU2_CONV,'x':x})