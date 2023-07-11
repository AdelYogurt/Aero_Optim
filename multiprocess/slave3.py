import os
import sys
import numpy as np
from scipy import io

if __name__ == "__main__":
    os.chdir('..')
    sys.path.append(os.path.dirname(sys.path[0]))
    
    from eval_model import evalX
    from model.util import MyLogger

    # initialize
    my_logger = MyLogger("multiprocess/slave3_out.log")
    partitions = 128

    data = io.loadmat('LHD.mat')
    X = data['X']
    dir_out = 'model/SU2_temp/slave3'
    eval_range=range(750,1000)
    evalX(X, partitions, dir_out, 2.5, 0.8, my_logger,eval_range)
    my_logger.logger.info('SU2 all done')
