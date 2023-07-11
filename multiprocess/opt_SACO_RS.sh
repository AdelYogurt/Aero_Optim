#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 128

# dos2unix *.py
# dos2unix *.m

# run command
# sbatch -p amd_512 -o opt_SACO_RS.log opt_SACO_RS.sh
source /public3/soft/modules/module.sh
module load python/3.7.6  mpi/intel/2022.1  anaconda/3-Python-3.8.3-phonopy-phono3py
export SU2_RUN=/public3/wshome/ws173/scg7758/zhangyao/program/SU2/bin/
export PATH=$SU2_RUN:$PATH
export PATH=/public3/wshome/ws173/scg7758/software-scg7758/matlab2022a/install22a/bin:$PATH

# opti by python
# /public3/wshome/ws173/scg7758/.conda/envs/zyHGV/bin/python3 opt_SACO_RS.py

# opti by matlab
matlab -nojvm -noFigureWindows -nosplash -softwareopengl < opt_SACO_RS.m
