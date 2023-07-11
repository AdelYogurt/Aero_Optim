#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 128

# dos2unix *.py
# dos2unix *.m

# run command
# sbatch -p amd_512 -o opt_SADE_RBF.log opt_SADE_RBF.sh
source /public3/soft/modules/module.sh
module load python/3.7.6  mpi/intel/2022.1  anaconda/3-Python-3.8.3-phonopy-phono3py
export SU2_RUN=/public3/wshome/ws173/scg7758/zhangyao/program/SU2/bin/
export PATH=$SU2_RUN:$PATH
export PATH=/public3/wshome/ws173/scg7758/software-scg7758/matlab2022a/install22a/bin:$PATH

# opti by python
# /public3/wshome/ws173/scg7758/.conda/envs/zyHGV/bin/python3 opt_SADE_RBF.py

# opti by matlab
matlab -nojvm -noFigureWindows -nosplash -softwareopengl < opt_SADE_RBF.m
