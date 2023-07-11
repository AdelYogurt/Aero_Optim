#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -t 240

# run command
# sbatch -p amd_512 -o config_CFD.log config_CFD.sh
source /public3/soft/modules/module.sh
module load gcc/7.3.0-wzm   python/3.7.6    swig/4.1.0  mpi/intel/2022.1  
export SU2_RUN=/public3/wshome/ws173/scg7758/zhangyao/program/SU2/bin/
export PATH=$SU2_RUN:$PATH
mpirun -np 128 SU2_CFD config_CFD.cfg
