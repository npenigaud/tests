#!/bin/bash
#SBATCH -N1
#SBATCH --time 00:05:00
##SBATCH --time 02:00:00
#SBATCH --exclusive
#SBATCH --export="NONE"
#SBATCH --partition ndl

module load nvhpc/22.11

set -x

#export DR_HOOK=1
#export DR_HOOK_OPT=prof
#export DR_HOOK_IGNORE_SIGNALS=-1
#export DR_NVTX=1

#blocking 
#export CUDA_LAUNCH_BLOCKING=1
#export NVCOMPILER_ACC_SYNCHRONOUS=1

ulimit -s unlimited
export OMP_STACK_SIZE=4G

cd $SLURM_SUBMIT_DIR

\rm *.grb

export SLURM_EXPORT_ENV=ALL
export MPIAUTOCONFIG=mpiauto.PGI.conf
#export MPIAUTOCONFIG=mpiauto.DDT.conf

srun nsys profile testmmv11b.x
#srun ncu --set full -o ncureport.ncu testmmv2
