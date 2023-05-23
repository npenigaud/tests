#!/bin/bash
#SBATCH -p ndl
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --time "00:05:00"
#SBATCH --switches=3

module load nvhpc/22.11

set -x

cd /home/gmap/mrpm/marguina/tmp/mpicudaaware


export SLURM_EXPORT_ENV='ALL'
export MPIAUTOCONFIG=mpiauto.PGI.conf

#mpirun -np 4 ./mpicudaaware.x / ajouter -nn 2 si -N 2
/opt/softs/mpiauto/mpiauto --verbose -openmp 1 --nouse-slurm-mpi --wrap --wrap-stdeo -np 2 -- ./mpc9.x

