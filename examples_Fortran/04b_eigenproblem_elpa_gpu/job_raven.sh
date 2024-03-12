#!/bin/bash -l
# Standard output and error:
##SBATCH -o ./job.out.%j
##SBATCH -e ./job.err.%j
# Initial working directory:
#SBATCH -D ./
# Job name
#SBATCH -J 04b_eigenproblem_elpa_gpu

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=18

#SBATCH --constraint="gpu"
#SBATCH --partition=gpudev

# Raven
#SBATCH --gres=gpu:a100:4
#SBATCH --mem=500000

#SBATCH --mail-type=all
#SBATCH --time=00:01:00

module purge
module load gcc/11 cuda/11.4 impi/2021.11 mkl/2023.1
module load elpa/mpi/standard/gpu/2023.11.001

export LD_LIBRARY_PATH=$ELPA_HOME/lib:$LD_LIBRARY_PATH

set -v

#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_PLACES=cores

FileName="eigenproblem_elpa_gpu"

srun $FileName 10000 10000 32
