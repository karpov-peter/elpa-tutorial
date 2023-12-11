#!/bin/bash -l
# Standard output and error:
##SBATCH -o ./job.out.%j
##SBATCH -e ./job.err.%j
# Initial working directory:
#SBATCH -D ./
# Job name
#SBATCH -J 00c_gemm_big_matrices

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=1

# LUMI-specific
#SBATCH --account=project_465000866

#SBATCH --partition=dev-g
##SBATCH --partition=standard-g

##SBATCH --mem=2G # 2G per core
#SBATCH --mem=12G
#SBATCH --mail-type=ALL
#SBATCH --time=00:02:00

source compile_lumi.sh

set -v

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_PLACES=cores

srun $HOME/pincheck/pincheck

srun $FileName 10000
