#!/bin/bash -l
# Standard output and error:
##SBATCH -o ./job.out.%j
##SBATCH -e ./job.err.%j
# Initial working directory:
#SBATCH -D ./
# Job name
#SBATCH -J 01_gemm_big_matrices

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

# LUMI-specific
#SBATCH --account=project_465000539

#SBATCH --partition=debug
##SBATCH --partition=small
##SBATCH --reservation=nomad_school_20231006

##SBATCH --mem=2G # 2G per core
#SBATCH --mem=12G
#SBATCH --mail-type=ALL
#SBATCH --time=00:02:00

source compile.sh

set -v

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_PLACES=cores

srun $FileName 10000
