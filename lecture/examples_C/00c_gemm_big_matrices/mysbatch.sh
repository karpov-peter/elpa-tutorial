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
#SBATCH --cpus-per-task=1

# --- default case: use a single GPU on a shared node ---

# !!! for correct running of Fortran tests we need 2 GPUs and 2 MPI tasks

# LUMI

#SBATCH --partition=debug
#SBATCH --account=project_465000539

#SBATCH --mem=2G
#SBATCH --time=00:01:00

source compile.sh

srun $FileName 1000

