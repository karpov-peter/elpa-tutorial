#!/bin/bash -l
# Standard output and error:
##SBATCH -o ./job.out.%j
##SBATCH -e ./job.err.%j
# Initial working directory:
#SBATCH -D ./
# Job name
#SBATCH -J 02_pdgemm

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 # MPI
#SBATCH --cpus-per-task=2   # OpenMP

# --- default case: use a single GPU on a shared node ---

# !!! for correct running of Fortran tests we need 2 GPUs and 2 MPI tasks

# LUMI
#SBATCH --account=project_465000539

#SBATCH --reservation=nomad_school_20231006
#SBATCH --partition=small
##SBATCH --partition=debug

#SBATCH --mem=2G
#SBATCH --time=00:01:00

source compile_lumi.sh

set -v

srun $FileName 10000 32

