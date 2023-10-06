#!/bin/bash -l
# Standard output and error:
##SBATCH -o ./job.out.%j
##SBATCH -e ./job.err.%j
# Initial working directory:
#SBATCH -D ./
# Job name
#SBATCH -J 02_pdgemm

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1

# --- default case: use a single GPU on a shared node ---

# !!! for correct running of Fortran tests we need 2 GPUs and 2 MPI tasks

# LUMI
#SBATCH --account=project_465000539

#SBATCH --reservation=nomad_school_20231006
#SBATCH --partition=small
##SBATCH --partition=debug

#SBATCH --mem=12G
#SBATCH --time=00:02:00

source compile.sh

set -v

srun $FileName 10000 32

