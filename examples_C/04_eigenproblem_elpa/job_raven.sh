#!/bin/bash -l
# Standard output and error:
##SBATCH -o ./job.out.%j
##SBATCH -e ./job.err.%j
# Initial working directory:
#SBATCH -D ./
# Job name
#SBATCH -J 04_eigenproblem_elpa

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1

#SBATCH --mem=12G
#SBATCH --mail-type=ALL
#SBATCH --time=00:03:00

source compile_raven.sh

export LD_LIBRARY_PATH=$ELPA_HOME/lib:$LD_LIBRARY_PATH

set -v

#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_PLACES=cores


srun $FileName 10000 0 32
