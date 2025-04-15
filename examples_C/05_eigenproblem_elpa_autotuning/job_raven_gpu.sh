#!/bin/bash -l
# Standard output and error:
##SBATCH -o ./job.out.%j
##SBATCH -e ./job.err.%j
# Initial working directory:
#SBATCH -D ./
# Job name
#SBATCH -J eigenproblem_elpa_autotuning

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

FileName="eigenproblem_elpa_autotuning"

# Run the program:
srun ./$FileName
