#!/bin/bash -l
# Standard output and error:
##SBATCH -o ./job.out.%j
##SBATCH -e ./job.err.%j
# Initial working directory:
#SBATCH -D ./
# Job name
#SBATCH -J iter_refinement
#
##SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --constraint="gpu"
#SBATCH --partition=gpudev

# Raven
#SBATCH --gres=gpu:a100:4
#SBATCH --mem=500000

#
# --- uncomment to use 1 GPUs on a shared node ---
# #SBATCH --gres=gpu:a100:1
# #SBATCH --cpus-per-task=18
# #SBATCH --mem=125000
#
# --- uncomment to use 2 GPUs on a shared node ---
# #SBATCH --gres=gpu:a100:2
# #SBATCH --cpus-per-task=36
# #SBATCH --mem=250000
#
# --- uncomment to use 4 GPUs on a full node ---
# #SBATCH --gres=gpu:a100:4
# #SBATCH --cpus-per-task=72
# #SBATCH --mem=500000
#
#SBATCH --mail-type=all
##SBATCH --mail-user=userid@example.mpg.de
#SBATCH --time=00:10:00

module purge
module load intel/21.2.0  
module load impi/2021.2
module load mkl/2021.2
module load cuda/11.2


#module load elpa/mpi/standard/2022.05.001

#FileName="main"

#mpiicc -I$ELPA_HOME/include/elpa-2022.05.001 $FileName.c -o $FileName $ELPA_LIBS -Wl,-rpath=$ELPA_HOME/lib


# Run the program:
srun ./main 1000
