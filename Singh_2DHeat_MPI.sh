#!/bin/bash

#Account and Email Information
###SBATCH -A gurpalsingh
#SBATCH --mail-type=end
#SBATCH --mail-user=gurpalsingh@u.boisestate.edu

# Specify parition (queue)
#SBATCH --partition=batch

# Join output and errors into output
#SBATCH -o output/MPI_SIM1_slurm.o%j
#SBATCH -e output/MPI_SIM1_slurm.e%j

# Specify job not to be rerunable
#SBATCH --no-requeue

# Job Name
#SBATCH --job-name="MPI_SIM1"

# Specify walltime
###SBATCH --time=48:00:00

# Specify number of requested nodes:
#SBATCH -N 3
# Specify the total number of requested procs:
#SBATCH -n 12
# Specify the procs per node:
#SBATCH --ntasks-per-node=4

module load openmpi/gcc-4.8.1/cuda75/1.10.1

cd $SLURM_SUBMIT_DIR

srun --mpi=pmi2 MPI_FD.exe 10000
