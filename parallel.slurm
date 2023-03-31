#!/bin/bash
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3700
#SBATCH --time=48:00:00

# NOTE: For MPI, we change ntasks-per-node and keep  cpus-per-task at 1

# Change into the directory where I submitted the job
cd $SLURM_SUBMIT_DIR 

# Load the modules used when compiling my code
module load GCC/8.3.0
module load intel/2019b

# Note we use srun rather than mpirun here
srun ./pluto

exit 0
