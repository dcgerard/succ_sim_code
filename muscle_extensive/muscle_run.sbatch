#!/bin/bash

#SBATCH --job-name=muscle
#SBATCH --nodes=1
#SBATCH --exclusive

module load R/3.2

# the openmpi module is not loaded by default with R
module load openmpi/1.8+intel-15.0

# Set R libraries to local directory
export R_LIBS_USER=$HOME/R_libs

# Always use -n 1 for the snow package. It uses Rmpi internally to spawn
# additional processes dynamically
mpirun -np 1 Rscript muscle.R
