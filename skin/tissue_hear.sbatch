#!/bin/bash

#SBATCH --job-name=my_first_snow_run
#SBATCH --nodes=4
#SBATCH --exclusive

module load R/3.2

# the openmpi module is not loaded by default with R
module load openmpi/1.8+intel-15.0

# Always use -n 1 for the snow package. It uses Rmpi internally to spawn
# additional processes dynamically
mpirun -np 1 Rscript two_step_succ_vrest.R