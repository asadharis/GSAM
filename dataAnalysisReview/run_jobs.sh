#!/bin/bash
#SBATCH --array=1-50
#SBATCH --time=24:59:00           # time (HH:MM:SS)
#SBATCH --mem=30G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --error=/dev/null
#SBATCH --output=/dev/null
#SBATCH --mail-user=asad.haris5862@gmail.com # Send email updates to you or someone else
#SBATCH --mail-type=ALL         # send an email in all cases (job started, job ended, job aborted)


## ARGUMENTS ARE:
## 1. seed
## 2. name: Name of dataset.
## 3. ncores: Number of cores for parallel process.

if [ -d "output$1" ]
then
	echo "Dir Exists"
else
	mkdir output$1
fi

module load nixpkgs/16.09 gcc/7.3.0 r/3.6.1


export R_LIBS=~/local/R_libs/
R CMD BATCH --no-save --no-restore "--args $SLURM_ARRAY_TASK_ID $1  $SLURM_CPUS_PER_TASK" simulations.R output$1/$SLURM_ARRAY_TASK_ID.Rout

#### SBATCH --account=def-rwplatt   # replace this with your own account
#### SBATCH --ntasks=4              # number of processes
#### SBATCH --mem-per-cpu=2048M      # memory; default unit is megabytes
#### SBATCH --time=0-01:15           # time (DD-HH:MM)

