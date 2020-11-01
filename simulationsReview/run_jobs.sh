#!/bin/bash
#SBATCH --array=1-100
#SBATCH --time=00:30:00           # time (HH:MM:SS)
#SBATCH --mem=500M
#SBATCH --error=/dev/null
#SBATCH --output=/dev/null
#SBATCH --mail-user=asad.haris5862@gmail.com # Send email updates to you or someone else
#SBATCH --mail-type=ALL         # send an email in all cases (job started, job ended, job aborted)


## ARGUMENTS ARE:
## 1. seed
## 2. n: Sample Size
## 3. p: Number of Variables
## 4. Variance of the noise
## 5. scen: Scenario number, we will use 3 default.


module load nixpkgs/16.09 gcc/7.3.0 r/4.0.2


export R_LIBS=~/local/R_libs/
R CMD BATCH --no-save --no-restore "--args $SLURM_ARRAY_TASK_ID $1  100 1 3" simulations.R bla.Rout

#### SBATCH --account=def-rwplatt   # replace this with your own account
#### SBATCH --ntasks=4              # number of processes
#### SBATCH --mem-per-cpu=2048M      # memory; default unit is megabytes
#### SBATCH --time=0-01:15           # time (DD-HH:MM)

