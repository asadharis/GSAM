#!/bin/bash
#SBATCH --array=1-100
#SBATCH --time=05:59:00           # time (HH:MM:SS)
#SBATCH --mem=1G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --error=err/%j.err
#SBATCH --output=err/%j.out
#SBATCH --mail-user=asad.haris5862@gmail.com # Send email updates to you or someone else
#SBATCH --mail-type=BEGIN         # send an email in all cases (job started, job ended, job aborted)


## ARGUMENTS ARE:
## 1. seed
## 2. n: Sample Size $1
## 3. p: Number of Variables $2
## 4. Variance of the noise
## 5. scen: Scenario number, we will use 3 default. $3
## 6. number of cores.


module load nixpkgs/16.09 gcc/7.3.0 r/3.6.1


export R_LIBS=~/local/R_libs/
R CMD BATCH --no-save --no-restore "--args $SLURM_ARRAY_TASK_ID $1 $2 1 $3 $SLURM_CPUS_PER_TASK" simulations.R bla.Rout

#### SBATCH --account=def-rwplatt   # replace this with your own account
#### SBATCH --ntasks=4              # number of processes
#### SBATCH --mem-per-cpu=2048M      # memory; default unit is megabytes
#### SBATCH --time=0-01:15           # time (DD-HH:MM)

