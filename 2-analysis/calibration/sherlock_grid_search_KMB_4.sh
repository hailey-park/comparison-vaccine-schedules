#!/bin/sh
#                     # lines starting with #SBATCH is an instruction to the job scheduler
#SBATCH --job-name=grid3   # Job name
#SBATCH --mail-type=ALL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kbubar@stanford.edu # Where to send mail  
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --mem-per-cpu=1gb           # Memory per processor
#SBATCH --time=10:00:00             # Time limit hrs:min:sec
#SBATCH --output=array_%A-%a.out    # Standard output and error log
#SBATCH --array=600-801:10           # Array range
#SBATCH -n 10                       # Array range



date
hostname

ml R/4.2.0

for i in {0..9}; do
    srun -n 1 Rscript run-grid-search.R $((SLURM_ARRAY_TASK_ID+i)) &
done

wait