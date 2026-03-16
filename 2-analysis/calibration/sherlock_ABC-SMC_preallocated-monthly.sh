#!/bin/sh
#                     # lines starting with #SBATCH is an instruction to the job scheduler
#SBATCH --job-name=ABCmonthly	        # Job name
#SBATCH --mail-type=ALL         	    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kbubar@stanford.edu # Where to send mail  
#SBATCH --nodes=1                   	# Use one node
#SBATCH --ntasks-per-node=16          # Number of cores per node
#SBATCH --time=48:00:00            	  # Time limit days-hrs:min:sec
#SBATCH --output=array_%A-%a.out    	# Standard output and error log
#SBATCH --array=1               	    # Array range
#SBATCH --mem=200G                    # before used 220G, OR can use mem-per-cpu = 8G per CPU core

date
hostname

ml R/4.4.2

Rscript run-ABC-SMC-preallocatedvax-monthly.R ${SLURM_ARRAY_TASK_ID}

wait