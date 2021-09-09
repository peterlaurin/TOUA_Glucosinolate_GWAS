#!/bin/bash
#SBATCH --job-name=ed_calculation_arrays
#SBATCH --array=1-32
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=4 
#SBATCH --mem=10GB 
#SBATCH --time=1:00:00 

# Peter Laurin
# Array script to parallelize distance calculations between accessions
# run by geo_boss.sh




module purge
module load intel/19.1.2 r/gcc/4.0.4

Rscript geo_fig_batch.R $SLURM_ARRAY_TASK_ID ed$SLURM_ARRAY_TASK_ID.txt

echo Task $SLURM_ARRAY_TASK_ID done!
