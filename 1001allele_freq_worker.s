#!/bin/bash
#SBATCH --job-name=allele_freq_sampling_array
#SBATCH --array=1-22
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=4 
#SBATCH --mem=10GB 
#SBATCH --time=1:00:00 

module purge
module load intel/19.1.2 r/intel/4.0.4

Rscript genAlleleFrequencies.R x$SLURM_ARRAY_TASK_ID.vcf test_allele_freqs$SLURM_ARRAY_TASK_ID.frq

echo Task $SLURM_ARRAY_TASK_ID done!
