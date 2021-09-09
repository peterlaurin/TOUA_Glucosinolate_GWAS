#!/bin/bash
#SBATCH --job-name=ADMIXTURE_RegMap
#SBATCH --array=7-8,10-15
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=12 
#SBATCH --mem=10GB 
#SBATCH --time=8:00:00 

module purge
module load intel/19.1.2 admixture/1.3.0

mkdir admixture_$SLURM_ARRAY_TASK_ID
cd admixture_$SLURM_ARRAY_TASK_ID

admixture ../../plink/1001_pca_03_8r.bed $SLURM_ARRAY_TASK_ID --cv -j12 --haploid="*"
echo Task $SLURM_ARRAY_TASK_ID done!
