#!/bin/bash
#SBATCH --job-name=TOUA_filtering
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=16 
#SBATCH --mem=10GB 
#SBATCH --time=03:00:00 
#SBATCH --output=TOUA_af_filter

#Peter Laurin
#Script to generate snp position subsets for comparison between two mapping
#populations, as well as generate new vcfs to be used in future analyses
#
#run on NYU Greene, or some hpc with SLURM-managed tasks 
#one node, 16 cpus, 10 gb memory, < 1 hr run time 
#could be run locally, with longer runtime, with no --threads argument
#
#TOUA_genotypes.vcf.gz from Frachon et al. 2017, available at https://lipm-browsers.toulouse.inra.fr/pub/Frachon2017-NEE/
#1001genomes_snp-short-indel_only_ACGTN.vcf.gz from the 1001 Genomes data, available at https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz

module purge
module load bcftools/intel/1.11


bcftools +fill-tags TOUA_genotypes.vcf.gz --threads 16 -Ob -o TOUA_with_af.bcf -- -t MAF,AF

bcftools query -f'%CHROM\t%POS\t%MAF\t%AF\n' TOUA_with_af.bcf > TOUA_pos_and_af.txt
bcftools query -f'%CHROM\t%POS\t%MAF\t%AF\n' -i'FILTER="PASS"' TOUA_with_af.bcf > TOUA_pass_pos_and_af.txt

bcftools +fill-tags 1001genomes_snp-short-indel_only_ACGTN.vcf.gz  --threads 16 -Ob -o 1001_with_af.bcf -- -t MAF,AF
bcftools query -f'%CHROM\t%POS\t%MAF\t%AF\n' 1001_with_af.bcf > 1001_pos_and_af.txt
bcftools query -f'%CHROM\t%POS\t%MAF\t%AF\n' -i'TYPE="snp"' 1001_with_af.bcf > 1001_snps_pos_and_af.txt


bcftools view -m 2 -M 2 --threads 16 -i 'TYPE="snp"' -Ob -o 1001_biallelic_only_with_af.bcf 1001_with_af.bcf
bcftools query -f'%CHROM\t%POS\t%MAF\t%AF\n' 1001_biallelic_only_with_af.bcf > 1001_biallelic_only_pos_and_af.txt


bcftools view --threads 16 -o 1001_biallelic_only_with_af.vcf 1001_biallelic_only_with_af.bcf 
bcftools view --threads 16 -o TOUA_pass.vcf -i'FILTER="PASS"' TOUA_with_af.bcf

