#!/bin/bash
#SBATCH --job-name=final_filtering1001
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=12 
#SBATCH --mem=10GB 
#SBATCH --time=8:00:00

#Peter Laurin
#
#Script to generate European subset (From Katz et al. 2021) 
#and prune resulting snps into marker set in linkage equilibrium 
#and 0.03 maf 
#
# 1001genomes_snp-short-indel_only_ACGTN.vcf from 1001 Genomes (https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz)
# KatzAccessions_1001_metadata.csv is 1001 Genomes metadata (https://1001genomes.org/accessions.html) subsetted to CS Numbers from Katz et al. 2021 sheet "D GSLs_data enmeans" (https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNjc3ODQvZWxpZmUtNjc3ODQtc3VwcDEtdjIueGxzeA--/elife-67784-supp1-v2.xlsx?_hash=DDyP7JsX8%2FBjy610ZvQD%2BgQAMe3BXtizrPYYHdByTmE%3D)  
#
#run on NYU Greene or some other slurm environment
#one node, 12 cpus, 10 gb memory,  2 hr run time 

module purge
module load intel/19.1.2 vcftools/intel/0.1.16 plink/1.90b6.21 


#prepare file to be 'plinkable'
tail -n +2 KatzAccessions_1001_metadata.csv | cut -f 1 > katz_accessions.txt 
vcftools --vcf ../../1001_genomes/1001genomes_snp-short-indel_only_ACGTN.vcf --keep katz_accessions.txt --recode --recode-INFO-all --out 1001_genomes_KATZ

pigz -p 12 1001_genomes_KATZ.recode.vcf 

#prune data set to 0.8 r^2, 0.03 maf
plink --allow-extra-chr --double-id --indep-pairwise 50 10 0.8 --maf 0.03 --out 1001_genomes_03_8r --set-missing-var-ids @:# --vcf 1001_genomes_KATZ.recode.vcf.gz 
plink --allow-extra-chr --double-id --extract 1001_genomes_03_8r.prune.in --make-bed --out 1001_pca_03_8r --pca --set-missing-var-ids @:# --vcf 1001_genomes_KATZ.recode.vcf.gz


#rewrite marker set to vcf to find 100KB window around snp
plink --bfile 1001_pca_03_8r --out 1001_genomes_LD_filter --recode vcf

#run admixture batch, from K=1-15
sbatch 1001_batch_admixture.s

sleep 60m 

#make folder to generate plots 
mkdir 1001
cd 1001
for i in {1..15}; do cp admixture$i/*.P ./


#find snps around window
cp ../1001_pop_strat_windows.sh ./
bash 1001_pop_strat_windows.sh

#generate positions for P files 
tail -n +11 1001_genomes_LD_filter.vcf | cut -f 3 > SNP_ids.txt

cd ..

