# TOUA_Glucosinolate_GWAS

Script (and partial data) repository for Gloss et al. 2021 "Genome-wide association mapping within a single Arabidopsis thaliana population reveals a richer genetic architecture for defensive metabolite diversity"

Files included in repo so far: 

vcf_subset.s - subset 1001G and TOUA for snp comparison, future analyses

snp_comp.R - short R script demonstrating SNP comparisons used in paper

1001_admixture.s - generate pruned marker set for ADMIXTURE (1001 Genomes)

	1001_batch_admixture.s - array job script to run ADMIXTURE (K=1-15)

TOUA_admixture.s - generate pruned marker set for ADMIXTURE (TOUA) 

	TOUA_batch_admixture.s - array job script to run ADMIXTURE (K=1-15)

1001allele_freq_boss.sh - runs several helper scripts to downsample alleles and calculate af (1001 Genomes)

	generate_1001_subset.R - R script which generates random sample of 1001G accessions in Europe (n = num in TOUA)
	
	1001allele_freq_worker.s - array job script to run genAlleleFrequencies.R on vcf part

	genAlleleFrequencies.r - R script which samples from smaller allele pool, calculates af

TOUA_allele_freq_boss.sh - runs several helper scripts to downsample alleles and calculate af (TOUA)

	TOUA_allele_freq_worker.s - array job script to run genAlleleFrequencies.R on vcf part

	genAlleleFrequencies.r - R script which samples from smaller allele pool, calculates af
