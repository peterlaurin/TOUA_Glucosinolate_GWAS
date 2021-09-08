# TOUA_Glucosinolate_GWAS

Script (and partial data) repository for Gloss et al. 2021 "Genome-wide association mapping within a single Arabidopsis thaliana population reveals a richer genetic architecture for defensive metabolite diversity"

Files included in repo so far:

#### Analysis 

1. vcf_subset.s - subset 1001G and TOUA for snp comparison, future analyses

2. snp_comp.R - short R script demonstrating SNP comparisons used in paper

3. 1001_admixture.s - generate pruned marker set, run ADMIXTURE, and extract snp windows (1001 Genomes)

	- 1001_batch_admixture.s - array job script to run ADMIXTURE (K=1-15)

	- 1001_pop_strat_windows.sh - extract pruned marker set for 100KB regions around GWAS-inferred loci

4. TOUA_admixture.s - generate pruned marker set, run ADMIXTURE, and extract snp windows (TOUA) 

	- TOUA_batch_admixture.s - array job script to run ADMIXTURE (K=1-15)

	- TOUA_pop_strat_windows.sh - extract pruned marker set for 100KB regions around GWAS-inferred loci

5. 1001allele_freq_boss.sh - runs several helper scripts to downsample alleles and calculate af (1001 Genomes)

	- generate_1001_subset.R - R script which generates random sample of 1001G accessions in Europe (n = num in TOUA)
	
	- 1001allele_freq_worker.s - array job script to run genAlleleFrequencies.R on vcf part

	- genAlleleFrequencies.r - R script which samples from smaller allele pool, calculates af

6. TOUA_allele_freq_boss.sh - runs several helper scripts to downsample alleles and calculate af (TOUA)

	- TOUA_allele_freq_worker.s - array job script to run genAlleleFrequencies.R on vcf part

	- genAlleleFrequencies.r - R script which samples from smaller allele pool, calculates af

7. calculate_TD.s - calculates Tajima's D for downsampled alleles in 1001 and TOUA

	- fake_genotypes_from_freqs.pl - creates vcf-like genotype data from downsampled allele frequency files

	- 1001_header.txt - simplified vcf header for 1001 Genomes vcf

	- TOUA_header.txt - simplified vcf header for TOUA vcf

8. geo_boss.sh - Script to run array, parallelize distance calculations between accessions in Katz et al.

	- ed_calc_worker.s - array job script which runs part of geo_fig_batch.R

	- geo_fig_batch.R - R script to calculate (part) of Haversine distances between accessions



#### Plots

1. geo_plot.R - R script to generate figure 1a - proportion of non-matching GSL chemotypes

2. alleleFreq_and_TD_plots.R - R script to generate figure 1b, 1c - allele frequency spectrum and Tajima's D

3. population_stratifcation_plots.R - R script to generate figure S7 - divergence in allele frequency...

