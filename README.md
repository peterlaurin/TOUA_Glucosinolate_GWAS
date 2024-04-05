# TOUA_Glucosinolate_GWAS

Script (and partial data) repository for Gloss et al. 2021 [https://doi.org/10.1098/rstb.2020.0512]("Genome-wide association mapping within a local *Arabidopsis thaliana* population more fully reveals the genetic architecture for defensive metabolite diversity")

Files included in repo:

### Population Genomic Analysis 

1. **vcf_subset.s** - subset 1001G and TOUA for snp comparison, future analyses

2. **snp_comp.R** - short R script demonstrating SNP comparisons used in paper

3. **1001_admixture.s** - generate pruned marker set, run ADMIXTURE, and extract snp windows (1001 Genomes)

	- **1001_batch_admixture.s** - array job script to run ADMIXTURE (K=1-15)

	- **1001_pop_strat_windows.sh** - extract pruned marker set for 100KB regions around GWAS-inferred loci

4. **TOUA_admixture.s** - generate pruned marker set, run ADMIXTURE, and extract snp windows (TOUA) 

	- **TOUA_batch_admixture.s** - array job script to run ADMIXTURE (K=1-15)

	- **TOUA_pop_strat_windows.sh** - extract pruned marker set for 100KB regions around GWAS-inferred loci

5. **1001allele_freq_boss.sh** - runs several helper scripts to downsample alleles and calculate af (1001 Genomes)

	- **generate_1001_subset.R** - R script which generates random sample of 1001G accessions in Europe (n = num in TOUA)
	
	- **1001allele_freq_worker.s** - array job script to run genAlleleFrequencies.R on vcf part

	- **genAlleleFrequencies.r** - R script which samples from smaller allele pool, calculates af

6. **TOUA_allele_freq_boss.sh** - runs several helper scripts to downsample alleles and calculate af (TOUA)

	- **TOUA_allele_freq_worker.s** - array job script to run genAlleleFrequencies.R on vcf part

	- **genAlleleFrequencies.r** - R script which samples from smaller allele pool, calculates af

7. **calculate_TD.s** - calculates Tajima's D for downsampled alleles in 1001 and TOUA

	- **fake_genotypes_from_freqs.pl** - creates vcf-like genotype data from downsampled allele frequency files

	- **1001_header.txt** - simplified vcf header for 1001 Genomes vcf

	- **TOUA_header.txt** - simplified vcf header for TOUA vcf

8. **geo_boss.sh** - Script to run array, parallelize distance calculations between accessions in Katz et al.

	- **ed_calc_worker.s** - array job script which runs part of geo_fig_batch.R

	- **geo_fig_batch.R** - R script to calculate (part) of Haversine distances between accessions



### Population Genomic Plots

1. **geo_plot.R** - R script to generate figure 1a - proportion of non-matching GSL chemotypes

2. **alleleFreq_and_TD_plots.R** - R script to generate figure 1b, 1c - allele frequency spectrum and Tajima's D

3. **population_stratifcation_plots.R** - R script to generate figure S7 - divergence in allele frequency...


### Glucosinolate Profiling

1. **functions_metaboliteQuantification.R** - Defines and describes custom R functions, which utilize various functions within xcms and related packages, for the peak quantification pipeline implemented in downstream scripts

2. **script_metaboliteQuantification.R** - Implementation of the peak quantification pipeline

3. **script_representativeEicPlot.R** - Production of a representative EIC plot, included as Fig. S2a

4. **script_compareSliceAndPeakPicking.R** - Production of plots comparing xcms peak integration and custom “slice” integration, included as Supplementary Note Fig. 2a

### Phenotypic Analyses and GWAS Prep (TOU-A)

1. **script_TOU_heritability.R** - Load raw phenotype data, fit LMMs with covariates to get heritability estimates and per-accession BLUPs for each glucosinolate molecule

2. **script_TOU_prepForGemma.R** - Reformat BLUPs to prep input phenotype file for downstream GWAS using GEMMA

3. **script_TOU_BLUPs2ratios4Gemma.R** - Transform BLUPs back to original (linear) scale phenotypes, then calculate log2(ratios) of precursor:product abundances and prep input phenotype file for downstream GWAS using GEMMA

### Phenotypic Analyses and GWAS Prep (other panels: Brachi, Katz, Wu datasets)

1. **script_\<DatasetName\>_parseGlucPhenotypes.R** - parse raw data from the published dataset indicated by \<DatasetName\> (for all datasets), fit LMMs with covariates to get heritability estimates and per-accession BLUPs for each glucosinolate molecule (for Brachi and Katz only)

2. **script_\<DatasetName\>_prepForGemma.R** - Reformat BLUPs (Brachi, Katz) or “means” (i.e., the single technical replicate reported per sample by Wu, which can be considered a mean of all the biological replicates pooled within in) to prep input phenotype file for downstream GWAS using GEMMA

3. **script_\<DatasetName\>_prepForGemmaRatios.R** - Transform BLUPs back to original (linear) scale phenotypes, then calculate log2(ratios) of precursor:product abundances and prep input phenotype file for downstream GWAS using GEMMA

### Running GWAS

1. **script_GWAS_kinship.txt**  - Kinship calculated for GWAS using GEMMA

2. **script_GWAS_lmm.txt** - Single-trait GWAS using GEMMA (LMM)

3. **script_GWAS_mvlmm.txt** - Multi-trait GWAS using GEMMA (mvLMM)

### Parsing / Analyzing GWAS 

1. **GS_OX1_plots.R** - R script to generate Figure S6 - Associations at a minor variant near GS-OX1 

### Phenotypic Comparisons

1. **script_corrPlots.R** - Creates plot showing correlation between per-accession BLUPs for each glucosinolate molecule within a given study (i.e., within TOU-A, Brachi, Katz, or Wu datasets).

2. **script_compareHeritability.R** - Creates plot comparing heritability values per molecule in TOU-A vs. Katz and Brachi.

### Post-GWAS Analysis

1. **script_quantifyGenomicInflation.R** - Calculate genomic inflation for each GWAS based on median P-values, and also extract observed P-values at different percentile rankings of the p-value distribution.

2. **script_defineLocusBoundaries.R** - Define locus boundaries around each glucosinolate biosynthetic gene. These extended windows will be considered the "gene regions" for downstream analysis of GWAS output, e.g. to associate significant SNPs with biosynthetic genes.

3. **script_getSnpsNearCandidateGenes.R** - Determine candidate SNPs by assigning SNPs to glucosinolate biosynthetic genes, if SNPs are within boundaries as specified in script_defineLocusBoundaries.R.

4. **script_combineGwasWideFormat.R** - Merges all GWAS output files within specified directories so that p-values for each GWAS are listed in the same table (rows = SNPs, columns = GWAS analyses)

5. **script_countSigSNPsPerLocus.R** - Tallies up the number of significant SNPs and the top p-value by candidate gene for each trait in a the merged p-value table created by script_combineGwasWideFormat.R.

6. **script_topSNPbarplots.R** - Creates barplot comparing top p-value per gene across each dataset (Tou-A, Brachi, Katz, Wu) as shown in supplementary figure.

7. **script_EffectSizeComparison_GWAS.R** - Creates heatmap of GWAS effect sizes in the TOU-A population for the leading SNP at each locus for each glucosinolate molecule.

8. **script_effectSizeComparison_bcat3.R** - Creates plot comparing effects of the BCAT3 mutant with GWAS effect sizes for the leading SNP at the the BCAT3 locus 

9. **script_ManhattansForFigures.R** - Create Manhattan plots showing top p-value per SNP across all combined traits (glucosinolate molecules) within a molecule class (e.g., aliphatic glucosinolates), highlighting SNPs inside candidate gene windows.

10. **script_qqPlots.R** - Creates qq-plots for selected GWAS shown in supplementary material.

11. **script_local_manhattan.R** - Creates local Manhattan plots for the candidate regions harboring significant SNPs as shown in supplementary material.

12. **script_PVE_Table.R** - Creates table of percent variance explained (PVE) for leading SNPs.

### Input/Output Files

Input and output files from some of the scripts that are not in the Dryad repository or straightforward to obtain in the proper format from the sources cited in the manuscript. Large files, like modified GWAS output, are either in the Dryad repository directly or can be regenerated from files in the repository.
