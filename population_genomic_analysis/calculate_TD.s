#!/bin/bash
# 
#SBATCH --job-name=CalulcateTD
#SBATCH --tasks-per-node=1 
#SBATCH --time=2:00:00
#SBATCH --mem=10GB


# Peter Laurin
#
# This script takes downsampled allele frequencies from 1001 Genomes
# and TOUA mapping populations and calculates Tajima's D values
# by making a 'vcf' of 1's and 0's (see script_FakeGenotypesFromFreqs.pl)
# and integer sample names but real positions and using vcftools
#
# 1001_subset.frq generated from 1001allele_freq_boss.sh
# 1001_subset.vcf generated from 1001allele_freq_boss.sh 
# TOUA_pass.frq generated from TOUA_allele_freq_boss.sh
# TOUA_pass.vcf generated from vcf_subset.s 
# 

module purge
module load perl/intel/5.32.0 vcftools/intel/0.1.16


#create fake vcf for 1001 genomes
perl fake_genotypes_from_freqs.pl 1001_subset.frq 1001_downsampled.vcf.tmp 100
tail -n +25 1001_subset.vcf | cut -f 1-7 > 1001_cols_start.tmp
line_count_1001=`wc -l 1001_cols_start.tmp | cut -d ' ' -f 1`
for i in $( seq 1 $line_count_1001 ); do echo "DP=15000" >> depth.tmp; echo "GT" >> geno.tmp; done;
paste 1001_cols_start.tmp depth.tmp geno.tmp > 1001_cols.tmp
paste 1001_cols.tmp 1001_downsampled.vcf.tmp > 1001_data.tmp 
cat 1001_header.txt 1001_data.tmp > 1001_downsampled.vcf
rm *.tmp

#create fake vcf for TOUA
perl fake_genotypes_from_freqs.pl TOUA_pass.frq TOUA_downsampled.vcf.tmp 100
tail -n +44 TOUA_pass.vcf | cut -f 1-7 > TOUA_cols_start.tmp
line_count_TOUA=`wc -l TOUA_cols_start.tmp | cut -d ' ' -f 1`
for i in $( seq 1 $line_count_TOUA ); do echo "DP=15000" >> depth.tmp; echo "GT" >> geno.tmp; done;
paste TOUA_cols_start.tmp depth.tmp geno.tmp > TOUA_cols.tmp
paste TOUA_cols.tmp TOUA_downsampled.vcf.tmp > TOUA_data.tmp
cat fake_vcf/TOUA_header.txt TOUA_data.tmp > TOUA_downsampled.vcf
rm *.tmp

#Calculate Tajima's D - windows of 50000 bp
vcftools --vcf 1001_downsampled.vcf --TajimaD 50000 --out 1001_TD_50000
vcftools --vcf TOUA_downsampled.vcf --TajimaD 50000 --out TOUA_TD_50000



