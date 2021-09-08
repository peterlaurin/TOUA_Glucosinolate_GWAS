#!/bin/bash
#SBATCH --job-name=subset_1001
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20GB
#SBATCH --time=2:00:00


#Peter Laurin
#
#Script to generate a subset of 1001 genomes to 
#a sample of European accessions (i.e. in Katz et al.)
# n = num samples in TOUA
#script to generate downsampled allele frequencies 
#runs other scripts, including  
#generate_1001_subset.R 
#genAlleleFrequencies.R and 
#batch array 1001allele_freq_worker.s
#
#run on NYU Greene, to submit batch array job 
#
#Need 1001_accession_meta_data.csv (https://1001genomes.org/accessions.html) and
#KatzAccessions.tsv (see 1001_admixture.s)


#generate subset of 1001G - sample using generate_1001_subset.R
module load bcftools/intel/1.11 r/gcc/4.0.3

Rscript --vanilla generate_1001_subset.R
tail -n +2 1001_genomes_subset.tsv | cut -f 1 > 1001_subset_accessions.txt
bcftools view --threads 20 -S 1001_subset_accessions.txt -o 1001_subset.vcf 1001_biallelic_only_with_af.vcf

#downsample alleles

#gen working directory, split file into 100,000bp chunks 
mkdir split_97
cd split_97/
split --numeric-suffixes=1 -l 500000 ../1001_subset.vcf
head -n 24 x01 > header.txt
for file in x*; do cat header.txt $file > $file.vcf; rm $file; done
tail -n +25 x01.vcf > tmp.vcf
mv tmp.vcf x01.vcf
for file in x0*; do REAL_NAME=`echo $file | tr -d 0`; mv $file $REAL_NAME; done;

#run array job to calculate parts of allele freqs using genAlleleFrequencies.R
cp ../1001allele_freq_worker.s ./
cp ../genAlleleFrequencies.R ./
sbatch 1001allele_freq_worker.s

#wait for batch jobs to finish, cat back together
sleep 30m
cat test_allele_freqs*.frq > 1001_subset.frq
cd ..
mv split_97/1001_subset.frq ./
rm -R split_97

