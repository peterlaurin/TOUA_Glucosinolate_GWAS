#!/bin/bash
#Peter Laurin
#
#script to generate downsampled allele frequencies 
#runs other scripts, including  genAlleleFrequencies.R and batch array TOUA_allele_freq_worker.s
#
#run on NYU Greene, to submit batch array job 
#
#TOUA_pass.vcf generated from vcf_subset.s 


#gen working directory, split file into 100,000bp chunks 
mkdir split_97
cd split_97/
split --numeric-suffixes=1 -l 100000 ../TOUA_pass.vcf
head -n 43 x01 > header.txt
tail -n +44 x01 > tmp
mv tmp x01
for file in x*; do cat header.txt $file > $file.vcf; rm $file; done
for file in x0*; do REAL_NAME=`echo $file | tr -d 0`; mv $file $REAL_NAME; done;

#run array job to calculate parts of allele freqs using genAlleleFrequencies.R
cp ../TOUA_allele_freq_worker.s ./
cp ../genAlleleFrequencies.R ./
sbatch TOUA_allele_freq_worker.s

#wait for batch jobs to finish
sleep 10m
cat test_allele_freqs*.frq > TOUA_pass.frq 
cd ..
mv split_97/TOUA_pass.frq ./
rm -R split_97


