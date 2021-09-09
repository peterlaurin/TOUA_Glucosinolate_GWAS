#!/bin/bash

# Peter Laurin
#
# Script to extract snps in 100 kb window around GWAS-inferred GSL loci 
# TOUA_LD_fil_new.vcf is rewritten marker set from TOUA_03_8r_pca (see TOUA_admixture.s) 
# Gene positions inferred from (https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff)

module purge
module load vcftools/intel/0.1.16


mkdir 100KB
cd 100KB

#GS-OX1
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 1 --to-bp 24550905 --freq --out var1 --from-bp 24450905
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 1 --to-bp 24550905 --out var1 --from-bp 24450905 --recode

#GS-OH
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 2 --to-bp 10880785 --freq --out var2 --from-bp 10780785
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 2 --to-bp 10880785 --out var2 --from-bp 10780785 --recode

#BCAT3
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 3 --to-bp 18474171 --freq --out var3 --from-bp 18374171
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 3 --to-bp 18474171 --out var3 --from-bp 18374171 --recode

#AOP2
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 4 --to-bp  1402917 --freq --out var4 --from-bp 1302917
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 4 --to-bp  1402917 --out var4 --from-bp 1302917 --recode

#MAM1
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 5 --to-bp 7754994 --freq --out var5 --from-bp 7654994
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 5 --to-bp 7754994 --out var5 --from-bp 7654994 --recode

#CYP79F2
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 1 --to-bp 5656316 --freq --out var6 --from-bp 5556316
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 1 --to-bp 5656316 --out var6 --from-bp 5556316 --recode

#CYP83A1
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 4 --to-bp 8041398 --freq --out var7 --from-bp 7941398
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 4 --to-bp 8041398 --out var7 --from-bp 7941398 --recode

#IGMT2
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 1 --to-bp 7446000 --freq --out var8 --from-bp 7436000
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 1 --to-bp 7446000 --out var8 --from-bp 7436000 --recode

#CYP81F4
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 4 --to-bp 17641846 --freq --out var9 --from-bp 17541846
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 4 --to-bp 17641846 --out var9 --from-bp 17541846 --recode

#CYP81F2
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 5 --to-bp 23238848 --freq --out var10 --from-bp 23138848
vcftools --vcf ../TOUA_LD_fil_new.vcf --chr 5 --to-bp 23238848 --out var10 --from-bp 23138848 --recode

#extract allele freq, position
for i in {1..10}; 
do 
   tail -n +11 var$i.recode.vcf | cut -f 3 > ID_raw_$i.tmp; 
   cat var$i.frq | cut -f 6 > af_raw_$i.tmp; 
   paste ID_raw_$i.tmp af_raw_$i.tmp > final_ID_af$i.txt; 
done;

rm *.tmp

cd ..


