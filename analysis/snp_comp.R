#Author: Peter Laurin
#Date: Aug 28, 2021
#SNP Comparison for glucosinulate paper

#all data read generated from script vcf_subset.s

library(tidyverse)


#setwd("/.../snp_comparison/")

####SNP_comparison####
TOUA_all <- read_tsv("TOUA_pos_and_af.txt")  %>% mutate(chromosome = substr(chromosome, 6, str_length(chromosome)), tog = paste0(chromosome,"_", position))
TOUA_pass <- read_tsv("TOUA_pass_pos_and_af.txt") %>%  mutate(chromosome = substr(chromosome, 6, str_length(chromosome)), tog = paste0(chromosome,"_", position))

TOUA_three <- TOUA_pass %>% filter(minor_allele_freq > 0.03)

all_1001 <- read_tsv("1001_pos_and_af.txt") %>% mutate(tog = paste0(chromosome,"_", position))
snps_w_multiallelic <- read_tsv("1001_snps_pos_and_af.txt") %>% mutate(tog = paste0(chromosome,"_", position))
biallelic_1001 <- read_tsv("1001_biallelic_only_pos_and_af.txt") %>% mutate(tog = paste0(chromosome,"_", position))

all_1001_three <- all_1001 %>% filter(minor_allele_freq > 0.03)
snps_w_multiallelic_three <- snps_w_multiallelic %>% filter(minor_allele_freq > 0.03)
bialleleic_1001_three <- biallelic_1001 %>% filter(minor_allele_freq > 0.03)

#number of variants in TOUA
nrow(TOUA_all)
nrow(TOUA_pass)
nrow(TOUA_three)


#all variants comparisons
intersect(TOUA_all$tog, all_1001$tog) %>% length()
intersect(TOUA_all$tog, snps_w_multiallelic$tog) %>% length()
intersect(TOUA_all$tog, biallelic_1001$tog) %>% length()
intersect(TOUA_pass$tog, all_1001$tog) %>% length()
intersect(TOUA_pass$tog, snps_w_multiallelic$tog) %>% length()
intersect(TOUA_pass$tog, biallelic_1001$tog) %>% length()

#variants with maf > 0.03
intersect(TOUA_three$tog, all_1001_three$tog) %>% length()
intersect(TOUA_three$tog, bialleleic_1001_three$tog) %>% length()








