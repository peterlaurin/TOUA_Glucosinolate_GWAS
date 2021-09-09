#PCA Analysis of RegMap, TOUA and 1001Genomes collections
#Author:Peter Laurin

# R script to plot variation in admixture groups' estimated ancestral alleles,
# as determined by ADMIXTURE's 'P' files
# all "P" files generated in 1001_admixture.s or TOUA_admixture.s
# all "final_ID_afX.txt" files generated in TOUA_pop_strat_windows.sh or 1001_pop_strat_windows.sh
# SNP_ids.txt files (for TOUA and 1001G) generated in 1001_admixture.s or TOUA_admixture.s


library(tidyverse)

####ancestral_allele_frequencies#####

#setwd(".../ancestral_allele_analysis")

load_snp_set <- function(path, window, K, type){
  wd <- str_split(path, pattern = '/', simplify = T)[,1]
  #choose file = which K
  K_file <- str_split(path, pattern = '/', simplify = T)[,2]
  setwd(wd)
  #load snps in 100KB window around locus, common name
  snps_1 <- read_tsv(paste0(window,"KB/final_ID_af1.txt"), col_types = 'cc') %>% mutate(snp = "GS-OX1")
  snps_2 <- read_tsv(paste0(window,"KB/final_ID_af2.txt"), col_types = 'cc') %>% mutate(snp = "GS-OH")
  snps_3 <- read_tsv(paste0(window,"KB/final_ID_af3.txt"), col_types = 'cc') %>% mutate(snp = "BCAT3")
  snps_4 <- read_tsv(paste0(window,"KB/final_ID_af4.txt"), col_types = 'cc') %>% mutate(snp = "AOP2")
  snps_5 <- read_tsv(paste0(window,"KB/final_ID_af5.txt"), col_types = 'cc') %>% mutate(snp = "MAM1")
  snps_6 <- read_tsv(paste0(window,"KB/final_ID_af6.txt"), col_types = 'cc') %>% mutate(snp = "CYP79F2")
  snps_7 <- read_tsv(paste0(window,"KB/final_ID_af7.txt"), col_types = 'cc') %>% mutate(snp = "CYP83A1")
  snps_8 <- read_tsv(paste0(window,"KB/final_ID_af8.txt"), col_types = 'cc') %>% mutate(snp = "IGMT2")
  snps_9 <- read_tsv(paste0(window,"KB/final_ID_af9.txt"), col_types = 'cc') %>% mutate(snp = "CYP81F4")
  snps_10 <- read_tsv(paste0(window,"KB/final_ID_af10.txt"), col_types = 'cc') %>% mutate(snp = "CYP81F2")
  
  big_window <- rbind(snps_1, snps_2, snps_3, snps_4, snps_5, snps_6, snps_7, snps_8, snps_9, snps_10)
  
  #read in all positions, to match with P files
  all_snps <- read_lines("SNP_ids.txt", skip = 1)
  
  KN_snps <- read_delim(K_file, delim = ' ', col_names = paste0('pop', 1:K))
  KN_snps$ID <- all_snps
  
  final_tb <- left_join(big_window, KN_snps, by = 'ID')
  final_tb <- final_tb %>%  rowwise() %>% mutate(allele_var = var(c_across(pop1:paste0('pop', K))), 
                                                 avg = mean(c_across(pop1:paste0('pop', K))))
  final_tb <- final_tb %>% mutate(coef_var = sqrt(allele_var) / avg)
  final_tb <- final_tb %>% rename(maf = X2)
  final_tb$maf <- str_split(final_tb$maf, pattern = ':', simplify = T)[,2] %>% as.double()
  setwd("..")
  final_tb$type = type
  return(final_tb)
}

load_full_snps <- function(K, window){
  TOUA_set <- load_snp_set(paste0("TOUA/TOUA_03_8r_pca.",K,".P"), window, K, "TOUA")
  T001_set <- load_snp_set(paste0("1001/1001_pca_03_8r.",K,".P"), window, K, "1001")
  full_set <- bind_rows(T001_set, TOUA_set)
  return(full_set)
}

#plot_snps_separately <- function( K, i, window){
#    part_df <- load_full_snps(K, window) %>% filter(snp == i)
#    ggplot(part_df, aes(maf, allele_var, color = type)) + 
#      geom_point(alpha = 0.3) + geom_smooth() + 
#      labs(title = paste("K =", K, "snp =", i, "window =", window))
#}

#lapply(1:5, FUN = plot_snps_separately, K = 5, window=30)

#five_two_hundred <- load_full_snps(5, 200)
five_one_hundred <- load_full_snps(5, 100)
#five_thirty <- load_full_snps(5, 30)
ggplot(five_one_hundred, aes(maf, allele_var, color = type)) + 
  geom_point(alpha = 0.3) + geom_smooth(se = T) + 
  facet_wrap(~snp) + theme_classic() + labs(title = "Significant loci, window = 100 KB") + 
  scale_color_manual(values = c("dodgerblue", "gold1")) + 
  guides(colour=guide_legend(title = "mapping panel", override.aes=list(shape=15, linetype = c(0,0), fill = NA, alpha = 1, size = 3))) + 
  xlab("minor allele frequency") + ylab("among-group variance in allele frequency")



  

  
  
