#Author:Peter Laurin
#date: Aug 23, 2021

# Script to plot minor associations near GS-OX1 - Figure S6
# requires associations files "lmm_maf03_ratio11.assoc.fixed.txt", 
# "lmm_maf03_PC2alilong.assoc.fixed.txt", and
# "lmm_maf03_BlupG4MSB.assoc.fixed.txt")

library(tidyverse)

get_GWAS_data <- function(path){
  full_gwas <- read_tsv(path)
  pval03 <- 0.05 / nrow(full_gwas)
  pval05 <- 0.05 / nrow(full_gwas %>% filter(af > 0.05))
  #100 KB window to either side 
  reduced_gwas <- full_gwas %>% filter(chr == 1 & ps < 24587900 & ps > 24387900)
  LD_data <- read_tsv("LD_comm_snp.list.geno.ld") %>% 
    rename(r2 = `R^2`, ps = POS2) %>% select(r2, ps)
  plot_data <- inner_join(reduced_gwas, LD_data, by = 'ps')
  
  #missing from LD analysis - it is the calculated position
  missing <- reduced_gwas %>%  filter(ps == 24487900) %>% mutate(r2 = NA)
  plot_data <- rbind(plot_data, missing)
  
  
  plot_data <- plot_data %>% mutate(log_p = -log10(p_wald))
  plot_data$log_p[plot_data$log_p > 20] <- 20
  
  
  plot_data <- plot_data %>% arrange(r2, af)
  plot_data <- plot_data %>% 
    mutate(af_bin=cut(af, breaks=c(0.03, 0.05, 0.1, 1), 
                      labels=c("0.03-0.05", "0.05-0.1", ">0.1")))
  plot_data$three <- pval03
  plot_data <- plot_data %>% select(ps, log_p, r2, af_bin, three)
  return(plot_data)
  
}

gwases <- c("lmm_maf03_ratio11.assoc.fixed.txt", 
              "lmm_maf03_PC2alilong.assoc.fixed.txt", 
              "lmm_maf03_BlupG4MSB.assoc.fixed.txt")

full_df <- tibble()

for(i in gwases){
  i_data <- get_GWAS_data(i)
  type <- substr(i, 11, nchar(i)) %>% str_split(pattern = '\\.', simplify = T)
  type <- type[,1]
  if(type == "BlupG4MSB"){
    type = "4mSOb"
  }
  if(type == "PC2alilong"){
    type = "long-chain GSLs (PC2)"
  }
  if(type == "ratio11"){
    type = "7mSh : 7mSOh"
  }
  i_data$type = type
  full_df <- rbind(full_df, i_data)
}

ggplot(full_df, aes(ps, log_p, fill = r2, size = af_bin)) + 
  geom_point(shape = 21) + 
  scale_fill_viridis_c() + 
  geom_hline(yintercept = -log10(full_df$three), color = 'indianred') + 
  scale_size_manual(values = c(2, 2.75, 3.5)) + 
  xlim(24434291,24534291) +
  theme_classic() +
  facet_wrap(~type) + 
  xlab("position on chromosome 1") + 
  ylab("-log10(P)") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))#, legend.position = 'none')


