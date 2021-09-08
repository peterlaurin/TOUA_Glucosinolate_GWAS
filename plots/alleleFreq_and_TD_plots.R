#Allele_Frequency_Plots_97
#Author: Peter Laurin
#Date: Aug 5, 2021

# R Script to generate allele frequency plot from downsampled allele frequencies 
# generated in 1001allele_freq_boss.sh and TOUA_allele_freq_boss.sh



library(tidyverse)

#####Allele Frequency Plots#####

#load allele freq. datasets
af_1001 <- read_tsv("1001_subset.frq", col_names = c("downsampled_freqs", "real_freqs"))
af_TOUA <- read_tsv("TOUA_pass.frq", col_names = c("downsampled_freqs", "real_freqs"))

af_TOUA$population = "TOUA"
af_1001$population = "1001"
total_freqs <- rbind(af_1001, af_TOUA)
total_freqs <- pivot_longer(total_freqs, 
                            cols = ends_with("freqs"), 
                            names_to = "freq_type", 
                            values_to = "frequency")

total_freqs <- drop_na(total_freqs)
non_zeros <- total_freqs %>% filter(frequency > 0)
downsampled <- non_zeros %>% filter(freq_type == "downsampled_freqs")
downsampled03 <- downsampled %>% filter(frequency > 0.03)
paper_theme <- theme_classic() + theme(panel.border = element_rect(fill = NA))

#proportion plots
ggplot(downsampled, aes(x = frequency, color = population, stat(density))) + 
  geom_freqpoly(binwidth = 1/100) + xlim(0.00000000001, 0.5) + 
  labs(title = "all minor allele frequencies") + 
  paper_theme +
  scale_color_brewer(palette = "Dark2") + xlab("minor allele frequency") + 
  ylab("percentage of total SNPs in mapping panel") + 
  guides(color=guide_legend(title="mapping panel"))

#maf > 0.03
ggplot(downsampled03, aes(x = frequency, color = population, stat(density))) + 
  geom_freqpoly(binwidth = 1/100) + xlim(0.030000000001, 0.5) + 
  labs(title = "minor alleles with frequency > 0.03") +
  paper_theme +
  scale_color_brewer(palette = "Dark2") + xlab("minor allele frequency") + 
  ylab("percentage of total SNPs in mapping panel") + 
  guides(color=guide_legend(title="mapping panel")) 


#####Tajima's D Plots#####

#50000 windoow
TD_TOUA <- read_tsv("TOUA_TD_50000.Tajima.D")
TD_1001 <- read_tsv("1001_TD_50000.Tajima.D")
TD_TOUA$type = "TOUA"
TD_1001$type = "1001"
full_TD <- rbind(TD_1001, TD_TOUA)

ggplot(data = full_TD, aes(type, TajimaD, fill = type)) + 
  geom_violin() + geom_boxplot(width = 0.2) + 
  labs(title = "Tajima's D, window size = 50000") + theme_classic() + 
  theme(legend.position = "none") + xlab("mapping panel") + ylab("Tajima's D")
