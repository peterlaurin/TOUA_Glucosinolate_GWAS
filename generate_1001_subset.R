#Peter Laurin 

#Script to generate list of 195 random European accessions 


library(tidyverse)

katz_accessions <- read_tsv("KatzAccessions.tsv")
thousand_accessions <- read_csv("1001_accession_meta_data.csv", col_types = 'dccccddcccfcf')

reduced_accessions <- thousand_accessions %>% filter(CS_Number %in% toupper(katz_accessions$CS))

set.seed(0221)
rand_accessions <- reduced_accessions %>% sample_n(199) 
#TOUA only has 195 not 199 samples - whoops 
write_tsv(rand_accessions[1:195,], file = "1001_genomes_subset.tsv")
