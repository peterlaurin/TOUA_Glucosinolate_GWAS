#Author: Peter Laurin

#Script to downsample alleles in mapping population, and generate new allele
#frequencies based on these subsets. Usually run as partition of larger vcf.




library(tidyverse)
library(vcfR)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1] #usually part of larger vcf
out_name <- args[2] #filename to write allele freqs

downsample_allele_freqs <- function(nAllele, ref_n, min_n, min_min_n, nDraws){
  ##to sample alleles from total allele pool, and generate allele frequencies
  #read error
  if(is.na(nAllele)){
    return(NA)
  }
  #not enough alleles to sample
  if((min_min_n >= 1 & min_n >= 1) | nAllele < nDraws){
    return(NA)
  }
  #return maf not alt af
  if(min_min_n > min_n){
    min_n <- min_min_n
  }
  minor_allele <- min(ref_n, min_n)
  major_allele <- nAllele - minor_allele
  allele_counts <- rhyper(nn = 1, m = minor_allele, n = major_allele, k = nDraws)
  allele_freq <- allele_counts / nDraws
  return(allele_freq)
}

generate_allele_freqs <- function(nAllele, ref_n, min_n, min_min_n, nDraws){
  #read error
  if(is.na(nAllele)){
    return(NA)
  }
  #not enough alleles to sample
  if((min_min_n >= 1 & min_n >= 1) | nAllele < nDraws){
    return(NA)
  }
  #return maf not alt af
  if(min_min_n > min_n){
    min_n <- min_min_n
  }
  freq <- min(ref_n, min_n) / nAllele
  return(freq)
}


main <- function(){
  vcf <- read.vcfR(filename)
  refs <- maf(vcf, element = 1)
  allele_counts <- tibble(nAllele = refs[,1] / 2, 
                         not_genotyped = refs[,2],
                         ref_n = refs[,3] / 2, 
                         min_n = maf(vcf, element = 2)[,3] / 2, 
                         min_min_n = maf(vcf, element = 3)[,3] / 2)
  
  downsampled_freqs <- mapply(FUN = downsample_allele_freqs, 
                              nAllele = allele_counts$nAllele, 
                              ref_n = allele_counts$ref_n,
                              min_n = allele_counts$min_n, 
                              min_min_n = allele_counts$min_min_n, 
                              nDraws = 100)

  
  true_freqs <- mapply(FUN = generate_allele_freqs,
                       nAllele = allele_counts$nAllele,
                       ref_n = allele_counts$ref_n,
                       min_n = allele_counts$min_n, 
                       min_min_n = allele_counts$min_min_n,
                       nDraws = 100)
  
  
  #write (part) of allele freqs - repress column names  
  write_tsv(tibble(downsampled_freqs,true_freqs), file = out_name, col_names = FALSE)
  
}

main()

