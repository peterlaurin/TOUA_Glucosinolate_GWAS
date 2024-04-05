# Purpose: Merges all GWAS output files within specified directories so that p-values for
#          each GWAS are listed in the same table (rows = SNPs, columns = GWAS analyses)

library("tidyverse") #str_subset, str_which
library("tools") # for file_path_sans_ext()

# note: before running this script, append "_omitted" in place of .txt the filenames of GWAS
#       that should be excluded due to genomic inflation

# define directory containing all subdirectories of GWAS output files
input.dir = "/output_gwas_assocs_lmm/"

# find all sub-directories containing "_assoc"
dirs.to.analyze = list.dirs(input.dir) [str_which(basename(list.dirs(input.dir)), pattern = "_assoc")]
dirs.to.analyze = dirs.to.analyze[dirs.to.analyze != input.dir] # exclude input path itself

# loop through all subdirectories, and create merged p-value table for each subdirectory
# through the steps below:
for (d in dirs.to.analyze){
  
  print( paste0("Begin directory ", d) )
  
  # set working directory
  setwd(d)
  
  # list all files in folder
  gwas.files = list.files(pattern = "assoc.fixed.txt")
  
  # use first loaded file to make data frame to store best p-values for each site,
  # including chromosome, position, allele freq
  print("loading first file...")
  gwas.all = read.delim(gwas.files[1], sep="\t")[,c("chr","ps","af","p_wald")]
  gwas.all$rs = paste(gwas.all$chr, gwas.all$ps, sep = "_")
  gwas.all = gwas.all[,c("rs","chr","ps","af","p_wald")]
  
  # apply minor allele frequency filter
  gwas.all = subset(gwas.all, af > 0.05)
  
  # rename p-value column to indicate trait
  molecule = gsub( "[:.*].*$", "", gsub("^.*_", "", file_path_sans_ext(gwas.files[1])) )
  colnames(gwas.all)[5] = paste0("p_",molecule)
  
  # merge in additional files
  for (i in 2:length(gwas.files)) {
    
    print(paste0("loading file ",i," of ",length(gwas.files),"..."))
    gwas.new = read.delim(gwas.files[i], sep="\t")[,c("chr","ps","p_wald")]
    gwas.new$rs = paste(gwas.new$chr, gwas.new$ps, sep = "_")
    
    # rename p-value column to indicate trait
    molecule = gsub( "[:.*].*$", "", gsub("^.*_", "", file_path_sans_ext(gwas.files[i])) )
    colnames(gwas.new)[3] = paste0("p_",molecule)
    
    # merge with all other trait p-values, assuming all "rs" in i'th file are in first file
    gwas.all = merge(gwas.all,gwas.new[,c("rs",paste0("p_",molecule))], by = "rs")
    
  }
  
  # get minimum observed p-value for each SNP
  gwas.all$p_wald_min = do.call(pmin, gwas.all[,5:ncol(gwas.all)])
  
  # output merged and best p-value SNP table into the same directory
  # housing the individual GWAS files being merged
  write.csv(gwas.all, "snp_table_all_pvals.csv", quote = F, row.names = F)
  
}
