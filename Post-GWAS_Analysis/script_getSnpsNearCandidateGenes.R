# Purpose: assign SNPs to glucosinolate biosynthetic genes, if SNPs are within boundaries as 
#          specified in script_defineLocusBoundaries.R

# load all SNPs interrogated in GWAS, requiring maf > 0.05, from analysis with each mapping panel
# (since all traits were phenotyped in the same individuals of a given mapping panel, the GEMMA GWAS output file for any trait is acceptable)
snps.toua        = subset(read.delim("/output_gwas_assocs_lmm/toua_FrachonGenoMiss10/blup_aliphatic_assoc/lmm_maf03_BlupG2H3B.assoc.fixed.txt", sep="\t")[,c("chr","ps","af")], af > 0.05)
snps.brachi.full = subset(read.delim("/output_gwas_assocs_lmm/brachi_full/blup_aliphatic_assoc/lmm_maf03_BlupG2P.assoc.fixed.txt", sep="\t")[,c("chr","ps","af")], af > 0.05)
snps.brachi.n192 = subset(read.delim("/output_gwas_assocs_lmm/brachi_N192/blup_aliphatic_assoc/lmm_maf03_N192_BlupG2P.assoc.fixed.txt", sep="\t")[,c("chr","ps","af")], af > 0.05)
snps.katz.full   = subset(read.delim("/output_gwas_assocs_lmm/katz_full/blup_aliphatic_assoc/lmm_maf03_BlupG2P.assoc.fixed.txt", sep="\t")[,c("chr","ps","af")], af > 0.05)
snps.katz.n192   = subset(read.delim("/output_gwas_assocs_lmm/katz_N192/blup_aliphatic_assoc/lmm_maf03_N192_BlupG2P.assoc.fixed.txt", sep="\t")[,c("chr","ps","af")], af > 0.05)
snps.wu.full     = subset(read.delim("/output_gwas_assocs_lmm/wu_full/blup_aliphatic_assoc/lmm_maf03_G2Ppm.assoc.fixed.txt", sep="\t")[,c("chr","ps","af")], af > 0.05)
snps.wu.n192     = subset(read.delim("/output_gwas_assocs_lmm/wu_N192/blup_aliphatic_assoc//lmm_maf03_N192_G2Ppm.assoc.fixed.txt", sep="\t")[,c("chr","ps","af")], af > 0.05)

# get full list of SNPs across all mapping panels
snps.all         = rbind(snps.brachi.full[,c("chr","ps")], snps.brachi.n192[,c("chr","ps")], 
                         snps.katz.full[,c("chr","ps")], snps.katz.n192[,c("chr","ps")],
                         snps.wu.full[,c("chr","ps")], snps.wu.n192[,c("chr","ps")], 
                         snps.toua[,c("chr","ps")])

# standardize SNP id ("rs") as chr_position for all SNPs
# (the rs from the vcf sometimes has chr first (tou), sometimes ps first (regmap/1001 imputed SNP set))
snps.all$rs      = paste(snps.all$chr, snps.all$ps, sep = "_")
snps.all         = unique(snps.all)

# load candidate gene list
candidate.genes  = read.delim("pathway_BiosyntheticGenes_HarunTable2.txt", sep = "\t")
candidate.genes  = subset(candidate.genes, aliphatic == "yes" | indolic == "yes")
nrow(candidate.genes) # 45 candidate genes

# load table of gene models (left and right boundaries)
gene.models      = read.delim("TAIR10_GFF3_geneBoundariesTableADG_extendedWindows.txt", sep = "\t")

# merge gene model info with candidate gene list
candidate.genes  = merge(candidate.genes, gene.models, by = "gene_id")
nrow(candidate.genes) # 45 genes still retained with gene model coordinates now included

### ALIPHATICS

# subset candidate genes to those affecting aliphatic GSL biosynthesis
candidate.genes.ali = subset(candidate.genes, aliphatic == "yes")
nrow(candidate.genes.ali) # 33 genes

candidate.snps.ali = vector()

# get SNPs within boundaries specified for each candidate gene
for (i in 1:nrow(candidate.genes.ali)){
  
  hits = subset(snps.all, chr == candidate.genes.ali[i,"chr"] & 
                          ps > candidate.genes.ali[i,"ps_left"] &
                          ps < candidate.genes.ali[i,"ps_right"]
               )
  
  candidate.snps.ali = c(candidate.snps.ali, as.vector(hits$rs))
  
  print(paste0("gene: ", candidate.genes.ali[i,"gene_id"], " - #", i, " of ", nrow(candidate.genes.ali)))
  print(paste0("hits: ", nrow(hits), " SNPs"))
  print(head(hits))

}

candidate.snps.ali = unique(candidate.snps.ali)

### INDOLICS

# subset candidate genes to those affecting indolic GSL biosynthesis
candidate.genes.ind = subset(candidate.genes, indolic == "yes")
nrow(candidate.genes.ind) # 33 genes

candidate.snps.ind = vector()

# get SNPs within boundaries specified for each candidate gene
for (i in 1:nrow(candidate.genes.ind)){
  
  hits = subset(snps.all, chr == candidate.genes.ind[i,"chr"] & 
                  ps > candidate.genes.ind[i,"ps_left"] &
                  ps < candidate.genes.ind[i,"ps_right"]
  )
  
  candidate.snps.ind = c(candidate.snps.ind, as.vector(hits$rs))
  
  print(paste0("gene: ", candidate.genes.ind[i,"gene_id"], " - #", i, " of ", nrow(candidate.genes.ind)))
  print(paste0("hits: ", nrow(hits), " SNPs"))
  print(head(hits))
  
}

candidate.snps.ind = unique(candidate.snps.ind)

# write candidate SNP lists to files
write(candidate.snps.ali, "candidate_snps_aliphatic_30kb_majorLociExtended_maf05.txt")
write(candidate.snps.ind, "candidate_snps_indolic_30kb_majorLociExtended_maf05.txt")


