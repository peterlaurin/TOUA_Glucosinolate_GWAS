# Purpose: Tallies up the number of significant SNPs and the top p-value by candidate gene
#          for each trait (column) in a the merged p-value table created by script_combineGwasWideFormat.R.
#          The output includes tables of the number of significant SNPs and top p-values per
#          glucosinolate candidate gene, and also a table summarizing per study (Brachi, Katz, WU, TOU-A)
#          rather than per trait.


# specify maf cutoff for the analysis
maf.cutoff = 0.05

##### DEFINE SNP LIST, CANDIDATE GENES #####

# load SNPs, maf > 0.05, TOU-A
snps.toua        = subset(read.delim("/output_gwas_assocs_lmm/tou_miss10/blup_aliphatic_assoc/lmm_maf03_BlupG2H3B.assoc.fixed.txt", sep="\t")[,c("chr","ps","af")], af > 0.05)
# load SNPs, maf > 0.05, 1001 Genomes + RegMap imputed
snps.brachi.full = subset(read.delim("/output_gwas_assocs_lmm/brachi/blup_aliphatic_assoc/lmm_maf03_BlupG2P.assoc.fixed.txt", sep="\t")[,c("chr","ps","af")], af > 0.05)
snps.brachi.n192 = subset(read.delim("/output_gwas_assocs_lmm/brachi_n192/blup_aliphatic_assoc/lmm_maf03_N192_BlupG2P.assoc.fixed.txt", sep="\t")[,c("chr","ps","af")], af > 0.05)
snps.katz.full   = subset(read.delim("/output_gwas_assocs_lmm/katz/blup_aliphatic_assoc/lmm_maf03_BlupG2P.assoc.fixed.txt", sep="\t")[,c("chr","ps","af")], af > 0.05)
snps.katz.n192   = subset(read.delim("/output_gwas_assocs_lmm/katz_n192/blup_aliphatic_assoc/lmm_maf03_N192_BlupG2P.assoc.fixed.txt", sep="\t")[,c("chr","ps","af")], af > 0.05)
snps.wu.full     = subset(read.delim("/output_gwas_assocs_lmm/wu/blup_aliphatic_assoc/lmm_maf03_G2Ppm.assoc.fixed.txt", sep="\t")[,c("chr","ps","af")], af > 0.05)
snps.wu.n192     = subset(read.delim("/output_gwas_assocs_lmm/wu_n192/blup_aliphatic_assoc//lmm_maf03_N192_G2Ppm.assoc.fixed.txt", sep="\t")[,c("chr","ps","af")], af > 0.05)

# get full list of SNPs
snps.all         = rbind(snps.brachi.full[,c("chr","ps")], snps.brachi.n192[,c("chr","ps")], 
                         snps.katz.full[,c("chr","ps")], snps.katz.n192[,c("chr","ps")],
                         snps.wu.full[,c("chr","ps")], snps.wu.n192[,c("chr","ps")], 
                         snps.toua[,c("chr","ps")])

# rs from vcf sometimes has chr first (tou), sometimes ps first (regmap/1001 imputed), so standardize
snps.all$rs      = paste(snps.all$chr, snps.all$ps, sep = "_")
snps.all         = unique(snps.all)

# # save SNPs list so it can be loaded directly later
# write.csv(snps.all, "full_SNP_list_maf05_wTouFrachonMiss10.csv")
# # to load:
# snps.all = read.csv("full_SNP_list_maf05_wTouFrachonMiss10.csv")

# load candidate gene list
candidate.genes  = read.delim("pathway_BiosyntheticGenes_HarunTable2.txt", sep = "\t")
candidate.genes  = subset(candidate.genes, aliphatic == "yes" | indolic == "yes")
nrow(candidate.genes) # 45 candidate genes

# load table of gene models (left and right boundaries)
gene.models      = read.delim("TAIR10_GFF3_geneBoundariesTableADG_extendedWindows.txt", sep = "\t")

# merge gene model info with candidate gene list
candidate.genes  = merge(candidate.genes, gene.models, by = "gene_id")
nrow(candidate.genes) # 45 genes still retained with gene model coordinates now included

##### ALIPHATICS - prep #####

    # Optional: subset candidate genes to the significant loci for aliphatic GSL traits
    # GS-OX1, GS-OH, BCAT3, AOP2, MAM1

    #candidate.genes.ali = subset(candidate.genes, gene_id %in% c("AT1G65860","AT2G25450","AT3G49680","AT4G03060","AT5G23010"))

candidate.genes.ali = subset(candidate.genes, aliphatic == "yes")

candidate.snps.per.gene = list()

# get list of candidate SNPs per gene
for (i in 1:nrow(candidate.genes.ali)){

  snps = subset(snps.all, chr == candidate.genes.ali[i,"chr"] &
                  ps > candidate.genes.ali[i,"ps_left"] &
                  ps < candidate.genes.ali[i,"ps_right"],
                select = rs
  )

  candidate.snps.per.gene[[ as.character(candidate.genes.ali[i,"gene_id"]) ]] = snps

}

##### ALIPHATICS - count sig SNPs #####

# specify which sets of GWAS to parse for top p-values and number of significant SNPs
blups.tou.m10  = subset(read.csv("/output_gwas_assocs_lmm/tou_miss10/blup_aliphatic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)
blups.bra.full = subset(read.csv("/output_gwas_assocs_lmm/brachi/blup_aliphatic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)
blups.bra.192  = subset(read.csv("/output_gwas_assocs_lmm/brachi_n192/blup_aliphatic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)
blups.kat.full = subset(read.csv("/output_gwas_assocs_lmm/katz/blup_aliphatic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)
blups.kat.192  = subset(read.csv("/output_gwas_assocs_lmm/katz_n192/blup_aliphatic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)
blups.wu.full  = subset(read.csv("/output_gwas_assocs_lmm/wu/blup_aliphatic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)
blups.wu.192   = subset(read.csv("/output_gwas_assocs_lmm/wu_n192/blup_aliphatic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)

ratios.tou.m10  = subset(read.csv("/output_gwas_assocs_lmm/tou_miss10/ratio_aliphatic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)
ratios.bra.full = subset(read.csv("/output_gwas_assocs_lmm/brachi/ratio_aliphatic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)
ratios.bra.192  = subset(read.csv("/output_gwas_assocs_lmm/brachi_n192/ratio_aliphatic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)
ratios.kat.full = subset(read.csv("/output_gwas_assocs_lmm/katz/ratio_aliphatic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)
ratios.kat.192  = subset(read.csv("/output_gwas_assocs_lmm/katz_n192/ratio_aliphatic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)
ratios.wu.full  = subset(read.csv("/output_gwas_assocs_lmm/wu/ratio_aliphatic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)
ratios.wu.192   = subset(read.csv("/output_gwas_assocs_lmm/wu_n192/ratio_aliphatic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)

gwas.to.parse = list(blups.tou.m10, blups.bra.full, blups.bra.192, blups.kat.full, blups.kat.192, blups.wu.full, blups.wu.192, 
                     ratios.tou.m10, ratios.bra.full, ratios.bra.192, ratios.kat.full, ratios.kat.192, ratios.wu.full, ratios.wu.192)
names(gwas.to.parse) = c("blups.tou.m10", "blups.bra.full", "blups.bra.192", "blups.kat.full", "blups.kat.192", "blups.wu.full", "blups.wu.192", 
                         "ratios.tou.m10", "ratios.bra.full", "ratios.bra.192", "ratios.kat.full", "ratios.kat.192", "ratios.wu.full", "ratios.wu.192")

# create data frame to hold top significant SNP pval
top.pvals.sig.best = data.frame(matrix(nrow = nrow(candidate.genes.ali), ncol = length(gwas.to.parse)+1))
colnames(top.pvals.sig.best) = c("gene_id", names(gwas.to.parse))
top.pvals.sig.best$gene_id = candidate.genes.ali$gene_id


for (i in 1:length(names(gwas.to.parse))){
  
  # print progress update
  print(paste0("Parsing: ", names(gwas.to.parse)[[i]], " ..."))
  
  # specify focal gwas
  focal.gwas = gwas.to.parse[[names(gwas.to.parse)[i]]]
  
  # create data frame to hold top SNP p-values
  # dimensions: rows = number of loci, columns = number of traits 
  top.pvals = data.frame(matrix(nrow = nrow(candidate.genes.ali), ncol = ncol(focal.gwas) - 3))
  colnames(top.pvals) = c("gene_id", colnames(focal.gwas[5:ncol(focal.gwas)]))
  top.pvals$gene_id = candidate.genes.ali$gene_id
  
  # create identical data frame, but for only considering significant SNPs
  top.pvals.sig = top.pvals
  
  # create identical data frame to hold number of sig. SNPs
  num.sig.snps = top.pvals
  
  for (j in top.pvals$gene_id) {
    
    #for (k in 4:ncol(focal.gwas)){
    for (k in colnames(focal.gwas)[5:ncol(focal.gwas)]){
      
      # subset to candidate SNPs for the trait
      sig.snps = focal.gwas[focal.gwas[,k] < maf.cutoff / nrow(focal.gwas) / (ncol(focal.gwas)-5)
                            & focal.gwas$rs %in% candidate.snps.per.gene[[j]]$rs, k]
      candidate.snps = focal.gwas[focal.gwas$rs %in% candidate.snps.per.gene[[j]]$rs, k]
      
      # top p value
      top.pvals[top.pvals$gene_id == j, k] = min(candidate.snps)
      
      # top p value (sig. SNPs only)
      top.pvals.sig[top.pvals$gene_id == j, k] = min(sig.snps)
      
      
      # number of significant SNPs
      num.sig.snps[num.sig.snps$gene_id == j, k] = length(sig.snps)
      
    }
    
  }
  
  top.pvals.sig.best[,names(gwas.to.parse)[i]] = top.pvals.sig[,"p_wald_min"]
  
  # individually name and save these as desired  
  write.csv(top.pvals, paste0("/comparativeGWAS/", "aliphatic.top.pvals.", names(gwas.to.parse)[i], ".csv"))
  write.csv(top.pvals.sig, paste0("/comparativeGWAS/", "aliphatic.top.pvals.sig.", names(gwas.to.parse)[i], ".csv"))
  write.csv(num.sig.snps, paste0("/comparativeGWAS/", "aliphatic.num.sig.snps.", names(gwas.to.parse)[i], ".csv"))
  
  # print data frame of top SNP p-values to file
  
}

write.csv(top.pvals.sig.best, "top_pvals_sig_best_aliphatic_all-wFrachonMiss10.csv")



##### INDOLICS - prep #####

candidate.genes.ind = subset(candidate.genes, indolic == "yes")

candidate.snps.per.gene = list()

# get list of candidate SNPs per gene
for (i in 1:nrow(candidate.genes.ind)){
  
  snps = subset(snps.all, chr == candidate.genes.ind[i,"chr"] & 
                  ps > candidate.genes.ind[i,"ps_left"] &
                  ps < candidate.genes.ind[i,"ps_right"],
                select = rs
  )
  
  candidate.snps.per.gene[[ as.character(candidate.genes.ind[i,"gene_id"]) ]] = snps
  
}

##### INDOLICS - count sig SNPs #####

# specify GWAS to parse for top p-values and number of significant SNPs

blups.tou.m10  = subset(read.csv("/output_gwas_assocs_lmm/tou_miss10/blup_indolic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)
blups.wu.full  = subset(read.csv("/output_gwas_assocs_lmm/wu/blup_indolic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)
blups.wu.192   = subset(read.csv("/output_gwas_assocs_lmm/wu_n192/blup_indolic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)

ratios.tou.m10  = subset(read.csv("/output_gwas_assocs_lmm/tou_miss10/ratio_indolic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)
ratios.wu.full  = subset(read.csv("/output_gwas_assocs_lmm/wu/ratio_indolic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)
ratios.wu.192   = subset(read.csv("/output_gwas_assocs_lmm/wu_n192/ratio_indolic_assoc/snp_table_all_pvals.csv"), af > maf.cutoff)

mvlmm.tou.m10   = subset(read.delim("mvlmmA_maf03.assoc.fixed.txt"), af > maf.cutoff)[,c("rs","chr","ps","af","p_wald")]
mvlmm.tou.m10$p_wald_min = mvlmm.tou.m10$p_wald

# 
gwas.to.parse = list(blups.tou.m10, blups.wu.full, blups.wu.192, ratios.tou.m10, ratios.wu.full, ratios.wu.192, mvlmm.tou.m10)
names(gwas.to.parse) = c("blups.tou.m10", "blups.wu.full", "blups.wu.192", "ratios.tou.m10", "ratios.wu.full", "ratios.wu.192", "mvlmm.tou.m10")

# create data frame to hold top significant SNP pval
top.pvals.sig.best = data.frame(matrix(nrow = nrow(candidate.genes.ind), ncol = length(gwas.to.parse)+1))
colnames(top.pvals.sig.best) = c("gene_id", names(gwas.to.parse))
top.pvals.sig.best$gene_id = candidate.genes.ind$gene_id

for (i in 1:length(names(gwas.to.parse))){
  
  # print progress update
  print(paste0("Parsing: ", names(gwas.to.parse)[[i]], " ..."))
  
  # specify focal gwas
  focal.gwas = gwas.to.parse[[names(gwas.to.parse)[i]]]
  
  # create data frame to hold top SNP p-values
  # dimensions: rows = number of loci, columns = number of traits 
  top.pvals = data.frame(matrix(nrow = nrow(candidate.genes.ind), ncol = ncol(focal.gwas) - 3))
  colnames(top.pvals) = c("gene_id", colnames(focal.gwas[5:ncol(focal.gwas)]))
  top.pvals$gene_id = candidate.genes.ind$gene_id
  
  # create identical data frame, but for only considering significant SNPs
  top.pvals.sig = top.pvals
  
  # create identical data frame to hold number of sig. SNPs
  num.sig.snps = top.pvals
  
  for (j in top.pvals$gene_id) {
    
    #for (k in 4:ncol(focal.gwas)){
    for (k in colnames(focal.gwas)[5:ncol(focal.gwas)]){
      
      # subset to candidate SNPs for the trait
      sig.snps = focal.gwas[focal.gwas[,k] < maf.cutoff / nrow(focal.gwas) / (ncol(focal.gwas)-5)
                            & focal.gwas$rs %in% candidate.snps.per.gene[[j]]$rs, k]
      candidate.snps = focal.gwas[focal.gwas$rs %in% candidate.snps.per.gene[[j]]$rs, k]
      
      # top p value
      top.pvals[top.pvals$gene_id == j, k] = min(candidate.snps)
      
      # top p value (sig. SNPs only)
      top.pvals.sig[top.pvals$gene_id == j, k] = min(sig.snps)
      
      # number of significant SNPs
      num.sig.snps[num.sig.snps$gene_id == j, k] = length(sig.snps)
      
    }
    
  }
  
  # # could individually name and save these as desired  
  write.csv(top.pvals, paste0("/comparativeGWAS/", "indolic.top.pvals.", names(gwas.to.parse)[i], ".csv"))
  write.csv(top.pvals.sig, paste0("/comparativeGWAS/", "indolic.top.pvals.sig.", names(gwas.to.parse)[i], ".csv"))
  write.csv(num.sig.snps, paste0("/comparativeGWAS/", "indolic.num.sig.snps.", names(gwas.to.parse)[i], ".csv"))
  
  # print data frame of top SNP p-values to file
  
  top.pvals.sig.best[,names(gwas.to.parse)[i]] = top.pvals.sig[,"p_wald_min"]
  
}

write.csv(top.pvals.sig.best, "top_pvals_sig_best_indolic_all-wFrachonMiss10.csv")



