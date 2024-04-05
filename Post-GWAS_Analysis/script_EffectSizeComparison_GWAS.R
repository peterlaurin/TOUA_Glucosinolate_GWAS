# Purpose: Creates heatmap of GWAS effect sizes in the TOU-A population for the leading SNP at each locus
#          on each glucosinolate molecule.

### aliphatic effect sizes

# load summarized results of GWAS
gwas.all.blup  = read.csv("/tou_miss10/blup_aliphatic_assoc/snp_table_all_pvals.csv")
gwas.all.ratio = read.csv("/tou_miss10/ratio_aliphatic_assoc/snp_table_all_pvals.csv")

# load candidate gene positions
gene.positions = read.delim("TAIR10_GFF3_geneBoundariesTableADG_extendedWindows.txt")

# subset to genes of interest that had been found to harbor a significant SNP in at least one GWAS
sig.genes.ali  = c("AT1G16400","AT1G65860","AT2G25450","AT3G49680","AT4G03060","AT4G13770","AT5G23010")
gene.positions.ali = subset(gene.positions, gene_id %in% sig.genes.ali)
gene.positions.ali$top.snp = "NA"

for (i in 1:length(sig.genes.ali)){
  
  gwas.sub.blup  = subset(gwas.all.blup,  chr == gene.positions.ali[i,"chr"] & ps > gene.positions.ali[i,"ps_left"] & ps < gene.positions.ali[i,"ps_right"])
  gwas.sub.ratio = subset(gwas.all.ratio, chr == gene.positions.ali[i,"chr"] & ps > gene.positions.ali[i,"ps_left"] & ps < gene.positions.ali[i,"ps_right"])

  gwas.sub       = merge(gwas.sub.blup[,c("rs","p_wald_min")], gwas.sub.ratio[,c("rs","p_wald_min")], by = "rs")
  
  gwas.sub$pmin  = do.call(pmin, gwas.sub[,c("p_wald_min.x", "p_wald_min.y")])
  
  top.snp        = subset(gwas.sub, pmin == min(gwas.sub$pmin))
  
  gene.positions.ali[i,"top.snp"] = as.vector(top.snp$rs)
  
}

gene.positions.ali[gene.positions.ali$gene_id == "AT1G65860","top.snp"] = "1_24516880" # GS-OX SNP with maf < 0.03


aliphatics.names = c("BlupG7MTH","BlupG8MTO",
                     "BlupG3MSP","BlupG4MSB","BlupG5MSP","BlupG6MSH","BlupG7MSH","BlupG8MSO",
                     "BlupG2P","BlupG3B","BlupG4P","BlupG2H3B","BlupG2H4P")

effect.sizes = as.data.frame( matrix(ncol = length(aliphatics.names), nrow=length(gene.positions.ali$top.snp), dimnames = list(gene.positions.ali$top.snp,aliphatics.names)) )
effect.pvals = as.data.frame( matrix(ncol = length(aliphatics.names), nrow=length(gene.positions.ali$top.snp), dimnames = list(gene.positions.ali$top.snp,aliphatics.names)) )

for (i in aliphatics.names){
  
  filename = paste0("/tou_miss10/blup_aliphatic_assoc/lmm_maf03_",i,".assoc.fixed.txt")
  
  print(paste0("reading GWAS: ", filename))
  
  gwas = read.delim(filename, sep="\t")
  
  for (j in gene.positions.ali$top.snp){
    
    gwas.focal.snp = subset(gwas, rs == j)
    effect.sizes[j,i] = gwas.focal.snp$beta
    effect.pvals[j,i] = gwas.focal.snp$p_wald
    
  }
  
}

library("heatmap3")

#write.csv(effect.sizes, "focalSNP_aliphatic_effectsizes_wFrachonGeno.csv", quote = F, row.names = T)
#write.csv(effect.pvals, "focalSNP_aliphatic_pvals_wFrachonGeno.csv", quote = F, row.names = T)

effect.sizes.capped = effect.sizes
effect.sizes.capped[effect.sizes.capped >  1] =  1
effect.sizes.capped[effect.sizes.capped < -1] = -1

# effect.pvals.capped = effect.pvals
# effect.pvals.capped[effect.pvals.capped > 0.01] = NA

effect.sizes.capped[effect.pvals > 0.01] = 0


pdf("heatmaps_ali.pdf", h = 6, w = 6)

effect.sizes[8,] = c(-1,0,0,0,0,0,0,0,0,0,0,0,1)
rownames(effect.sizes)[8] = "dummy"

heatmap3(effect.sizes.capped, Rowv = NA, Colv = NA, scale = "none", balanceColor = T)

dev.off()





### indolic effect sizes

# load summarized results of GWAS
gwas.all.blup  = read.csv("/tou_miss10/blup_indolic_assoc/snp_table_all_pvals.csv")
gwas.all.ratio = read.csv("/tou_miss10/blup_indolic_assoc/snp_table_all_pvals.csv")
mvlmm.tou.m10  = subset(read.delim("mvlmmA_maf03.assoc.fixed.txt"), af > 0.05)[,c("rs","chr","ps","af","p_wald")]

# incorporate mvlmm into ratios for TOU-A
gwas.all.ratio  = merge(subset(gwas.all.ratio, select=-c(p_wald_min)), mvlmm.tou.m10[,c("rs","p_wald")], by = "rs")
gwas.all.ratio$p_wald_min = do.call(pmin, gwas.all.ratio[,c(5:ncol(gwas.all.ratio))])

# load candidate gene positions
gene.positions = read.delim("TAIR10_GFF3_geneBoundariesTableADG_extendedWindows.txt")

# subset to genes of interest that had been found to harbor a significant SNP in at least one GWAS
sig.genes.ind  = c("AT1G21100","AT4G37410","AT5G57220")
gene.positions.ind = subset(gene.positions, gene_id %in% sig.genes.ind)
gene.positions.ind$top.snp = "NA"

for (i in 1:length(sig.genes.ind)){
  
  gwas.sub.blup  = subset(gwas.all.blup,  chr == gene.positions.ind[i,"chr"] & ps > gene.positions.ind[i,"ps_left"] & ps < gene.positions.ind[i,"ps_right"])
#  gwas.sub.ratio = subset(gwas.all.ratio, chr == gene.positions.ind[i,"chr"] & ps > gene.positions.ind[i,"ps_left"] & ps < gene.positions.ind[i,"ps_right"])
  
#  gwas.sub       = merge(gwas.sub.blup[,c("rs","p_wald_min")], gwas.sub.ratio[,c("rs","p_wald_min")], by = "rs")
  
#  gwas.sub$pmin  = do.call(pmin, gwas.sub[,c("p_wald_min.x", "p_wald_min.y")])
  
      gwas.sub = gwas.sub.blup
      gwas.sub$pmin = gwas.sub$p_wald_min
  
  top.snp        = subset(gwas.sub, pmin == min(gwas.sub$pmin))
  
  gene.positions.ind[i,"top.snp"] = as.vector(top.snp$rs)
  
}

indolics.names = c("BlupG1MOIM","BlupG4HIM","BlupG4MOIM","BlupGIM")

effect.sizes = as.data.frame( matrix(ncol = length(indolics.names), nrow=length(gene.positions.ind$top.snp), dimnames = list(gene.positions.ind$top.snp,indolics.names)) )
effect.pvals = as.data.frame( matrix(ncol = length(indolics.names), nrow=length(gene.positions.ind$top.snp), dimnames = list(gene.positions.ind$top.snp,indolics.names)) )

for (i in indolics.names){
  
  filename = paste0("/tou_miss10/blup_indolic_assoc/lmm_maf03_",i,".assoc.fixed.txt")
  
  print(paste0("reading GWAS: ", filename))
  
  gwas = read.delim(filename, sep="\t")
  
  for (j in gene.positions.ind$top.snp){
    
    gwas.focal.snp = subset(gwas, rs == j)
    effect.sizes[j,i] = gwas.focal.snp$beta
    effect.pvals[j,i] = gwas.focal.snp$p_wald
    
  }
  
}

library("heatmap3")

# write.csv(effect.sizes, "focalSNP_indolic_effectsizes_wFrachonGeno.csv", quote = F, row.names = T)
# write.csv(effect.pvals, "focalSNP_indolic_pvals_wFrachonGeno.csv", quote = F, row.names = T)

effect.sizes.capped = effect.sizes
effect.sizes.capped[effect.sizes.capped >  1] =  1
effect.sizes.capped[effect.sizes.capped < -1] = -1

# effect.pvals.capped = effect.pvals
# effect.pvals.capped[effect.pvals.capped > 0.01] = NA

effect.sizes.capped[effect.pvals > 0.01] = 0


pdf("heatmaps_ind.pdf", h = 6, w = 6)

effect.sizes[4,] = c(-.5,0,0,.5)
rownames(effect.sizes)[4] = "dummy"

heatmap3(effect.sizes.capped, Rowv = NA, Colv = NA, scale = "none", balanceColor = T)

dev.off()
