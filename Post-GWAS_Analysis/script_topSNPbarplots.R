# Purpose: Creates barplot comparing top p-value per gene across each dataset (Tou-A, Brachi, Katz, Wu)
#          as shown in supplementary figure.

# load table of top p-values by candidate gene across each dataset
top.pvals.sig.best.all = read.csv("top_pvals_sig_best_aliphatic_all-wFrachonMiss10.csv")
top.pvals.sig.best.all = subset(top.pvals.sig.best.all, select = -X)

# reformat to only include genes that were significant in at least one analysis using BLUPs as traits for GWAS
top.pvals.sig.best.blups.matrix = as.matrix(subset(top.pvals.sig.best.all, select = c(blups.tou.m10,blups.bra.192,blups.wu.192,blups.kat.192), gene_id %in% c("AT1G16400","AT1G65860","AT2G25450","AT3G49680","AT4G03060","AT4G13770","AT5G23010")))
rownames(top.pvals.sig.best.blups.matrix) = top.pvals.sig.best.all$gene_id [top.pvals.sig.best.all$gene_id %in% c("AT1G16400","AT1G65860","AT2G25450","AT3G49680","AT4G03060","AT4G13770","AT5G23010")]

# reformat to only include genes that were significant in at least one analysis using precursor/product ratios as traits for GWAS
top.pvals.sig.best.ratios.matrix = as.matrix(subset(top.pvals.sig.best.all, select = c(ratios.tou.m10,ratios.bra.192,ratios.wu.192,ratios.kat.192), gene_id %in% c("AT1G16400","AT1G65860","AT2G25450","AT3G49680","AT4G03060","AT4G13770","AT5G23010")))
rownames(top.pvals.sig.best.ratios.matrix) = top.pvals.sig.best.all$gene_id [top.pvals.sig.best.all$gene_id %in% c("AT1G16400","AT1G65860","AT2G25450","AT3G49680","AT4G03060","AT4G13770","AT5G23010")]

# create PDF for plot
pdf("comparison_pvals_n192_InflatedTraitsRmd_FrachonGenoMiss10.pdf", h = 5, w = 6.5)
par(mfrow = c(2,7))

# make barplot
for (i in 1:nrow(top.pvals.sig.best.blups.matrix)){
  
  barplot(-log10(top.pvals.sig.best.blups.matrix[i,]), main = paste0("Blups ", rownames(top.pvals.sig.best.blups.matrix)[i]),
          names.arg = colnames(top.pvals.sig.best.blups.matrix), horiz = T, xlim = c(0,max(-log10(top.pvals.sig.best.all.types.matrix))))
  
}

for (i in 1:nrow(top.pvals.sig.best.ratios.matrix)){
  
  barplot(-log10(top.pvals.sig.best.ratios.matrix[i,]), main = paste0("Blups ", rownames(top.pvals.sig.best.ratios.matrix)[i]),
          names.arg = colnames(top.pvals.sig.best.ratios.matrix), horiz = T, xlim = c(0,max(-log10(top.pvals.sig.best.all.types.matrix))))
  
}

dev.off()
