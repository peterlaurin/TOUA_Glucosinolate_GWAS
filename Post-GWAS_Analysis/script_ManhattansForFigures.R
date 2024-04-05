# Purpose: Create Manhattan plots showing top p-value across combined traits (glucosinolate molecules) within a
#          molecule class (e.g., aliphatic glucosinolates). Indicate p-value threshold and color SNPs in
#          candidate gene windows as blue.
#
# Note:    Aesthetics and labels of the outputted plots were manually refined 
#          in Adobe Illustrator to create the plots in the main figures of the paper.

# Define a custom function for Manhattan plots
manhattan_single = function(rs, chr, ps, pval, max.y, outfile, snps2highlight, sig.thresh, sig.thresh.multi, point.scale = 1, gap.size = 1e6){
  
  ## create object storing positions to plot
  
  ps.plot = ps
  
  # loop over chromosomes, assuming they have consecutive numeric names starting at 1
  for (i in 2:length(unique(chr))){
    
    # adjust positions of each SNP on x-axis to add "gap" between each chromosome
    ps.plot[chr >= i] = ps[chr >= i] + max(ps.plot[chr == i-1]) + gap.size
    
  }
  
  # # qq Plot (not used; ignore)
  # 
  # CMplot(data.frame(rs, chr, ps, pval), plot.type = "q",
  #        main = file_path_sans_ext(infilename), memo = file_path_sans_ext(infilename))
  
  ## Manhattan Plot
  
  # convert between -log10 and unscaled p values
  pval[-log10(pval) > max.y] = 10^-max.y
  
  # set filename for output Manhattan plot
  png(filename = outfile, width = 4800, height = 2400)
  
  # plot points for all SNPs
  print( plot(-log10(pval) ~ ps.plot, ylim = c(0,max.y+1), axes = F, xlab = "", ylab = "",  
              pch = 16, col = "gray90", cex = 2.2 * point.scale) )

    # plot points for significant SNPs  
  print( points(-log10(pval[pval<sig.thresh]) ~ ps.plot[pval<sig.thresh], 
                pch = 16, col = "gray75", cex = 2.2 * point.scale) )
  
  # plot points for non-significant SNPs in candidate genes
  print( points(-log10(pval[rs %in% snps2highlight & pval > sig.thresh.multi]) ~ 
                  ps.plot[rs %in% snps2highlight & pval > sig.thresh.multi], 
                pch = 16, col = "dodgerblue", cex = 3.5 * point.scale) ) 
  
  # plot points for significant SNPs in candidate genes
  print( points(-log10(pval[rs %in% snps2highlight & pval < sig.thresh.multi & pval != 10^-max.y]) ~ 
                  ps.plot[rs %in% snps2highlight & pval < sig.thresh.multi & pval != 10^-max.y], 
                pch = 16, col = "dodgerblue", cex = 5.0 * point.scale) )
  
  # optionally, if p-values are equal to the y axis limit, make these points larger
  # (can be useful if setting p-values above the limit to equal the limit to avoid stretched axis)
  print( points(-log10(pval[rs %in% snps2highlight & pval == 10^-max.y]) ~ 
                  ps.plot[rs %in% snps2highlight & pval == 10^-max.y], 
                pch = 17, col = "dodgerblue", cex = 6.5 * point.scale) ) 
  
  # draw lines for single GWAS and and multi-GWAS significance thresholds
  print( abline(h = -log10(sig.thresh), col = "red", lty = 2, lwd = 4) )
  print( abline(h = -log10(sig.thresh.multi), col = "black", lty = 2, lwd = 4) )
  
  # plot axes
  print( axis(side = 2, cex.axis = 5, cex.lab = 7, lwd = 5))
  
  dev.off()
  
}

########## Aliphatic GWAS

# Load candidate SNP files
candidate.snps = list()
candidate.snps[["aliphatic"]] = read.table("candidate_snps_aliphatic_30kb_majorLociExtended_maf05_wTouFrachonGenoMiss10.txt", h = F) 

# Load table with GWAS p-values combined across all traits for Tou-A 
tou.ali.blup.bestp    = read.csv("/blup_aliphatic_assoc/snp_table_all_pvals.csv")
tou.ali.ratio.bestp   = read.csv("/ratio_aliphatic_assoc/snp_table_all_pvals.csv")

# Run function above to make Manhattan Plots
# (user: manually modify gwas.all and gwas.name to use names from directly above)
gwas.all  = tou.ali.ratio.bestp
gwas.name = "tou.ali.ratio.bestp"
gsl.type  = "aliphatic"
manhattan_single(rs               = gwas.all$rs, # position of each SNP
                 chr              = gwas.all$chr, # chromosome of each SNP
                 ps               = gwas.all$ps, # unique ID of each SNP
                 pval             = gwas.all$p_wald, # p value for each SNP
                 max.y            = 25, # specify max -log10(p) value
                 outfile          = paste0("~/Desktop/Figures/Manhattan_",gwas.name,".png"), # specify location for 
                 snps2highlight   = subset(gwas.all, rs %in% candidate.snps[[gsl.type]]$V1)[,"rs"], # specify candidate SNPs
                 sig.thresh       = .05 / nrow(gwas.all), # single GWAS significance threshold
                 sig.thresh.multi = .05 / nrow(gwas.all)/ (ncol(gwas.all) - 4 - 1), # multi-GWAS combined significance threshold
                 point.scale      = 1.3, # size of points
                 gap.size         = 2e6 # size of gap (number of positions) between each chromosome on x-axis of plot
)

########## Indolic GWAS

# Load candidate SNP files
candidate.snps = list()
candidate.snps[["indolic"]]   = read.table("candidate_snps_indolic_30kb_majorLociExtended_maf05_wTouFrachonGenoMiss10.txt", h = F)

# Load table with GWAS p-values combined across all traits for Tou-A 
tou.ind.blup.bestp    = read.csv("/blup_indolic_assoc/snp_table_all_pvals.csv")
tou.ind.ratio.bestp   = read.csv("/ratio_indolic_assoc/snp_table_all_pvals.csv")
tou.mvlmm             = subset(read.delim("mvlmmA_maf03.assoc.fixed.txt"), af > 0.05)

# Run function above to make Manhattan Plots
# (user: manually modify gwas.all and gwas.name to use names from directly above)
gwas.all  = tou.ind.ratio.bestp
gwas.name = "tou.ind.ratio.bestp"
gsl.type  = "indolic"
manhattan_single(rs               = gwas.all$rs,
                 chr              = gwas.all$chr,
                 ps               = gwas.all$ps,
                 pval             = gwas.all$p_wald,
                 max.y            = 15,
                 outfile          = paste0("~/Desktop/Figures/Manhattan_",gwas.name,".png"),
                 snps2highlight   = subset(gwas.all, rs %in% candidate.snps[[gsl.type]]$V1)[,"rs"],
                 sig.thresh       = .05 / nrow(gwas.all),
                 sig.thresh.multi = .05 / nrow(gwas.all)/ (ncol(gwas.all) - 4 - 1),
                 point.scale      = 1.3,
                 gap.size         = 2e6
)

gwas.all  = tou.mvlmm
gwas.name = "tou.mvlmm"
gsl.type  = "indolic"
manhattan_single(rs               = gwas.all$rs,
                 chr              = gwas.all$chr,
                 ps               = gwas.all$ps,
                 pval             = gwas.all$p_wald,
                 max.y            = 15,
                 outfile          = paste0("~/Desktop/Figures/Manhattan_",gwas.name,".png"),
                 snps2highlight   = subset(gwas.all, rs %in% candidate.snps[[gsl.type]]$V1)[,"rs"],
                 sig.thresh       = .05 / nrow(gwas.all),
                 sig.thresh.multi = 1e-999,
                 point.scale      = 1.8,
                 gap.size         = 2e6
)
