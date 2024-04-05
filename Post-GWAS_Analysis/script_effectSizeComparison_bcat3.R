# Purpose: Creates plot comparing effects of the BCAT3 mutant with
#          GWAS effect sizes for the leading SNP at the the BCAT3 locus 

##### Prepare data frame to store results #####

XmSOo.names = c("Blup3mSOp","Blup4mSOb","Blup5mSOp","Blup6mSOh","Blup7mSOh","Blup8mSOo")
XmSOo.names.short = c("3C","4C","5C","6C","7C","8C")

focal.loci = c("MAM","BCAT3")
focal.snps = c("5_7694175","3_18425806")

effect.sizes = data.frame(locus              = rep(focal.loci, each = length(XmSOo.names)),
                          snp                = rep(focal.snps, each = length(XmSOo.names)),
                          method             = "GWAS",
                          population         = "TOU-A",
                          molecule           = rep(XmSOo.names, length(focal.snps)),
                          molecule.short     = c("3C","4C","5C","6C","7C","8C"),
                          zscore.log10       = NA,
                          beta.zscore.log10  = NA,
                          se.zscore.log10    = NA,
                          beta.scaled.log2fc = NA,
                          p.wald             = NA,
                          se.scaled.log2fc.upper = NA,
                          se.scaled.log2fc.lower = NA
                          )

##### Get effect size (z score + SE, p-value) from each GWAS #####

for (i in XmSOo.names){
  
  filename = paste0("/blup_aliphatic_assoc/lmm_maf03_",i,".assoc.fixed.txt")
  print(paste0("reading file: ", filename, "..."))
  
  gwas = read.delim(filename, sep="\t")
  
  for (j in focal.snps){
    
    gwas.focal.snp = subset(gwas, rs == j)
    
    effect.sizes[effect.sizes$snp == j & effect.sizes$molecule == i,]$beta.zscore.log10 = gwas.focal.snp$beta
    effect.sizes[effect.sizes$snp == j & effect.sizes$molecule == i,]$se.zscore.log10 = gwas.focal.snp$se
    effect.sizes[effect.sizes$snp == j & effect.sizes$molecule == i,]$p.wald = gwas.focal.snp$p_wald
    
  }
  
}

##### Use distribution of BLUP values to convert effect size from z-score to log2(fold-change) #####

dat.blups     = read.csv("phen_glucosinolatesBLUPs_2021-08-02.csv")
dat.blups$set = stringr::str_extract(dat.blups$ID, "[^_]+")
dat.blups     = subset(dat.blups, set != 5)

XmSOo.blup.names = c("gsl.3mSOp.blup","gsl.4mSOb.blup","gsl.5mSOp.blup","gsl.6mSOh.blup","gsl.7mSOh.blup","gsl.8mSOo.blup")

for (i in 1:length(XmSOo.blup.names)){
  
  blup.id = XmSOo.blup.names[i]
  mol.id  = XmSOo.names[i]
  
  effect.sizes[effect.sizes$molecule == mol.id,]$zscore.log10 = sd(dat.blups[,blup.id])
  
}

effect.sizes$beta.log10 = effect.sizes$beta.zscore.log10 * effect.sizes$zscore.log10 # SE <---
effect.sizes$se.log10   = effect.sizes$se.zscore.log10   * effect.sizes$zscore.log10 # SE <---

effect.sizes$fc = 10^(1+effect.sizes$beta.log10) / 10^(1)
effect.sizes$log2fc = log2(effect.sizes$fc)

effect.sizes$se.fc.upper = 10^(1 + effect.sizes$beta.log10 + effect.sizes$se.log10) / 10^(1) # SE <---
effect.sizes$se.fc.lower = 10^(1 + effect.sizes$beta.log10 - effect.sizes$se.log10) / 10^(1) # SE <---
effect.sizes$ci.fc.upper = 10^(1 + effect.sizes$beta.log10 + effect.sizes$se.log10*2) / 10^(1) # CI <---
effect.sizes$ci.fc.lower = 10^(1 + effect.sizes$beta.log10 - effect.sizes$se.log10*2) / 10^(1) # CI <---
effect.sizes$se.log2fc.upper = log2(effect.sizes$se.fc.upper) # SE <---
effect.sizes$se.log2fc.lower = log2(effect.sizes$se.fc.lower) # SE <---
effect.sizes$ci.log2fc.upper = log2(effect.sizes$ci.fc.upper) # CI <---
effect.sizes$ci.log2fc.lower = log2(effect.sizes$ci.fc.lower) # CI <---

##### Make plots #####

library("ggplot2")

level_order = rev(c("3C","4C","5C","6C","7C","8C"))

# BCAT3 GWAS, z-score + 95%CI
temp = subset(effect.sizes, snp == "3_18425806")
g = ggplot( temp, aes(x=beta.zscore.log10, y=factor(molecule.short, level = level_order)) ) + 
  geom_pointrange(aes(xmin=beta.zscore.log10-se.zscore.log10*2, xmax=beta.zscore.log10+se.zscore.log10*2), alpha = 0.99) +
  xlim(-1,1)
ggsave("/effect_size_dotPlots_MAM-vs-BCAT3/effect_GWAS-z_BCAT3.pdf", height = 2, width = 2)

# MAM GWAS, z-score + 95%CI
temp = subset(effect.sizes, snp == "5_7694175")
g = ggplot( temp, aes(x=beta.zscore.log10, y=factor(molecule.short, level = level_order)) ) + 
  geom_pointrange(aes(xmin=beta.zscore.log10-se.zscore.log10*2, xmax=beta.zscore.log10+se.zscore.log10*2), alpha = 0.99) +
  xlim(-1,1)
ggsave("~/effect_size_dotPlots_MAM-vs-BCAT3/effect_GWAS-z_MAM.pdf", height = 2, width = 2)

# BCAT3 GWAS, log2fc + 95%CI
temp = subset(effect.sizes, snp == "3_18425806")
g = ggplot( temp, aes(x=log2fc, y=factor(molecule.short, level = level_order)) ) + 
  geom_pointrange(aes(xmax=ci.log2fc.upper, xmin=ci.log2fc.lower), alpha = 0.99) +
  xlim(-1,1)
ggsave("/effect_size_dotPlots_MAM-vs-BCAT3/effect_GWAS-log2fc_BCAT3.pdf", height = 2, width = 2)

# MAM GWAS, log2fc + 95%CI
temp = subset(effect.sizes, snp == "5_7694175")
g = ggplot( temp, aes(x=log2fc, y=factor(molecule.short, level = level_order)) ) + 
  geom_pointrange(aes(xmax=ci.log2fc.upper, xmin=ci.log2fc.lower), alpha = 0.99) +
  xlim(-1,1)
ggsave("/effect_size_dotPlots_MAM-vs-BCAT3/effect_GWAS-log2fc_MAM.pdf", height = 2, width = 2)


###

mutants = read.csv("effectSizes_bcat3-vs-MAM_mutants.csv")

# BCAT3 mutant, log2(fc)
temp = subset(mutants, gene == "BCAT3")
temp$log2fc = log2(temp$KO / temp$WT)
g = ggplot( temp, aes(x=log2fc, y=factor(molecule.short, level = level_order)) ) + 
    geom_point(size = 2.75, alpha = 0.99) + xlim(-5,5)
ggsave("/effect_size_dotPlots_MAM-vs-BCAT3/effect_mutant-log2fc_BCAT3.pdf", height = 2, width = 2)

# MAM1 mutant, log2(fc)
temp = subset(mutants, gene == "MAM1")
temp$log2fc = log2(temp$KO / temp$WT)
g = ggplot( temp, aes(x=log2fc, y=factor(molecule.short, level = level_order)) ) + 
    geom_point(size = 2.75, alpha = 0.99) + xlim(-5,5)
ggsave("/effect_size_dotPlots_MAM-vs-BCAT3/effect_mutant-log2fc_MAM1.pdf", height = 2, width = 2)

# MAM3 mutant, log2(fc)
temp = subset(mutants, mutant == "gsm2-avg")
temp$log2fc = log2(temp$KO / temp$WT)
g = ggplot( temp, aes(x=log2fc, y=factor(molecule.short, level = level_order)) ) + 
    geom_point(size = 2.75, alpha = 0.99) + xlim(-5,5)
ggsave("/effect_size_dotPlots_MAM-vs-BCAT3/effect_mutant-log2fc_MAM3.pdf", height = 2, width = 2)

# pvals

temp = subset(effect.sizes, snp == "3_18425806")
temp$logp = -log10(temp$p.wald)
g = ggplot( temp, aes(x = factor(molecule.short, level = level_order), y = logp)) +
  geom_bar(stat="identity") + ylim(0,20) + coord_flip()
ggsave("/effect_size_dotPlots_MAM-vs-BCAT3/effect_GWAS_pval_BCAT3.pdf", height = 2, width = 2)

temp = subset(effect.sizes, snp == "5_7694175")
temp$logp = -log10(temp$p.wald)
g = ggplot( temp, aes(x = factor(molecule.short, level = level_order), y = logp)) +
  geom_bar(stat="identity") + ylim(0,20) + coord_flip()
ggsave("/effect_size_dotPlots_MAM-vs-BCAT3/effect_GWAS_pval_MAM.pdf", height = 2, width = 2)


