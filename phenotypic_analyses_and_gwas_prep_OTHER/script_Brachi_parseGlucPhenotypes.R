#library("plyr")
library("lme4")
library("HLMdiag") # for varcomp.mer

##### 1. LOAD DATA #####

# include accessions mapped in PNAS paper (in/near Europe, but also some elsewhere), using
# data reported in Brachi et al.

phen = read.delim("/brachi/glucs_data_2repmin.csv", h=T, sep=",", na.string="NA")
phen$ID = as.factor(phen$ID)

# ensure ecotype IDs match those used in the imputed genotype file
fam = read.delim("arabidopsis_2029_Maf001_filter95.fam", sep = " ", header = F)[,1]
table(phen$ID %in% fam)

##### 2. RUN LMMs FOR EACH GSL MOLECULE #####

# create table to store some model output when initiating a loop to run the LMM on each molecule

model.info = data.frame(molecule = colnames(phen)[c(-1:-4)],
                        intercept = NA,
                        p = NA,
                        pve_id = NA)

# create a data frame to store BLUPs (effect of each plant ID) from these models

gls.blups = data.frame(ID = unique(phen$ID))

# create lists to store the model objects created by the LMMs with and without plant ID as a factor

gls.models.full = list()
gls.models.null = list()

# loop over the glucosinolate molecule names

for (i in (1:nrow(model.info))){
  
  # run LMM with plant ID as a factor (alternative model)
  
  f1a = paste("log10(", model.info$molecule[i], "+ 1) ~ (1|ID) + (1|batch)")
  m1a = lmer(f1a, data = phen)
  gls.models.full[[i]] = m1a
  
  # run LMM without plant ID as a factor (null model)
  f1n = paste("log10(", model.info$molecule[i], "+ 1) ~ (1|batch)")
  m1n = lmer(f1n, data = phen)
  gls.models.null[[i]] = m1n
  
  # get intercept of model 
  # (useful if converting BLUPs in log10 scale back to predicted phenotypes in linear scale)
  
  model.info$intercept[i] = fixef(m1a)[[1]]
  
  # use likelihood ratio test to assess significance of plant ID term, save in model.info table
  
  lrt = as.numeric(-2*logLik(m1n)+2*logLik(m1a))
  df  = df.residual(m1n) - df.residual(m1a)
  model.info$p[i] = pchisq(lrt, df, lower=F)
  
  # calculate variance explained by plant ID (genotype effect),
  # save in model.info table
  
  model.info$pve_id[i] = varcomp.mer(m1a)[["D11"]] / (varcomp.mer(m1a)[["D11"]] + varcomp.mer(m1a)[["sigma2"]])
  
  # save BLUPs for plant ID
  
  blups = data.frame( rownames(ranef(m1a)$ID), ranef(m1a)$ID[,1] )
  colnames(blups) = c( "ID", paste0(model.info$molecule[i],".blup") )
  gls.blups = merge(gls.blups, blups, by = "ID", all = T)
  
}

write.table(model.info, "/brachi/model_info.txt", sep = "\t", quote = F, eol = "\n", row.names = F)
write.table(gls.blups, "/brachi/phen_glucosinolatesBLUPs_2021-07-28.csv", sep = ",", quote = F, eol = "\n", row.names = F)










