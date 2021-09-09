
##### 1. LOAD DATA #####

phen = read.delim("/katz/glucs_data_katz.csv", h=T, sep=",", na.string="NA")

# include accessions mapped in eLife paper (only in/near Europe)
acc.eur = read.delim("/katz/emmeans_eLife_supTableD.csv", h=T, sep=",", na.string="NA")[,"CS"]
phen = subset(phen, CS %in% acc.eur)

id.key = read.delim("/katz/1001genomes_accesions_IDs_key.csv", h=T, sep=",", na.string="NA")
phen = merge(id.key[,c("ID","CS")], phen, by = "CS")
phen$round_plate = paste0(phen$Round, "_", phen$Plate)

phen$ID = as.factor(phen$ID)
phen$Round = as.factor(phen$Round)

# ensure ecotype IDs match those used in the imputed genotype file
fam = read.delim("arabidopsis_2029_Maf001_filter95.fam", sep = " ", header = F)[,1]
table(phen$ID %in% fam)

# inspect phenotype distributions
# add smallest value observed in dataset for log transformation
# sort(unique(as.vector(as.matrix(phen[,8:27]))))[1:5] # --> smallest value is 5.86898e-06, so add 5e-6
for(i in 8:(ncol(phen)-1)) {hist(log10(phen[,i]+5e-6), main = colnames(phen)[i])}
# note: log transform is warranted; very right skewed dataset with little separation between 0 and non-zero values
#       in original distribution, and normal distribution of non-zero values with separation from 0s in transformed
#       distribution. (I investigated this for G3HP).


##### 2. RUN LMMs FOR EACH GSL MOLECULE #####

# create table to store some model output when initiating a loop to run the LMM on each molecule

model.info = data.frame(molecule = colnames(phen)[c(-1:-7,-28)],
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
  
  f1a = paste("log10(", model.info$molecule[i], "+ 5e-6) ~ (1|ID) + (1|round_plate)")
  m1a = lmer(f1a, data = phen)
  gls.models.full[[i]] = m1a
  
  # run LMM without plant ID as a factor (null model)
  f1n = paste("log10(", model.info$molecule[i], "+ 5e-6) ~ (1|round_plate)")
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

# remove non-significant molecules (p-value_ID > 0.05)
gls.blups = subset(gls.blups, select = -c(G2H4P.blup,G4MOI3M.blup,G4HI3M.blup))

# write.table(model.info, "/katz/model_info.txt_significantOnly.csv", sep = "\t", quote = F, eol = "\n", row.names = F)
# write.table(gls.blups, "/katz/phen_glucosinolatesBLUPs_2021-07-30.csv", sep = ",", quote = F, eol = "\n", row.names = F)



#####

# Note: it's unsurprising that G2H4P is not significantly heritable in the LMM
# gplots::venn(list(unique(subset(phen,G2H4P>0))$ID, unique(subset(phen,G2H4P==0)$ID)))
# --> 27 accessions had it in 1 of 2 reps, 1 accession had it both times, and the rest lacked it

# Note: null model failed to converge for i=1 (molecule: G3HP), but this is fine since at worst
# it biases the p value, but ID already explains the majority of the variance in the 
# alternate model

