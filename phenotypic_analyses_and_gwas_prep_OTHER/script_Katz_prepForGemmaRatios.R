# RATIOS

# load in model info (need: intercepts)
model.info = read.csv("/katz/model_info_significantOnly.txt", sep = "\t")

# load in katz BLUPs
dat = read.csv("/katz/phen_glucosinolatesBLUPs_2021-07-30.csv")

# merge with list of genotypes, ordered as in the .fam file
# (from the set of PLINK binary files used as genotypes in GWAS with GEMMA)
fam.file.IDs = read.csv("arabidopsis_2029_Maf001_filter95.fam", sep = " ", header = F)[,1]

dat = dplyr::left_join(data.frame(ID=fam.file.IDs), dat, by = "ID")

# verify order is same as in the .fam genotype file
id.check = data.frame(fam.IDs = fam.file.IDs, blup.ids = dat$ID)
id.check$dif = id.check$fam.IDs - id.check$blup.ids
table(id.check$dif)

# convert from 0-centered BLUPs to predicted values (intercept + BLUP, then log10 --> normal scale)

dat.predict = data.frame(matrix(nrow = nrow(dat), ncol = ncol(dat)))

colnames(dat.predict)[1] = "ID"
dat.predict$ID = dat$ID

# could adjust this and preceding lines to be ncol = nrow in model info, not ncol in dat, due to n.s. tests w/ no BLUPs
colnames(dat.predict)[2:ncol(dat.predict)] = as.character(model.info$molecule)


for (i in model.info$molecule){
  
  dat.predict[,i] = 10^ ( dat[,paste0(i,".blup")] + model.info[model.info$molecule == i,"intercept"] )
  
}

par(mfrow = c(2,2))
for (i in colnames(dat.predict)[2:ncol(dat.predict)]){
  temp = dat.predict[,c("ID",i)]
  hist(temp[,2], main=i)
}

# check if smallest value = 0
min(dat.predict[,2:ncol(dat.predict)], na.rm = T) # no, so don't need to add constant

ratios.2.calculate = read.csv("/katz/BLUPs/ratios_to_calculate.csv", h = T)

blups.ratios = data.frame(matrix(nrow = nrow(dat), ncol = nrow(ratios.2.calculate)+1))
colnames(blups.ratios)[1] = "ID"
blups.ratios$ID = dat.predict$ID

for (i in 1:nrow(ratios.2.calculate)){
  
  mol1 = as.character( ratios.2.calculate[i,"molecule_1a"] )
  mol2 = as.character( ratios.2.calculate[i,"molecule_2a"] )
  
  blups.ratios[,i+1] = log2( dat.predict[,mol1] / dat.predict[,mol2] )
  
  colnames(blups.ratios)[i+1] = paste0("ratio",i)
  
}

par(mfrow = c(2,2))
for (i in colnames(blups.ratios)[2:ncol(blups.ratios)]){
  temp = blups.ratios[,c("ID",i)]
  hist(temp[,2], main=i)
}

# write .csv output file for GEMMA (phenotypes, no headers)
write.table(blups.ratios[,2:ncol(blups.ratios)], "/katz/BLUPs/dat_tot_ratios_forGemma.csv",
            sep = ",", row.names = F, col.names = F, quote = F, eol = "\n")

# write .csv output file with BLUPs converted to original scale (not log-transformed, scaled or centered, etc.)
write.table(dat.predict, "/katz/BLUPs/dat_tot_original_scale.csv",
            sep = ",", row.names = F, col.names = T, quote = F, eol = "\n")



