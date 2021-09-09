
library("tidyr")

# load in Tou-A BLUPs
dat = read.csv("/toua_wFrachonGenotypes/phen_glucosinolatesBLUPs_2021-08-02.csv")

# load in model info (need: intercepts)
model.info = read.csv("/toua_wFrachonGenotypes/model_info.txt", sep = "\t")

# convert from 0-centered BLUPs to predicted values (intercept + BLUP, then log10 --> normal scale)
dat.predict = data.frame(matrix(nrow = nrow(dat), ncol = ncol(dat)))

colnames(dat.predict)[1] = "ID"
dat.predict$ID = dat$ID

colnames(dat.predict)[2:ncol(dat.predict)] = as.character(model.info$molecule)

# log-transform
for (i in model.info$molecule){
  
  dat.predict[,i] = 10^ ( dat[,paste0(i,".blup")] + model.info[model.info$molecule == i,"intercept"] )

}

# note to self: for loop line is written weirdly (but functional), be more concise
par(mfrow = c(2,2))
for (i in colnames(dat.predict)[2:ncol(dat.predict)]){
  temp = dat.predict[,c("ID",i)]
  hist(temp[,2], main=i)
}

# check if smallest value = 0
min(dat.predict[,2:ncol(dat.predict)]) # no, so don't need to add constant

# log2 ratios to calculate (after adding small value to every observation avoid 0's)
#
#   ALIPHATIC:
#   G3MSP:G4MSB   G3MSP:G5MSP   G3MSP:G8MSO   G4MSB:G5MSP   G4MSB:G8MSO   G5MSP:G8MSO
#   G3B:G4P       G2H3B:G2H4P
#   G3MTP:G3MSP   G4MTB:G4MSB   G7MTH:G7MSH   G8MTH:G8MSO
#   G3MSP:G2P     G4MSB:G3B
#   G3B:G2H3B     G4P:G2H4P
#
#   INDOLIC:
#   GI3M:G1MOI3M
#   GI3M:G4MOI3M
#   GI3M:G4HI3M
#   G4HI3M:G4MOI3M

ratios.2.calculate = read.csv("/toua_wFrachonGenotypes/ratios_to_calculate.csv", h = T)

blups.ratios = data.frame(matrix(nrow = nrow(dat), ncol = nrow(ratios.2.calculate)+1))
colnames(blups.ratios)[1] = "ID"
blups.ratios$ID = dat.predict$ID

# calculate ratios
for (i in 1:nrow(ratios.2.calculate)){
  
  mol1 = as.character( ratios.2.calculate[i,"molecule_1"] )
  mol2 = as.character( ratios.2.calculate[i,"molecule_2"] )

  blups.ratios[,i+1] = log2( dat.predict[,mol1] / dat.predict[,mol2] )
  
  colnames(blups.ratios)[i+1] = paste0("ratio",i)
  
}

par(mfrow = c(2,2))
for (i in colnames(blups.ratios)[2:ncol(blups.ratios)]){
  temp = blups.ratios[,c("ID",i)]
  hist(temp[,2], main=i)
}

# write output for GEMMA


acc = read.csv("/toua_wFrachonGenotypes/TOU_accessions.csv")

# merge with list of genotypes, ordered as in the .fam file
# (from the set of PLINK binary files used as genotypes in GWAS with GEMMA)
fam.file.IDs = read.csv("/TOUA_glucs/toua_wFrachonGenotypes/genotypes.fam", sep = " ", header = F)[,1]

dat = merge(acc, blups.ratios, by = "ID", all.x = T)

dat = dplyr::left_join(data.frame(ecotype_id=fam.file.IDs), dat, by = "ecotype_id")

table( dat$ecotype_id - fam.file.IDs )

dat = subset(dat, select = c(-ecotype_id, -ID, -SAMPLE, -pop_code))

# write .csv output file for GEMMA (phenotypes, no headers)
write.table(dat, "/toua_wFrachonGenotypes/dat_tot_ratios_forGemma.csv",
            sep = ",", row.names = F, col.names = F, quote = F, eol = "\n")

# write .csv output file with BLUPs converted to original scale (not log-transformed, scaled or centered, etc.)
write.table(dat.predict, "/toua_wFrachonGenotypes/dat_tot_original_scale.csv",
            sep = ",", row.names = F, col.names = T, quote = F, eol = "\n")


