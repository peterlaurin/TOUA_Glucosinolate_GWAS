# RATIOS

# load in Wu means
dat = read.csv("/wu/phen_glucosinolatesMeans_2021-07-30.csv")

# merge with list of genotypes, ordered as in the .fam file
# (from the set of PLINK binary files used as genotypes in GWAS with GEMMA)
fam.file.IDs = read.csv("arabidopsis_2029_Maf001_filter95.fam", sep = " ", header = F)[,1]

dat = dplyr::left_join(data.frame(ID=fam.file.IDs), dat, by = "ID")

# verify order is same as in the .fam genotype file
id.check = data.frame(fam.IDs = fam.file.IDs, blup.ids = dat$ID)
id.check$dif = id.check$fam.IDs - id.check$blup.ids
table(id.check$dif)



# check if smallest value = 0
min(dat[,2:ncol(dat)], na.rm = T) # yes, add constant
sort(unique(as.vector(as.matrix(dat[,2:ncol(dat)]))))[1:5] # smallest value (non-zero) is 7.632979
dat[,2:ncol(dat)] = dat[,2:ncol(dat)] + 1 # add 1

ratios.2.calculate = read.csv("/wu/means/ratios_to_calculate.csv", h = T)

blups.ratios = data.frame(matrix(nrow = nrow(dat), ncol = nrow(ratios.2.calculate)+1))
colnames(blups.ratios)[1] = "ID"
blups.ratios$ID = dat$ID

for (i in 1:nrow(ratios.2.calculate)){
  
  mol1 = as.character( ratios.2.calculate[i,"molecule_1a"] )
  mol2 = as.character( ratios.2.calculate[i,"molecule_2a"] )
  
  blups.ratios[,i+1] = log2( dat[,mol1] / dat[,mol2] )
  
  colnames(blups.ratios)[i+1] = paste0("ratio",i)
  
}

par(mfrow = c(2,2))
for (i in colnames(blups.ratios)[2:ncol(blups.ratios)]){
  temp = blups.ratios[,c("ID",i)]
  hist(temp[,2], main=i)
}

# write .csv output file for GEMMA (phenotypes, no headers)
write.table(blups.ratios[,2:ncol(blups.ratios)], "/wu/means/dat_tot_ratios_forGemma.csv",
            sep = ",", row.names = F, col.names = F, quote = F, eol = "\n")

# write .csv output file with BLUPs converted to original scale (not log-transformed, scaled or centered, etc.)
write.table(dat.predict, "/wu/means/dat_tot_original_scale.csv",
            sep = ",", row.names = F, col.names = T, quote = F, eol = "\n")



