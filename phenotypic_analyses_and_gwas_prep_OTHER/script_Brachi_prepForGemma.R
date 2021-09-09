
# load in Brachi BLUPs
dat = read.csv("/brachi/phen_glucosinolatesBLUPs_2021-07-28.csv")

# merge with list of genotypes, ordered as in the .fam file
# (from the set of PLINK binary files used as genotypes in GWAS with GEMMA)
fam.file.IDs = read.csv("arabidopsis_2029_Maf001_filter95.fam", sep = " ", header = F)[,1]

dat = dplyr::left_join(data.frame(ID=fam.file.IDs), dat, by = "ID")

# verify order is same as in the .fam genotype file
id.check = data.frame(fam.IDs = fam.file.IDs, blup.ids = dat$ID)
id.check$dif = id.check$fam.IDs - id.check$blup.ids
table(id.check$dif)

# inspect phenotype distributions
par(mfrow = c(2,2))
for (i in colnames(dat)[2:ncol(dat)]){
  temp = dat[,c("ID",i)]
  hist(temp[,2], main=i)
}

# scale and center phenotypes
for (i in colnames(dat)[2:ncol(dat)]){
  dat[,i] = scale(dat[,i])
}

# write .csv output file for GEMMA (phenotypes, no headers)
write.table(dat[,c(2:ncol(dat))], "/brachi/BLUPs/dat_tot_forGemma.csv",
            sep = ",", row.names = F, col.names = F, quote = F, eol = "\n")



