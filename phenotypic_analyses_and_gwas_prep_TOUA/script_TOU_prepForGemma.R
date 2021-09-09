
# load in TOU-A BLUPs
dat = read.csv("/toua_wFrachonGenotypes/phen_glucosinolatesBLUPs_2021-08-02.csv")
acc = read.csv("/toua_wFrachonGenotypes/TOU_accessions.csv")

# merge with list of genotypes, ordered as in the .fam file
# (from the set of PLINK binary files used as genotypes in GWAS with GEMMA)
fam.file.IDs = read.csv("/toua_wFrachonGenotypes/genotypes.fam", sep = " ", header = F)[,1]

dat = merge(acc, dat, by = "ID", all.x = T)

dat = dplyr::left_join(data.frame(ecotype_id=fam.file.IDs), dat, by = "ecotype_id")

# glucosinolates to retain
dat = dat[,c("ecotype_id","gsl.7mSh.blup","gsl.8MTO.blup","gsl.3mSOp.blup","gsl.4mSOb.blup",
               "gsl.5mSOp.blup","gsl.6mSOh.blup","gsl.7mSOh.blup","gsl.8mSOo.blup","gsl.Pren.blup",
               "gsl.Buen.blup","gsl.Peen.blup","gsl.S2hBuen.blup","gsl.2hPeen.blup","gsl.IM.blup",
               "gsl.1moIM.blup","gsl.1hIM.blup","gsl.4moIM.blup")]

colnames(dat)[1] = "ID"

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
write.table(dat[,c(2:ncol(dat))], "/toua_wFrachonGenotypes/dat_tot_forGemma.csv",
            sep = ",", row.names = F, col.names = F, quote = F, eol = "\n")

# calculate PCA (note: only followed up on long-chain PC2, reflecting mSx vs. mSOx abundances)
gsl.names.ali.short = c("gsl.3mSOp.blup","gsl.4mSOb.blup","gsl.5mSOp.blup","gsl.Pren.blup","gsl.Buen.blup","gsl.Peen.blup","gsl.S2hBuen.blup","gsl.2hPeen.blup")
gsl.names.ali.long  = c("gsl.7mSh.blup","gsl.8MTO.blup","gsl.6mSOh.blup","gsl.7mSOh.blup","gsl.8mSOo.blup")
gsl.names.ali.all   = c(gsl.names.ali.short, gsl.names.ali.long)
gsl.names.ind       = c("gsl.IM.blup","gsl.1moIM.blup","gsl.1hIM.blup","gsl.4moIM.blup")




par(mfrow = c(3,2))

dat.pca = data.frame(ID = na.omit(dat)$ID)

pca = prcomp(na.omit(dat[,gsl.names.ali.short]), center = F, scale = F)
summary(pca)
barplot(pca$rot[,1], las=2); barplot(pca$rot[,2], las=2); barplot(pca$rot[,3], las=2); barplot(pca$rot[,4], las=2); barplot(pca$rot[,5], las=2)
dat.pca = cbind(dat.pca, pca$x[,c(1:5)])
colnames(dat.pca)[2:6] = paste(colnames(dat.pca)[2:6], "alishort", sep=".")

pca = prcomp(na.omit(dat[,gsl.names.ali.long]), center = F, scale = F)
summary(pca)
barplot(pca$rot[,1], las=2); barplot(pca$rot[,2], las=2); barplot(pca$rot[,3], las=2); barplot(pca$rot[,4], las=2); barplot(pca$rot[,5], las=2)
dat.pca = cbind(dat.pca, pca$x[,c(1:5)])
colnames(dat.pca)[7:11] = paste(colnames(dat.pca)[7:11], "alilong", sep=".")

pca = prcomp(na.omit(dat[,gsl.names.ali.all]), center = F, scale = F)
summary(pca)
barplot(pca$rot[,1], las=2); barplot(pca$rot[,2], las=2); barplot(pca$rot[,3], las=2); barplot(pca$rot[,4], las=2); barplot(pca$rot[,5], las=2)
dat.pca = cbind(dat.pca, pca$x[,c(1:5)])
colnames(dat.pca)[12:16] = paste(colnames(dat.pca)[12:16], "aliall", sep=".")

pca = prcomp(na.omit(dat[,gsl.names.ind]), center = F, scale = F)
summary(pca)
barplot(pca$rot[,1], las=2); barplot(pca$rot[,2], las=2); barplot(pca$rot[,3], las=2); barplot(pca$rot[,4], las=2);
dat.pca = cbind(dat.pca, pca$x[,c(1:4)])
colnames(dat.pca)[17:20] = paste(colnames(dat.pca)[17:20], "ind", sep=".")

dat.pca = dplyr::left_join(data.frame(ID=fam.file.IDs), dat.pca, by = "ID")

# write .csv output file for GEMMA (phenotypes, no headers)
write.table(dat.pca[,c(2:20)], "/toua_wFrachonGenotypes/dat_tot_pca_forGemma.csv",
            sep = ",", row.names = F, col.names = F, quote = F, eol = "\n")





