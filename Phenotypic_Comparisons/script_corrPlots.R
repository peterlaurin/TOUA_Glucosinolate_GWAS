# Purpose: Create plot showing correlation between per-accession BLUPs for each glucosinolate molecule within a given study

# Define list of molecules shared across studies
shared.mols = c("7mSh","8mSo","3mSOp","4mSOb","5mSOp","6mSOh","7mSOh","8mSOo","Pren","Buen","Peen","2hBuen")

library("corrplot")

# Note: in-line comments only shown for TOU-A section of this script but are the same across all sections

### TOU-A
# load table of BLUPs
dat = read.csv("phen_glucosinolatesBLUPs_2021-08-02.csv")

# define sub-categories of molecules
mol.mt   = c("gsl.7mSh.blup","gsl.8MTO.blup")
mol.ms   = c("gsl.3mSOp.blup","gsl.4mSOb.blup","gsl.5mSOp.blup","gsl.6mSOh.blup","gsl.7mSOh.blup","gsl.8mSOo.blup")
mol.alk  = c("gsl.Pren.blup","gsl.Buen.blup","gsl.Peen.blup")
mol.halk = c("gsl.S2hBuen.blup","gsl.2hPeen.blup")
mol.ind  = c("gsl.IM.blup","gsl.1moIM.blup","gsl.1hIM.blup","gsl.4moIM.blup")

# order molecules by subcategory
dat = dat[,c(mol.mt,mol.ms,mol.alk,mol.halk,mol.ind)]

# rename some molecules for consistency with other datasets
mol.mt   = c("7mSh","8mSo")
mol.ms   = c("3mSOp","4mSOb","5mSOp","6mSOh","7mSOh","8mSOo")
mol.alk  = c("Pren","Buen","Peen")
mol.halk = c("2hBuen","2hPeen")
mol.ind  = c("IM","1moIM","1hIM","4moIM")

# name columns
colnames(dat) = c(mol.mt,mol.ms,mol.alk,mol.halk,mol.ind)

# calculate correlations
dat.cor = cor(dat, method = "spearman")

# create corrPlot
pdf("gsl_correlogram_toua.pdf", h = 4.75, w = 5.25)
corrplot(dat.cor, method="circle")
dev.off()

# optionally, plot only the molecules shared across all studies
shared.tou = dat.cor[shared.mols,shared.mols]

### Brachi

dat = read.csv("/brachi/phen_glucosinolatesBLUPs_2021-07-28.csv")

mol.mt   = c("G3MTP.blup","G4MTB.blup","G7MTH.blup","G8MTO.blup")
mol.ms   = c("G3MSP.blup","G4MSB.blup","G5MSP.blup","G6MSH.blup","G7MSH.blup","G8MSO.blup")
mol.halkyl = c("G3HP.blup","G4HB.blup")
mol.alkenyl  = c("G2P.blup","G3B.blup","G4P.blup")
mol.halkenyl = c("G2H3B.blup","G2H4P.blup")

dat = dat[,c(mol.mt,mol.ms,mol.halkyl,mol.alkenyl,mol.halkenyl)]

# new names
mol.mt   = c("3mSp","4mSb","7mSh","8mSo")
mol.ms   = c("3mSOp","4mSOb","5mSOp","6mSOh","7mSOh","8mSOo")
mol.halkyl = c("3HP","4HP")
mol.alkenyl  = c("Pren","Buen","Peen")
mol.halkenyl = c("2hBuen","2hPeen")

colnames(dat) = c(mol.mt,mol.ms,mol.halkyl,mol.alkenyl,mol.halkenyl)

dat.cor = cor(dat, method = "spearman")

pdf("gsl_correlogram_brachi.pdf", h = 4.75, w = 5.25)
corrplot(dat.cor, method="circle")
dev.off()

shared.bra = dat.cor[shared.mols,shared.mols]


### Katz

dat = read.csv("/katz/phen_glucosinolatesBLUPs_2021-07-30.csv")

mol.mt   = c("G3MTP.blup","G4MTB.blup","G7MTH.blup","G8MTO.blup")
mol.ms   = c("G3MSP.blup","G4MSB.blup","G5MSP.blup","G6MSH.blup","G7MSH.blup","G8MSO.blup")
mol.halkyl = c("G3HP.blup","G4HB.blup")
mol.alkenyl  = c("G2P.blup","G3B.blup","G4P.blup")
mol.halkenyl = c("G2H3B.blup")
mol.ind  = c("GI3M.blup")

dat = dat[,c(mol.mt,mol.ms,mol.halkyl,mol.alkenyl,mol.halkenyl,mol.ind)]

# new names
mol.mt   = c("3mSp","4mSb","7mSh","8mSo")
mol.ms   = c("3mSOp","4mSOb","5mSOp","6mSOh","7mSOh","8mSOo")
mol.halkyl = c("3HP","4HP")
mol.alkenyl  = c("Pren","Buen","Peen")
mol.halkenyl = c("2hBuen")
mol.ind  = c("IM")

colnames(dat) = c(mol.mt,mol.ms,mol.halkyl,mol.alkenyl,mol.halkenyl,mol.ind)

dat.cor = cor(dat, method = "spearman")

pdf("gsl_correlogram_katz.pdf", h = 4.75, w = 5.25)
corrplot(dat.cor, method="circle")
dev.off()

shared.kat = dat.cor[shared.mols,shared.mols]


### Wu

dat = read.csv("/wu/phen_glucosinolatesMeans_2021-07-30.csv")

mol.mt   = c("G3MTPnm","G4MTBnm","G7MTHpm","G8MTOnm")
mol.ms   = c("G3MSPpm","G4MSBpm","G5MSPnm","G6MSHnm","G7MSHpm","G8MSOpm")
mol.halkyl = c("G3HPpm","G4HBnm")
mol.alkenyl  = c("G2Ppm","G3Bpm","G4Pnm")
mol.halkenyl = c("G2H3Bpm")
mol.ind  = c("GI3Mpm","G1MOI3Mpm","G4HI3Mnm","G4MOI3Mpm")

dat = dat[,c(mol.mt,mol.ms,mol.halkyl,mol.alkenyl,mol.halkenyl,mol.ind)]

# new names
mol.mt   = c("3mSp","4mSb","7mSh","8mSo")
mol.ms   = c("3mSOp","4mSOb","5mSOp","6mSOh","7mSOh","8mSOo")
mol.halkyl = c("3HP","4HP")
mol.alkenyl  = c("Pren","Buen","Peen")
mol.halkenyl = c("2hBuen")
mol.ind  = c("GI3M","G1MOI3M","G4HI3M","G4MOI3M")

colnames(dat) = c(mol.mt,mol.ms,mol.halkyl,mol.alkenyl,mol.halkenyl,mol.ind)

dat.cor = cor(dat, method = "spearman")

pdf("/gsl_correlogram_wu.pdf", h = 4.75, w = 5.25)
corrplot(dat.cor, method="circle")
dev.off()

shared.wu = dat.cor[shared.mols,shared.mols]






