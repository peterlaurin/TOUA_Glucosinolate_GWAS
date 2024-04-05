# Purpose: Create table of percent variance explained (PVE) for leading SNPs.

##### 1. GET KINSHIP MATRIX FOR TOU-A ACCESSIONS #####

# get order of genotypes in kinship matrix
fam.file.IDs = read.csv("/pve/genotypes.fam", sep = " ", header = F)[,1]

# load kinship matrix
kmat = read.delim("/pve/dat_tot_forGemma_kmat.cXX.txt", sep = "\t", header = F)

# assign accession IDs as row/column names in kinship matrix

kmat = as.matrix(kmat)
rownames(kmat) = fam.file.IDs
colnames(kmat) = fam.file.IDs



##### 2. GET TOU-A GENOTYPES FOR LEADING SNP AT 10 SIGNIFICANT LOCI #####

### A. Get VCF file of genotypes.

#install.packages("vcfR")
library("vcfR")

vcf = read.vcfR(file = "/pve/Andy_positions.recode.vcf")

### B. Convert to table of genotypes

#install.packages("adegenet")
library("adegenet")

genind = data.frame( vcfR2genind(vcf)@tab )
genind$SAMPLE = rownames(genind)

### C. Re-order genotype file to match kinship matrix

tou.accessions = read.csv("/pve/TOU_accessions.csv")
tou.accessions = subset(tou.accessions, ecotype_id %in% fam.file.IDs)

geno = dplyr::left_join(tou.accessions[,c("SAMPLE","ecotype_id")], genind, by = "SAMPLE")
geno = dplyr::left_join(data.frame(ecotype_id = fam.file.IDs), geno, by = "ecotype_id" )

# View(data.frame(geno = geno$ecotype_id,
#                 fam  = fam.file.IDs,
#                 kin  = rownames(kmat)))



##### 3. CALCULATE % VARIANCE EXPLAINED (PVE) PER SNP WITH LME4QTL ######

# install.packages("remotes")
# remotes::install_github("variani/lme4qtl")
library("lme4qtl")

# load phenotypes

pheno.blup  = read.csv("/pve/phen_glucosinolatesBLUPs_2021-08-02.csv")
pheno.ratio = read.csv("/pve/phen_ratios.csv")


# merge phenotypes
phen = dplyr::left_join(pheno.blup, pheno.ratio, by = "ID")

# merge phenotype with genotypes
dat = dplyr::left_join(geno, phen, by = "ecotype_id")

# View(data.frame(geno = geno$ecotype_id,
#                 fam  = fam.file.IDs,
#                 kin  = rownames(kmat),
#                 dat  = dat$ecotype_id))

dat$ecotype_id = as.character(dat$ecotype_id)


# model


m1a = relmatLmer(ratio2 ~ X5_7698197.0 + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
# HLMdiag::varcomp.mer(m1a)[2] / sum(HLMdiag::varcomp.mer(m1a)) # explained by kmat
# car::Anova(m1a)
# lme4::fixef(m1a)
# MuMIn::r.squaredGLMM(m1a)[1]



# CYP79F1/2
m1a = relmatLmer(ratio13 ~ (1|X1_5609258.0) + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
HLMdiag::varcomp.mer(m1a)[3] / sum(HLMdiag::varcomp.mer(m1a)) # explained by SNP 0.8260421
subset ( read.delim("lmm_maf03_ratio13.assoc.fixed.txt"),
         rs == "1_5609258")
# chr        rs      ps n_miss allele1 allele0    af     beta        se   logl_H1  l_remle       p_wald
# 21754   1 1_5609258 5609258     19       T       A 0.075 1.556184 0.1547353 -119.9017 23.05345 2.391523e-19
m1a.f = relmatLmer(ratio13 ~ X1_5609258.0 + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
MuMIn::r.squaredGLMM(m1a.f)[1] # 0.399822

# IGMT2
m1a = relmatLmer(ratio21 ~ (1|X1_7398643.0) + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
HLMdiag::varcomp.mer(m1a)[3] / sum(HLMdiag::varcomp.mer(m1a)) # explained by SNP 0.3331985
subset ( read.delim("lmm_maf03_ratio21.assoc.fixed.txt"),
         rs == "1_7398643")
# chr        rs      ps n_miss allele1 allele0    af       beta         se   logl_H1  l_remle         p_wald
# 29000   1 1_7398643 7398643     14       A       G 0.483 -0.3560667 0.05340043 -114.2267 3.294033   2.746939e-10
m1a.f = relmatLmer(ratio21 ~ X1_7398643.0 + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
MuMIn::r.squaredGLMM(m1a.f)[1] # 0.2041138

# GS-OX1 <--- chose ratio11 (7MTH:7MSH) to reflect known GS-OX1 function (7MTH:7MSH had better p-value than 8MTO:8MSO)
m1a = relmatLmer(ratio11 ~ (1|X1_24516880.0) + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
HLMdiag::varcomp.mer(m1a)[3] / sum(HLMdiag::varcomp.mer(m1a)) # explained by SNP 0.8150414
subset ( read.delim("lmm_maf03_ratio11.assoc.fixed.txt"),
         rs == "1_24516880")
# chr         rs       ps n_miss allele1 allele0    af      beta         se   logl_H1   l_remle       p_wald
# 125729   1 1_24516880 24516880     11       T       C 0.033 0.5411118 0.07454482 -67.10372 0.1462106 9.711327e-12
m1a.f = relmatLmer(ratio11 ~ X1_24516880.0 + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
MuMIn::r.squaredGLMM(m1a.f)[1] # 0.2238156

# GS-OH
m1a = relmatLmer(ratio17 ~ (1|X2_10795336.0) + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
HLMdiag::varcomp.mer(m1a)[3] / sum(HLMdiag::varcomp.mer(m1a)) # explained by SNP 0.5841938
subset ( read.delim("lmm_maf03_ratio17.assoc.fixed.txt"),
         rs == "2_10795336")
# chr         rs       ps n_miss allele1 allele0    af      beta         se   logl_H1  l_remle       p_wald
# 218355   2 2_10795336 10795336     16       T       C 0.062 0.4930957 0.06917878 -64.11902 3.745845 2.065129e-11
m1a.f = relmatLmer(ratio17 ~ X2_10795336.0 + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
MuMIn::r.squaredGLMM(m1a.f)[1] # 0.1442538

# BCAT3
m1a = relmatLmer(ratio10 ~ (1|X3_18442070.0) + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
HLMdiag::varcomp.mer(m1a)[3] / sum(HLMdiag::varcomp.mer(m1a)) # explained by SNP 0.5038035
subset ( read.delim("lmm_maf03_ratio10.assoc.fixed.txt"),
         rs == "3_18442070")
# chr         rs       ps n_miss allele1 allele0   af      beta         se   logl_H1  l_remle       p_wald
# 346596   3 3_18442070 18442070     13       T       C 0.24 0.3920634 0.04420245 -29.68121 6.602609 5.356244e-16
m1a.f = relmatLmer(ratio10 ~ X3_18442070.0 + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
MuMIn::r.squaredGLMM(m1a.f)[1] # 0.2740293

# AOP2
m1a = relmatLmer(ratio15 ~ (1|X4_1248628.0) + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
HLMdiag::varcomp.mer(m1a)[3] / sum(HLMdiag::varcomp.mer(m1a)) # explained by SNP 0.4155394
subset ( read.delim("lmm_maf03_ratio15.assoc.fixed.txt"),
         rs == "4_1248628")
# chr        rs      ps n_miss allele1 allele0    af      beta         se   logl_H1  l_remle       p_wald
# 375897   4 4_1248628 1248628     17       C       T 0.103 0.4197446 0.06755447 -106.9232 4.875777 3.205191e-09
m1a.f = relmatLmer(ratio15 ~ X4_1248628.0 + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
MuMIn::r.squaredGLMM(m1a.f)[1] # 0.1193737

# CYP83A1
m1a = relmatLmer(gsl.2hPeen.blup ~ (1|X4_8000492.0) + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
HLMdiag::varcomp.mer(m1a)[3] / sum(HLMdiag::varcomp.mer(m1a)) # explained by SNP 0.3893483
subset ( read.delim("lmm_maf03_BlupG2H4P.assoc.fixed.txt"),
         rs == "4_8000492")
# chr        rs      ps n_miss allele1 allele0    af      beta        se   logl_H1  l_remle       p_wald
# 412815   4 4_8000492 8000492     10       C       A 0.055 -1.152461 0.1744192 -230.8767 32.86944 3.831281e-10
m1a.f = relmatLmer(gsl.2hPeen.blup ~ X4_8000492.0 + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
MuMIn::r.squaredGLMM(m1a.f)[1] # 0.06386437

# CYP81F1/3/4
m1a = relmatLmer(ratio18 ~ (1|X4_17567085.0) + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
HLMdiag::varcomp.mer(m1a)[3] / sum(HLMdiag::varcomp.mer(m1a)) # explained by SNP 0.2719625
subset ( read.delim("lmm_maf03_ratio18.assoc.fixed.txt"),
         rs == "4_17567085")
# chr         rs       ps n_miss allele1 allele0    af       beta         se   logl_H1  l_remle       p_wald
# 465192   4 4_17567085 17567085     10       A       C 0.379 -0.2578379 0.04038493 -45.27507 8.421293 1.286664e-09
m1a.f = relmatLmer(ratio18 ~ X4_17567085.0 + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
MuMIn::r.squaredGLMM(m1a.f)[1] # 0.1534442

# MAM1/2
m1a = relmatLmer(ratio2 ~ (1|X5_7698197.0) + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
HLMdiag::varcomp.mer(m1a)[3] / sum(HLMdiag::varcomp.mer(m1a)) # explained by SNP 0.5982064
subset ( read.delim("lmm_maf03_ratio2.assoc.fixed.txt"),
         rs == "5_7698197")
# chr        rs      ps n_miss allele1 allele0    af       beta         se   logl_H1  l_remle       p_wald
# 508831   5 5_7698197 7698197      8       G       T 0.505 -0.8717641 0.07423788 -173.8975 3.497276 2.715677e-24
m1a.f = relmatLmer(ratio2 ~ X5_7698197.0 + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
MuMIn::r.squaredGLMM(m1a.f)[1] # 0.4296903

# CYP81F2
m1a = relmatLmer(gsl.4moIM.blup ~ (1|X5_23187786.0) + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
HLMdiag::varcomp.mer(m1a)[3] / sum(HLMdiag::varcomp.mer(m1a)) # explained by SNP 0.1837962
subset ( read.delim("lmm_maf03_BlupG4MOIM.assoc.fixed.txt"),
         rs == "5_23187786")
# chr         rs       ps n_miss allele1 allele0   af       beta         se   logl_H1  l_remle       p_wald
# 600825   5 5_23187786 23187786     16       T       A 0.29 -0.4984922 0.09623073 -210.1594 9.492969 5.628478e-07
m1a.f = relmatLmer(gsl.4moIM.blup ~ X5_23187786.0 + (1|ecotype_id), data = dat, relmat = list(ecotype_id = kmat), REML = T)
MuMIn::r.squaredGLMM(m1a.f)[1] # 0.08857372




