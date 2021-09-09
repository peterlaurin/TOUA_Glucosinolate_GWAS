library("lme4")
library("tidyr")
library("plyr")
library("HLMdiag")

# load integrated peak areas from HPLC-MS/MS
temp = list.files(path = "/glucFiles_mzML/integratedOutput/", full.names = T, pattern = "*.csv")
myfiles = lapply(temp, read.delim, sep = ",")
datMS = rbind.fill(myfiles)[,c("sample","peakname","slice_intb")]

# get sample metadata for metabolite profiling runs, filter samples and columns
datMS = separate(datMS, sample, into = c("Set", "Plate"), sep = "_plate")
datMS = separate(datMS, Plate, into = c("Plate", "SampleName"), sep = "_P1_") # NA's from external standards, don't have sample names
                                                                              # Note: there are some P2's that may have been excluded due to this!
                                                                              # e.g., "set2_7_P2_A01_MS1converted.mzML"
datMS = separate(datMS, Plate, into = c("Plate", "StandardName"), sep = "_r") # NA's from samples, don't have external standard names
datMS$TrayMS = paste(datMS$Set, datMS$Plate, sep = "_")
datMS = separate(datMS, SampleName, into = c("PositionMS", "ext"), sep = "_MS1converted")
datMS[c(1:3, 5)] <- NULL # remove un-needed columns
datMSsubset = subset(datMS, PositionMS != F) # removes standards

# merge with additional metadata about plant growth
dat_harvest = rbind(
  read.csv("ExtractionCD_ADG.csv"),
  read.csv("ExtractionEF_ADG.csv")
)[,c("ID_short","Day","Block","TrayNum","Weight","Plaque","Ligne","Colonne")]

dat_harvest$Colonne = formatC(dat_harvest$Colonne, width = 2, flag = "0")

dat_harvest$TrayMS = paste0("set",dat_harvest$Day,"_",dat_harvest$Plaque)
dat_harvest$PositionMS = paste0(dat_harvest$Ligne,dat_harvest$Colonne)

dat_harvest$TrayPlanting = paste0("Block",dat_harvest$Block,"_Tray",dat_harvest$TrayNum)

dat = merge(dat_harvest, datMSsubset, by = c("TrayMS", "PositionMS"))

# adjust for biomass per sample
dat$slice_intb_adj = dat$slice_intb / dat$Weight






##### 2. RUN LMMs FOR EACH GSL MOLECULE #####

# create table to store some model output when initiating a loop to run the LMM on each molecule

model.info = data.frame(molecule = unique(dat$peakname),
                        intercept = NA,
                        p = NA,
                        pve_id = NA)

# create a data frame to store BLUPs (effect of each plant ID) from these models

gls.blups = data.frame(ID = unique(dat$ID))

# create lists to store the model objects created by the LMMs with and without plant ID as a factor

gls.models.full = list()
gls.models.null = list()

# loop over the glucosinolate molecule names

for (i in (1:nrow(model.info))){
  
  # subset to focal GSL molecule
  temp = subset(dat, peakname == model.info$molecule[i])
  
  # run LMM with plant ID as a factor (alternative model)
  m1a = lmer(log10(slice_intb_adj -44) ~ (1|ID_short) + (1|TrayMS) + (1|TrayPlanting), data = temp)
  gls.models.full[[i]] = m1a
  
  # run LMM without plant ID as a factor (null model)
  m1n = lmer(log10(slice_intb_adj -44) ~ (1|TrayMS) + (1|TrayPlanting), data = temp)
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

write.table(model.info, "/toua_wFrachonGenotypes/model_info.txt", sep = "\t", quote = F, eol = "\n", row.names = F)
write.table(gls.blups, "/toua_wFrachonGenotypes/phen_glucosinolatesBLUPs_2021-08-02.csv", sep = ",", quote = F, eol = "\n", row.names = F)











