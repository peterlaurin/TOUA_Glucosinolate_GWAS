##### 1. LOAD DATA #####

library(dplyr)

# load data from Wu et al. supplement, control condition
phen = read.delim("/wu/metabolites_data_normalized_control_wu.csv", h=T, sep=",", na.string="NA")

gluc.key = read.delim("/wu/metabolites_names_control_wu.csv", h=T, sep=",", na.string="NA")
gluc.key = subset(gluc.key, gluc_id != "" & gluc_id %in% colnames(phen))

phen = cbind( phen$ID, phen[, colnames(phen) %in% gluc.key$gluc_id] )

colnames(phen)[1] = "ID"

phen = phen %>% rename_at(vars(as.character(gluc.key$gluc_id)), ~ as.character(gluc.key$gluc_abbreviation))

phen[is.na(phen)] = 0

write.table(phen, "/wu/phen_glucosinolatesMeans_2021-07-30.csv", sep = ",", quote = F, eol = "\n", row.names = F)



