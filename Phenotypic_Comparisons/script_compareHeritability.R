# Purpose: Create plot comparing heritability values per molecule between TOU-A and Katz, Brachi
#
# Note:    The outputted figure was manually refined in Adobe Illustrator.

# load files with H^2 from LMMs for each mapping population
h2.tou = read.delim("model_info.txt")
h2.bra = read.delim("/brachi/model_info.txt")
h2.kat = read.delim("/katz/model_info.txt")

# rename molecules in TOU-A file to match other files
h2.tou$molecule = c("G6MSH","G1MOI3M","G2H3B","G4MSB","G4HI3M","G8MTO","exclude","exclude","G4P",
                    "G3MSP","exclude","exclude","exclude","exclude","exclude","G4MOI3M","exclude",
                    "exclude","G7MTH","GI3M","G5MSP","G2H4P","G7MSH","exclude","G3B","G8MSO","G2P")

h2.tou = subset(h2.tou, molecule != "exclude")

# combine in wide format
h2.all = merge(h2.tou[,c("molecule","p","pve_id")], h2.bra[,c("molecule","p","pve_id")], by = "molecule", all = T)
h2.all = merge(h2.all, h2.kat[,c("molecule","p","pve_id")], by = "molecule", all = T)
colnames(h2.all) = c( "molecule", paste0( c("p","pve"), ".tou"), paste0( c("p","pve"), ".bra"), paste0( c("p","pve"), ".kat"))

# If mol. isn't measured in Brachi or Katz, set PVE = 1.1
h2.all[is.na(h2.all$p.bra) & is.na(h2.all$p.kat),"pve.bra"] = 1.1

# plot H^2 values in TOU-A (y) vs. Brachi, Katz (x)
plot  (h2.all$pve.tou ~ h2.all$pve.bra, xlim = c(0,1.2), ylim = c(0,1), pch = 16, col = "dodgerblue") # brachi = blue points
points(h2.all$pve.tou ~ h2.all$pve.kat, pch = 16, col = "gold") # katz = gold points
abline(0,1)

h2.all[is.na(h2.all$p.bra) & is.na(h2.all$p.kat),]

# make histograms of H^2 distributions to align next to x and y axes with Illustrator
hist(subset(h2.all, pve.tou != "NA")$pve.bra, breaks = 5, xlim = c(0,1))
hist(subset(h2.all, pve.tou != "NA")$pve.kat, breaks = 10, xlim = c(0,1))
hist(subset(h2.all, pve.tou != "NA")$pve.tou, breaks = 10, xlim = c(0,1))

# sign test to ask if heritability significantly differs between TOU-A and Brachi, Katz
library("BSDA")
SIGN.test(x = h2.all$pve.tou, y = h2.all$pve.bra)
SIGN.test(x = h2.all$pve.tou, y = h2.all$pve.kat)

sort(na.omit(h2.all$pve.tou - h2.all$pve.bra))
sort(na.omit(h2.all$pve.tou - h2.all$pve.kat))







