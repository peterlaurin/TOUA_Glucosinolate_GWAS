# Purpose: Create local Manhattan plots for the candidate regions harboring significant SNPs
#          as shown in supplementary material.

library("ggplot2")
library("gggenes")

# gggenes links
# https://cran.r-project.org/web/packages/gggenes/vignettes/introduction-to-gggenes.html
# https://cran.r-project.org/web/packages/gggenes/gggenes.pdf
# https://rdrr.io/cran/gggenes/man/example_genes.html

# test example from "gggenes" package
ggplot(example_genes, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")


##### 1. Modify GFF file to format for "gggenes" package #####

gff.tair10 = read.delim("TAIR10_GFF3_genes.gff", sep = "\t", h = F)

# get genes only file 
# andy@10-17-12-132 PhilTrans_Revisions % grep "TAIR10    gene" TAIR10_GFF3_genes.gff > TAIR10_GFF3_genes_onlyGenes.gff  

# modify genes only file
# find: (Chr.)\tTAIR10\tgene\t(\d+)\t(\d+)\t\.\t(.)\t\.\tID=(AT.G\d+);\S+\n
# replace: \1\t\5\t\2\t\3\t\4\tforward\n
# header: molecule	gene  start	end orientation strand

# get CDS only file
# andy@10-17-12-132 PhilTrans_Revisions % grep "TAIR10    CDS" TAIR10_GFF3_genes.gff > TAIR10_GFF3_genes_onlyCDS.gff   

# modify CDS only file
# find: (Chr.)\tTAIR10\tCDS\t(\d+)\t(\d+)\t\.\t(.)\t\d\tParent=(AT.G\d+)\.\d\,(AT.G\d+)\.(\d+)\-Protein;\n
# replace: \1\t\5\t\6\-\7\t\2\t\3\n
# find: \+
# replaxe: 1
# find: \-
# replaxe: -1
# header: molecule	gene	subgene	from	to

gff.genes = read.delim("TAIR10_GFF3_genes_onlyGenes.tsv", sep = "\t", h = T)
gff.cds   = read.delim("TAIR10_GFF3_genes_onlyCDS.tsv", sep = "\t", h = T)
gff.all   = merge(gff.genes, gff.cds, by = c("molecule","gene"))

gff.all = subset(gff.all, select = c(molecule,gene,start,end,strand,subgene,from,to,orientation))
gff.all$from = as.numeric(gff.all$from)

##### 2. Get gggenes-style coordinates for gene regions of interest #####

gff.gene.pos = read.delim("TAIR10_GFF3_geneBoundariesTableADG.txt")
candidate.genes = read.delim("pathway_BiosyntheticGenes_HarunTable2.txt")
gff.gene.pos.cand = merge(gff.gene.pos, candidate.genes, by = "gene_id")
gff.gene.pos.cand.sig = subset(gff.gene.pos.cand, gene_abbreviation %in% c("MAM1","BCAT3","CYP79F2","CYP83A1","FMOGS-OX1","AOP2","GSL-OH","IGMT2","CYP81F4","CYP81F2"))
gff.gene.pos.cand.sig$pos = floor((gff.gene.pos.cand.sig$ps_left + gff.gene.pos.cand.sig$ps_right) / 2)

gff.all.candidates = setNames(data.frame(matrix(ncol=length(colnames(gff.all)),nrow=0)), colnames(gff.all))

window.size = 2e5

for (i in 1:nrow(gff.gene.pos.cand.sig)){
  
  temp = subset(gff.all, molecule == paste0("Chr",gff.gene.pos.cand.sig[i,"chr"]) & 
                  start > gff.gene.pos.cand.sig[i,"pos"] - window.size/2 & 
                  end < gff.gene.pos.cand.sig[i,"pos"] + window.size/2)
  
  temp$molecule = gff.gene.pos.cand.sig[i,"gene_abbreviation"]
  
  # add 
  temp.dummy = data.frame(molecule = temp[1,"molecule"],
                          gene     = "FAKE",
                          start    = gff.gene.pos.cand.sig[i,"pos"] - window.size/2,
                          end      = gff.gene.pos.cand.sig[i,"pos"] + window.size/2,
                          strand   = "forward",
                          subgene  = "FAKE-1",
                          from     = gff.gene.pos.cand.sig[i,"pos"] - window.size/2,
                          to       = gff.gene.pos.cand.sig[i,"pos"] + window.size/2,
                          orientation = 1
                          )
  temp = rbind(temp.dummy,temp)
  
  # iteratively append to models for all candidate gene regions
  gff.all.candidates = rbind(gff.all.candidates, temp)
  
}

##### 3. Plot using "gggenes" #####

temp = gff.all.candidates

# shaded genes
ggplot(temp, aes(xmin = start, xmax = end, y = molecule, fill = gene, forward = orientation)) +
     geom_gene_arrow(show.legend = F, arrowhead_width = grid::unit(1, "mm")) +
     facet_wrap(~ molecule, scales = "free", ncol = 1) +
     theme_genes()

# # shaded cds regions
# ggplot(temp, aes(xmin = start, xmax = end, y = molecule, forward = orientation)) +
#   facet_wrap(~ molecule, scales = "free", ncol = 1) +
#   geom_gene_arrow(fill = "white", show.legend = F, arrowhead_width = grid::unit(1, "mm")) +
#   geom_subgene_arrow(data = temp, show.legend = F, arrowhead_width = grid::unit(1, "mm"), 
#                      aes(xmin = start, xmax = end, y = molecule, fill = gene,
#                          xsubmin = from, xsubmax = to), color="black", alpha=.7) +
#   theme_genes()

##### 4. Find and plot locus x phenotype with top-scoring GWAS association for each significant gene #####

# Locus              TOU-A       Brachi      Katz        Wu         
# GS-OH (AT2G25450)  ratio17     BlupG2H4P   ratio17     G2H3Bpm
# AOP2  (AT4G03060)  ratio15     BlupG3HP    BlupG2P     G3HPpm
# MAM1  (AT5G23010)  ratio2      ratio7      BlupG4MTB   G5MSPnm

# Locus                 TOU-A           Other
# CYP79F2 (AT1G16400)   ratio13         Katz (G3HP)
# GS-OX1  (AT1G65860)   BlupG2H4P       x
# BCAT3   (AT3G49680)   ratio10         x
# CYP83A1 (AT4G13770)   BlupG2H4P       x
# IGMT2   (AT1G21120)   ratio21         x
# CYP81F4 (AT4G37410)   mvlmm3          x
# CYP81F2 (AT5G57220)   mvlmm3          Wu (ratio19)

setwd("/TOUA_glucs/")

files.to.plot = read.csv("individual_gwas_to_plot.csv")
#files.to.plot = read.csv("individual_gwas_to_plot_N192.csv")

files.to.plot$order = 1:nrow(files.to.plot)
files.to.plot = merge(files.to.plot, gff.gene.pos.cand.sig[,c("gene_abbreviation","chr","pos")])
files.to.plot = files.to.plot[order(files.to.plot$order), ]

for (i in 1:nrow(files.to.plot)){
  
  assoc.file = read.delim(paste0(files.to.plot[i,"directory"], files.to.plot[i,"file"]))
  
  cutoff.bon = 0.05 / nrow(subset(assoc.file, af > 0.05))
  
  assoc.file = subset(assoc.file, af > 0.05 & chr == files.to.plot[i,"chr"] &
                        ps > files.to.plot[i,"pos"] - window.size/2 &
                        ps < files.to.plot[i,"pos"] + window.size/2 &
                        p_wald < 0.1, ylab = "-log10p"
                      )
  
  # adjust labels to units of x10^6 (Mbp)
  assoc.file$ps = assoc.file$ps / 1e6
  
  pdf(file = paste0("manhat_local_",
                    files.to.plot[i,"mapping_panel"],"_",
                    files.to.plot[i,"gene_abbreviation"],"_",
                    files.to.plot[i,"file"],".pdf"),
      height = 4, width = 8)
  
  ### ADD XLIM based on pos +/- window.size/2
  
  plot(-log10(assoc.file$p_wald) ~ assoc.file$ps, pch = 16, col = "dodgerblue2",
       ylim = c(0, max(-log10(assoc.file$p_wald)*1.2)),
       xlim = c((files.to.plot[i,"pos"] - window.size/2) / 1e6, (files.to.plot[i,"pos"] + window.size/2) / 1e6)
       )
  
  abline(h = -log10(cutoff.bon), col = "red", lty = 2, lwd = 2)
  abline(v = (files.to.plot[i,"pos"] - window.size/2) / 1e6, col = "gold")
  abline(v = (files.to.plot[i,"pos"] + window.size/2) / 1e6, col = "gold")
  
  dev.off()
  
}
