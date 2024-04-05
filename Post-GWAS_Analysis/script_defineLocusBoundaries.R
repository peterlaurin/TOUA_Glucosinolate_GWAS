# PURPOSE: Define locus boundaries around each glucosinolate biosynthetic gene. These extended windows will
#          be considered the "gene regions" for downstream analysis of GWAS output, e.g. to associate significant 
#          SNPs with biosynthetic genes. Boundaries are defined as a 30kb extension of the left/right boundaries
#          of the TAIR10 gene model for most genes; for AOP, GSOH, and MAM, boundaries are further extended until
#          to capture 90% of the significant SNPs (Brachi/Katz/Wu, all univariate GWAS for single molecule BLUPs)
#          in a 1 Mb window centered on the gene.

# define minor allele frequency (maf) cutoff for downstream analyses
maf.thresh = 0.05

##### 1. DEFINE BOUNDARIES FOR GSOH/AOP/MAM #####

# specify folders with GWAS outputs for single molecule BLUPs for Brachi/Katz/Wu datasets
my.directories = c("/output_gwas_assocs_lmm/brachi_full/blup_aliphatic_assoc/",
                   "/output_gwas_assocs_lmm/katz_full/blup_aliphatic_assoc/",
                   "/output_gwas_assocs_lmm/wu_full/blup_aliphatic_assoc/"
                   )

# list filenames in specified folders
my.filenames = list.files(path = my.directories, pattern = "assoc.fixed.txt", full.names = T)

# create dataframe that will hold significant SNPs from Brachi/Katz/Wu
gwas.sig.snps = data.frame(matrix(nrow = 0, ncol = 5))
colnames(gwas.sig.snps) = c("analysis","chr","ps","af","p_wald")                          

# extract Bonferonni-corrected significant SNPs from each GWAS
for (i in 1:length(my.filenames)){

  print( paste0("Processing file ", i, " of ", length(my.filenames), "...") )
  
  gwas.i = read.delim(my.filenames[i])
  
  gwas.i = subset(gwas.i, p_wald < .05/nrow(gwas.i) & af > maf.thresh)
  
  if (nrow(gwas.i) > 0){
    
    gwas.i$analysis = my.filenames[i]
    
    gwas.sig.snps = rbind(gwas.sig.snps, gwas.i[,c("analysis","chr","ps","af","p_wald")])
    
  }
  
}

# # optional: write output to save for later
# write.table(gwas.sig.snps, "sig_snps_maf05_blupsAliphatic_brachiKatzWuFull.csv", sep = ",", row.names = F, quote = F)


library("plyr")

# # optional: load significant SNP list if starting new R session
# gwas.sig.snps = read.delim("sig_snps_maf05_blupsAliphatic_brachiKatzWuFull.csv", sep = ",")

# bin significant SNPs into 10kb windows
gwas.sig.snps$ps.bin = cut(gwas.sig.snps$ps, breaks = seq(0,3e7, by = 1e4), labels = seq(1e4,3e7, by = 1e4))
gwas.sig.snps$ps.bin = as.numeric(as.character(gwas.sig.snps$ps.bin))

# define positions of GSOH/AOP/MAM loci; then run code below for each locus separately
gene.chr = 2; gene.pos = 10.83e6 #GSOH
gene.chr = 4; gene.pos = 1.35e6  #AOP2
gene.chr = 5; gene.pos = 7.70e6  #MAM1

# run this separately for each locus

    # get significant within 1 Mbp of the locus 
    gwas.sig.snps.chr.i = subset(gwas.sig.snps, chr == gene.chr & ps > gene.pos - 5e5 & ps < gene.pos + 5e5)
    gwas.sig.snps.chr.i.per.bin = ddply(gwas.sig.snps.chr.i, .(ps.bin), summarize, count = length(chr))
    
    # determine proportion of sig. SNPs falling within each 10 kbp bin
    gwas.sig.snps.chr.i.per.bin$proportion = gwas.sig.snps.chr.i.per.bin$count / sum(gwas.sig.snps.chr.i.per.bin$count)
    
    # find window centered on gene that captures 90% of sig. SNPs
    # try increments of 0.01e6 (=10kb)
    for (n in seq(0.01e6, 0.5e6, by = 0.01e6)){
      prop.captured = sum( subset(gwas.sig.snps.chr.i.per.bin, ps.bin %in% seq(gene.pos-n+1e4, gene.pos+n, by = 1e4))$proportion ) # left window +1e4 b/c bins named by right border
      print(paste0("window: ",n ,"     proportion sig. SNPs captured: ", prop.captured))
    }
    
    # Locus   90% captured    95% captured   [<-- report half window sizes]
    # GSOH    80,000          90,000
    # AOP2    120,000         150,000
    # MAM1    110,000         170,000


##### 2. PERFORM GENE BOUNDARY ADJUSTMENTS #####


# note: "major" loci are...

    # number	gene_id	  gene_abbreviation	  gene_name
    # 1	      AT5G23010	MAM1	              methylthioalkylmalate synthase 1
    # 2	      AT5G23020	MAM3	              methylthioalkylmalate synthase 3
    # 46	    AT4G03060	AOP2              	2-oxoglutarate-dependent dioxygenase
    # 47    	AT4G03050	AOP3	              2-oxoglutarate-dependent dioxygenase
    # 48	    AT2G25450	GSL-OH	            2-oxoacid-dependent dioxygenase

### load files

gene.boundaries = read.delim("TAIR10_GFF3_geneBoundariesTableADG.txt")

gsl.genes = read.delim("pathway_BiosyntheticGenes_HarunTable2.txt")
gsl.genes = subset(gsl.genes, aliphatic == "yes" | indolic == "yes")

major.loci  = c("AT2G25450","AT4G03050","AT4G03060","AT5G23010","AT5G23020")
other.genes = subset(gsl.genes, !(gene_id %in% major.loci))$gene_id

### adjust major loci

# initial boundaries
subset(gene.boundaries, gene_id %in% major.loci)

    #       chr  ps_left ps_right   gene_id
    # 1       4  1351980  1353854 AT4G03060
    # 9430    2 10829916 10831655 AT2G25450
    # 17997   4  1344182  1346501 AT4G03050
    # 24201   5  7703092  7706896 AT5G23010
    # 24202   5  7718121  7721839 AT5G23020

# adjustment

gene.boundaries[gene.boundaries$gene_id == "AT2G25450",]$ps_left  = gene.boundaries[gene.boundaries$gene_id == "AT2G25450",]$ps_left  - 8e4
gene.boundaries[gene.boundaries$gene_id == "AT2G25450",]$ps_right = gene.boundaries[gene.boundaries$gene_id == "AT2G25450",]$ps_right + 8e4

gene.boundaries[gene.boundaries$gene_id == "AT4G03050",]$ps_left  = gene.boundaries[gene.boundaries$gene_id == "AT4G03050",]$ps_left  - 12e4
gene.boundaries[gene.boundaries$gene_id == "AT4G03050",]$ps_right = gene.boundaries[gene.boundaries$gene_id == "AT4G03050",]$ps_right + 12e4

gene.boundaries[gene.boundaries$gene_id == "AT4G03060",]$ps_left  = gene.boundaries[gene.boundaries$gene_id == "AT4G03060",]$ps_left  - 12e4
gene.boundaries[gene.boundaries$gene_id == "AT4G03060",]$ps_right = gene.boundaries[gene.boundaries$gene_id == "AT4G03060",]$ps_right + 12e4

gene.boundaries[gene.boundaries$gene_id == "AT5G23010",]$ps_left  = gene.boundaries[gene.boundaries$gene_id == "AT5G23010",]$ps_left  - 11e4
gene.boundaries[gene.boundaries$gene_id == "AT5G23010",]$ps_right = gene.boundaries[gene.boundaries$gene_id == "AT5G23010",]$ps_right + 11e4

gene.boundaries[gene.boundaries$gene_id == "AT5G23020",]$ps_left  = gene.boundaries[gene.boundaries$gene_id == "AT5G23020",]$ps_left  - 11e4
gene.boundaries[gene.boundaries$gene_id == "AT5G23020",]$ps_right = gene.boundaries[gene.boundaries$gene_id == "AT5G23020",]$ps_right + 11e4

# updated boundaries
subset(gene.boundaries, gene_id %in% major.loci)

    #       chr  ps_left ps_right   gene_id
    # 1       4  1231980  1473854 AT4G03060
    # 9430    2 10749916 10911655 AT2G25450
    # 17997   4  1224182  1466501 AT4G03050
    # 24201   5  7593092  7816896 AT5G23010
    # 24202   5  7608121  7831839 AT5G23020

# double check window sizes
subset(gene.boundaries, gene_id %in% major.loci)$ps_right - subset(gene.boundaries, gene_id %in% major.loci)$ps_left
    # [1] 241874 161739 242319 223804 223718

### adjust other loci

# initial boundaries
head(gene.boundaries[gene.boundaries$gene_id %in% other.genes,],4)

    #      chr ps_left ps_right   gene_id
    # 1256   1 4118594  4120871 AT1G12130
    # 1257   1 4121369  4124595 AT1G12140
    # 1259   1 4126066  4128310 AT1G12160
    # 1727   1 5605159  5607473 AT1G16400
    # ...

# adjust boundaries
gene.boundaries[gene.boundaries$gene_id %in% other.genes,]$ps_left  = gene.boundaries[gene.boundaries$gene_id %in% other.genes,]$ps_left  - 3e4
gene.boundaries[gene.boundaries$gene_id %in% other.genes,]$ps_right = gene.boundaries[gene.boundaries$gene_id %in% other.genes,]$ps_right + 3e4

# updated boundaries
head(gene.boundaries[gene.boundaries$gene_id %in% other.genes,],4)

    #      chr ps_left ps_right   gene_id
    # 1256   1 4088594  4150871 AT1G12130
    # 1257   1 4091369  4154595 AT1G12140
    # 1259   1 4096066  4158310 AT1G12160
    # 1727   1 5575159  5637473 AT1G16400

# double check window sizes
summary( gene.boundaries[gene.boundaries$gene_id %in% other.genes,]$ps_right - gene.boundaries[gene.boundaries$gene_id %in% other.genes,]$ps_left )
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 60959   61530   62176   62126   62453   63938 

### save file with updated boundaries
write.table(gene.boundaries, "TAIR10_GFF3_geneBoundariesTableADG_extendedWindows.txt",
            sep = "\t", row.names = F, quote = F)






