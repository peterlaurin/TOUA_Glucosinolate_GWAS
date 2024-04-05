# Purpose: calculate genomic inflation for each GWAS based on median P-values, and also extract observed P-values
#          at different percentile rankings of the p-value distribution (which can be used to identify genomic
#          inflation that is not evident from median-based approaches)

# calculation of lambda from p-value distribution is described here:
#     described in https://bioinformaticsngs.wordpress.com/2016/03/08/genomic-inflation-factor-calculation/
#              and http://genometoolbox.blogspot.com/2014/08/how-to-calculate-genomic-inflation.html
#              and https://stats.stackexchange.com/questions/110755/calculate-inflation-observed-and-expected-p-values-from-uniform-distribution-in


# note: before downstream analyses with this and other scripts, GWAS output files were sorted into folders
#       based on trait type (BLUP or ratio), glucosinolate type (aliphatic or indole), and dataset (TOU-A, Brachi, Katz, Wu)
my.directories = c("/output_gwas_assocs_lmm/toua_N192/blup_aliphatic_assoc/",
                "/output_gwas_assocs_lmm/toua_N192/blup_indolic_assoc/",
                "/output_gwas_assocs_lmm/toua_N192/ratio_aliphatic_assoc/",
                "/output_gwas_assocs_lmm/toua_N192/ratio_indolic_assoc/",
                "/output_gwas_assocs_lmm/brachi_full/blup_aliphatic_assoc/",
                "/output_gwas_assocs_lmm/brachi_full/ratio_aliphatic_assoc/",
                "/output_gwas_assocs_lmm/brachi_N192/blup_aliphatic_assoc/",
                "/output_gwas_assocs_lmm/brachi_N192/ratio_aliphatic_assoc/",
                "/output_gwas_assocs_lmm/katz_full/blup_aliphatic_assoc/",
                "/output_gwas_assocs_lmm/katz_full/ratio_aliphatic_assoc/",
                "/output_gwas_assocs_lmm/katz_N192/blup_aliphatic_assoc/",
                "/output_gwas_assocs_lmm/katz_N192/ratio_aliphatic_assoc/",
                "/_katz_full_moved/blup_indolic_assoc/",
                "/_katz_N192_moved/blup_indolic_assoc/",
                "/output_gwas_assocs_lmm/wu_full/blup_aliphatic_assoc/",
                "/output_gwas_assocs_lmm/wu_full/blup_indolic_assoc/",
                "/output_gwas_assocs_lmm/wu_full/ratio_aliphatic_assoc/",
                "/output_gwas_assocs_lmm/wu_full/ratio_indolic_assoc/",
                "/output_gwas_assocs_lmm/wu_N192/blup_aliphatic_assoc/",
                "/output_gwas_assocs_lmm/wu_N192/blup_indolic_assoc/",
                "/output_gwas_assocs_lmm/wu_N192/ratio_aliphatic_assoc/",
                "/output_gwas_assocs_lmm/wu_N192/ratio_indolic_assoc/"
                )

# determine full list of files within specified directories of GWAS output
my.filenames = list.files(path = my.directories, pattern = "assoc.fixed.txt", full.names = T)

# create dataframe to store lambda and percentile p-values 
# (e.g., p50 is the observed 50th percentile p-value, p005 is the 0.5 percentile p-value, etc.)
gwas.summary = data.frame(matrix(nrow = length(my.filenames), ncol = 32))
colnames(gwas.summary) = c("gwas.filename","lambda","p0","p005","p010","p015","p020","p025","p030","p035","p040","p045","p050",
                           "p10","p15","p20","p25","p30","p35","p40","p45","p50","p55","p60","p65","p70","p75","p80","p85","p90","p95","p100")
gwas.summary$gwas.filename = my.filenames

# Parse each file, filling out dataframe defined above
for (i in 1:length(my.filenames)){
    
  print( paste0("Processing file ", i, " of ", length(my.filenames), "...") )
  
  gwas.i = read.delim(my.filenames[i])
  
  gwas.i = subset(gwas.i, af > 0.05)
  
  # quantile p-values
  gwas.summary[i,3:13] = quantile(gwas.i$p_wald, probs = seq(0, 0.05, 0.005))
  gwas.summary[i,14:32] = quantile(gwas.i$p_wald, probs = seq(0.1, 1, 0.05))
  
  # genomic inflation factor (labmda)
  gwas.i$chisq = qchisq(gwas.i$p_wald, 1, lower.tail=FALSE)
  lambda = median(gwas.i$chisq) / qchisq(0.5, 1)
  gwas.summary[i,"lambda"] = lambda

}

# write output
write.table(gwas.summary, "/pvalue_inflation_summary.txt", sep = "\t", row.names = F, quote = F)








