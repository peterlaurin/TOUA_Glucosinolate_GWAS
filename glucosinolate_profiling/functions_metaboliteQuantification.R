# Overview: These scripts provide wrappers for using functions in xcms and related packages for
#           a pipeline to quantify parent molecule abundances from the output of a QQQ HPLC-MS/MS
#           run in scan mode for a single product molecule. Specifically, this pipeline was used
#           to quantify putative glucosinolates (parent molecules, yielding a sulfate moiety) as
#           described in:
#
#           Gloss et al. (submitted), "Genome-wide association mapping within a single Arabidopsis 
#           thaliana population reveals a richer genetic architecture for defensive metabolite
#           diversity.
#
# Author:   Andrew Gloss, New York University, ag8612 [at] nyu.edu
#
# Note #1:  These functions were created in August 2019 on a macbook with the following R environment:
#           - R       version 3.5.1 (2018-07-02)
#           - MSnbase version 2.8.3
#           - mzR     version 2.16.2
#           - xcms    version 3.4.4 
#
# Note #2:  Anecdotally, there appears to be some issues with running these scripts on a different
#           system; this may be a mac vs. windows issue, or an issue with the version of R or
#           or the packages being used. The source of these issues have not yet been diagnosed.

##########################################################
#####     function #1: fileConverter_qqqMS2toMS1     #####
##########################################################

# purpose:  Takes .mzML files from QQQ HPLC-MS/MS run in scan mode, and creates modified files in which
#           MS2 scan data is coded as MS1 data to enable easy functionality with xcms and similar packages.
# 
# requires: library("mzR")
#           library("tools") # for file_path_sans_ext
#           library("progress") # for progress bar
#
# usage:    fileConverter_qqqMS2toMS1(filenames_vector = files, output_directory_path = outdir)
#           where ** files is a vector of file names (full paths)
#                 ** outdir is the path to the PRE-EXISTING directory where converted files will be written.
#
# example:  filepath = "/your/directory/of/mzMLs/andSubFolders/"
#           files    = list.files(filepath, pattern = "*.mzML", recursive = TRUE, full.names = TRUE)
#           outdir   = "/your/directory/of/mzMLs/convertedOutput/"
#           fileConverter_qqqMS2toMS1(files, outdir)

fileConverter_qqqMS2toMS1 = function(filenames_vector, output_directory_path) {

  # create progress bar
  progress.bar = progress_bar$new(total = length(filenames_vector))

  for (i in filenames_vector){
    # update progress bar
    progress.bar$tick()
    # open file
    dat.temp <- openMSfile(i)
    # get file header
    hdr <- header(dat.temp)
    # get spectra data
    pks <- spectra(dat.temp)
    # remove all scans but ms2 scans
    pks <- pks[hdr$msLevel==2]
    hdr <- hdr[hdr$msLevel==2,]
    # provide new scan/acquisition numbers
    hdr$seqNum<-hdr$acquisitionNum<-seq(nrow(hdr))
    # overwrite msLevel to 'pretend' to be MS1 data
    hdr$msLevel <- 1
    # get file name without extension
    file_baseName_noExt = file_path_sans_ext(basename(i))
    # write out the new 'MS1' data
    writeMSData(object = pks, file = paste(outdir,file_baseName_noExt,"_MS1converted.mzML",sep=""), header=hdr)
  }
}

##########################################################
#####     function #2: eicPlotter                    #####
##########################################################

# purpose:  Creates EIC Plots for one or more chromatograms, filtered for a user-specified range of
#           m/z and retention time values, for a list of files. This is particularly useful for
#           comparing across chromatograms to determine the typical range of observed retention times
#           for a moleclue of interest.
#
# requires: library("xcms")
#           library("magrittr") # for pipe function
#
# usage:    eicPlotter(msObjs = MSnExpObj, mzR = c(num,num))
#           where ** msObjs is an object of class MSnExp read by xcms
#                 ** mzR and rtR are the m/z and retention time ranges for chromatogram filtering
#
# example:  filepath = "/your/directory/of/mzMLs/convertedOutput/"
#           files <- list.files(filepath, pattern = "*_MS1converted.mzML", recursive = TRUE, full.names = TRUE)
#           MSnExpObj = readMSData(files, msLevel. = 1, mode = "onDisk")
#           eicPlotter(MSnExpObj, mzR = c(447,448), rtR = c(0,12000))
#
# note:     Give a nonsensically high value (like 1e6) for the upper range of rtR
#           if you don't want to filter by retention time.

eicPlotter = function(msObjs, mzR, rtR){
  MSnExpObj %>%
    filterRt(rt = rtR) %>%
    filterMz(mz = mzR) %>%
    chromatogram(aggregationFun = "max") %>%
    plot(col = "dodgerblue") 
}

##########################################################
#####     function #3: peakIntegrator                #####
##########################################################

# purpose:  Implements various approaches for integrating peak abundances within a specified m/z
#           and retention time range.
#
# output:   - A table of integrated peak areas for each sample, with some metadata about the
#             approach(es) used to quantify peak areas.
#           - A PDF with a single plot of all samples overlaid, and individual samples with
#             shading to indicate the peak area that was integrated.
#
# requires: library("xcms")
#
# usage:    peakIntegrator(msObjs = MSnExpObj, method = c("integrate","slice","summit"), 
#                          peak_name = "whatever", mzR = c(num,num), rtR = c(num,num),
#                          peakwidth_param = c(num,num), snthresh_param = num, prefilter_param = c(num,num), 
#                          summitPoints = num, sliceProp = num)
#           where ** msObjs is an object of class MSnExp read by xcms
#                 ** method is a vector of one, two, or three ways to analyze peaks
#                    - "integrate" uses xcms::findChromPeaks()
#                    - "slice" sums the distances from the baseline intensity to actual intensity for every
#                       point in the upper nth (user specified proportion) percentile of highest intensity
#                       readings in the chromatogram
#                    - "summit" sums the distances from the baseline intensity to actual intensity for the
#                       n points (specified with summitPoints = num) with the highest intensity in the chromatogram
#                 ** peak_name is the name the peak will be called in the output table (not required)
#                 ** mzR and rtR are the m/z and retention time ranges for chromatogram filtering, applied
#                    prior to peak integration
#                 ** peakwidthparam, snthresparam, prefilterparam are settings for "peakwidth",
#                    "snthresparam", and "prefilterparam" for the xcms CentWaveParam() function;
#                    [only used if the "integrate" method is chosen, in which case it's required]
#                 ** summitPoints specifies the n highest intensity points to use
#                    [only used if the "summit" method is chosen, in which case it's required]
#                 ** sliceProp specifies the proportion (upper percentile) of points to use,
#                    e.g., sliceProp = 0.25 only sums slices created using the 25% of points
#                    with the highest intensity
#                    [only used if the "slice" method is chosen, in which case it's required]
#
# example:  filepath = "/your/directory/of/mzMLs/convertedOutput/"
#           files <- list.files(filepath, pattern = "*_MS1converted.mzML", recursive = TRUE, full.names = TRUE)
#           MSnExpObj = readMSData(files, msLevel. = 1, mode = "onDisk")
#           peakOutput = peakIntegrator(MSnExpObj, method = "integrate", mzR = c(447,448), rtR = c(750,900), 
#                                       peakwidth_param = c(5,10), snthresh_param = 4, prefilter_param = c(1,20))
#
# note:     findChromPeaks should return an object of class "XChromatogram". However, it is simply
#           a vector of numbers in the environment used to develop this script, which suggests that 
#           xcms may not be giving the full output. Thus, this function may need some slight modifications
#           to work in other environments. Specifically, if findChromPeaks correctly returns an XChromatogram object, 
#           what is called "peaks.focal" in the script needs to be called "peaks.focal@chromPeaks" instead. This can 
#           be done by changing the following line of code:
#                  peaks.focal = as.data.frame( findChromPeaks(data_chr[[i]], param = cwp) )
#           ... to ...
#                  peaks.focal = as.data.frame( findChromPeaks(data_chr[[i]], param = cwp)@chromPeaks )
#           However, I'm not sure what, if anything else, may also need to be changed to accomodate the full
#           XChromatogram object. Please contact the author of this script (Andrew Gloss) if you desire to run
#           it in your local computing environment, are encountering this issue, and would like to troubleshoot.

peakIntegrator = function(msObjs, method, peak_name, mzR, rtR, peakwidth_param, snthresh_param, prefilter_param, summitPoints, sliceProp) {
  
  print(paste( "Progress: began filtering chromatograms... (step 1 of 3)", date(), sep =".   "))
  
  data_chr = chromatogram(msObjs, rt = rtR, mz = mzR, msLevel = 1, aggregationFun = "max")

  print(paste( "Progress: began plotting combined chromatogram... (step 2 of 3)", date(), sep =".   "))

  par(mfrow = c(5,5))
  colors1 = rainbow(nrow(msObjs@phenoData))
  colors2 = c("dodgerblue","red","gold","darkorchid1","tan1","green3",rep("gray50",100))
  plot(data_chr, col = colors1[1:length(data_chr)])
  
  print(paste( "Progress: began quantifying peaks... (step 3 of 3)", date(), sep =".   "))
  
  colnames_vector = c("sample","peakname","mzr","rtr","numPts","baseline")
  if ("integrate" %in% method) {colnames_vector = c(colnames_vector,"rt","rtmin","rtmax","into","intb","maxo","sn",
                                                    "numPeaks","peakwidth","snthresh","prefilter")}
  if ("slice"     %in% method) {colnames_vector = c(colnames_vector,"slice_intb","slice_prop")}
  if ("summit"    %in% method) {colnames_vector = c(colnames_vector,"summit_intb","summit_numPts")}
  
  peaks = as.data.frame(matrix(NA, ncol = length(colnames_vector), nrow = nrow(msObjs@phenoData)))
  colnames(peaks) = colnames_vector
  peaks$mzr = paste(mzR[1],mzR[2],sep="_")
  peaks$rtr = paste(rtR[1],rtR[2],sep="_")
  if (!(missing(peak_name))) { peaks$peakname = peak_name }

  if ("integrate" %in% method){
    peaks$numPeaks  = 0
    peaks$peakwidth = paste(peakwidth_param[1],peakwidth_param[2],sep="_")
    peaks$snthresh  = snthresh_param
    peaks$prefilter = paste(prefilter_param[1],prefilter_param[2],sep="_")
  }
  if ("slice" %in% method) { peaks$slice_prop = sliceProp }
  if ("summit" %in% method) {peaks$summit_numPts = summitPoints}
  
  # for each file
  for (i in c(1:length(data_chr))){
    
    peaks[i,"sample"]   = rownames(data_chr@phenoData)[i]
    
    data_chr_df = as.data.frame(data_chr[[i]])
    
    plot(data_chr[[i]], main = paste(i,substring( rownames(data_chr@phenoData)[i], 1,25),sep=". "))
    
    # can also be "data_chr@.Data[[i]]@intensity" ...???
    peaks[i,"baseline"]   = median(data_chr_df$intensity, na.rm = T)
    peaks[i,"numPts"] = nrow(data_chr_df)
    
    # "slice"
    if ("slice" %in% method) {
      
      lastSliceNum          = as.integer( length(data_chr_df$intensity) * sliceProp )
      dat_chr_df_kept       = data_chr_df[order(data_chr_df$intensity,decreasing=T),][1:lastSliceNum,]
      peaks[i,"slice_intb"] = sum(dat_chr_df_kept$intensity - peaks[i,"baseline"])
      
      segments(x0 = dat_chr_df_kept$rtime, y0 = peaks[i,"baseline"], y1 = dat_chr_df_kept$intensity, col = "gold", lwd = 0.2)
      
    }

    # "summit"
    if ("summit" %in% method){
      
      dat_chr_df_kept = data_chr_df[order(data_chr_df$intensity,decreasing=T),][1:summitPoints,]
      peaks[i,"summit_intb"] = sum(dat_chr_df_kept$intensity - peaks[i,"baseline"])
      
      segments(x0 = dat_chr_df_kept$rtime, y0 = peaks[i,"baseline"], y1 = dat_chr_df_kept$intensity, col = "black", lwd = 0.3)
      
    }
    
    # integrate
    if ("integrate" %in% method){
      
      cwp = CentWaveParam()
      cwp = CentWaveParam(peakwidth = peakwidth_param, snthresh = snthresh_param, prefilter = prefilter_param)
      
      peaks.focal.catch = try( findChromPeaks(data_chr[[i]], param = cwp), silent = T )
      
      if (class(peaks.focal.catch) == "try-error") {
        peaks[i,"numPeaks"] = "integration_error"
        next
      }
      
      if (class(peaks.focal.catch) != "try-error") {
        peaks.focal = as.data.frame(peaks.focal.catch)
      }
      
      peaks[i,"numPeaks"] = nrow( peaks.focal )
      
      if (peaks[i,"numPeaks"] > 0) {
        
        peaks[i,"rt"]       = mean( min(peaks.focal$rtmin), max(peaks.focal$rtmax) )
        peaks[i,"rtmin"]    = min( peaks.focal$rtmin )
        peaks[i,"rtmax"]    = max( peaks.focal$rtmax )
        peaks[i,"into"]     = sum( peaks.focal$into )
        peaks[i,"intb"]     = sum( peaks.focal$intb )
        peaks[i,"maxo"]     = max( peaks.focal$maxo )
        peaks[i,"sn"]       = mean( peaks.focal$sn )
        
        # for each peak
        for (j in c(1:nrow(peaks.focal))) {
          p.polygon = NULL; p.polygon.kept = NULL
          p.polygon = data.frame(x = data_chr_df$rtime, y = data_chr_df$intensity)
          p.polygon.kept = subset(p.polygon, x >= peaks.focal[j,"rtmin"] & x <= peaks.focal[j,"rtmax"])
          polygon(x = c(min(p.polygon.kept$x), p.polygon.kept$x, max(p.polygon.kept$x)),
                  y = c(0, p.polygon.kept$y, 0), col = adjustcolor(colors2[j],alpha.f=0.5)) # old: col = "#1b98e0"
          
        }
        
      }
      
    }
    
  }
  
  return(peaks)
  
}

