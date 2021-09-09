# This R script was used to quantify peak areas, using functions in the file
# "functions_metaboliteQuantification.R"

################################################################################
##### 1. convert .mzML files so that MS2 data is encoded as MS1 data       #####
################################################################################

library("xcms"); library("mzR"); library("tools"); library("progress")

filepath = "/glucFiles_mzML/Originals/"
files    = list.files(filepath, pattern = "*.mzML", recursive = TRUE, full.names = TRUE)
outdir   = "/glucFiles_mzML/Converted4xcms/"

fileConverter_qqqMS2toMS1(files, outdir)

### sometimes returns an error (related to stochastic I/O problems on MacOS): 
###   Error in object@backend$getPeakList(x) : 
###   [MSData::Spectrum::getMZIntensityPairs()] Sizes do not match. 
### when this occurs: restart on the next file after the last successful one
#fileConverter_qqqMS2toMS1(files[ (length(list.files(outdir))+1):length(files) ], outdir)

### note #1: Encountered errors on the "blank" files run to flush the LCMS after each plate:
###   Error in `$<-.data.frame`(`*tmp*`, "acquisitionNum", value = 1:0) : 
###   replacement has 2 rows, data has 0
### Removed these files because we don't need to analyze them...

### note #2: Agilent .d files were converted to .mzML files with ProteoWizard 
### prior to running this script, since it requires .mzML input files.

################################################################################
##### 2. explore EICs for each Arabidopsis GSL m/z values of interest      #####
################################################################################

library("xcms"); library("magrittr")

## Inspect EICs to determine which molecules are present, and their approximate retention
## time ranges.
## Focus on molecules of interest in chromatogram sets 2 and 3 only, since set 1 had
## inconsistent retention times (exluded from the published analysis).

## example:
filepath = "/glucFiles_mzML/Converted4xcms/sets23_only/"
files <- list.files(filepath, pattern = "*_MS1converted.mzML", recursive = TRUE, full.names = TRUE)
MSnExpObj = readMSData(files[seq(1,length(files),30)], msLevel. = 1, mode = "onDisk") # every 30th file
  MSnExpObj = readMSData(files[1170], msLevel. = 1, mode = "onDisk") # individual file (#1170)
eicPlotter(MSnExpObj, mzR = c(476.6,477.6), rtR = c(800,1200))
abline(v=c(920,970),col="red")

################################################################################
##### 3. test integration parameters for each Arabidopsis GSL m/z value    #####
################################################################################

library("xcms"); library("progress")

## At this step, integrated peak areas for reach molecule were inspected under different parameters
## to determine optimal values for the full pipeline to quantify all peaks in the next step.

## load data
filepath = "/glucFiles_mzML/Converted4xcms/sets23_only/"
files <- list.files(filepath, pattern = "*_MS1converted.mzML", recursive = TRUE, full.names = TRUE)
#MSnExpObj2 = readMSData(files, msLevel. = 1, mode = "onDisk") # all files
MSnExpObj2 = readMSData(files[seq(1,1224,53)], msLevel. = 1, mode = "onDisk") # subset of files

## inspect EIC plots for many samples combined, example:
eicPlotter(MSnExpObj2, mzR = c(437.5,438.5), rtR = c(0,1200))

## inspect individual peak integration results, example:
mzr = c(437.5,438.5)
rtr = c(500,600)
pdf("/Users/andy/Desktop/temp_gluc_plot.pdf", height = 12, width = 18)
# try various pararmeter settings individually
peakOutput = peakIntegrator(MSnExpObj2, mzR = mzr, rtR = rtr, method = "integrate",
                           peakwidth_param = c(5,10), snthresh_param = 4, prefilter_param = c(1,20))
# peakOutput = peakIntegrator(MSnExpObj2, mzR = mzr, rtR = rtr, method = "integrate",
#                             peakwidth_param = c(5,14), snthresh_param = 4, prefilter_param = c(1,20))
# peakOutput = peakIntegrator(MSnExpObj2, mzR = mzr, rtR = rtr, method = "integrate",
#                            peakwidth_param = c(7,20), snthresh_param = 4, prefilter_param = c(1,20))
# peakOutput = peakIntegrator(MSnExpObj2, mzR = mzr, rtR = rtr, method = "integrate",
#                             peakwidth_param = c(8,16), snthresh_param = 4, prefilter_param = c(1,20))
# peakOutput = peakIntegrator(MSnExpObj2, mzR = mzr, rtR = rtr, method = "integrate",
#                             peakwidth_param = c(8,20), snthresh_param = 4, prefilter_param = c(1,20))
dev.off()

################################################################################
##### 4. integrate peaks                                                   #####
################################################################################

library("xcms")

filepath = "/glucFiles_mzML/Converted4xcms/sets23_only/"
files     = list.files(filepath, pattern = "*_MS1converted.mzML", recursive = TRUE, full.names = TRUE)
peak.int.input = read.csv("~/glucMetadata/GSL_List_corrected_ADG2019Aug09.csv", h=T) # molecule names, RT and m/z ranges

# for full run
MSnExpObj = readMSData(files, msLevel. = 1, mode = "onDisk")

# for subsets of files,
# change first line below to different chunks of 100s files to limit memory usage, 
# leave everything else the same, e.g.
#read.range = c(1101:length(files)) # <--------------------- CHANGE ME!!!!!!!!!!!!!!
#MSnExpObj2 = readMSData(files[read.range], msLevel. = 1, mode = "onDisk")

for (i in c(1:nrow(peak.int.input))){

  out.plot.name = paste("/glucFiles_mzML/integratedOutput/",peak.int.input[i,"name"],"_samples_",min(read.range),"_",max(read.range),"_plot.pdf", sep="")
  pdf(out.plot.name, height = 12, width = 18)

  peakOutput = peakIntegrator(MSnExpObj2, method = c("integrate","slice","summit"),
                              peak_name = peak.int.input[i,"name"],
                              mzR = c(peak.int.input[i,"mz.min"],peak.int.input[i,"mz.max"]),
                              rtR = c(peak.int.input[i,"rt.min"],peak.int.input[i,"rt.max"]), 
                              peakwidth_param = c(peak.int.input[i,"peakwidth.1"],peak.int.input[i,"peakwidth.2"]), 
                              snthresh_param = 4, prefilter_param = c(1,20), 
                              sliceProp = 0.33, summitPoints = 3)
  
  dev.off()
  
  out.table.name = paste("/glucFiles_mzML/integratedOutput/",peak.int.input[i,"name"],"_samples_",min(read.range),"_",max(read.range),"_table.csv", sep="")
  write.csv(peakOutput, out.table.name, quote=F)
  
}


