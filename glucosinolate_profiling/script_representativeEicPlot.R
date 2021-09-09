library("xcms"); library("magrittr")

eicPlotter = function(msObjs, mzR, rtR, plot.title){
  MSnExpObj %>%
    filterRt(rt = rtR) %>%
    filterMz(mz = mzR) %>%
    chromatogram(aggregationFun = "max") %>%
    plot(col = "dodgerblue", main = plot.title) 
}

peak.int.input = read.csv("~/glucMetadata/GSL_List_corrected_ADG2019Aug09.csv", h=T) # molecule names, RT and m/z ranges

## representative samples
filepath = "/glucFiles_mzML/Converted4xcms/sets23_only/"
files     = list.files(filepath, pattern = "*_MS1converted.mzML", recursive = TRUE, full.names = TRUE)
MSnExpObj = readMSData(files[1150], msLevel. = 1, mode = "onDisk") # pick one sample, #1150 in this instance

mols = c("gsl.R2hBuen", "gsl.Pren", "gsl.S2hBuen","gsl.4mSOb","gsl.5mSOp","gsl.2hPeen","gsl.Buen",   
         "gsl.6mSOh","gsl.1hIM","gsl.7mSOh","gsl.Peen","gsl.8mSOo","gsl.IM","gsl.4moIM","gsl.1moIM",
         "gsl.7mSh","gsl.8MTO","gsl.3mSOp")

# all peaks together (full m/z and retention time range for the entire QQQ run)
pdf("/Figures/eic_rep_sample.pdf", h = 4, w = 6.5)
eicPlotter(MSnExpObj, mzR = c(300,520), rtR = c(200,1150), plot.title = "representative sample")
dev.off()

# individual molecule peaks separately (plotting the full retention time range for each m/z range)
pdf("/Figures/eic_rep_sample_indiv.pdf", h = 30, w = 6.5)
par(mfrow = c(9,1))
for(i in mols[1:9]){
  focal.mol = peak.int.input[peak.int.input$name == i,]
  mzR1 = focal.mol$mz.min
  mzR2 = focal.mol$mz.max
  rt1  = focal.mol$rt.min
  rt2  = focal.mol$rt.max
  eicPlotter(MSnExpObj, rtR = c(200,1200), mzR = c(mzR1, mzR2), plot.title = paste0(focal.mol$name,"--",rt1,"-",rt2))
}
dev.off()

pdf("/Figures/eic_rep_sample_indiv2.pdf", h = 30, w = 6.5)
par(mfrow = c(9,1))
for(i in mols[10:18]){
  focal.mol = peak.int.input[peak.int.input$name == i,]
  mzR1 = focal.mol$mz.min
  mzR2 = focal.mol$mz.max
  rt1  = focal.mol$rt.min
  rt2  = focal.mol$rt.max
  eicPlotter(MSnExpObj, rtR = c(200,1200), mzR = c(mzR1, mzR2), plot.title = paste0(focal.mol$name,"--",rt1,"-",rt2))
}
dev.off()
