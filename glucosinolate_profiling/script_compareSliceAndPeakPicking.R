
### merge output tables from peak picking runs

gsl.filenames = list.files("/glucFiles_mzML/integratedOutput/", pattern="*.csv", full.names = T)
gsl.files     = lapply(gsl.filenames, read.csv)

dat = do.call(rbind, gsl.files)
dat[ is.na(dat$intb), ]$intb = 0

### compare peak picking and slicing methods for 3-Butenyl GSL

pdf("/Figures/Buen_plot.pdf", h = 7.5, w = 2.2)
par(mfrow = c(3,1))
dat.buen = subset(dat, peakname == "gsl.Buen")
plot(log10((dat.buen$slice_intb+2)/2) ~ log10(dat.buen$intb+1), xlim = c(0,5), ylim = c(0,5), col = scales::alpha("dodgerblue", 0.3), pch = 16)
abline(0,1)

hist(log10((dat.buen$slice_intb+2)/2), xlim = c(0,5), breaks = 10)
hist(log10(dat.buen$intb+1), xlim = c(0,5), breaks = 20)
dev.off()

### compare peak picking and slicing methods for 4-Pentenyl GSL

pdf("/Figures/Peen_plot.pdf", h = 7.5, w = 2.2)
par(mfrow = c(3,1))
dat.peen = subset(dat, peakname == "gsl.Peen")
plot(log10((dat.peen$slice_intb+2)/2) ~ log10(dat.peen$intb+1), xlim = c(0,5), ylim = c(0,5), col = scales::alpha("dodgerblue", 0.3), pch = 16)
abline(0,1)

hist(log10((dat.peen$slice_intb+2)/2), xlim = c(0,5), breaks = 10)
hist(log10(dat.peen$intb+1), xlim = c(0,5), breaks = 20)
dev.off()

### compare peak picking and slicing methods for 4-hydroxyindol-3-ylmethyl GSL

pdf("/Figures/4hIM_plot.pdf", h = 7.5, w = 2.2)
par(mfrow = c(3,1))
dat.him = subset(dat, peakname == "gsl.1hIM")
plot(log10((dat.him$slice_intb+2)/2) ~ log10(dat.him$intb+1), xlim = c(0,5), ylim = c(0,5), col = scales::alpha("dodgerblue", 0.3), pch = 16)
abline(0,1)

hist(log10((dat.him$slice_intb+2)/2), xlim = c(0,5), breaks = 10)
hist(log10(dat.him$intb+1), xlim = c(0,5), breaks = 20)
dev.off()
