# This script summarizes the x (e.g. 100) most intense peaks/point of each
# spectrum in one single spectrum
# The increasing number indicates the number of processed spectra
# The final spectrum is exported into a .mzML file and printed to a .tiff
# A density distribution map is printed to a publication quality .tiff
# Further, peaks are picked and exported to a .csv text file and added to a zoom spectrum
# Please adjust the paramenters for peak picking (e.g. signal/noise) and spectrum zoom below

# You can change the number of selected peaks with (this would be the setting for 100 signals):
#if (length(newspec)>"200"){        //2x the number of peaks
#  newspec <- newspec[1:100,]}      //1:x strongest signals

# Clear Workspace
rm(list=ls())

# install MALDIquant
#install.packages(c("MALDIquantForeign", "MALDIquant"))
library("MALDIquant")
library("MALDIquantForeign")

# Load the imzML data
filename <- file.path("..", "data", "ltpms-chilli.imzML")

# Define a new function to find and return the highest n peaks in each spectrum.
topN <- function(peaklist, n=100) {
  lapply(peaklist, function(x) {
    x[sort.int(intensity(x), decreasing=TRUE, index.return=TRUE)$ix[1:n]]
  })
}

imagespectra <- importImzMl(filename, centroided=TRUE)

imagespectra <- topN(imagespectra, n=100)

mzdensity <- density(sapply(imagespectra, mass), bw=0.01)

## as needs MALDIquant 1.11.4 from github: install_github("sgibb/MALDIquant")
## otherwise you could use:
## peaks <- mergeMassPeaks(imagespectra.sg, method="sum")
## masterspectrum <- createMassSpectrum(mass(peaks), intensity(peaks))
masterspectrum <- as(mergeMassPeaks(imagespectra, method="sum"),
                     "MassSpectrum")

# peak picking and export to a peaklist.csv file
peaks <- detectPeaks(masterspectrum, SNR=3)


# printing spectrum. format options could be e.g. png, eps, tiff, pdf.
# You can control the resolution with res
tiff(filename="spectrumdensity.tiff", res=1200, compression="lzw",
     height=200, width=200, units="mm")
plot(mzdensity, xlab="m/z", ylab="density", main="m/z density distribution")
polygon(mzdensity, col="grey")
dev.off()

# print master spectrum
tiff(filename="masterspectrum.tiff", res=1200, compression="lzw",
     height=200, width=200,units="mm")
plot(masterspectrum, xlab="m/z")
dev.off()

tiff(filename="spectrum-with-peaks.tiff", res=1200, compression="lzw",
     height=200, width=200, units="mm")
plot(masterspectrum, xlim=c(40, 550), xlab="m/z")
points(peaks, pch=18)
dev.off()


# export a single mzML spectrum
exportMzMl(masterspectrum, force=TRUE, file="masterspectrum.mzML")

# export of peaks
exportCsv(peaks, file="peaklist.csv")

message("Files generated: Pseudo master spectrum (tiff and mzML), ",
        "spectrum density (tiff), ",
        "spectrum-zoom with peaks (tiff, peaklist.csv)")
