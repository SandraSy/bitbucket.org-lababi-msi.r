
library("msir")

masterspecFunction<-function(filename){
imagespectra <- importImzMl(filename, centroided=TRUE)
imagespectra <- topN(imagespectra, n=100)
mzdensity <- density(unlist(lapply(imagespectra, mass)), bw=0.01)
masterspectrum <- masterSpectrum(imagespectra, method="sum")
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
        "spectrum-zoom with peaks (tiff, peaklist.csv)")}


#GUI
require(gWidgets)
options("guiToolkit"="RGtk2")
win <- gwindow("MSI.R GUI", visible=TRUE)
group <- ggroup(horizontal = FALSE, container=win)
obj <- gbutton("Choose .imzML file",container=group, handler = function(h,...) {filename<<-gfile()})
obj <- gbutton("Create master spectrum",container=group, handler = function(h,...) masterspecFunction(filename))
               