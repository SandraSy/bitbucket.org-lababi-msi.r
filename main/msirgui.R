
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
        
featurescanningFunction<-function(filename){
imagespectra <- importImzMl(filename, centroided=TRUE)
s <- slides(imagespectra, range=c(83, 85), step=0.2, tolerance=0.4)
center <- attr(s, "center")
tolerance <- attr(s, "tolerance")
filenames <- sprintf("%08.1f_%.1f.png", center, tolerance)
titles <- sprintf("mz %.1f+/-%.1f", center, tolerance)
col <- rainbow(100L, start=1L/5L)
pb <- txtProgressBar(0L, length(center), style=3)
## loop through intervals
for (i in seq(along=center)) {
  png(filenames[i])
  print(levelplot(s[,,i],
                  main=titles[i], xlab="x/ pixel", ylab="y/ pixel",
                  scales=list(draw=FALSE), contour=TRUE, pretty=TRUE,
                  col.regions=col))
  ## useRaster=TRUE, decreases run time dramatically but the visual
  ## representation changes as well
  #print(levelplot(s[,,i],
  #                main=titles[i], xlab="x/ pixel", ylab="y/ pixel",
  #                scales=list(draw=FALSE), contour=TRUE, pretty=TRUE,
  #                col.regions=col, useRaster=TRUE, interpolate=TRUE))
  setTxtProgressBar(pb, i)
  dev.off()
}
close(pb)
print("Scanning finished.")
}

ionimageFunction<-function(filename){
imagespectra <- importImzMl(filename, centroided=TRUE)
#plotImsSlice(imagespectra)
elementsimzML<-length(imagespectra)
#determine size and position of image
imagepos <- sapply(imagespectra, function(x)metaData(x)$imaging$pos)
minx<-min(imagepos[1,])
maxx<-max(imagepos[1,])
miny<-min(imagepos[2,])
maxy<-max(imagepos[2,])
xrange<-maxx-minx+1
yrange<-maxy-miny+1
intmatrix <- matrix(0,yrange,xrange)
colnames(intmatrix) <- seq(minx,maxx)
rownames(intmatrix) <- seq(miny,maxy)
specno <-1
centermz<-84.1 #in m/z, mass trace of interest
scantolerance<-0.4 #in m/z
lowermass<-centermz-scantolerance
highermass<-centermz+scantolerance
while (specno<=elementsimzML)
{
mzs <- mass(imagespectra[[specno]])
counts <- intensity(imagespectra[[specno]])
specmatrix <- cbind(mzs, counts)
mzrangeindex <- (specmatrix[,1] > lowermass & specmatrix[,1] < highermass)
decisionvalue<- sum(mzrangeindex)
if (decisionvalue=="0"){mzpixintensity<-0}
if (decisionvalue=="1"){mzrangeints<-specmatrix[mzrangeindex, ]
mzpixintensity<-mzrangeints['counts']}
if (decisionvalue>"1"){
mzrangeints<-specmatrix[mzrangeindex, ]
mzpixintensity<-sum(mzrangeints[,2])}
specposition<-metaData(imagespectra[[specno]])$imaging$pos
specpositionx<-paste(specposition[1])
specpositiony<-paste(specposition[2])
intmatrix[specpositiony,specpositionx]<-mzpixintensity
print(specno)
specno=specno+1
}
# printing spectrum. format options could be e.g. png, eps, tiff, pdf. You can control the resolution with res
prefix<-""
if (centermz<1000000){prefix<-"0"}
if (centermz<100000){prefix<-"00"}
if (centermz<10000){prefix<-"000"}
if (centermz<1000){prefix<-"0000"}
if (centermz<100){prefix<-"00000"}
if (centermz<10){prefix<-"000000"}
tifffilename<-paste(prefix,toString(centermz),"_",toString(scantolerance),".tiff",sep="")
tifffiletitle<-paste("m/z ",toString(centermz),"+/-",toString(scantolerance),sep="")
tiff(filename=tifffilename,res=1200,compression="lzw",height=200,width=200,units="mm")
#if you don't want contour lines, change to contour=FALSE
#HCL sequential colors
print(levelplot(intmatrix,main=tifffiletitle,xlab="x/ pixel",ylab="y/ pixel",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,col.regions = heat_hcl(100, h = c(0, 250), c. = c(30, 150), l = c(100, 25), power = c(1/5, 3), gamma = NULL, fixup = TRUE, alpha = 1)))
#terrain (map-like) colors
#print(levelplot(intmatrix,main=tifffiletitle,xlab="x/ pixel",ylab="y/ pixel",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,col.regions = terrain.colors(100)))
#greyscale
#print(levelplot(intmatrix,main=tifffiletitle,xlab="x/ pixel",ylab="y/ pixel",scales = list(draw = TRUE),contour=FALSE,pretty=TRUE,col.regions = grey.colors(100,start = 1, end = 0, gamma = 2.2)))
dev.off()
# search for spectrum with highest intensity for this m/z
maxintind = which(intmatrix == max(intmatrix), arr.ind = TRUE)
#max(intmatrix)
ymaxint=rownames(intmatrix)[maxintind[,1]]
xmaxint=colnames(intmatrix)[maxintind[,2]]
specposition['y']<-as.numeric(ymaxint)
specposition['x']<-as.numeric(xmaxint)
xs<-imagepos[1,]==specposition['x']
ys<-imagepos[2,]==specposition['y']
specno<-which(xs&ys)
#extract spectrum data
mzs <- mass(imagespectra[[specno]])
counts <- intensity(imagespectra[[specno]])
#plot to tiff
# OUTPUT
tiffspecname<-paste(toString(specposition['x']),"x_",toString(specposition['y']),"y_maxInt_",prefix,toString(centermz),"_",toString(scantolerance),".tiff",sep="")
tiffspectitle<-paste("x=",toString(specposition['x']),", y=",toString(specposition['y'])," m/z ",toString(centermz),"+/-",toString(scantolerance),sep="")
mzMLspecname<-paste(toString(specposition['x']),"x_",toString(specposition['y']),"y_maxInt_",prefix,toString(centermz),"_",toString(scantolerance),".mzML",sep="")
tiff(filename=tiffspecname,res=1200,compression="lzw",height=200,width=200,units="mm")
plot(mzs,counts,xlab="m/z",ylab="intensity",main=tiffspectitle,"h")
abline(h=0)
dev.off()
#export .mzML data
s <- list(createMassSpectrum(mass=mzs, intensity=counts))
exportMzMl(s[[1]], force=TRUE, file=mzMLspecname)
print("mz image and spectra saved to files.")
}

threeimagesFunction<-function(filename){
imagespectra <- importImzMl(filename, centroided=TRUE)
#plotImsSlice(imagespectra)
elementsimzML<-length(imagespectra)
#determine size and position of image
imagepos <- sapply(imagespectra, function(x)metaData(x)$imaging$pos)
minx<-min(imagepos[1,])
maxx<-max(imagepos[1,])
miny<-min(imagepos[2,])
maxy<-max(imagepos[2,])
xrange<-maxx-minx+1
yrange<-maxy-miny+1
intmatrix <- matrix(0,yrange,xrange)
colnames(intmatrix) <- seq(minx,maxx)
rownames(intmatrix) <- seq(miny,maxy)
# Image R
specno <-1
centermz<-306.4 #in m/z, mass trace of interest
scantolerance<-0.4 #in m/z
lowermass<-centermz-scantolerance
highermass<-centermz+scantolerance
while (specno<=elementsimzML)
{
  mzs <- mass(imagespectra[[specno]])
  counts <- intensity(imagespectra[[specno]])
  specmatrix <- cbind(mzs, counts)
  mzrangeindex <- (specmatrix[,1] > lowermass & specmatrix[,1] < highermass)
  decisionvalue<- sum(mzrangeindex)
  if (decisionvalue=="0"){mzpixintensity<-0}
  if (decisionvalue=="1"){mzrangeints<-specmatrix[mzrangeindex, ]
                          mzpixintensity<-mzrangeints['counts']}
  if (decisionvalue>"1"){
    mzrangeints<-specmatrix[mzrangeindex, ]
    mzpixintensity<-sum(mzrangeints[,2])}
  specposition<-metaData(imagespectra[[specno]])$imaging$pos
  specpositionx<-paste(specposition[1])
  specpositiony<-paste(specposition[2])
  intmatrix[specpositiony,specpositionx]<-mzpixintensity
  print(specno)
  specno=specno+1
}
R<-intmatrix
minR<-0.1*max(R)
R[ R < minR ] <- 0
col.R <- colorRampPalette(c('white','red'))(100) 
tiff(filename="R.tiff",res=1200,compression="lzw",height=200,width=200,units="mm")
print(levelplot(R,main="",xlab="",ylab="",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,colorkey=FALSE,regions=TRUE,col.regions = col.R))
dev.off()
# Image G
specno <-1
centermz<-84.1 #in m/z, mass trace of interest
scantolerance<-0.4 #in m/z
lowermass<-centermz-scantolerance
highermass<-centermz+scantolerance
while (specno<=elementsimzML)
{
  mzs <- mass(imagespectra[[specno]])
  counts <- intensity(imagespectra[[specno]])
  specmatrix <- cbind(mzs, counts)
  mzrangeindex <- (specmatrix[,1] > lowermass & specmatrix[,1] < highermass)
  decisionvalue<- sum(mzrangeindex)
  if (decisionvalue=="0"){mzpixintensity<-0}
  if (decisionvalue=="1"){mzrangeints<-specmatrix[mzrangeindex, ]
                          mzpixintensity<-mzrangeints['counts']}
  if (decisionvalue>"1"){
    mzrangeints<-specmatrix[mzrangeindex, ]
    mzpixintensity<-sum(mzrangeints[,2])}
  specposition<-metaData(imagespectra[[specno]])$imaging$pos
  specpositionx<-paste(specposition[1])
  specpositiony<-paste(specposition[2])
  intmatrix[specpositiony,specpositionx]<-mzpixintensity
  print(specno)
  specno=specno+1
}
G<-intmatrix
minG<-0.08*max(G)
G[ G < minG ] <- 0
col.G <- colorRampPalette(c('white','green'))(100) 
tiff(filename="G.tiff",res=1200,compression="lzw",height=200,width=200,units="mm")
print(levelplot(G,main="",xlab="",ylab="",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,colorkey=FALSE,regions=TRUE,col.regions = col.G))
dev.off()
# Image B
specno <-1
centermz<-62.1 #in m/z, mass trace of interest
scantolerance<-0.4 #in m/z
lowermass<-centermz-scantolerance
highermass<-centermz+scantolerance
while (specno<=elementsimzML)
{
  mzs <- mass(imagespectra[[specno]])
  counts <- intensity(imagespectra[[specno]])
  specmatrix <- cbind(mzs, counts)
  mzrangeindex <- (specmatrix[,1] > lowermass & specmatrix[,1] < highermass)
  decisionvalue<- sum(mzrangeindex)
  if (decisionvalue=="0"){mzpixintensity<-0}
  if (decisionvalue=="1"){mzrangeints<-specmatrix[mzrangeindex, ]
           mzpixintensity<-mzrangeints['counts']}
  if (decisionvalue>"1"){
    mzrangeints<-specmatrix[mzrangeindex, ]
    mzpixintensity<-sum(mzrangeints[,2])}
  specposition<-metaData(imagespectra[[specno]])$imaging$pos
  specpositionx<-paste(specposition[1])
  specpositiony<-paste(specposition[2])
  intmatrix[specpositiony,specpositionx]<-mzpixintensity
  print(specno)
  specno=specno+1
}
B<-intmatrix
minB<-0.1*max(B)
B[ B < minB ] <- 0
col.B <- colorRampPalette(c('white','blue'))(100) 
tiff(filename="B.tiff",res=1200,compression="lzw",height=200,width=200,units="mm")
print(levelplot(B,main="",xlab="",ylab="",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,colorkey=FALSE,regions=TRUE,col.regions = col.B))
dev.off()
print("Images for RGB generated. You can combine them e.g. with ImageMagick $convert -combine R.tiff G.tiff B.tiff RGB.tiff")
}

#GUI
require(gWidgets)
options("guiToolkit"="RGtk2")
win <- gwindow("MSI.R GUI", visible=TRUE)
group <- ggroup(horizontal = FALSE, container=win)
obj <- gbutton("Choose .imzML file",container=group, handler = function(h,...) {filename<<-gfile()})
obj <- gbutton("Create master spectrum",container=group, handler = function(h,...) masterspecFunction(filename))
obj <- gbutton("Scan and create images",container=group, handler = function(h,...) featurescanningFunction(filename))
obj <- gbutton("Individual ion analysis",container=group, handler = function(h,...) ionimageFunction(filename))
obj <- gbutton("Create RGB ion images",container=group, handler = function(h,...) threeimagesFunction(filename))
               
