# This script generates mass images for a selected m/z range and m/z tolerance
library("msir")

# Load the imzML data
filename <- file.path("..", "data", "ltpmsi-chilli.imzML")

imagespectra <- importImzMl(filename, centroided=TRUE)

s <- slides(imagespectra, range=c(40, 540), step=0.2, tolerance=0.4)
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
