# This script generates mass images for a selected m/z range and m/z tolerance

# Clear Workspace
rm(list=ls())

# install MALDIquant
#install.packages(c("MALDIquantForeign", "MALDIquant"))
library("MALDIquant")
library("MALDIquantForeign")

# library for parallel computing and creation of cluster. 
library(iterators)
library(doParallel)
# create R cluster. Please register the number of CPUs you want to use.
cl<-makeCluster(2)
registerDoParallel(cl)

# Library for 2D maps
#install.packages("latticeExtra")
library(latticeExtra)

# Load the imzML data
filename <- "ftms_example.imzML"

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

#startmz<-40 #in m/z, start of scanning
#centermz<-startmz #in m/z, start of scanning
#scanmaxmz<-540 #in m/z, end of scanning
#scansteps<-0.2 #in m/z, should be smaller than scantolerance
scantolerance<-0.001 #in m/z
#scanrangevec<-seq(startmz, scanmaxmz, by = scansteps)

scanrangevec<-c(499.0274047852,
  500.0307312012,
  501.042755127,
  514.9743041992,
  534.2953491211,
  542.9702758789,
  544.3393554688,
  551.0149536133,
  558.9431152344,
  560.3111572266,
  566.3214111328,
  566.9888305664,
  567.9919433594,
  570.3551635742,
  571.3583374023,
  577.5189819336,
  580.9262084961,
  581.9546508789,
  582.2954101563,
  582.9628295898,
  583.2998046875,
  583.966003418,
  584.9602661133,
  591.0534057617,
  592.3377075195,
  596.9206542969,
  598.9364624023,
  599.5026855469,
  599.940246582,
  600.9345092773,
  603.9372558594,
  604.9451293945,
  608.3110351563,
  609.314453125,
  611.5391845703,
  614.9105224609,
  619.0488891602,
  619.9096679688,
  620.9185180664,
  625.5186767578,
  626.5220336914,
  627.5344238281,
  637.05859375,
  653.0536499023,
  681.4859619141,
  682.4572143555,
  703.5746459961,
  704.5775756836,
  706.9585571289,
  718.978515625,
  725.5568847656,
  726.560546875,
  731.6062011719,
  732.5540771484,
  734.5690917969,
  734.9529418945,
  735.573425293,
  739.4666137695,
  740.5222167969,
  741.5305786133,
  742.5343017578,
  742.9962158203,
  744.4931640625,
  753.5880737305,
  756.5512084961,
  756.9340209961,
  757.5554199219,
  758.5693969727,
  758.9714355469,
  759.5718383789,
  759.9739990234,
  760.5850219727,
  760.9691162109,
  761.4522094727,
  761.5881958008,
  766.5379638672,
  769.5618896484,
  770.5095214844,
  770.5648193359,
  771.4922485352,
  771.5130615234,
  772.4971923828,
  772.5253295898,
  772.9296875,
  773.5290527344,
  773.9371948242,
  774.9449462891,
  775.9486694336,
  776.9429931641,
  778.4780273438,
  780.5512084961,
  780.9538574219,
  781.5543823242,
  782.5661010742,
  782.5697021484,
  783.5696411133,
  783.5733032227,
  784.5616455078,
  784.576171875,
  784.5850830078,
  786.6008911133,
  787.4676513672,
  787.6040039063,
  788.4694213867,
  788.616027832,
  790.5145263672,
  790.9190673828,
  791.9218139648,
  792.9175415039,
  794.6060791016,
  795.6083984375,
  795.9190063477,
  796.5255126953,
  796.9260253906,
  797.5288696289,
  798.5409545898,
  799.5444335938,
  800.4611206055,
  800.5382080078,
  800.5481567383,
  803.9844360352,
  804.4934082031,
  804.5509033203,
  805.5544433594,
  806.5699462891,
  806.892578125,
  808.5852661133,
  809.588684082,
  810.592956543,
  810.6010742188,
  811.6046142578,
  811.8938598633,
  812.9014282227,
  816.4344482422,
  816.5877685547,
  817.5910644531,
  820.5255737305,
  821.528137207,
  822.5232543945,
  822.5324707031,
  822.541015625,
  823.5435180664,
  824.5565185547,
  825.5596923828,
  826.4548950195,
  826.5726318359,
  827.5745849609,
  828.4689941406,
  828.5512695313,
  830.5463867188,
  830.5668945313,
  831.548828125,
  831.5707397461,
  832.5621337891,
  832.5744628906,
  832.5828857422,
  833.5642089844,
  836.6162719727,
  837.619140625,
  842.4497070313,
  844.5250854492,
  845.5282592773,
  846.5407104492,
  847.5444946289,
  848.5390014648,
  848.5477294922,
  848.5569458008,
  849.5418701172,
  849.5598144531,
  857.4822387695,
  873.4552001953,
  874.5724487305,
  876.5862426758,
  905.7106323242,
  917.6510009766,
  919.7152099609,
  920.7183227539,
  934.9808349609,
  943.6658325195,
  944.6690673828,
  945.7321166992,
  948.9381103516,
  950.9542236328,
  951.9573364258,
  966.9281005859,
  967.9309082031,
  968.9258422852,
  968.9314575195,
  972.936340332,
  982.9026489258,
  993.7306518555,
  1004.885559082,
  1022.5497436523,
  1022.7663574219,
  1104.2752685547,
  1104.3623046875,
  1104.4207763672,
  1104.4510498047,
  1104.5007324219,
  1104.6979980469,
  1104.7890625,
  1104.9027099609,
  1105.1485595703,
  1105.1734619141,
  1105.21484375,
  1105.2414550781,
  1105.361328125,
  1105.4573974609,
  1105.5971679688,
  1105.6979980469,
  1105.7359619141,
  1105.8089599609,
  1105.9031982422,
  1106.0244140625,
  1106.1110839844,
  1106.2154541016,
  1106.7215576172,
  1126.9619140625,
  1142.9378662109,
  1158.9094238281,
  1163.1281738281,
  1182.8903808594,
  1280.0397949219,
  1287.1346435547,
  1287.6241455078,
  1287.6772460938,
  1287.775390625,
  1287.9870605469,
  1288.0289306641,
  1288.1899414063,
  1288.2648925781,
  1288.453125,
  1288.5529785156,
  1288.673828125,
  1288.7591552734,
  1288.8188476563,
  1288.8615722656,
  1288.890625,
  1288.9309082031,
  1288.9680175781,
  1289.0573730469,
  1289.154296875,
  1289.2155761719,
  1289.2882080078,
  1289.330078125,
  1289.4541015625,
  1289.4995117188,
  1289.7509765625,
  1289.9914550781,
  1290.1176757813,
  1290.2113037109,
  1311.5435791016,
  1318.9449462891,
  1334.9196777344,
  1366.2967529297,
  1464.5776367188,
  1508.6885986328,
  1570.1096191406,
  1638.3363037109,
  1704.41796875,
  1743.6278076172,
  1763.447265625,
  1767.9663085938,
  1770.9979248047,
  1879.3156738281,
  1880.6632080078,
  1994.4298095703,
  2012.0270996094,
  2033.6062011719,
  2125.3942871094,
  2133.7761230469,
  2222.5529785156,
  2231.9675292969,
  2243.880859375,
  2257.3754882813,
  2362.8701171875,
  2385.3837890625,
  2388.1352539063,
  2441.6171875,
  2457.5505371094,
  2468.4624023438,
  2490.5610351563,
  2494.4787597656,
  2510.8857421875,
  2559.3410644531,
  2565.5009765625,
  2584.2314453125,
  2631.779296875,
  2659.3815917969,
  2662.7065429688,
  2677.4252929688,
  2685.3212890625,
  2704.7272949219,
  2720.7841796875,
  2751.3012695313,
  2795.1118164063,
  2823.7194824219,
  2841.5068359375,
  2853.9692382813,
  2881.4467773438,
  2899.8405761719,
  2940.1489257813,
  2951.7175292969,
  2955.2043457031,
  2972.6552734375,
  2997.1909179688)
iteratorvec<-iter(scanrangevec)
  
# For running in sequential mode, change %dopar% to %do%

foreach(centermz=iteratorvec, .packages=c('MALDIquant','latticeExtra','iterators')) %dopar%
{
  # This will be name and title of the .png file
  # To sort the files according to the m/z mass, a prefix is used
  # The actual pngfiletitle is printed on the console for monitoring the progress of the script
  prefix<-""
  if (centermz<1000000){prefix<-"0"}
  if (centermz<100000){prefix<-"00"}
  if (centermz<10000){prefix<-"000"}
  if (centermz<1000){prefix<-"0000"}
  if (centermz<100){prefix<-"00000"}
  if (centermz<10){prefix<-"000000"}
  pngfilename<-paste(prefix,toString(centermz),"_",toString(scantolerance),".png",sep="")
  pngfiletitle<-paste("m/z ",toString(centermz),"+/-",toString(scantolerance),sep="")
  print(pngfiletitle)

specno <-1

while (specno<=elementsimzML)
{
  
  mzs <- mass(imagespectra[[specno]])
  counts <- intensity(imagespectra[[specno]])
  specmatrix <- cbind(mzs, counts)
  
  lowermass<-centermz-scantolerance
  highermass<-centermz+scantolerance
  
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

png(pngfilename)
#if you don't want contour lines, change to contour=FALSE
#rainbow colors
print(levelplot(intmatrix,main=pngfiletitle,xlab="x/ pixel",ylab="y/ pixel",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,col.regions = rainbow(100,start=1/5)))
#terrain (map-like) colors
#print(levelplot(intmatrix,main=pngfiletitle,xlab="x/ pixel",ylab="y/ pixel",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,col.regions = terrain.colors(100)))
#greyscale
#print(levelplot(intmatrix,main=pngfiletitle,xlab="x/ pixel",ylab="y/ pixel",scales = list(draw = FALSE),contour=TRUE,pretty=TRUE,col.regions = grey.colors(100)))
dev.off()

}

# Information about employed CPU workers and closing the cluster
workers <- getDoParWorkers()
print('CPU workers employed: ')
print(workers)
stopCluster(cl)

