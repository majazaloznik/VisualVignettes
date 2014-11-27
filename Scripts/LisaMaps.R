library(spdep)
library(RColorBrewer)
library(classInt)
library(rgdal)
library(maptools)
gpclibPermit()
library(dplyr)
library(vcd)
oldpar <- par(no.readonly=TRUE)
#########################################
## Load data and Shape file
## For this example:
## LA.shape - E&W Local authority SpatialPolygon with LA names as ID
## LA.admin - dataframe with codes and names of administrative areas, just in case
## cars.data - 2001 census, proportion households in LA with two or more cars
## lw.foq - first order queen list of weights
########################################
load("Data/EnglandWales.R")


########################################
## TO CHANGE using own data:
## You need a SpatialPolygons object (with correctly named polygons - area ID)
## a verctor(or dataframe) of data of same length (also named with area ID)
## a listw object created from the SpatialPolygons object
########################################
shape <- shape
list.of.weights <- lw.foq
mydata <- data.frame(variable=cars.data)
mydata <- data.frame(variable=distwork)

rm(cars.data, LA.shape, lw.foq, distwork)
#########################################
## Data frame setup
#########################################
## This calculates
## - lagged variable values
## - z-scores and lagged z-scores
## - local moran stats
## - lagged values
## - quadrant
#########################################
lisa.data <-  mutate(mydata, 
                     lagged = lag.listw(list.of.weights,variable),
                     z = variable-mean(variable),
                     z = z/sqrt((sum(z^2)/n())),
                     laggedz = lag.listw(list.of.weights, z),
                     lisa = localmoran(variable, list.of.weights, alternative="two.sided")[,1],
                     sig = localmoran(variable, list.of.weights, alternative="two.sided")[,5],
                     quadrant = ifelse(z >= 0, ifelse(laggedz > 0,1,2), ifelse(laggedz >= 0,4,3)))


#########################################
## Regular Colour palette 
#########################################
FunPalette <- function (x=lisa.data, bin.n=6, nians="YlOrRd", style="sd", z.score=FALSE){
  if (z.score == TRUE) {x <- lisa.data$z} else {x <- lisa.data$variable} 
  bins <- classIntervals(x, bin.n, style=style)
  pal <- brewer.pal(length(bins$brks)-1, nians)
  col <- findColours(bins, pal)
  pal.names <- paste(round(bins$brks[-length(bins$brks)],2), 
                     round(bins$brks[-1], 2), sep= " - ")
  pal.n.from <- round(bins$brks[-length(bins$brks)],2)
  pal.n.to <- round(bins$brks[-1],2)
  return(list(col, pal, pal.names, pal.n.from, pal.n.to, bins$brks, z.score))
}
plot(shape, col=FunPalette()[[1]])

#########################################
## Significance: Function for Symbol Size 
#########################################
## takes vector of significances 
## and scaling factor that can be adjusted
## returns list of things used to plot symbol sizes:
## [1] cexez - sizes of symbols
## [2] cex.legend - unique sizes for legend
## [3] cex.names - labels for legend
## [4] pchz - shapes of symbols
## [5] pch.legend - symbols for legend
## [6] bins  - return bins fyi
## [7] sigz - return significances fyi
##########################################

FunSigSize <- function(df=lisa.data, scale=1) {
  sigz <- df[,6]
  bins <- c(0, 0.001, 0.005, 0.01, 0.05,  1 )
  cex.legend <- c(2.60, 2.34, 1.95, 1.30, 0.91)*scale
  pch.legend <- c(19, 19, 19, 19, 4)       # change this for different symbols
  cexez <- cex.legend[findInterval(sigz ,bins, rightmost.closed=T)]
  cex.names <- paste(round(bins[-length(bins)],4),"-",round(bins[-1],4))
  pchz <- pch.legend[findInterval(sigz ,bins, rightmost.closed=T)]
  return(list(cexez, cex.legend, cex.names, pchz, pch.legend,bins))
}
  
FunSigSize()

#########################################
## Colour palette centered on zero - 2-way or 4-way
#########################################
## Takes dataframe "lisa.data" created above
## Either two palettes - to distinguish positive and negative spat.autocorr only
## or four palettes - HH, HL, LL, LH
## steps.pos and steps.neg determines number of categories in each
## Returns list of things used to plot symbol sizes:
## [1] cexez - sizes of symbols
## [2] cex.legend - unique sizes for legend
## [3] cex.names - labels for legend
## [4] pchz - shapes of symbols
## [5] pch.legend - symbols for legend
## [6] bins  - return bins fyi
## [7] sigz - return significances fyi
########################################

FunPaletteMoran <- function (df=lisa.data,
                            steps = 3,
                            four.way = FALSE,
                            pal.hh = "PuBu", 
                            pal.hl = "Oranges",
                            pal.ll = "RdPu",
                            pal.lh = "Greens"){
  if (steps > 7) stop("Sorry, can't have more than seven steps !")
  if (four.way == FALSE) {pal.ll <- pal.hh 
                          pal.lh <- pal.hl}
  bins.pos <- seq(0, max(df$lisa), length=steps+1)
  bins.neg <- seq(min(df$lisa),0, length=steps+1)

  pal.hh <- brewer.pal(steps+2, pal.hh)[-c(1:2)]
  pal.hl <- rev(brewer.pal(steps+2, pal.hl)[-c(1:2)])
  pal.ll <- brewer.pal(steps+2, pal.ll)[-c(1:2)]
  pal.lh <- rev(brewer.pal(steps+2, pal.lh)[-c(1:2)])
  df <- df %>%
    mutate (int = ifelse(lisa >= 0, findInterval(lisa, bins.pos, rightmost.closed=TRUE), 
                         findInterval(lisa, bins.neg, rightmost.closed=TRUE)),
            colz = ifelse(quadrant==1, pal.hh[int], 
                          ifelse(quadrant == 2, pal.hl[int],
                                 ifelse(quadrant == 3, pal.ll[int], pal.lh[int]))))
  pal.n.from <- c(round(bins.neg[-length(bins.neg)],2),round(bins.pos[-length(bins.pos)],2) )
  pal.n.to <- c(round( bins.neg[-1],2), round( bins.pos[-1],2))
  pal.names <- paste( pal.n.from, pal.n.to,sep= " - ")                        
  return(list(df$colz, pal.hh, pal.hl, pal.ll, pal.lh, pal.names))
}

###############################################
## Legends
##############################################

FunLegend <- function(palette=FunPalette){
  plot(c(0,1), c(0,1),  col = NA, ann = FALSE, axes = FALSE)
  points(rep(0.3,length(palette[[2]])), 1/length(palette[[2]])*
           seq(1, length(palette[[2]])),pch=15, col=palette[[2]], cex=2)   
  points(rep(0.3,length(palette[[2]])),1/length(palette[[2]])*
           seq(1, length(palette[[2]])),pch=0, lwd=1, cex=2)
  text(rep(0.6,length(palette[[2]])), 1/length(palette[[2]])*
         seq(1, length(palette[[2]])), palette[[3]])
}

FunLegendSig <- function(sizes=FunSigSize()){
  par(mar=c(0,0,2,0)) 
  plot(c(0,1), c(0,1),  col = NA, ann = FALSE, axes = FALSE)
  title(main="Significance Levels", font.main=1)
  points(rep(0.3,length(sizes[[2]])), 1/length(sizes[[2]])*
           seq(1, length(sizes[[2]])),pch=sizes[[5]], cex=sizes[[2]], col="gray", lwd=2)   
  points(rep(0.3,length(sizes[[2]])), 1/length(sizes[[2]])*
           seq(1, length(sizes[[2]])),  pch=c(rep(1,length(sizes[[5]])-1), NA),
         cex=sizes[[2]][sizes[[5]]!=4], lwd=0.5)
  text(rep(0.7,length(sizes[[2]])), 1/length(sizes[[2]])*
         seq(1, length(sizes[[2]])), sizes[[3]])}

FunLegendMoran <- function(palette=FunPaletteMoran()){
  par(mar=c(0,0,2,0)) 
  plot(c(0,1), c(0,1),  col = NA, ann = FALSE, axes = FALSE)
  title(main="Local Moran's I", font.main=1)
  points(rep(0.3,2*length(palette[[2]])), 1/(2*length(palette[[2]]))*
           seq(1, 2*length(palette[[2]])),pch=19,  col=c(palette[[3]], palette[[2]]), cex=1.5)
  points(rep(0.3,2*length(palette[[2]])), 1/(2*length(palette[[2]]))*
           seq(1, 2*length(palette[[2]])),  cex=1.5)
  if (!identical(palette[[3]], palette[[5]])){
  points(rep(0.2,2*length(palette[[2]])), 1/(2*length(palette[[2]]))*
           seq(1, 2*length(palette[[2]])),pch=19,  col=c(palette[[5]], palette[[4]]), cex=1.5)
  points(rep(0.2,2*length(palette[[2]])), 1/(2*length(palette[[2]]))*
           seq(1, 2*length(palette[[2]])),  cex=1.5)}
  text(rep(0.7,2*length(palette[[2]])), 1/(2*length(palette[[2]]))*
         seq(1, 2*length(palette[[2]])), palette[[6]])
}

########################################################
## Moran Scatterplot
########################################################
## Can be used independently, or called from FunLisaMap
## Takes dataframe "lisa.data" created above
## and logical whether to use z.scores or raw variable values
## extra parameters get passed on to plot
#######################################################

FunMoranScatterplot <- function(lisa.data, z.score=FALSE, 
                                palette= FunPaletteMoran(),
                                sizes = FunSigSize(), ...) {  
  if (z.score == TRUE) {x <- lisa.data$z
                        lx <- lisa.data$laggedz
                        xlab <- "z-score"
                        ylab <- "lagged z-score"} else {
                          x <- lisa.data$variable
                          lx <- lisa.data$lagged
                          xlab <- "variable"
                          ylab <- "lagged variable"} 
  plot(x, lx, col=palette[[1]], cex=sizes[[1]], pch=sizes[[4]], xlab=xlab, ylab=ylab, ...)
  points(x[sizes[[4]]!=4], lx[sizes[[4]]!=4], cex=sizes[[1]][sizes[[4]]!=4])
  abline(v= mean(x), lty=2, lwd=0.5)
  ## option1? x=lx and mean(x) on both
  abline(a=0, b=1, lwd=0.5) 
  abline(h= mean(x), lty=2, lwd=0.5)  
  ## option2? linear model and mean of lagged
  #abline(lm(lx~x), lwd=0.5) 
  #abline(h= mean(lx), lty=2, lwd=0.5)
}
FunMoranScatterplot(lisa.data)

################################################
## Regular cloropleth map of variable (or z-scores) with London Inset
################################################
## Takes dataframe "lisa.data" created above
## and shape SpatialPolygon
## and  palette from FunPalette()
################################################
## CHANGE logical in FunPalette() to switch to z.scores
## Default inset centered on London City, otherwise change
## Depending on your map shape and device size, CHANGE also
## position of Inset
## position variable density distribution
## position legend above dens.
#################################################

FunUKMap <- function(lisa.data, shape, palette=FunPalette()){
  oldpar <- par(no.readonly=TRUE)
  X11(width=10, height=7)  
  par(mar=c(0,0,0,0)) 
  plot(shape, col=palette[[1]])
  par("usr")
  # London inset, change CoL and w1, and h1 to suit
  CoL <- coordinates(shape)[32,] # center on London City (00BK)
  w1 = abs(par("usr")[1] - par("usr")[2])/30
  h1 = abs(par("usr")[3] - par("usr")[4])/30
  rect( CoL[1]-w1, CoL[2]-h1,CoL[1]+w1, CoL[2]+h1, lwd=2)
  leftX <- 0.78 # position London Inset
  rightX <- 0.99
  bottomY <- 0.2
  par(fig=c(leftX,rightX,bottomY, bottomY+rightX-leftX), new =TRUE)
  plot(shape, col=palette[[1]], xlim=c(CoL[1]-w1, CoL[1]+w1),
       ylim=c(CoL[2]-h1, CoL[2]+h1))
  box(lwd=2)  
  leftX <- 0 # position variable density distribution
  rightX <- 0.25
  bottomY <- 0.35
  par(fig=c(leftX,rightX,bottomY, bottomY+rightX-leftX), new =TRUE,mar=c(2,1,0,0))
  if (palette[[7]] == TRUE) {x <- lisa.data$z} else {x <- lisa.data$variable} 
  plot(density(x, na.rm=TRUE), axes=FALSE, ann=FALSE, lwd=2)
  axis(1)
  abline(v=palette[[6]], col="red", lty=c(2,2,2,1,2,2,2))  
  bottomY <- bottomY+rightX-leftX # position legend above dens.
  par(fig=c(leftX,rightX,bottomY, bottomY+rightX-leftX),mar=c(0,0,0,0), new =TRUE)
  FunLegend(palette)
  par(oldpar)
  
}

FunUKMap(lisa.data, shape, palette=FunPalette(z.score=TRUE))

FunUKMoranMap <- function(df=lisa.data, shape, palette=FunPaletteMoran(four.way=TRUE),
                          sizes=FunSigSize()){
  oldpar <- par(no.readonly=TRUE)
  #X11(width=10, height=7)  
  par(mar=c(0,0,0,0)) 
  plot(shape,  lwd=0.5)
  points(coordinates(shape), col=palette[[1]], 
         cex=sizes[[1]],pch=sizes[[4]], lwd=2 )
  points(coordinates(shape)[sizes[[4]]!=4,], cex=sizes[[1]][sizes[[4]]!=4] )
  par("usr")
  # London inset, change CoL and w1, and h1 to suit
  CoL <- coordinates(shape)[32,] # center on London City (00BK)
  w1 = abs(par("usr")[1] - par("usr")[2])/30
  h1 = abs(par("usr")[3] - par("usr")[4])/30
  rect( CoL[1]-w1, CoL[2]-h1,CoL[1]+w1, CoL[2]+h1, lwd=2)
  leftX <- 0.78 # position London Inset
  rightX <- 0.99
  bottomY <- 0.2
  par(fig=c(leftX,rightX,bottomY, bottomY+rightX-leftX), new =TRUE)
  plot(shape,  lwd=0.5,xlim=c(CoL[1]-w1, CoL[1]+w1),
       ylim=c(CoL[2]-h1, CoL[2]+h1))
  points(coordinates(shape), col=palette[[1]], 
         cex=sizes[[1]],pch=sizes[[4]], lwd=2 )
  points(coordinates(shape)[sizes[[4]]!=4,], cex=sizes[[1]][sizes[[4]]!=4] )
  box(lwd=2)  
  leftX <- 0 # position legend significances
  rightX <- 0.30
  bottomY <- 0.3
  par(fig=c(leftX,rightX,bottomY, bottomY+rightX-leftX), new =TRUE,mar=c(2,1,0,0))
  FunLegendSig(sizes)
  bottomY <- bottomY+rightX-leftX# position legend Morans I
  par(fig=c(leftX,rightX,bottomY, bottomY+rightX-leftX), new =TRUE,mar=c(2,1,0,0))
  FunLegendMoran(palette)
  par(oldpar)
  leftX <- 0.70 # position moran scatterplot
  rightX <- 0.99
  bottomY <- 0.65
  par(fig=c(leftX,rightX,bottomY, bottomY+rightX-leftX), new =TRUE,mar=c(0,0,0,0))
  FunMoranScatterplot(df, palette=palette, sizes=sizes)
  par(oldpar)
  
}
  
FunUKMoranMap(lisa.data, shape, palette=FunPaletteMoran(four.way=TRUE))

