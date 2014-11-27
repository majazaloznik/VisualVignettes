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
                            pal.hl = "RdPu",
                            pal.ll = "Greens",
                            pal.lh = "Oranges"){
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

# PLOT 3
# plot function for regular UK map, points an dpoint sizes! 
FunUKMapPoints <- function(shape, palette, cex.palette){
  par(oldpar)
  X11(width=10, height=7)
  par(mar=c(0,0,0,0)) 
  xx<-cex.palette[[8]]
  plot(shape,  lwd=0.5, xlim=c(100000, 900000))
  points(coordinates(shape), col=palette[[1]], 
         cex=cex.palette[[1]],pch=cex.palette[[8]], lwd=2 )
  points(coordinates(shape)[xx!=4,], cex=cex.palette[[1]][xx!=4], lwd=0.5 )
  a=55000
  b=45000
  polygon(x=c(505000,505000+a,505000+a,505000), 
          y=c(155000,155000,155000+b,155000+b), lwd=2)
  lines(x=c(505000, 700000), y=c(155000, 100000), lwd=1)
  box.ratio <- a/b
  plot.ratio <-(par()$usr[2] - par()$usr[1]) /(par()$usr[4] - par()$usr[3])
  x1 <- (-par()$usr[1]+700000)/(par()$usr[2] - par()$usr[1])
  y1 <- (-par()$usr[3]+100000)/(par()$usr[4] - par()$usr[3])
  par(fig=c(x1,x1+0.25*box.ratio/plot.ratio, y1, y1+0.25), new =TRUE)
  plot(shape,  xlim=c(505000,560000), ylim=c(155000,200000), 
       lwd=0.5)
  points(coordinates(shape), col=palette[[1]], cex=cex.palette[[1]]*1.5, 
         pch=cex.palette[[8]] , lwd=2)
  points(coordinates(shape)[xx!=4,], cex=cex.palette[[1]][xx!=4]*1.5, lwd=0.5 )
  plot(shape.regions[7], xlim=c(505000,560000), ylim=c(155000,200000), 
       lwd=2, add=TRUE)
  box(lwd=2)
  par(fig=c(0,0.27, 0.35, 0.75), new =TRUE,mar=c(0.3,0,0,0))
  plot(c(0,1), c(0,1),  col = NA, ann = FALSE, axes = FALSE)
  points(rep(0.05,length(cex.palette[[2]])), 
         0.08*seq(1, length(cex.palette[[2]])), 
         col="gray", pch=cex.palette[[9]], cex=cex.palette[[2]], lwd=2) 
  points(rep(0.05,length(cex.palette[[2]])), 
         0.08*seq(1, length(cex.palette[[2]])), 
         pch=c(rep(1,length(cex.palette[[9]])-1), NA),
         cex=cex.palette[[2]][cex.palette[[9]]!=4], lwd=0.5)
  text(rep(0.2,length(cex.palette[[2]])), 0.08*seq(1, length(cex.palette[[2]])), cex.palette[[4]] )
  text(rep(0.35,length(cex.palette[[2]])), 0.08*seq(1, length(cex.palette[[2]])), rep(" - ", length(palette[[2]]))  )
  text(rep(0.5,length(cex.palette[[2]])), 0.08*seq(1, length(cex.palette[[2]])), cex.palette[[5]]  )  
  text(0.35, 0.57, "Significance level")
  par(fig=c(0,0.27, 0.6, 1), new =TRUE,mar=c(0.3,0,0,0) )
  plot(c(0,1), c(0,1),  col = NA, ann = FALSE, axes = FALSE)
  points(rep(0.05,length(palette[[2]])), 0.08*seq(1, length(palette[[2]])), 
         col=palette[[2]], pch=19, cex=2)   
  points(rep(0.05,length(palette[[2]])), 0.08*seq(1, length(palette[[2]])),
         pch=1, lwd=0.5, cex=2)
  text(rep(0.2,length(palette[[2]])), 0.08*seq(1, length(palette[[2]])), palette[[4]] )
  text(rep(0.35,length(palette[[2]])), 0.08*seq(1, length(palette[[2]])), rep(" - ", length(palette[[2]]))  )
  text(rep(0.5,length(palette[[2]])), 0.08*seq(1, length(palette[[2]])), palette[[5]]  )
  text(0.35, 0.65, "Local Moran's I")
  par(oldpar)
}

# Plot MoransI plot using size and col palette
FunPlotMoran <- function(x, lw.foq, ppZ, ppS){
  X11(width=7, height=5)
  plot(x, lag(lw.foq, x, NAOK=TRUE), col=ppZ[[1]], cex=ppS[[1]], pch=ppS[[8]],
       ylab="lagged x", axes=FALSE) 
  points(x[ppS[[8]]!=4], lag(lw.foq, x, NAOK=TRUE)[ppS[[8]]!=4], 
         cex=ppS[[1]][ppS[[8]]!=4], lwd=0.5) 
  abline(h=mean(lag(lw.foq, x, NAOK=TRUE), na.rm=TRUE), lty=2, lwd=0.5)
  abline(v=mean(x, na.rm=TRUE), lty=2,lwd=0.5)
  abline(lm(lag(lw.foq, x, NAOK=TRUE)~x), lwd=0.5)
  axis(1)
  axis(2)
}

# Plot lm b~a which also returns palette of residuals
FunPlotLM <- function(N, b.neg=3, b.poz=3){
  X11(width=7, height=5)
  a<-as.vector(collapsed.la[[N]][,5])
  b<-as.vector(collapsed.la[[N]][,6])
  res <- lm(b~a)$residual
  ppR<- FunPaletteZeroLO(res, nians.neg="Oranges", 
                         nians.poz="Greens",
                         b.neg=b.neg, b.poz=b.poz)
  plot(a,b, xlab=colnames(collapsed.la[[N]])[5],col=ppR[[1]],
       ylab=colnames(collapsed.la[[N]])[6], axes=FALSE, pch=19, cex=1.5) 
  points(a,b, cex=1.5, lwd=0.5)
  abline(lm(b~a), lwd=0.5)
  axis(1)
  axis(2)
  return(list(ppR, res))
}

FunMosaic <- function(x, split_vertical=c(TRUE, FALSE), lwd=2,...){
  X11(height=7, width=7)
  mosaic(x, spacing=spacing_equal(unit(1, "lines")),  labeling=NULL,...)
}

###
###END OF FUNCTIONS
###

# find out which tables have no missing logodds!
x<-sapply(collapsed.la, function(x)  sum(is.na(as.vector(x[,8]))))
good <- x==0
select.stats[,9]<- as.character(select.stats[,9])
select.stats[,12]<- seq(1, dim(select.stats)[1])
good.select.stats<-select.stats[good,]

nn<- order(good.select.stats[,10])[3300:3316]
good.select.stats[nn,9]
good.select.stats[10,12]
select.stats[15,]
N=15
select.stats[N,9]
#two
good.select.stats[174,12]
select.stats[403,]
N=403
select.stats[N,9]
# map of a 
a<-as.vector(collapsed.la[[N]][,5])
ppA<- FunPalette(a, bin.n=6)
FunUKMap(shape, ppA)
dev.copy2eps(file="FPt01.eps", family="ComputerModern")

# map of a lisa
lisa<-localmoran(a, lw.foq, alternative="two.sided", na.action=na.omit,
                 zero.policy=TRUE)
ppZ <- FunPaletteZero(lisa[,1])
ppS <- FunPaletteSigSize(lisa[,5])
FunUKMapPoints(shape, ppZ, ppS)
dev.copy2eps(file="FPt02.eps", family="ComputerModern")
# moranplot
FunPlotMoran(a, lw.foq, ppZ, ppS)
dev.copy2eps(file="FPt021.eps", family="ComputerModern")

# map of b
b<-as.vector(collapsed.la[[N]][,6])
ppB<- FunPalette(b, bin.n=6)
FunUKMap(shape, ppB)
dev.copy2eps(file="FPt03.eps", family="ComputerModern")

# map of b lisa
lisa<-localmoran(b, lw.foq, alternative="two.sided", na.action=na.omit,
                 zero.policy=TRUE)
ppZ <- FunPaletteZero(lisa[,1])
ppS <- FunPaletteSigSize(lisa[,5])
FunUKMapPoints(shape, ppZ, ppS)
dev.copy2eps(file="FPt04.eps", family="ComputerModern")
# moranplot
FunPlotMoran(b, lw.foq, ppZ, ppS)
dev.copy2eps(file="FPt041.eps", family="ComputerModern")

#  residual  and scatterplot
lm.stuff <- FunPlotLM(N, b.neg=5, b.poz=3)
dev.copy2eps(file="FPt05.eps", family="ComputerModern")

#  residual  map
FunUKMap(shape, lm.stuff[[1]])
dev.copy2eps(file="FPt06.eps", family="ComputerModern")

# map of residual lisa
lisa<-localmoran(lm.stuff[[2]], lw.foq, alternative="two.sided", na.action=na.omit,
                 zero.policy=TRUE)
ppZ <- FunPaletteZero(lisa[,1])
ppS <- FunPaletteSigSize(lisa[,5])
FunUKMapPoints(shape, ppZ, ppS)
dev.copy2eps(file="FPt07.eps", family="ComputerModern")
# moranplot
FunPlotMoran(lm.stuff[[2]], lw.foq, ppZ, ppS)
dev.copy2eps(file="FPt071.eps", family="ComputerModern")
log(1.4)
# logodds 
N=403
x<-as.vector(collapsed.la[[N]][,8])
pp<- FunPaletteZeroLO(x, b.neg=5, b.poz=2)
FunUKMap(shape, pp)
dev.copy2eps(file="FPt08.eps", family="ComputerModern")
# moran of logods
lisa<-localmoran(x, lw.foq, alternative="two.sided", na.action=na.omit,
                 zero.policy=TRUE)
# removed <- as.vector(attributes(lisa)$na.action)
# LI<-rep(NA, length(x))
# LI[-removed] <- lisa[,1]
# ppZ <- FunPaletteZero(LI)
# LP <-rep(NA, length(x))
# LP[-removed] <- lisa[,5]
# ppS <- FunPaletteSigSize(LP)
which.max(x)
cbind(x, lisa[,1], lisa[,5])
colnames(lisa)
ppZ <- FunPaletteZero(lisa[,1])
ppS <- FunPaletteSigSize(lisa[,5])
FunUKMapPoints(shape, ppZ, ppS)
dev.copy2eps(file="FPt09.eps", family="ComputerModern")
FunPlotMoran(x, lw.foq, ppZ, ppS)
dev.copy2eps(file="FPt10.eps", family="ComputerModern")

N=15
pp[[3]]
nat.15 <- matrix(apply(collapsed.la[[15]][,1:4],2, sum), c(2,2), byrow=TRUE)
which.max(collapsed.la[[15]][,7])
max.15 <- matrix(collapsed.la[[15]][19,1:4], c(2,2), byrow=TRUE)
which.min(collapsed.la[[15]][,7])
min.15 <- matrix(collapsed.la[[15]][262,1:4], c(2,2), byrow=TRUE)
FunOdds(max.15)
FunOdds(min.15)
FunOdds(nat.15)
odz <- c(FunOdds(max.15), FunOdds(min.15), FunOdds(nat.15))
log(odz)
col <- pp[[2]][findInterval(log(odz), pp[[6]], rightmost.closed=TRUE)]

FunMosaic(max.15, gp=gpar(fill=c(col[1])))
dev.copy2eps(file="FP11.eps", family="ComputerModern")
FunMosaic(min.15, gp=gpar(fill=c(col[2])))
dev.copy2eps(file="FP12.eps", family="ComputerModern")
FunMosaic(nat.15, gp=gpar(fill=c(col[3])))
dev.copy2eps(file="FP13.eps", family="ComputerModern")


N=403
pp[[2]]
nat.403 <- matrix(apply(collapsed.la[[403]][,1:4],2, sum), c(2,2), byrow=TRUE)
which.max(collapsed.la[[403]][,7])
max.403 <- matrix(collapsed.la[[403]][199,1:4], c(2,2), byrow=TRUE)
which.min(collapsed.la[[403]][,7])
min.403 <- matrix(collapsed.la[[403]][31,1:4], c(2,2), byrow=TRUE)
FunOdds(max.403)
FunOdds(min.403)
FunOdds(nat.403)
odz <- c(FunOdds(max.403), FunOdds(min.403), FunOdds(nat.403))
log(odz)
col <- pp[[2]][findInterval(log(odz), pp[[6]], rightmost.closed=TRUE)]

FunMosaic(max.403, gp=gpar(fill=c(col[1])))
dev.copy2eps(file="FPt11.eps", family="ComputerModern")
FunMosaic(min.403, gp=gpar(fill=c(col[2])))
dev.copy2eps(file="FPt12.eps", family="ComputerModern")
FunMosaic(nat.403, gp=gpar(fill=c(col[3])))
dev.copy2eps(file="FPt13.eps", family="ComputerModern")


### Plots of Geographical hierarchies in UK
rat <- dist(shape@bbox[2,])/dist(shape@bbox[1,])
X11(width=4, height=4*rat)
par(mar=c(0,0,0,0))
plot(shape.UK,  col="palegoldenrod", lwd=2)
dev.copy2eps(file="UK01.eps", family="ComputerModern")

X11(width=4, height=4*rat)
par(mar=c(0,0,0,0))
plot(shape, col="palegoldenrod", lwd=0.5)
plot(shape.UK, lwd=2, add=TRUE)
dev.copy2eps(file="UK02.eps", family="ComputerModern")

X11(width=4, height=4*rat)
par(mar=c(0,0,0,0))
plot(shape, col="palegoldenrod", lwd=0.5, border="gray50")
plot(shape.counties, lwd=2, add=TRUE)
dev.copy2eps(file="UK03.eps", family="ComputerModern")

X11(width=4, height=4*rat)
par(mar=c(0,0,0,0))
plot(shape, col="palegoldenrod", lwd=0.5, border="gray50")
plot(shape.regions, lwd=2, add=TRUE)
dev.copy2eps(file="UK04.eps", family="ComputerModern")
par(oldpar)


2624857*58

# Plots of FOQ might be useful
X11(width=4, height=4*rat)
par(mar=c(0,0,0,0))
plot(shape.UK)
plot(nb.foq, coordinates(shape), add=TRUE)
dev.copy2eps(file="UK05.eps", family="ComputerModern")

plot(shape.counties)
plot(nb.foq.cnty, coordinates(shape.counties),add=TRUE)
plot(shape.regions)
plot(nb.foq.gor, coordinates(shape.regions),add=TRUE)
par(oldpar)

var.cats1 <- sapply(tables.gor, function(x) dim(x)[2])
var.cats2 <- sapply(tables.gor, function(x) dim(x)[2])
sum(var.cats1*var.cats2)
var.cats <- c(4, var.cats[1:56])
var.cats <- as.vector(var.cats)
prod(var.cats)
rm(tables.gor)

# seems 4 by 13 becomes 28,665 tables!?
mos01 <- tables.la[[1]]
x <- margin.table(mos01, c(1,2))
x <-x[,c(1,8,2,3,4,5,6,7,9,10,11,12,13)]
FunMosaic(x, labeling=TRUE, labeling_args=dimnames(x))
dev.copy2eps(file="MOS01.eps", family="ComputerModern")
y <- matrix(apply(collapsed.la[[1]][,1:4],2, sum), c(2,2), byrow=TRUE)
dimnames(x)[1][1]
dimnames(y)
ny <- dimnames(x)
ny[[1]] <- ny[[1]][1]
ny[[1]] <- c(ny[[1]], "!1")
ny[[2]] <- ny[[2]][1]
ny[[2]][1] <- "0-5"
ny[[2]][2] <- "!0-5"
dimnames(y) <- ny
FunMosaic(y)
dev.copy2eps(file="MOS02.eps", family="ComputerModern")

mar01 <- matrix(c(1,1), 1,2)
FunMosaic(mar01)
dev.copy2eps(file="MOS03.eps", family="ComputerModern")

mar02 <- matrix(c(1,1), 2,1)
FunMosaic(mar02)
dev.copy2eps(file="MOS04.eps", family="ComputerModern")

tt <- matrix(c(1,1,1,1), 2,2)
start <- matrix(c(4,1,1,1), 2,2)
odds01 <-loglin(tt, c(1,2), start=start, fit=TRUE)$fit
FunMosaic(odds01)
dev.copy2eps(file="MOS05.eps", family="ComputerModern")
start <- matrix(c(10,1,1,1), 2,2)
odds02 <-loglin(tt, c(1,2), start=start, fit=TRUE)$fit
FunMosaic(odds02)
dev.copy2eps(file="MOS06.eps", family="ComputerModern")
start <- matrix(c(0.5,1,1,1), 2,2)
odds03 <-loglin(tt, c(1,2), start=start, fit=TRUE)$fit
FunMosaic(odds03)
dev.copy2eps(file="MOS07.eps", family="ComputerModern")
start <- matrix(c(1,1,1,1), 2,2)
odds05 <-loglin(tt, c(1,2), start=start, fit=TRUE)$fit
FunMosaic(odds05)
dev.copy2eps(file="MOS08.eps", family="ComputerModern")


colnames(collapsed.la[[15]])
y <- matrix(apply(collapsed.la[[15]][,1:4],2, sum), c(2,2), byrow=TRUE)
sum(y[,1])/sum(y)
sum(y[1,])/sum(y)
select.stats[15,]

sort(collapsed.la[[403]][,8])[370:374]

y <- matrix(apply(collapsed.la[[403]][,1:4],2, sum), c(2,2), byrow=TRUE)
sum(y[,1])/sum(y)
sum(y[1,])/sum(y)
select.stats[15,]


### OK, canI do a hierarchical one quick!

N=403
x<-as.vector(collapsed.la[[N]][,8])
pp<- FunPaletteZeroLO(x, b.neg=5, b.poz=2)
FunUKMap(shape, pp)
dev.copy2eps(file="H01.eps", family="ComputerModern")
x<-as.vector(collapsed.cnty[[N]][,8])
pp2 <- FunPaletteStolen(x, pp)
FunUKMapCnty(shape.counties, pp2)
dev.copy2eps(file="H02.eps", family="ComputerModern")
x<-as.vector(collapsed.gor[[N]][,8])
pp2 <- FunPaletteStolen(x, pp)
FunUKMapCnty(shape.regions, pp2)
dev.copy2eps(file="H03.eps", family="ComputerModern")
rbind(select.stats[403,], select.stats.cnty[403,],select.stats.gor[403,])



N=15
x<-as.vector(collapsed.la[[N]][,8])
pp<- FunPaletteZeroLO(x, b.neg=5, b.poz=2)
FunUKMap(shape, pp)
dev.copy2eps(file="H04.eps", family="ComputerModern")
x<-as.vector(collapsed.cnty[[N]][,8])
pp2 <- FunPaletteStolen(x, pp)
FunUKMapCnty(shape.counties, pp2)
dev.copy2eps(file="H05.eps", family="ComputerModern")
x<-as.vector(collapsed.gor[[N]][,8])
pp2 <- FunPaletteStolen(x, pp)
FunUKMapCnty(shape.regions, pp2)
dev.copy2eps(file="H06.eps", family="ComputerModern")
rbind(select.stats[15,], select.stats.cnty[15,],select.stats.gor[15,])


xx <- lapply(collapsed.la[101:200],function(x) cbind(localmoran(x[,6], lw.foq, alternative="two.sided")[,c(1,5)],
                                                   x[,6]-mean(x[,6])))

lapply(xx, function(x)  filter(data.frame(x),Pr.z....0. <= 0.05 , Ii<0))

head(data.frame(xx[[1]]))
lapply(collapsed.la[1:10],function(x) dimnames(x)[[2]][5])
data.frame(localmoran(collapsed.la[[4]][,5], lw.foq))

FunPalette <- function (x, bin.n=6, nians="YlOrRd", style="sd"){
  x.nr <- x[!is.na(x)]
  bins <- classIntervals(x.nr, bin.n, style=style)
  pal <- brewer.pal(length(bins$brks)-1, nians)
  col <- findColours(bins, pal)
  pal.names <- paste(round(bins$brks[-length(bins$brks)],2), 
                     round(bins$brks[-1], 2), sep= " - ")
  pal.n.from <- round(bins$brks[-length(bins$brks)],2)
  pal.n.to <- round(bins$brks[-1],2)
  return(list(col, pal, pal.names, pal.n.from, pal.n.to, bins$brks, x))
}

pp <- FunPalette(collapsed.la[[4]][,5])
plot(shape, col=pp[[1]])
FunUKMap(shape,FunPalette(collapsed.la[[4]][,5]) )

distwrk <- data.frame(variable=collapsed.la[[37]][,6])
test <- data.frame(variable=collapsed.la[[146]][,6])
