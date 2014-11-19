library(spdep)
library(RColorBrewer)
library(classInt)
library(rgdal)
library(maptools)
gpclibPermit()
library(dplyr)
library(vcd)


#########################################
## Load data and Shape file
## For this example:
## LAshape - E&W Local authority SpatialPolygon with administrative names
## LAdata - 2001 census proportion of persons living in detached hosuing
########################################
## TO CHANGE:
## You need 
oldpar <- par(no.readonly=TRUE)
load("Data/EnglandWales.R")



#########################################
## Data frame setup
#########################################
## Start with single column dataframe with local statistic only. 
## And appropriate list of weights (e.g. from nb2listw)
## This calculates
## - local moran stats
## - lagged values
## - quadrant
#########################################


FunDataSetup <- function(data, list.of.weights){
  data <- cbind(data, localmoran(data[,1], list.of.weights, alternative="two.sided", na.action=na.omit,
                                 zero.policy=TRUE))
  data$lagged <- (lag.listw(list.of.weights,data[,1]))
  mean.true <- mean(data[,1])
  mean.lagged <- mean(data$lagged)
  data$quadrant <- apply(data, 1, function(r) 
    if(r[1] >= mean.true){
      if(r[7] >= mean.lagged){1} else{2}} else {
        if(r[7] >= mean.lagged){3} else{4}})
  return(data)
}

LISA.data <- FunDataSetup(start.data, lw.foq)
  

plot(LISA.data[,1], lag(lw.foq, LISA.data[,1], NAOK=TRUE), col=data[,8], pch= data[,2]/abs(data[,2])+2)

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




# # neighbourhoods and lists
# shape.UK <-unionSpatialPolygons(shape.regions, rep(1,10), threshold=400000)
# nb.foq <- poly2nb(shape, snap=1)
# # check for no neighbours ans attach them back
# which(card(nb.foq)==0)
# nb.foq[[115]] <- as.integer(116)
# nb.foq[116] <- list(sort(c(nb.foq[[116]], as.integer(115))))
# nb.foq[[114]] <- as.integer(226)
# nb.foq[226] <- list(sort(c(nb.foq[[226]], as.integer(114))))
# which(card(nb.foq)==0)
# # create list of weights
# lw.foq <- nb2listw(nb.foq, style='W')
# # Then we need regions too
# nb.foq.gor <- poly2nb(shape.regions)
# lw.foq.gor <- nb2listw(nb.foq.gor, style='W')
# # and for counties
# nb.foq.cnty <- poly2nb(shape.counties)
# which(card(nb.foq.cnty)==0)
# nb.foq.cnty[[41]] <- as.integer(16)
# nb.foq.cnty[16] <- list(sort(as.integer(c(unlist(nb.foq.cnty[16]), 41))))
# lw.foq.gor <- nb2listw(nb.foq.cnty, style='W')

###
###BEGINING OF FUNCTIONS
###

# ODDS RATIO
FunOdds <-function(x) {
  x[1,1]*x[2,2]/(x[1,2]*x[2,1])
}

# NORMAL PALETTE
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

# new variable but stolen pallete
FunPaletteStolen <- function(x, pallette) {
  bins <- pallette[[6]]
  pal <- pallette[[2]]
  col <- pal[findInterval(x, bins, rightmost.closed=TRUE)]
  pal.n.from <- pallette[[4]]
  pal.n.to <- pallette[[5]]
  pal.names <- pallette[[3]]
  old <-  pallette[[7]]
  return(list(col, pal, pal.names, pal.n.from, pal.n.to, bins, x,old))
}
x<-lisa[,1]
# palete for significances and sizes 
FunPaletteSigSize <- function(x) {
  bins <- c(0, 0.001, 0.005, 0.01, 0.05,  1 )
  cex.pal <- c(2, 1.8, 1.5, 1, 0.7)*1.3
  cexez <- cex.pal[findInterval(x ,bins, rightmost.closed=T)]
  cex.names <- paste(round(bins[-length(bins)],4),"-",round(bins[-1],4))
  cex.n.from <- round(bins[-length(bins)],4)
  cex.n.to <- round(bins[-1],4)
  pch.pal <- c(19, 19, 19, 19, 4)
  pchz <- pch.pal[findInterval(x ,bins, rightmost.closed=T)]
  return(list(cexez, cex.pal, cex.names, cex.n.from, cex.n.to, bins, x, pchz,
              pch.pal))
}

#  PALETTE centered on ZERO
FunPaletteZero <- function (x, b.neg=3, b.poz=3,nians.poz="PuBuGn", nians.neg="RdPu"){
  x.nr <- x[!is.na(x)]
  bins <- c(seq(min(x.nr),0, length=b.neg+1), seq(0, max(x.nr), length=b.poz+1))
  bins <- bins[-(b.neg+1)]
  pal <- c(rev(brewer.pal(b.neg+2, nians.neg)), brewer.pal(b.poz+2, nians.poz))
  pal <- pal[-c((b.neg+1):(b.neg+4))]
  col <- pal[findInterval(x, bins, rightmost.closed=TRUE)]
  pal.n.from <- round(bins[-length(bins)],2)
  pal.n.to <- round(bins[-1],2)
  pal.names <- paste(round(bins[-length(bins)],2), round(bins[-1], 2), 
                     sep= " - ")
  return(list(col,  pal, pal.names, pal.n.from, pal.n.to,bins, x))
}

#  PALETTE centered on ZERO with NAs for logodds
FunPaletteZeroLO <- function (x, b.neg=3, b.poz=3, nians.neg="Reds", 
                              nians.poz="Blues" ){
  x.nr <- x[!is.na(x)]
  bins <- c(seq(min(x.nr),0, length=b.neg+1), seq(0, max(x.nr), length=b.poz+1))
  bins <- bins[-(b.neg+1)]
  pal <- c(rev(brewer.pal(b.neg+1, nians.neg)), brewer.pal(b.poz+1, nians.poz))
  pal <- pal[-c((b.neg+1):(b.neg+2))]
  col <- pal[findInterval(x, bins, rightmost.closed=TRUE)]
  pal.n.from <- round(bins[-length(bins)],2)
  pal.n.to <- round(bins[-1],2)
  pal.names <- paste(round(bins[-length(bins)],2), round(bins[-1], 2), 
                     sep= " - ")
  return(list(col,  pal, pal.names, pal.n.from, pal.n.to,bins, x))
}

# PLOT 1
# Plot function for County or Region map - no inset - cloropleth
FunUKMapCnty <- function(shape, palette){
  par(oldpar)
  X11(width=10, height=7)
  par(mar=c(0,0,0,0)) 
  plot(shape, col=palette[[1]], lwd=0.5, xlim=c(100000, 900000))
  par(fig=c(0,0.17, 0.35, 0.6), new =TRUE,mar=c(4,0,0,0))
  d1 <- density(palette[[8]], na.rm=TRUE) 
  d2 <- density(palette[[7]], na.rm=TRUE) 
  plot(range(d1$x), range( d2$y), type = "n", xlab = "x", ylab = "Density",
       axes=FALSE, ann=FALSE) 
  lines(d2, lwd=2)
  axis(1)
  abline(v=palette[[6]], col="red", lty=2)
  par(fig=c(0,0.27, 0.6, 1), new =TRUE,mar=c(0.3,0,0,0) )
  plot(c(0,1), c(0,1),  col = NA, ann = FALSE, axes = FALSE)
  points(rep(0.05,length(palette[[2]])), 0.08*seq(1, length(palette[[2]])),pch=15, 
         col=palette[[2]], cex=2)   
  points(rep(0.05,length(palette[[2]])), 0.08*seq(1, length(palette[[2]])),pch=0, lwd=1,
         cex=2)
  text(rep(0.2,length(palette[[2]])), 0.08*seq(1, length(palette[[2]])), palette[[4]] )
  text(rep(0.35,length(palette[[2]])), 0.08*seq(1, length(palette[[2]])), rep(" - ", length(palette[[2]]))  )
  text(rep(0.5,length(palette[[2]])), 0.08*seq(1, length(palette[[2]])), palette[[5]]  )
  par(oldpar)
}

# PLOT 2
# plot function for regular UK map, cloropleth 
FunUKMap <- function(shape, palette){
  par(oldpar)
  X11(width=10, height=7)
  par(mar=c(0,0,0,0)) 
  plot(shape, col=palette[[1]], lwd=0.5, xlim=c(100000, 900000))
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
  plot(shape, col=palette[[1]], xlim=c(505000,560000), ylim=c(155000,200000), lwd=0.5)
  plot(shape.regions[7], xlim=c(505000,560000), ylim=c(155000,200000), lwd=1.5, add=TRUE)
  box(lwd=2)
  par(fig=c(0,0.17, 0.35, 0.6), new =TRUE,mar=c(4,0,0,0))
  plot(density(palette[[7]], na.rm=TRUE), axes=FALSE, ann=FALSE, lwd=2)
  axis(1)
  abline(v=palette[[6]], col="red", lty=2)
  par(fig=c(0,0.27, 0.6, 1), new =TRUE,mar=c(0.3,0,0,0) )
  plot(c(0,1), c(0,1),  col = NA, ann = FALSE, axes = FALSE)
  points(rep(0.05,length(palette[[2]])), 0.08*seq(1, length(palette[[2]])),pch=15, 
         col=palette[[2]], cex=2)   
  points(rep(0.05,length(palette[[2]])), 0.08*seq(1, length(palette[[2]])),pch=0, lwd=1,
         cex=2)
  text(rep(0.2,length(palette[[2]])), 0.08*seq(1, length(palette[[2]])), palette[[4]] )
  text(rep(0.35,length(palette[[2]])), 0.08*seq(1, length(palette[[2]])), rep(" - ", length(palette[[2]]))  )
  text(rep(0.5,length(palette[[2]])), 0.08*seq(1, length(palette[[2]])), palette[[5]]  )
  par(oldpar)
}

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