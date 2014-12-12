###################################################
## Four plots of lexis graph to explain how it works:
## 1. normal
## 2. period effect
## 3. age effect
## 4. cohort effect
##################################################
require(Epi)

# A small bogus cohort (reworked from Lexis() example)
xcoh <- structure( list( id = c("A", "B", "C"),
                         birth = c("14/07/1922", "01/04/1937", "10/06/1947"),
                         entry = c("04/08/1925", "08/09/1972", "23/12/1991"),
                         exit = c("27/06/1997", "23/05/1995", "24/07/1998"),
                         fail = c(1, 0, 1) ),
                   .Names = c("id", "birth", "entry", "exit", "fail"),
                   row.names = c("1", "2", "3"),
                   class = "data.frame" )

# Convert the character dates into numerical variables (fractional years)
xcoh$bt <- cal.yr( xcoh$birth, format="%d/%m/%Y" )
xcoh$en <- cal.yr( xcoh$entry, format="%d/%m/%Y" )
xcoh$ex <- cal.yr( xcoh$exit , format="%d/%m/%Y" )

# Define as Lexis object with timescales calendar time and age
Lcoh <- Lexis( entry = list( CT=en ),
               exit = list( CT=ex, A=ex-bt ),
               exit.status = fail,
               data = xcoh )

# lexis1
plot( Lcoh, grid=list(seq(1920,2000,10), seq(0,80,10)), col="black", lwd=2,
      xlab="Calendar time", ylab="Age")
dev.copy2eps(file="Figures/lexis01.eps")

# lexis2
plot( Lcoh, grid=list(seq(1920,2000,10), seq(0,80,10)), col="black", lwd=2,
      xlab="Calendar time", ylab="Age")
rect(1980,0,1990,par("usr")[4], col="gray90", border="black" )
lines(Lcoh,col="black", lwd=2)
dev.copy2eps(file="Figures/lexis02.eps")

# lexis3
plot( Lcoh, grid=list(seq(1920,2000,10), seq(0,80,10)), col="black", lwd=2,
      xlab="Calendar time", ylab="Age")
rect(par("usr")[1],50,par("usr")[2],60, col="gray90", border="black" )
lines(Lcoh,col="black", lwd=2)
dev.copy2eps(file="Figures/lexis03.eps")

# lexis4
plot( Lcoh, grid=list(seq(1920,2000,10), seq(0,80,10)), col="black", lwd=2,
      xlab="Calendar time", ylab="Age")
polygon(c(1940,1950,par("usr")[2], par("usr")[2]), 
        c(0,0,par("usr")[2]-1950,par("usr")[2]-1940),
        col="gray90")
lines(Lcoh,col="black", lwd=2)
dev.copy2eps(file="Figures/lexis04.eps")
