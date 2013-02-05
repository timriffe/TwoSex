
# Author: triffe
###############################################################################
# check against HFD sums:
HFD.all <- read.table("/home/triffe/DATA/HFD/birthsRR.txt", 
                            header = TRUE, skip = 2, as.is = TRUE, stringsAsFactors = FALSE)
USA <- HFD.all[HFD.all$Code == "USA", ]
all.years <- tapply(USA$Total, USA$Year, sum)
my.years <- all.years[as.character(1969:2010)]
plot(1969:2010, my.years, type = 'l') # the HFD totals


# get my totals:
All <- paste0("Bxy",1969:2010,".Rdata")
#All.paths <- file.path("/home/triffe/git/Dissertation/DISSERTATION/DATA/USbirths", All)
All.paths <- file.path("/home/triffe/DATA/CDC/BIRTHS/Bxy", All)

All.my <- lapply(All.paths, function(x){
            local(get(load(x)))
        })
myB <- unlist(lapply(All.my, function(x){
            sum(x$BIRTHS)
        }))
# I have somewhat higher totals. not sure why
plot(1969:2010, myB, type = 'l')
lines(1969:2010, my.years, col = "red", lty = 2)
# probably HFD removed non-resident births...
#--------------------------------------------------------------------
# Try Adrien's decomposition of mean age? for males and females? 2d?
# it had a 'quantum' component I think. Maybe that's worth meting out
# load in Bxy
B                   <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy10_65.Rdata")))

#Bcube               <- abind::abind(B, along = 3)
#array_names         <- dimnames(Bcube)
#names(array_names)  <- c("Males", "Females", "Year")
#dimnames(Bcube)     <- array_names
E                   <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))

# some helper functions 
dimnames(B[[1]]) # Males in rows and Females in columns
# loads wmean() et al for MAC
source("/home/triffe/git/DISS/R/MeanFunctions.R")

Fxm     <- Fxf  <- list()
ages    <- 10:65
years   <- 1969:2009
# male and female rate matrices kept oriented same way
for (yr in 1:length(years)){
    Fxm[[yr]] <- B[[yr]] / E[E$Year == years[yr] & E$Age %in% ages, "Male"]
    Fxf[[yr]] <- t(t(B[[yr]]) / E[E$Year == years[yr] & E$Age %in% ages, "Female"])
}

MACm    <- TFRm     <- MACf    <- TFRf      <- vector(length = length(years))
MACim   <- TFRim    <- MACif   <- TFRif     <- matrix(ncol = length(years), nrow = length(ages), 
                                                    dimnames = list(Age = ages, Year = years))
age.mids    <- ages + .5
plot(rowSums(Fxm[[yr]], na.rm = TRUE) )
plot(rowSums(Fxm[[yr]] * age.mids, na.rm = TRUE) / rowSums(Fxm[[yr]], na.rm = TRUE))
plot(age.mids)

# literal, where MACi is mean age contrib for males age 20, bzw females age 20
for (yr in 1:length(years)){
    # males
    TFRim[,yr]      <-  rowSums(Fxm[[yr]], na.rm = TRUE) # male age in rows!
    TFRm[yr]        <-  sum(Fxm[[yr]], na.rm = TRUE)
    MACim[,yr]      <-  rowSums(Fxm[[yr]] * age.mids, na.rm = TRUE) / TFRim[, yr]
    MACm[yr]        <-  wmean(age.mids, rowSums(Fxm[[yr]], na.rm = TRUE))
    
    # females
    TFRif[,yr]      <-  colSums(Fxf[[yr]], na.rm = TRUE) # female age in columns!
    TFRf[yr]        <-  sum(Fxf[[yr]], na.rm = TRUE)
    MACif[,yr]      <-  colSums(t(t(Fxf[[yr]]) * age.mids), na.rm = TRUE) / TFRif[, yr]
    MACf[yr]        <-  wmean(age.mids, colSums(Fxf[[yr]], na.rm = TRUE))  
}
# MACi is pointless on second thought- the contribution of males aged 20 is 20.5 no matter how the females
# they reproduce with are aged.

# or do we wish to say for males age 20 what the contribution of females from each age was?
for (yr in 1:length(years)){
    # males
    TFRim[,yr]      <-  rowSums(Fxm[[yr]], na.rm = TRUE) # male age in rows!
    TFRm[yr]        <-  sum(Fxm[[yr]], na.rm = TRUE)
    MACim[,yr]      <-  rowSums(Fxm[[yr]] * age.mids, na.rm = TRUE) / TFRim[, yr]
    MACm[yr]        <-  wmean(age.mids, rowSums(Fxm[[yr]], na.rm = TRUE))
    
    # females
    TFRif[,yr]      <-  colSums(Fxf[[yr]], na.rm = TRUE) # female age in columns!
    TFRf[yr]        <-  sum(Fxf[[yr]], na.rm = TRUE)
    MACif[,yr]      <-  colSums(t(t(Fxf[[yr]]) * age.mids), na.rm = TRUE) / TFRif[, yr]
    MACf[yr]        <-  wmean(age.mids, colSums(Fxf[[yr]], na.rm = TRUE))  
}


# explore.
plot(years, MACm,type='l',col = "blue",ylim=c(25,32))
lines(years,MACf,col = "red")

# first differences, with male-female difference years marked
plot(years[-1], diff(MACm), type = 'l',col = "blue", ylim = c(-.15, .15))
lines(years[-1], diff(MACf), col = "red")
rects <- sign(diff(MACm)) != sign(diff(MACf))
abline(h=0)
rect(years[-1][rects] - .5, -.16,years[-1][rects] + .5, .16, col = "#CCCCCC40", border = NA)

# same thing but relative differences Rdiff:
plot(years[-1], diff(MACm) / ((MACm[-1] + MACm[-length(MACm)]) / 2), type = 'l',col = "blue", ylim = c(-.006, .006))
lines(years[-1], diff(MACf) / ((MACf[-1] + MACf[-length(MACf)]) / 2), col = "red")
abline(h=0)
rect(years[-1][rects] - .5, -.007, years[-1][rects] + .5, .007, col = "#CCCCCC40", border = NA)

# TFR plots
# good, crossover too
plot(years, TFRm, type='l', col = "blue", ylim = c(1.6, 2.8))
lines(years, TFRf, col = "red")

# first differences for TFR (oppisite sides of 0) - literally in childer born per female
plot(years[-1], diff(TFRm), type = 'l',col = "blue", ylim = c(-.3, .3))
lines(years[-1], diff(TFRf), col = "red")
rects <- sign(diff(TFRm)) != sign(diff(TFRf))
abline(h=0)
rect(years[-1][rects] - .5, -.35, years[-1][rects] + .5, .35, col = "#CCCCCC40", border = NA)

# same thing but relative differences Rdiff:
plot(years[-1], diff(TFRm) / ((TFRm[-1] + TFRm[-length(TFRm)]) / 2), type = 'l', col = "blue", ylim = c(-.12, .12))
lines(years[-1], diff(TFRf)/ ((TFRf[-1] + TFRf[-length(TFRf)]) / 2), col = "red")
abline(h = 0)
rect(years[-1][rects] - .5, -.13, years[-1][rects] + .5, .13, col = "#CCCCCC40", border = NA)

# now the decomposition
Dfertm  <- Dschedulem <- Dfertf  <- Dschedulef <- matrix(nrow = length(ages), ncol = length(years) - 1,
                                            dimnames = list(Age = ages, Year = years[-1]))
for (yr in 1:(length(years) - 1)){ # uses yr and yr + 1
    # males
    Dschedulem[, yr] <- (TFRim[, yr + 1] / TFRm[yr + 1] + TFRim[, yr] / TFRm[yr]) / 2 * (MACim[, yr + 1] - MACim[, yr])
    Dfertm[, yr]  <- (TFRim[, yr + 1] / TFRm[yr + 1] - TFRim[, yr] / TFRm[yr]) * ((MACim[, yr + 1] + MACim[, yr]) / 2)
    # females
    Dschedulef[, yr] <- (TFRif[, yr + 1] / TFRf[yr + 1] + TFRif[, yr] / TFRf[yr]) / 2 * (MACif[, yr + 1] - MACif[, yr])
    Dfertf[, yr]  <- (TFRif[, yr + 1] / TFRf[yr + 1] - TFRif[, yr] / TFRf[yr]) * ((MACif[, yr + 1] + MACif[, yr]) / 2)
}
sum( TFRim[, yr] / TFRm[yr])

plot(MACim[, yr + 1] - MACim[, yr])

range(c(Dschedulef, Dschedulem), na.rm = TRUE)
range(c(Dfertf, Dfertm), na.rm = TRUE)

brks <- pretty(rep(max(abs(range(c(Dschedulef, Dschedulem), na.rm = TRUE))),2) * c(-1,1), n = 50)
colfun <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"), space = "Lab")

dev.new(height = 5, width = 9)
par(mfrow = c(1, 2))
image(x = years + .5, y = ages + .5, t(Dschedulef), main = "Schedule Females", asp = 1, 
        xlim = range(years) + c(0, 1), ylim = range(ages) + c(0, 1),
        breaks = brks, col = rev(colfun(length(brks) - 1)))
fields::image.plot(x = years + .5, y = ages + .5, t(Dschedulem), main = "Schedule Males", asp = 1, 
        xlim = range(years) + c(0, 1), ylim = range(ages) + c(0, 1),
        breaks = brks, col = rev(colfun(length(brks) - 1)), zlim = range(brks))

par(mfrow = c(1, 2))
image(t(Dfertm), main = "Fert Females")
image(t(Dfertf), main = "Fert Females")

plot(1970:2009, colSums(Dweights, na.rm = TRUE),type = 'l', ylim = c(-.15,.15), col = gray(.5), lwd = 2)
lines(1970:2009, colSums(Dschedule, na.rm = TRUE), col = gray(.2),lwd=1.5, lty=2)

# so what does the MAC have to do with this anyway?


