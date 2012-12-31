
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
B                   <- local(get(load("/home/triffe/git/Dissertation/DISSERTATION/DATA/USbirths/USBxy10_65.Rdata")))

Bcube               <- abind::abind(B, along = 3)
array_names         <- dimnames(Bcube)
names(array_names)  <- c("Males", "Females", "Year")
dimnames(Bcube)     <- array_names
E                   <- local(get(load("/home/triffe/git/Dissertation/DISSERTATION/DATA/Exposures/USexp.Rdata")))

#dimnames(Bcube)
# Bcube: males in rows, females in columns:

# MAC is a weighted average of age where fx are the weights.
# as 
wmean <- function(x,w){
    sum(w * x) / sum(w)
}
wvar  <- function(x,w){
    sum(w * ((x -  wmean(x, w)) ^2)) / (sum(w) - 1)    
}
#bt <- B[[1]] 
#et <- E[E$Year == 1969 & E$Age %in% (10:65), "Male"]
#
#ft          <- rowSums(bt / et)
#TFRmt       <- sum(ft)
#age.mids    <- 10.5 : 65.5
#(MACmt       <- wmean(age.mids,ft))
#(VACmt       <- wvar(age.mids,ft))
#sqrt(VACmt)
#
#decompDeltaMean <- function(x, w1, w2){
#    c(w=(sum(w1*x) + sum(w2*x)) / 2 * (sum(w2) - sum(w1)),
#    f=(sum(w2*x) - sum(w1*x))  * ((sum(w2) + sum(w1)) / 2))
#}
#wmean(age.mids,rowSums(B[["1970"]]  / E[E$Year == 1970 & E$Age %in% (10:65), "Male"])) - MACmt
#decompDeltaMean(age.mids, w1 = rowSums(B[["1969"]]  / E[E$Year == 1969 & E$Age %in% (10:65), "Male"]),
#        w2 = rowSums(B[["1970"]]  / E[E$Year == 1970 & E$Age %in% (10:65), "Male"]))
#b1 <- B[["1969"]]
#b2 <- B[["1970"]]
#e1 <- E[E$Year == 1969 & E$Age %in% (10:65), "Male"]
#e2 <- E[E$Year == 1970 & E$Age %in% (10:65), "Male"]
#
#TFR1i <- colSums(b1 / e1)
#TFR2i <- colSums(b2 / e2)
#TFR1  <- sum(TFR1i)
#TFR2  <- sum(TFR2i)
#
#MAC1i <- colSums((b1 / e1) * age.mids) / colSums(b1 / e1)
#MAC2i <- colSums((b2 / e2) * age.mids) / colSums(b2 / e2)
#MAC1  <- wmean(age.mids,rowSums(b1) / e1)
#MAC2  <- wmean(age.mids,rowSums(b2) / e2)
#
#plot(10:65, MAC1i, type = 'l')
#lines(10:65, MAC2i, col = "red")
#Dschedule <- (TFR2i / TFR2 + TFR1i / TFR1) / 2 * (MAC2i - MAC1i)
#Dweights  <- (TFR2i / TFR2 - TFR1i / TFR1)  * ((MAC2i + MAC1i) / 2)
#
#sum(Dschedule, na.rm = TRUE) + sum(Dweights, na.rm = TRUE)
#MAC2-MAC1 # equal, nice
#plot(10:65, Dweights, type = 'l')
#lines(10:65, Dschedule, col = "red")

Fxm <- Fxf <- list()
ages <- 10:65
years <- 1969:2009
for (yr in 1:length(years)){
    Fxm[[yr]] <- B[[yr]] / E[E$Year == years[yr] & E$Age %in% ages, "Male"]
    Fxf[[yr]] <- t(t(B[[yr]]) / E[E$Year == years[yr] & E$Age %in% ages, "Female"])
}

MACm    <- TFRm     <- MACf    <- TFRf      <- vector(length = length(years))
MACim   <- TFRim    <- MACif   <- TFRif     <- matrix(ncol = length(years), nrow = length(ages), 
                                                    dimnames = list(Age = ages, Year = years))
age.mids    <- ages + .5
for (yr in 1:length(years)){
    # males
    TFRim[,yr]      <-  colSums(Fxm[[yr]], na.rm = TRUE) 
    TFRm[yr]        <-  sum(Fxm[[yr]], na.rm = TRUE)
    MACim[,yr]      <-  colSums(Fxm[[yr]] * age.mids, na.rm = TRUE) / TFRim[, yr]
    MACm[yr]        <-  wmean(age.mids, rowSums(Fxm[[yr]], na.rm = TRUE))
    
    # females
    TFRif[,yr]      <-  colSums(t(Fxf[[yr]]), na.rm = TRUE) 
    TFRf[yr]        <-  sum(t(Fxf[[yr]]), na.rm = TRUE)
    MACif[,yr]      <-  colSums(t(Fxf[[yr]]) * age.mids, na.rm = TRUE) / TFRif[, yr]
    MACf[yr]        <-  wmean(age.mids, rowSums(t(Fxf[[yr]]), na.rm = TRUE))  
}

# 
plot(years, MACm,type='l',col = "blue",ylim=c(25,32))
lines(years,MACf,col = "red")

# good, crossover too
plot(years, TFRm, type='l', col = "blue", ylim = c(1.5,2.8))
lines(years, TFRf, col = "red")

Dfertm  <- Dschedulem <- Dfertf  <- Dschedulef <- matrix(nrow = length(ages), ncol = length(years) - 1,
                                            dimnames = list(Age = ages, Year = years[-1]))
for (yr in 1:(length(years)-1)){
    Dschedulem[, yr] <- (TFRim[, yr+1] / TFRm[yr+1] + TFRim[,yr] / TFRm[yr]) / 2 * (MACim[, yr+1] - MACim[, yr])
    Dfertm[, yr]  <- (TFRim[, yr+1] / TFRm[yr+1] - TFRim[,yr] / TFRm[yr])  * ((MACim[, yr+1] + MACim[,yr]) / 2)
    Dschedulef[, yr] <- (TFRif[, yr+1] / TFRf[yr+1] + TFRif[,yr] / TFRf[yr]) / 2 * (MACif[, yr+1] - MACif[, yr])
    Dfertf[, yr]  <- (TFRif[, yr+1] / TFRf[yr+1] - TFRif[,yr] / TFRf[yr])  * ((MACif[, yr+1] + MACif[,yr]) / 2)
}

par(mfrow = c(1, 2))
image(t(Dschedulef), main = "Schedule Females")
image(t(Dschedulem), main = "Schedule Males")

par(mfrow = c(1, 2))
image(t(Dfertm), main = "Fert Females")
image(t(Dfertf), main = "Fert Females")

plot(1970:2009, colSums(Dweights, na.rm = TRUE),type = 'l', ylim = c(-.15,.15), col = gray(.5), lwd = 2)
lines(1970:2009, colSums(Dschedule, na.rm = TRUE), col = gray(.2),lwd=1.5, lty=2)

# so what does the MAC have to do with this anyway?


