setwd("/home/triffe/git/DISS/")
# scratch code, unused officially
# Author: triffe
###############################################################################

source("R/MeanFunctions.R")
#stolarsky.mean(x=c(50,100),r=.5,s=.5)
#stolarsky.mean(x=c(50,100),r=.01,s=5)


B <- local(get(load("DATA/ESbirths/ESbirths.Rdata")))
E <- local(get(load("DATA/Exposures/ESexp.Rdata")))

years <- 1975:2009
years.char <- as.character(years)
Fxlist <- list()
library(reshape2)

for (yr in 1:length(years)){
    Ef              <- with(E, Female[Year == years[yr]])
    Em              <- with(E, Male[Year == years[yr]])
    Bmat4males      <- acast(B[[years.char[yr]]], FAGE~MAGE, sum, value.var = "Births")
    Bmat4females    <- acast(B[[years.char[yr]]], MAGE~FAGE, sum, value.var = "Births")
    
    Fxmales         <- Bmat4males / Em
    Fxfemales       <- Bmat4females / Ef
    
    Fxlist[[years.char[yr]]] <- list(Males_mf =  Fxmales, 
            Females_fm =  Fxfemales, 
            Fxm = rowSums(Fxmales),
            Fxf = rowSums(Fxfemales))
}

# so the question is what kind of mean to take, given the rates from year t-i, given the exposures of year t
# start with i = 1


i <- 1
Bpred1 <- list()
for (yri in (1+i):length(years)){
    Fxm.lag          <- Fxlist[[years.char[yri - i ]]]$Fxm 
    Fxf.lag          <- Fxlist[[years.char[yri - i ]]]$Fxf 
    Eft              <- with(E, Female[Year == years[yri]])
    Emt              <- with(E, Male[Year == years[yri]])
    Bm.pred          <- Fxm.lag * Emt
    Bf.pred          <- Fxf.lag * Eft
    Bpred1[[years.char[yri]]] <- list(Bm.pred = Bm.pred, Bf.pred = Bf.pred)
}
#plot(0:110, Bpred1[[1]]$Bm.pred, type ='l', col = "blue")
#lines(0:110, Bpred1[[1]]$Bf.pred, col = "red")

Bm.pred <- do.call(cbind,lapply(Bpred1, function(x){
            x$Bm.pred
        }))
Bf.pred <- do.call(cbind,lapply(Bpred1, function(x){
                    x$Bf.pred
                }))
plot(1976:2009, colSums(Bm.pred, na.rm = TRUE) - colSums(Bf.pred, na.rm = TRUE), 
        type = 'l')


TFRf <- colSums(do.call(cbind,lapply(Fxlist, function(x){
                    x$Fxf
                })), na.rm = TRUE)
TFRm <- colSums(do.call(cbind,lapply(Fxlist, function(x){
                            x$Fxm
                        })), na.rm = TRUE)





plot(years, TFRf, type = 'l', col = "red")
lines(years, TFRm, col = "blue")
#par(new = TRUE)
#plot(years, c(NA,colSums(Bm.pred, na.rm = TRUE) - 
#                colSums(Bf.pred, na.rm = TRUE)), type = 'l')
# observed total births:
Bobs        <- unlist(lapply(B, function(x){
                        sum(x$Births)
                    }))[years.char]
Bmp         <- colSums(Bm.pred, na.rm = TRUE)
Bfp         <- colSums(Bf.pred, na.rm = TRUE)
chg.scale   <- Bobs[2:length(Bobs)] / Bobs[1:(length(Bobs)-1)]
chg.rate.f    <- TFRf[2:length(TFRf)] / TFRf[1:(length(TFRf)-1)]
chg.rate.m    <- TFRm[2:length(TFRm)] / TFRm[1:(length(TFRm)-1)]
# TODO: I think the thing to do is scale rates according to the age-specific change in the rate..
#       but then I'm not sure how to make sense of that, since it's sex-specific. Further thinking

plot(1975:2009, Bobs, type = 'l', main = "Difference not fully corrected")
lines(1976:2009, Bmp * chg.scale * chg.rate.m, col = "blue")
lines(1976:2009, Bfp * chg.scale * chg.rate.f, col = "red")

plot(1976:2009, chg.rate.m, type ='l', col = "blue")
lines(1976:2009, chg.rate.f, col = "red")

# PROBLEM, observed Births don't fall between male and female predicted.
colfun <- grDevices:::colorRampPalette(RColorBrewer:::brewer.pal(9,"RdBu"), space = "Lab")

plot(1975:2009, Bobs, type = 'l')
lines(1976:2009, Bmp + diff(Bobs), col = "blue")
lines(1976:2009, Bfp + diff(Bobs), col = "red")

lines(1975:2009, cumsum(c(Bobs[1],  diff(Bobs))), col = "green")


# -----------------------------------------------------------------------------------
# take a look for USA. problem with ES- overall trend drowns out sex differences...

B <- local(get(load("DATA/USbirths/USbirths.Rdata")))
E <- local(get(load("DATA/Exposures/USexp.Rdata")))
years <- 1969:2009

years.char <- as.character(years)
Fxlist <- list()

for (yr in 1:length(years)){
    Ef              <- with(E, Female[Year == years[yr]])
    Em              <- with(E, Male[Year == years[yr]])
    Bmat4males      <- B[[years.char[yr]]]
    Bmat4females    <- t(B[[years.char[yr]]])
    
    Fxmales         <- Bmat4males / Em
    Fxfemales       <- Bmat4females / Ef
    
    Fxlist[[years.char[yr]]] <- list(Males_mf =  Fxmales, 
            Females_fm =  Fxfemales, 
            Fxm = rowSums(Fxmales, na.rm = TRUE),
            Fxf = rowSums(Fxfemales, na.rm = TRUE))
}

# so the question is what kind of mean to take, given the rates from year t-i, given the exposures of year t
# start with i = 1


i <- 1
Bpred1 <- list()
for (yri in (1+i):length(years)){
    Fxm.lag          <- Fxlist[[years.char[yri - i ]]]$Fxm
    Fxf.lag          <- Fxlist[[years.char[yri - i ]]]$Fxf
    Eft              <- with(E, Female[Year == years[yri]])
    Emt              <- with(E, Male[Year == years[yri]])
    Bm.pred          <- Fxm.lag * Emt
    Bf.pred          <- Fxf.lag * Eft
    Bpred1[[years.char[yri]]] <- list(Bm.pred = Bm.pred, Bf.pred = Bf.pred)
}


Bm.pred <- do.call(cbind,lapply(Bpred1, function(x){
                    x$Bm.pred
                }))
Bf.pred <- do.call(cbind,lapply(Bpred1, function(x){
                    x$Bf.pred
                }))

plot(1970:2009, colSums(Bm.pred, na.rm = TRUE) - 
                colSums(Bf.pred, na.rm = TRUE), type = 'l')


TFRf <- colSums(do.call(cbind,lapply(Fxlist, function(x){
                            x$Fxf
                        })), na.rm = TRUE)
TFRm <- colSums(do.call(cbind,lapply(Fxlist, function(x){
                            x$Fxm
                        })), na.rm = TRUE)

plot(years, TFRf, type = 'l', col = "red")
lines(x=years, y=TFRm, col = "blue")

par(new = TRUE)
plot(years, c(NA,colSums(Bm.pred, na.rm = TRUE) - 
                        colSums(Bf.pred, na.rm = TRUE)), type = 'l')
# observed total births:
Bobs <- unlist(lapply(B, function(x){
                    sum(x$Births)
                }))
plot(1975:2010, Bobs, type = 'l')
lines(1975:2008, colSums(Bm.pred, na.rm = TRUE), col = "blue")
lines(1975:2008, colSums(Bf.pred, na.rm = TRUE), col = "red")

# 1979, 2004, 2005, incomplete USA births?

# ------------------------------------------------------------------------------------

source("R/MeanFunctions.R")
Bxy <- local(get(load("DATA/ESbirths/ESBxy.Rdata")))
E <- local(get(load("DATA/Exposures/ESexp.Rdata")))

years <- 1975:2009
yr <- 1

# year 1
Ef1             <- with(E, Female[Year == years[yr]])
Em1             <- with(E, Male[Year == years[yr]])
Bm1             <- Bxy[[yr]]
Bf1             <- t(Bxy[[yr]])
Fxm1            <- Bm1 / Em1 # note, males and females oriented differently! can't compare without t()!!!
Fxf1            <- Bf1 / Ef1
# year 2
Ef2             <- with(E, Female[Year == years[yr+1]])
Em2             <- with(E, Male[Year == years[yr+1]])
Bm2             <- Bxy[[yr+1]]
Bf2             <- t(Bxy[[yr+1]])
Fxm2            <- Bm2 / Em2
Fxf2            <- Bf2 / Ef2
# hypothetical Fx2 - gives same total births under prior year's Fxy distribution (rescaled)
Fxm2.hat        <- (Bm1 * (sum(Bxy[[yr+1]]) / sum(Bxy[[yr]]))) / Em2
Fxf2.hat        <- (Bf1 * (sum(Bxy[[yr+1]]) / sum(Bxy[[yr]]))) / Ef2


sum(Fxm2.hat * Em2)
sum(Fxf2.hat * Ef2)
sum(Fxm2 * Em2)
sum(Fxf2 * Ef2)
sum(abs(Fxm2.hat - Fxm2))
sum(Fxm2)
sum(Fxf2.hat)
sum(Fxf2)




sum(Bm1);sum(Bf1)
sum(Bm2);sum(Bf2)

DecomposeBirthsChange <- function(Fx1, Fx2, Ex1, Ex2, return.all = FALSE){
            RateComponent       <- (Fx2 - Fx1) * (Ex2 + Ex1) / 2
            ExposureComponent   <- (Fx2 + Fx1) / 2 * (Ex2 - Ex1)
            if (return.all){
                return(list(FxCmat = RateComponent,           # check out plots of these matrices to see how they differ
                                ExCmat = ExposureComponent, 
                                FxC = sum(RateComponent, na.rm = TRUE),      # these two should add to the change in B
                                ExC = sum(ExposureComponent, na.rm = TRUE)))
            } 
            out <- c(sum(RateComponent, na.rm = TRUE), sum(ExposureComponent, na.rm = TRUE))
            names(out) <- c("FxC","ExC")
            out
        }
        

# Dec 16th, 2012:
# ideas: 1) rescale rates in year t - 1 so that the total births in year t-1 equals
#        total births in year t. 
#        2) rescale male F(t - 1) so that when applied to E(t), you get B(t), for males and females
#        separately. What do we see?
Bxy <- local(get(load("DATA/ESbirths/ESBxy.Rdata")))
E   <- local(get(load("DATA/Exposures/ESexp.Rdata")))
        
is.list(E) 
image(t(Bxy[[1]])[12:65,12:65])
unique(Bxy[[1]][,14]) # columns females
unique(Bxy[[1]][14,]) # males rows

# females:
Fxf1   <- t(Bxy[[1]]) / with(E, Female[Year == 1975])
Bxyf.1  <- Fxf1 *  with(E, Female[Year == 1976])
Fxf.2  <- Fxf1 * (sum(Bxy[[2]])/sum(Bxyf.1)) # females t-1 rescaled
sum(Fxf.2 * with(E, Female[Year == 1976])) == sum(Bxy[[2]])
# males:
Fxm1    <- Bxy[[1]] / with(E, Male[Year == 1975])
Bxym.1  <- Fxm1 *  with(E, Male[Year == 1976])
Fxm.2   <- Fxm1 * (sum(Bxy[[2]])/sum(Bxym.1)) # females t-1 rescaled
sum(Fxm.2 * with(E, Male[Year == 1976])) == sum(Bxy[[2]])
sum(Fxm1)
sum(Fxm.2)

sum(Fxf1)
sum(Fxf.2)
plot(c(1,1),c(sum(Fxm1),sum(Fxf1)), col = c("blue", "red"), pch = "1",ylim = c(2.7,2.9))
        points(c(1,1), c(sum(Fxm.2),sum(Fxf.2)), col = c("blue","red"), pch = "2")
        
sum(Fxm.2) / sum(Fxm1)
sum(Fxf.2) / sum(Fxf1) # females increased more than males?

# compare to actual year t rates?
Fxm2   <- Bxy[[2]] / with(E, Male[Year == 1976])
Fxf2   <- t(Bxy[[2]]) / with(E, Female[Year == 1976])

sum(Fxm2) / sum(Fxm.2)
sum(Fxf2) / sum(Fxf.2) 

sum(abs(Fxm2 - Fxm.2))
sum(abs(Fxf2 - Fxf.2))

# OK, step 1 
Fxf1B    <- Fxm1B  <- Fxf1A    <- Fxm1A <- list()
years    <- 1975:2009
for (i in 1:(length(years) - 1)){
    # males
    Em1             <- with(E, Male[Year == years[i]])
    Em2             <- with(E, Male[Year == years[i + 1]])
    Bm1             <- Bxy[[i]] 
    Bm2             <- Bxy[[i+1]]
    Fxm1            <- Bm1 / Em1
    Fxm1B[[as.character(years[i+1])]]    <- Fxm1 * sum(Bm2) / sum(Fxm1 *  Em2, na.rm = TRUE)
    Fxm1A[[as.character(years[i+1])]]    <- Bm2 / Em2
    # females
    Ef1             <- with(E, Female[Year == years[i]])
    Ef2             <- with(E, Female[Year == years[i + 1]])
    Bf1             <- t(Bxy[[i]] )
    Bf2             <- t(Bxy[[i+1]])
    Fxf1            <- Bf1 / Ef1
    Fxf1B[[as.character(years[i+1])]]    <- Fxf1 * sum(Bf2) / sum(Fxf1 *  Ef2, na.rm = TRUE)
    Fxf1A[[as.character(years[i+1])]]    <- Bf2 / Ef2
}

plot(1976:2009, unlist(lapply(Fxm1A,sum,na.rm=TRUE)), 
        type ='l', col = "blue", ylim = c(1.2, 2.9), main = "WTF does this mean?")
lines(1976:2009, unlist(lapply(Fxf1A,sum,na.rm=TRUE)), col = "red")
lines(1976:2009,unlist(lapply(Fxf1B,sum,na.rm=TRUE)), col = "red", lty = 2)
lines(1976:2009,unlist(lapply(Fxm1B,sum,na.rm=TRUE)), col = "blue", lty = 2)

FMAD <- mapply(Fxf1A,Fxf1B,FUN = function(x,y){sum(abs(x-y), na.rm = TRUE)})
MMAD <- mapply(Fxm1A,Fxm1B,FUN = function(x,y){sum(abs(x-y), na.rm = TRUE)})

plot(1976:2009, FMAD, type = 'l', col = "red", main = "sum absolute deviation")
lines(1976:2009, MMAD, col = "blue")

plot(1976:2009, FMAD / unlist(lapply(Fxf1A, sum, na.rm = TRUE)), type = 'l', 
        col = "red", main = "sum relative deviations")
lines(1976:2009, MMAD / unlist(lapply(Fxm1A, sum, na.rm = TRUE)), col = "blue")




