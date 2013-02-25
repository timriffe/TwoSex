source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
# BxES is 0:110, years 1975:2009
BxUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
BxES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

#BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf10_65.Rdata")))
#BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf10_65.Rdata")))

# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
dxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfES.Rdata"))) 

dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)

yearsUS <- 1969:2009
yearsES <- 1975:2009

ages <- 0:110
# example for 1975 US:

Mex75US <- rowSums(ExpectedDx(Px = with(ExUS, Male[Year == 1975]), dx = dxmUS[,"1975"]))
Fex75US <- rowSums(ExpectedDx(Px = with(ExUS, Female[Year == 1975]), dx = dxfUS[,"1975"]))


ExBxy <- ExpectedDxMxFmatrix( BxUS[["1975"]], dxmUS[,"1975"], dxfUS[,"1975"])


Fxym <- ExBxy %row% Mex75US
Fxyf <- ExBxy %col% Fex75US

par(mfrow = c(1,2))
image(x = ages + .5, y = ages + .5, t(Fxym), ylim = c(0,111),xlim = c(0,111))
image(x = ages + .5, y = ages + .5, t(Fxyf), ylim = c(0,111),xlim = c(0,111))

# observed vs expected:

expected    <- outer(rowSums(ExBxy), colSums(ExBxy), "*") / sum(ExBxy)

#image(x = ages + .5, y = ages + .5, t(expected), ylim = c(0,111),xlim = c(0,111))
#image(x = ages + .5, y = ages + .5, t( expected / sum(expected)) - t(ExBxy / sum(ExBxy)) , ylim = c(0,111),xlim = c(0,111))
# total variation distance:
(ToalVar <- sum(abs(ExBxy / sum(ExBxy) - expected / sum(expected))) / 2 )

# overlap 
(coefOverlap     <- sum(pmin(ExBxy / sum(ExBxy), expected / sum(expected))))

# much more perfect distribution.


















