

# well-drawn diagram required for exPopulation Pyramid.

# highlight an age, say 0-5?
# take 1975 US

source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxES is 0:110, years 1975:2009
BxUS  <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))

# exposures, as such, straight from HMD, all ages 0-110, long form
ExUS  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
# get Lx estimates for R0, r
dxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata"))) 

dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)

ExmUS1975 <- rowSums(ExpectedDx( with(ExUS, Male[Year == 1975]), dxmUS[, "1975"]))
BxmUS1975 <- rowSums(ExpectedDx( rowSums(BxUS[["1975"]]), dxmUS[, "1975"]))
ExfUS1975 <- rowSums(ExpectedDx( with(ExUS, Female[Year == 1975]), dxfUS[, "1975"]))
BxfUS1975 <- rowSums(ExpectedDx( colSums(BxUS[["1975"]]), dxfUS[, "1975"]))

FxmUS1975 <- BxmUS1975 / ExmUS1975
FxfUS1975 <- BxfUS1975 / ExfUS1975

































