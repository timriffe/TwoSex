
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxES is 0:110, years 1975:2009
BxUS  <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
BxES  <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

yearsUS <- 1969:2009
yearsES <- 1975:2009

BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))
names(BxymfES) <- yearsES
names(BxymfUS) <- yearsUS

# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
dxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfES.Rdata"))) 

dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)

#---------------------------------------------
Bxy <- BxUS[["1969"]]
Exm1 <- with(ExUS, Male[Year == 1969])
Exf1 <- with(ExUS, Female[Year == 1969])
Exm2 <- with(ExUS, Male[Year == 1970])
Exf2 <- with(ExUS, Female[Year == 1970])
# eq 9
# estimator? from eq 15
m0 <- Exm1 - rowSums(Bxy)
f0 <- Exf1 - colSums(Bxy)
Pij <- log(Bxy / sqrt(outer(m0, f0, "*")))



function(Bxy, Exm1, Exf1, Exm2, Exf2){
    m0  <- Exm1 - rowSums(Bxy)
    f0  <- Exf1 - colSums(Bxy)
   
    ((pijt - pi0 - p0j) / 2) - log(Bxy / sqrt(outer(m0, f0, "*")))    
}





