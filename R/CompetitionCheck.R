source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

yearsUS <- 1969:2009
yearsES <- 1975:2009

# Bxy totals  (not used()
BxUS  <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
BxES  <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

# repeated below for better SRB assumptions
BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))

# exposures, as such, straight from HMD, all ages 0-110, long form
ExUS  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
dxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfES.Rdata"))) 
# make sum to 1
dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)


# 1) make a harmonic mean rate a x a'
yr <- "1975"
Bxy <- BxymfUS[[yr]][[1]] + BxymfUS[[yr]][[2]]
Bxyex <- ExpectedDxMxFmatrix(Bxy, dxm = dxmUS[,yr], dxf = dxfUS[,yr])
PMa <- with(ExUS, Male[Year == as.integer(yr)])
PFa <- with(ExUS, Female[Year == as.integer(yr)])
MexM <- ExpectedDx(PMa, dxmUS[,yr])
FexM <- ExpectedDx(PFa, dxfUS[,yr])
HM <- function(x,y){
   Minf0(Mna0((2 * x * y) / (x + y)))
}
ExpHM <- outer(rowSums(MexM),rowSums(FexM),HM)

mHM <-  Minf0(Mna0(Bxyex / ExpHM))
sum(mHM)

# 1) increase age 25 males by 2x

PMa2 <- PMa
PMa2[26] <- PMa2[26] * 2
MexM2 <- ExpectedDx(PMa2, dxmUS[,yr])


ExpHM2 <- outer(rowSums(MexM2),rowSums(FexM),HM)
#plot(rowSums(MexM2),type='l')
#lines(rowSums(MexM),col = "blue")
Bxyex2 <- mHM * ExpHM2

fields::image.plot(Bxyex2 - Bxyex)

# collapse Births to males, and distribute back out to ages and RYL
Bxm1 <- rowSums(Bxyex) * Minf0(Mna0(MexM / rowSums(MexM)))
Bredm <- ExpectedDx(rowSums(Bxy), dxmUS[,yr])

Bxm2 <- rowSums(Bxyex2) * Minf0(Mna0(MexM2 / rowSums(MexM2)))

# new age dist
colSums(Bxm2)
plot(0:110, colSums(Bxm2), type= 'l')
lines(0:110, colSums(Bxm1), col = "red", lty = 2)

# likewise, check before and after for females
Bxf1 <- colSums(Bxyex) * Minf0(Mna0(t(t(FexM) / colSums(FexM))))
Bxf2 <- colSums(Bxyex2) * Minf0(Mna0(t(t(FexM) / colSums(FexM))))

plot(0:110, colSums(Bxf2), type= 'l')
lines(0:110, colSums(Bxf1), col = "red", lty = 2)

plot(0:110, colSums(Bxf2) / colSums(Bxf1))

plot(0:110, colSums(Bxm2) / colSums(Bxm1))

sum(colSums(Bxm2) - colSums(Bxm1))









