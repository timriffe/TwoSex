setwd("/home/triffe/git/DISS/")
source("R/UtilityFunctions.R")
source("R/MeanFunctions.R")

yearsUS <- 1969:2009
yearsES <- 1975:2009

# Bxy totals  (not used()
BxUS  <- local(get(load("Data/USbirths/USBxy0_110.Rdata")))
BxES  <- local(get(load("Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

# repeated below for better SRB assumptions
BxymfES <- local(get(load("Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("Data/USbirths/USBxymf0_110.Rdata")))

# exposures, as such, straight from HMD, all ages 0-110, long form
ExUS  <- local(get(load("Data/Exposures/USexp.Rdata")))
ExES  <- local(get(load("Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
dxmUS <- local(get(load("Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("Data/HMD_dx/dxfES.Rdata"))) 
# make sum to 1
dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)

# Harmonic mean function
HM <- function(x,y){
    Minf0(Mna0((2 * x * y) / (x + y)))
}
# 1) make a harmonic mean rate a x a'
yr <- "1975"
Bxy <- BxymfUS[[yr]][[1]] + BxymfUS[[yr]][[2]]
Bxyex <- ExpectedDxMxFmatrix(Bxy, dxm = dxmUS[,yr], dxf = dxfUS[,yr])
PMa <- with(ExUS, Male[Year == as.integer(yr)])
PFa <- with(ExUS, Female[Year == as.integer(yr)])
MexM <- ExpectedDx(PMa, dxmUS[,yr])
FexM <- ExpectedDx(PFa, dxfUS[,yr])

ExpHM <- outer(rowSums(MexM),rowSums(FexM),HM)

mHM <-  Minf0(Mna0(Bxyex / ExpHM))
sum(mHM)

# 1) increase age 25 males by 2x

PMa2 <- PMa
PMa2[26] <- PMa2[26] * 2
MexM2 <- ExpectedDx(PMa2, dxmUS[,yr])
BexM2 <- ExpectedDx(rowSums(Bxy), dxmUS[,yr])

MexM3 <- MexM2
MexM3[,rowSums(Bxy) == 0] <- 0

# new harmonic avg:
ExpHM2 <- outer(rowSums(MexM2),rowSums(FexM),HM)
# new predicted births:
Bxyex2 <- mHM * ExpHM2
# problem: these are cross-classified by remaining years and not age.
# need to get back to male and female age marginals.

Male_eyMarg <- rowSums(Bxyex2)
# smooth as expected
# plot(Male_eyMarg)
MexM2


#  now do it the hard way...
Exm2array <- array(0,dim=c(111,111,111))
for (i in 1:111){
    Exm2array[i,,] <- outer(MexM2[,i],rowSums(FexM),HM) 
}

#plot(rowSums(MexM2),type='l')
#lines(rowSums(MexM),col = "blue")


plot(colSums(Minf0(Mna0(rowSums(Bxyex2) / MexM2)), na.rm=TRUE))


#fields::image.plot(Bxyex2 - Bxyex)

# collapse Births to males, and distribute back out to ages and RYL
Bxm1 <- rowSums(Bxyex) * Minf0(Mna0(MexM / rowSums(MexM)))
Bredm <- ExpectedDx(rowSums(Bxy), dxmUS[,yr])

Bxm2 <- rowSums(Bxyex2) * Minf0(Mna0(MexM2 / rowSums(MexM2)))

# new age dist
colSums(Bxm2)
plot(0:110, colSums(Bxm2), type= 'l')
lines(0:110, colSums(Bxm1), col = "red", lty = 2)
ratio <- colSums(Bxm2) / colSums(Bxm1)
ratio[26] <- NA
plot(ratio[15:40],type='l') # hm, need to work on proper redist to ages.

# likewise, check before and after for females
Bxf1 <- colSums(Bxyex) * Minf0(Mna0(t(t(FexM) / colSums(FexM))))
Bxf2 <- colSums(Bxyex2) * Minf0(Mna0(t(t(FexM) / colSums(FexM))))

plot(0:110, colSums(Bxf2), type= 'l')
lines(0:110, colSums(Bxf1), col = "red", lty = 2)

plot(0:110, colSums(Bxf2) / colSums(Bxf1))

plot(0:110, colSums(Bxm2) / colSums(Bxm1))

sum(colSums(Bxm2) - colSums(Bxm1))

# 2nd try - only repro ages?

PMa3 <- PMa
PMa3[rowSums(Bxy)==0] <- 0
PFa3 <- PFa
PFa3[colSums(Bxy)==0] <- 0

MexM3 <- ExpectedDx(PMa3, dxmUS[,yr])
FexM3 <- ExpectedDx(PFa3, dxfUS[,yr])

# colSums of these will give back age...
FemaleRatesA <-MaleRatesA <- array(dim=c(111,111,111))
for (i in 1:111){
    MaleRatesA[,,i] <- ExpectedDx(Bxy[,i] / PMa, dxmUS[,yr])
    FemaleRatesA[,,i] <- ExpectedDx(Bxy[i,] / PFa, dxfUS[,yr])
}
# each stack is key

rowSums(Bxy) / PMa

image(MaleRates)
ExpHM3 <- outer(rowSums(MexM3),rowSums(FexM3),HM)

mHM3 <-  Minf0(Mna0(Bxyex / ExpHM3))

# now chg exp
PMa4 <- PMa3
PMa4[26] <- PMa4[26] * 2
MexM4 <- ExpectedDx(PMa4, dxmUS[,yr])
BexM4 <- ExpectedDx(rowSums(Bxy), dxmUS[,yr])
# new harmonic avg:
ExpHM4 <- outer(rowSums(MexM4),rowSums(FexM3),HM)
# new predicted births:
Bxyex4 <- mHM3 * ExpHM4

sum(Bxyex4)
plot(colSums(Mna0(Minf0(rowSums(Bxyex4) / MexM4)))) # nope

sum(Mna0(Minf0(BexM4 / MexM4)))

sum(Mna0(Minf0(colSums(BexM4) / colSums(MexM4))))

seq(-100,100,by=10)
p <- seq(-1,1,length.out=101)
plot(p,stolarsky.mean.v(0,4,p),type = 'l')
ls()
stolarsky.mean.v(0,4,0)


