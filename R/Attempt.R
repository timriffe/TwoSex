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

ExBxyAallUS <- lapply(as.character(yearsUS), function(yr, .BxUS, .dxmUS, .dxfUS){
            ExpectedDxMxFmatrix( .BxUS[[yr]], .dxmUS[,yr], .dxfUS[,yr])
        }, .BxUS = BxUS, .dxmUS = dxmUS, .dxfUS = dxfUS)
ExBxyAallES <- lapply(as.character(yearsES), function(yr, .BxES, .dxmES, .dxfES){
            ExpectedDxMxFmatrix( .BxES[[yr]], .dxmES[,yr], .dxfES[,yr])
        }, .BxES = BxES, .dxmES = dxmES, .dxfES = dxfES)


Ex <- with(ExUS, Male[Year == 2000])
A <- rowSums(ExpectedDx(Ex, dx = dxmUS[,"2000"]))
Ex <- with(ExUS, Male[Year == 2000])
dxi <- dxmUS[,"2000"]
N <- length(dxi)
EDx2      <-EDx      <- matrix(0, nrow = N, ncol = N, dimnames = list(Ex = 0:(N-1), Age =  0:(N-1)))
# Population age loop
for (i in 1:N){
    # distribute each age of Populatin over death times
    EDx[1:length(dxi), i]    <- dxi
    # remove firs element and rescale
    dxi                      <- dxi[2:length(dxi)] / sum(dxi[2:length(dxi)], na.rm = TRUE)
}
dxi <- dxmUS[,"2000"]
for (i in 1:N){
    # distribute each age of Populatin over death times
    EDx2[1:length(dxi), i]    <- rev(dxi)
    # remove firs element and rescale
    dxi                       <- dxi[2:length(dxi)] / sum(dxi[2:length(dxi)], na.rm = TRUE)
}
#EDx2 <- EDx2 / rowSums(EDx2, na.rm = TRUE)

colSums(EDx2 * A)

A <- ExpectedDx(Ex, dx = dxmUS[,"2000"])
RedistMat <- A / rowSums(A)

age<- seq(.5,110.5)

wmean(age, RedistMat[1,])

plot(rowMeans(RedistMat %row% 1/seq(.5,110.5)))

colSums(rowSums(A) * RedistMat) - Ex

plot(RedistMat[1,])
image(RedistMat)

rowSums(ExpectedDx(Ex, dx = dxmUS[,"2000"]))

barplot(A)




