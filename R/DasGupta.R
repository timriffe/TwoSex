# preamble:
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
# BxES is 0:110, years 1975:2009
BxUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy10_65.Rdata")))
BxES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65
# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))

#--------------------------
# Das Gupta 1972:








# implementing Das Gupta 1978:


Bxy         <- BxUS[["1970"]]

Pxm <- with(ExUS, Male[Age >= 10 & Age <= 65 & Year == 1970])
Pxf <- with(ExUS, Female[Age >= 10 & Age <= 65 & Year == 1970])        

Bxy <- BxUS[["1970"]]
# eq 5.4
Uxy <- Mna0(Bxy %row% rowSums(Bxy))
Vxy <- Mna0(Bxy %col% colSums(Bxy))

Uxy <- Uxy / sum(Uxy)
Vxy <- Vxy / sum(Vxy)
thetam <- 1.05 / 2.05
thetaf <- 1 / 2.05

Mxy <- Uxy %row% (Pxm / 1)
Fxy <- Uxy %col% (Pxf / 1)


mxy <- Mna0(Bxy / (Uxy %row% (Pxm / 1) + Vxy %col% (Pxf / 1)))

mx <- rowSums(Bxy) / Pxm
my <- colSums(Bxy) / Pxf
mxm  <- thetam * mx
myf  <- thetaf * my
mxy <- (thetam / mxm + thetaf / myf) ^ -1
sum(mxy) 
sum(mxy)



plot(10:65, mxy, col = "red", type = "l")
polygon(c(10:65, 65:10), c(mx, rev(my)), col = "#55555550")

Bt <- sum((thetam * Uxy %row% (Pxm/1) + thetaf * Vxy %col% (Pxf/1)) * )

sum(apply(cbind(mxm / (2 * thetam), myf / (2 * thetaf)), 1, harmonic.mean))
sum(apply(cbind(mxm / (thetam), myf / ( thetaf)), 1, harmonic.mean))
sum(apply(cbind(mx, my), 1, harmonic.mean))

(TFRm <- sum(mxm / thetam))
(TFRf <- sum(myf / thetaf))
harmonic.mean(c(TFRm, TFRf)) #hmmm

# I've come to the conclusion that Das Gupta essentially banks on a *harmonic mean*, but with no
# lengthy justification.



