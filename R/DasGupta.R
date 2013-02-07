source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

BxUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy10_65.Rdata")))
# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))



Bxy         <- BxUS[["1970"]]
rowTotals   <- rowSums(Bxy)
colTotals   <- colSums(Bxy)
nOfCases    <- sum(rowTotals)
expected    <- outer(rowSums(Bxy), colSums(Bxy), "*") / nOfCases



par(mfrow=c(1,2))
image(x=10:65,y=10:65,log(t(Bxy)),asp=1, zlim = c(0,12))
image(x=10:65,y=10:65,log(t(expected)),asp=1, zlim = c(0,12))
log(60000)
max(expected)

Diff <- Bxy - expected
image(x=10:65,y=10:65,t(Diff),asp=1)
## Fisher's exact test gives p = 0.4857 ...
# not sure what I'm supposed to see in equation 4.1 everything cancels!


Pxm <- with(ExUS, Male[Age >= 10 & Age <= 65 & Year == 1970])
Pxf <- with(ExUS, Female[Age >= 10 & Age <= 65 & Year == 1970])        

Bxy <- BxUS[["1970"]]
# eq 5.4
Uxy <- Mna0(Bxy %row% (1 / rowSums(Bxy)))
Vxy <- Mna0(Bxy %col% (1 / colSums(Bxy)))

Uxy <- Uxy / sum(Uxy)
Vxy <- Vxy / sum(Vxy)
thetam <- 1.05 / 2.05
thetaf <- 1 / 2.05

Mxy <- Uxy %row% Pxm 
Fxy <- Uxy %col% Pxf 


mxy <- Mna0(Bxy / (Uxy %row% Pxm + Vxy %col% Pxf))

mx <- rowSums(Bxy) / Pxm
my <- colSums(Bxy) / Pxf
mxm  <- thetam * mx
myf  <- thetaf * my
mxy <- (thetam / mxm + thetaf / myf) ^ -1
sum(mxy) 
sum(mxy)



plot(10:65, mxy, col = "red", type = "l")
polygon(c(10:65, 65:10), c(mx, rev(my)), col = "#55555550")

Bt <- sum((thetam * Uxy %row% Pxm + thetaf * Vxy %col% Pxf) * )

sum(apply(cbind(mxm / (2 * thetam), myf / (2 * thetaf)), 1, harmonic.mean))
sum(apply(cbind(mxm / (thetam), myf / ( thetaf)), 1, harmonic.mean))
sum(apply(cbind(mx, my), 1, harmonic.mean))

(TFRm <- sum(mxm / thetam))
(TFRf <- sum(myf / thetaf))
harmonic.mean(c(TFRm, TFRf)) #hmmm

# I've come to the conclusion that Das Gupta essentially banks on a harmonic mean, but with no
# lengthy justification.



