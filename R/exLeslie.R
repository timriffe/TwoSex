source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxES is 0:110, years 1975:2009
BxUS  <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
BxES  <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

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

yearsUS <- 1969:2009
yearsES <- 1975:2009


# 1975 US
ExmUS1975 <- rowSums(ExpectedDx( with(ExUS, Male[Year == 1975]), dxmUS[, "1975"]))
BxmUS1975 <- rowSums(ExpectedDx( rowSums(BxUS[["1975"]]), dxmUS[, "1975"]))
ExfUS1975 <- rowSums(ExpectedDx( with(ExUS, Female[Year == 1975]), dxmUS[, "1975"]))
BxfUS1975 <- rowSums(ExpectedDx( colSums(BxUS[["1975"]]), dxfUS[, "1975"]))

FxmUS1975 <- BxmUS1975 / ExmUS1975
FxfUS1975 <- BxfUS1975 / ExfUS1975

plot(0:110, FxmUS1975, type = 'l', col = "blue", ylim = c(0,.07))
lines(0:110, FxfUS1975, col = "red")
Fx <- FxmUS1975
dx <- dxmUS[, "1975"]
MakeLeslie <- function(Fx, dx){
    rowSums(outer(dx,Fx[2:N], "*"))
    rowSums(ExpectedDx(Fx[2:N],dx))
   
    N <- length(Fx)
    rbind(cbind(0,diag(Fx[2:N] + 1)),0)
}
MakeLeslie <- function(Fx, dx){
    N <- length(Fx)
    cbind(0,outer(dx, Fx[2:N], "*")) + 
          rbind(cbind(0,diag(N -1)),0)
}
#install.packages("popbio")
library(popbio)
EigAn <- eigen.analysis(Lex, zero=TRUE)
names(EigAn)
plot(0:110,EigAn$stable.stage)
stage.vector.plot(pop.projection(Lex, FxmUS1975, 100)$stage.vectors)

Lex <- MakeLeslie(FxmUS1975, dxmUS[, "1975"])
eigen.analysis(Lex, zero=TRUE)[1]
eigen.analysis( MakeLeslie(FxfUS1975, dxfUS[, "1975"]), zero=TRUE)[1]

ExmUS1975i <- ExmUS1975
plot(0:110, ExmUS1975i / sum(ExmUS1975i), type = 'l', ylim = c(0,.05))
cols <- paste0(gray(.5),"50")
for (i in 400){
    ExmUS1975i <- c(Lex %*% ExmUS1975i)
    lines(0:110,  ExmUS1975i / sum(ExmUS1975i), col = cols[i])
}
lines()
ExmUS19752 <- c(Lex %*% ExmUS1975i)
log(sum(ExmUS19752)/sum(ExmUS1975i))
EigAn[[1]]


sum(ExmUS1975i) - sum(ExmUS1975)

(B <- sum(ExmUS1975 * FxmUS1975))
(Deaths <- ExmUS1975[1])

(B - Deaths) - (sum(ExmUS1975i) - sum(ExmUS1975))
(ExmUS1975 * FxmUS1975)[1] # infant mortality
ExmUS1975i[1] - ExmUS1975[2]

plot(0:110, ExmUS1975i / sum(ExmUS1975i), type = 'l')
lines(0:110, ExmUS1975i / sum(ExmUS1975i), col = "red")








