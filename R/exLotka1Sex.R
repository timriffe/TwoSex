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


ExmUS1975 <- rowSums(ExpectedDx( with(ExUS, Male[Year == 1975]), dxmUS[, "1975"]))
BxmUS1975 <- rowSums(ExpectedDx( rowSums(BxUS[["1975"]]), dxmUS[, "1975"]))
ExfUS1975 <- rowSums(ExpectedDx( with(ExUS, Female[Year == 1975]), dxfUS[, "1975"]))
BxfUS1975 <- rowSums(ExpectedDx( colSums(BxUS[["1975"]]), dxfUS[, "1975"]))

FxmUS1975 <- BxmUS1975 / ExmUS1975
FxfUS1975 <- BxfUS1975 / ExfUS1975

# minimizer function for 1 sex ex-perspective renewal function:

exOneSexMin <- function(r, dx, Fex, .a = .5:110.5){
    # get the overlapped / staggered dx structure
    dxM  <- matrix(0,ncol=111,nrow=111)
    dxi  <- rev(dx)
    for (i in 1:111){
        dxM[i:111, i] <- dxi 
        dxi <- dxi[1:(length(dxi)-1)]
    }     
    (1 - sum(rowSums(dxM %col% (1 /  exp(-r * .a))) * Fex)) ^ 2
}
optimize(exOneSexMin, interval = c(-.2,.2), dx = dxfUS[, "1975"], Fex = FxfUS1975 * (1/2.05))$minimum
optimize(exOneSexMin, interval = c(-.2,.2), dx = dxmUS[, "1975"], Fex = FxmUS1975 * (1.05/2.05))$minimum


exOneSexCoaleR <- function(Fex, dx, .a = .5:110.5, maxit = 1e2, tol = 1e-15, T.guess, r.start = .01){  
    N               <- length(Fex)
    dxM             <- matrix(0,ncol=N,nrow=N)
    dxi             <- rev(dx)
    for (i in 1:N){
        dxM[i:N, i] <- dxi 
        dxi         <- dxi[1:(length(dxi)-1)]
    } 
    
    if(missing(T.guess)) {
        T.guess     <- wmean(.a, rowSums(dxM) * Fex)
    }
    ri              <- r.start
    # be careful to discount Fex by SRB appropriately for males / females
    # prior to specification
    # Based on Coale (1957)
    for (i in 1:maxit){ # 15 is more than enough!
        deltai <- 1 - sum(rowSums(dxM %col% (1 /  exp(-ri * .a))) * Fex)
        # the mean generation time self-corrects 
        # according to the error produced by the Lotka equation
        ri <- ri + (deltai / (T.guess - (deltai / ri)))     
    }
    return(ri)  
}


exOneSexCoaleR(FxfUS1975 * (1/2.05),dxfUS[, "1975"], T.guess = 30, r.start = .001)






