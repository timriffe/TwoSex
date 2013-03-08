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


MakeLExOneSexProjMatrix <- function(Fx, dx, lambda){
    N       <- length(Fx)
    # discount for part of infant mortality not surviving until end of year
    dx[1]   <- dx[1] * lambda
    
    # NxN matrix
               # fertility component
    Y       <- outer(dx, Fx, "*") + 
               # add survival element-wise
               rbind(cbind(0,diag(N - 1)),0)
       
    # reduce e0 fertility by 1/2, as only exposed for part of year
    Y[, 1]  <- Y[, 1] / 2
    # do not allow for Inf or NA values: impute 0
    Y       <- Mna0(Minf0(Y))
    # return projection matrix
    Y
}

# 1) sigma, male weight

MakeLExTwoSexProjMatrix <- function(dxm, dxf, FexFF, FexFM, FexMM, FexMF, lambdaM, lambdaF, sigma = .5){
    N <- length(dxm)
    dxm[1] <- dxm[1] * lambdaM
    dxf[1] <- dxf[1] * lambdaF
    
    # define matrix, then discount first column for year t mortality
    Y <-
    cbind(
    # --------------------------------------------
      # Male side (left half)
      rbind(
        
        # upper left block
        # male-male fert
        sigma * outer(dxm, Fexm, "*") + 
        # add element-wise
        # male survival 
        rbind(cbind(0, diag(N - 1)),0),
        
        # lower left block
        # male- female fert
        sigma * outer(dxf, Fexm, "*") 
      ),
    # --------------------------------------------
      # Female side (right half)
        rbind(
        # upper right block
        # female- male fert
        (1 - sigma) * outer(dxm, Fexf, "*"), 
    
        # lower right block
        # female survival and female-female fert
        (1 - sigma) * outer(dxf, Fexf, "*") + 
        rbind(cbind(0, diag(N - 1)), 0)
      )
    )
    # discount fertility of those dying in year t (column 1) by half.
    Y[,1] <- Y[,1]/2
    
    # do not allow for Inf or NA values: impute 0
    Y     <- Mna0(Minf0(Y))
    # return projection matrix
    Y
}

#Lex2 <- MakeLExTwoSexProjMatrix()

