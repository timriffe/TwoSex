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

PxUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxUS.Rdata"))) 
PxES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxUS.Rdata"))) 

# get male and female lambdas
DTLTUUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_TLTU/DTLTUUS.Rdata"))) 
DTLTUES <- local(get(load("/home/triffe/git/DISS/Data/HMD_TLTU/DTLTUES.Rdata"))) 


getLambda <- function(DTLTU, sex = "Female"){
    D0 <- reshape2::acast(DTLTU[DTLTU$Age == 0,],Year~Cohort,sum,value.var=sex)
    D0[D0==0] <- NA
    D0 <- do.call(cbind,apply(D0,2,function(x){
                        if (sum(!is.na(x))==2){
                            x[!is.na(x)]
                        }
                    }))
    lambda        <- D0[1,] / colSums(D0)
    years         <- as.integer(names(lambda))
    years         <-c(years,years[length(years)]+1)
    lambda        <- c(lambda,lambda[length(lambda)])
    names(lambda) <- years
    lambda
}

lambdamUS <- getLambda(DTLTUUS, "Male")
lambdafUS <- getLambda(DTLTUUS, "Female")
lambdamES <- getLambda(DTLTUES, "Male")
lambdafES <- getLambda(DTLTUES, "Female")

# take a look at lambda:
# pretty stable, but clear pattern.
# don't bother including in thesis, just mention
# range observed.
plot(yearsUS, lambdamUS,  col = "blue", ylim = c(.84,.91), pch = 19)
points(yearsUS, lambdafUS, col = "red", pch = 19)
points(yearsES, lambdamES, col = "blue", pch = 17)
points(yearsES, lambdafES, col = "red", pch = 17)
lines(lowess(lambdafES~yearsES), col = "red", lty = 2)
lines(lowess(lambdamES~yearsES), col = "blue", lty = 2)
lines(lowess(lambdafUS~yearsUS), col = "red")
lines(lowess(lambdamUS~yearsUS), col = "blue")




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

# ------------------------------------------------------------------
# make 1-sex Leslie matrices for males and females
Y_1sexUS <- lapply(as.character(yearsUS), function(yr, .dxm, .dxf, .Ex, .Bxymf, .lambdaf, .lambdam){
            yri     <- as.integer(yr)
            .dxm.   <- .dxm[, yr]
            .dxf.   <- .dxf[, yr]
            ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
            BxMM    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
            
            ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
            BxFF    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
         
            # sex-sex-ex- specific rates:
            FxMM    <- BxMM / ExM
            FxFF    <- BxFF / ExF
  
            
            list(Yf = MakeLExOneSexProjMatrix(FxFF, .dxf., .lambdaf[yr]),
            Ym = MakeLExOneSexProjMatrix(FxMM, .dxm., .lambdam[yr]))
        }, .dxm = dxmUS, .dxf = dxfUS, 
        .Ex = ExUS, .Bxymf = BxymfUS, 
        .lambdaf = lambdafUS, .lambdam = lambdafES)

names(Y_1sexUS) <- yearsUS
# repeat for Spain
Y_1sexES <- lapply(as.character(yearsES), function(yr, .dxm, .dxf, .Ex, .Bxymf, .lambdaf, .lambdam){
            yri     <- as.integer(yr)
            .dxm.   <- .dxm[, yr]
            .dxf.   <- .dxf[, yr]
            ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
            BxMM    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
            
            ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
            BxFF    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
            
            # sex-sex-ex- specific rates:
            FxMM    <- BxMM / ExM
            FxFF    <- BxFF / ExF
            
            
            list(Yf = MakeLExOneSexProjMatrix(FxFF, .dxf., .lambdaf[yr]),
                    Ym = MakeLExOneSexProjMatrix(FxMM, .dxm., .lambdam[yr]))
        }, .dxm = dxmES, .dxf = dxfES, 
        .Ex = ExES, .Bxymf = BxymfES, 
        .lambdaf = lambdafES, .lambdam = lambdafES)

names(Y_1sexES) <- yearsES

library(popbio)

ev <- eigen(Y_1sexES[[1]][[1]])
ev$values
Mod(ev$values)
lmax<-which.max(Re(ev$values))
lmax
Re(ev$values)[lmax]
log(max(Re(ev$values)))
rES1sex <- do.call(rbind,lapply(Y_1sexES, function(x){
            c(r.f = log(max(Re(eigen(x[["Yf"]])$values))),
            r.m = log(max(Re(eigen(x[["Ym"]])$values))))
        }))
rUS1sex <- do.call(rbind,lapply(Y_1sexUS, function(x){
                    c(r.f = log(max(Re(eigen(x[["Yf"]])$values))),
                            r.m = log(max(Re(eigen(x[["Ym"]])$values))))
                }))

plot(yearsUS, rUS1sex[,"r.f"], type = 'l', col = "red", ylim = c(-.02,.02))
lines(yearsUS, rUS1sex[,"r.m"], col = "blue")
lines(yearsES, rES1sex[,"r.f"], col = "red",lty=2)
lines(yearsES, rES1sex[,"r.m"], col = "blue",lty=2)



