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


ExM     <- rowSums(ExpectedDx( with(ExUS, Male[Year == 1975]), dxmUS[, "1975"]))
BxMM    <- rowSums(ExpectedDx( rowSums(BxymfUS[["1975"]][["Bxym"]], na.rm = TRUE), dxmUS[, "1975"]))
BxMF    <- rowSums(ExpectedDx( rowSums(BxymfUS[["1975"]][["Bxyf"]], na.rm = TRUE), dxmUS[, "1975"]))

BxM    <- rowSums(ExpectedDx( rowSums(BxUS[["1975"]], na.rm = TRUE), dxmUS[, "1975"]))
BxF    <- rowSums(ExpectedDx( colSums(BxUS[["1975"]], na.rm = TRUE), dxfUS[, "1975"]))
ExF     <- rowSums(ExpectedDx( with(ExUS, Female[Year == 1975]), dxfUS[, "1975"]))
BxFF    <- rowSums(ExpectedDx( colSums(BxymfUS[["1975"]][["Bxyf"]], na.rm = TRUE), dxfUS[, "1975"]))
BxFM    <- rowSums(ExpectedDx( colSums(BxymfUS[["1975"]][["Bxym"]], na.rm = TRUE), dxfUS[, "1975"]))
# some rates:
FxMM    <- BxMM / ExM
FxFF    <- BxFF / ExF
FxMF    <- BxMF / ExM
FxFM    <- BxFM / ExF
# minimizer function for 1 sex ex-perspective renewal function:


exTwoSexLinearMin <- compiler::cmpfun(function(r, dxm, dxf, FexFF, FexFM, FexMM, FexMF, .a = .5:110.5, sigma = .5){
    # all input vectors (except r) must be same length!
    # get the overlapped / staggered dx structure
    N               <- length(FexFF)
    dxM  <- dxF   <- matrix(0, ncol = N, nrow = N)
    # remaining years go down rows. ages over columns
    dxmi  <- dxm
    dxfi  <- dxf
    for (i in 1:N){
        dxM[i, 1:length(dxmi)  ] <- dxmi 
        dxmi <- dxmi[2:length(dxmi) ]
        
        dxF[i, 1:length(dxfi)  ] <- dxfi 
        dxfi <- dxfi[2:length(dxfi) ]
    }    
    
    
    (1 - sum(
              sigma * rowSums(dxM %col% (1 / exp(-r * .a))) * FexMM + # male - male
              sigma * rowSums(dxF %col% (1 / exp(-r * .a))) * FexMF + # male - female
          (1-sigma) * rowSums(dxM %col% (1 / exp(-r * .a))) * FexFM + # female - male
          (1-sigma) * rowSums(dxF %col% (1 / exp(-r * .a))) * FexFF   # female - female
                )
                ) ^ 2
})

#optimize(exTwoSexLinearMin, interval = c(-.2,.2), 
#        dxm = dxmUS[, "1975"], dxf = dxfUS[, "1975"],
#        FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF, sigma = .5, tol = 1e-15)$minimum

exTwoSexLinearCoaleR <- compiler::cmpfun(function(dxm, dxf, FexFF, FexFM, FexMM, FexMF, 
        .a = .5:110.5, sigma = .5, maxit = 2e2, tol = 1e-15, r.start = 0.001){  
    
    N               <- length(FexFF)
    dxM    <- dxF   <- matrix(0, ncol = N, nrow = N)
    # remaining years go down rows. ages over columns
    dxmi            <- dxm
    dxfi            <- dxf
    for (i in 1:N){
        dxM[i, 1:length(dxmi)  ] <- dxmi 
        dxmi                     <- dxmi[2:length(dxmi) ]
        
        dxF[i, 1:length(dxfi)  ] <- dxfi 
        dxfi                     <- dxfi[2:length(dxfi) ]
    }  
    
    R0      <- (sigma * sum(dxM*FexMM) +
               sigma * sum(dxF*FexMF) + 
               (1 - sigma) * sum(dxM*FexFM) +
               (1 - sigma) * sum(dxF*FexFF) ) / 2
       
    T.guess <- wmean(.a,
            sigma * rowSums(dxM) * FexMM + 
            sigma * rowSums(dxF) * FexMF +
            (1 - sigma) * rowSums(dxM) * FexFM +
            (1 - sigma) * rowSums(dxF) * FexFF 
    ) # assuming r = 0
    
    r2      <- log(R0) / T.guess
    
    # be careful to discount Fex by SRB appropriately for males / females
    # prior to specification
    # Based on Coale (1957)
    for (i in 1:maxit){ # 15 is more than enough!
        cat(r2,i,"\n")
        r1 <- r2
        deltai <- 2 - sum(
                sigma * rowSums(dxM %col% (1 / exp(-r1 * .a))) * FexMM + # male - male
                        sigma * rowSums(dxF %col% (1 / exp(-r1 * .a))) * FexMF + # male - female
                        (1-sigma) * rowSums(dxM %col% (1 / exp(-r1 * .a))) * FexFM + # female - male
                        (1-sigma) * rowSums(dxF %col% (1 / exp(-r1 * .a))) * FexFF   # female - female
        )
        # the mean generation time self-corrects 
        # according to the error produced by the Lotka equation
        r2 <- r1 - (deltai / (T.guess - (deltai / r1))) 
        if (abs(r2 - r1) <= tol | zapsmall(abs(deltai)) <= tol){
            break
        }
        
    }
    if (i == maxit){
        cat("WARNING: max iterations reached, r may not be solution")
    }
    return(r2)  
})

exTwoSexLinearTy <- compiler::cmpfun(function(r, dxm, dxf, FexFF, FexFM, FexMM, FexMF, .a = .5:110.5, sigma = .5){
    N               <- length(FexFF)
    dxM    <- dxF   <- matrix(0, ncol = N, nrow = N)
    # remaining years go down rows. ages over columns
    dxmi            <- dxm
    dxfi            <- dxf
    for (i in 1:N){
        dxM[i, 1:length(dxmi)  ] <- dxmi 
        dxmi                     <- dxmi[2:length(dxmi) ]
        
        dxF[i, 1:length(dxfi)  ] <- dxfi 
        dxfi                     <- dxfi[2:length(dxfi) ]
    }     
    
    wmean(.a,
            sigma * rowSums(dxM %col% (1 / exp(-r * .a))) * FexMM + 
                    sigma * rowSums(dxF %col% (1 / exp(-r * .a))) * FexMF +
                    (1 - sigma) * rowSums(dxM %col% (1 / exp(-r * .a))) * FexFM +
                    (1 - sigma) * rowSums(dxF %col% (1 / exp(-r * .a))) * FexFF 
    )
})

#r <- exTwoSexLinearCoaleR(dxm = dxmUS[, "1975"], dxf = dxfUS[, "1975"],
#       FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF, sigma = .5)


rUS.5 <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxymf, .dxm, .dxf, .Ex, .sigma){
            yri     <- as.integer(yr)
            .dxm. <- .dxm[, yr]
            .dxf. <- .dxf[, yr]
            ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
            BxMM    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
            BxMF    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxm.))
            
            ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
            BxFF    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
            BxFM    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxf.))
            # sex-sex-ex- specific rates:
            FxMM    <- BxMM / ExM
            FxFF    <- BxFF / ExF
            FxMF    <- BxMF / ExM
            FxFM    <- BxFM / ExF
            
            r = exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                    sigma = .sigma)
            c(r = r, Ty = exTwoSexLinearTy(r, dxm = .dxm., dxf = .dxf., 
                            FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                            sigma = .sigma))
            
        }, .Bxymf = BxymfUS, .dxm = dxmUS, .dxf = dxfUS, .Ex = ExUS, .sigma = .5))

rES.5 <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxymf, .dxm, .dxf, .Ex, .sigma){
                    yri     <- as.integer(yr)
                    .dxm. <- .dxm[, yr]
                    .dxf. <- .dxf[, yr]
                    ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
                    BxMM    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
                    BxMF    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxm.))
                    
                    ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
                    BxFF    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
                    BxFM    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxf.))
                    # sex-sex-ex- specific rates:
                    FxMM    <- Minf0(Mna0(BxMM / ExM))
                    FxFF    <- Minf0(Mna0(BxFF / ExF))
                    FxMF    <- Minf0(Mna0(BxMF / ExM))
                    FxFM    <- Minf0(Mna0(BxFM / ExF))
                    
                    r = exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
                            FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                            sigma = .sigma)
                    c(r = r, Ty = exTwoSexLinearTy(r, dxm = .dxm., dxf = .dxf., 
                                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                                    sigma = .sigma))
                    
                    
                }, .Bxymf = BxymfES, .dxm = dxmES, .dxf = dxfES, .Ex = ExES, .sigma = .5))

plot(yearsUS, rUS[,1], type = 'l', ylim = c(-.015,.01))
lines(yearsES, rES[,1], col = "red")

plot(yearsUS, rUS[,2], type = 'l', ylim = c(44,51))
lines(yearsES, rES[,2], col = "red")

exRepUS <- cbind(rUS, R0 = exp(rUS[,1]*rUS[,2]) )
exRepES <- cbind(rES, R0 = exp(rES[,1]*rES[,2]) )

colnames(exRepES) <- colnames(exRepUS) <- c("$r^\upsilon$","$T^\upsilon$","$R_0^\upsilon$")
expression(r^Upsilon1)


rES1 <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxymf, .dxm, .dxf, .Ex, .sigma, .maxit, .tol){
                    yri     <- as.integer(yr)
                    .dxm. <- .dxm[, yr]
                    .dxf. <- .dxf[, yr]
                    ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
                    BxMM    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
                    BxMF    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxm.))
                    
                    ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
                    BxFF    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
                    BxFM    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxf.))
                    # sex-sex-ex- specific rates:
                    FxMM    <- Minf0(Mna0(BxMM / ExM))
                    FxFF    <- Minf0(Mna0(BxFF / ExF))
                    FxMF    <- Minf0(Mna0(BxMF / ExM))
                    FxFM    <- Minf0(Mna0(BxFM / ExF))
                    
                    r = exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
                            FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                            sigma = .sigma, maxit = .maxit, tol = .tol)
                    c(r = r, Ty = exTwoSexLinearTy(r, dxm = .dxm., dxf = .dxf., 
                                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                                    sigma = .sigma))
                    
                    
                }, .Bxymf = BxymfES, .dxm = dxmES, .dxf = dxfES, .Ex = ExES, .sigma = 1, .max = 5e2, .tol = 1e-15))

rUS1 <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxymf, .dxm, .dxf, .Ex, .sigma, .maxit, .tol){
                    yri     <- as.integer(yr)
                    .dxm. <- .dxm[, yr]
                    .dxf. <- .dxf[, yr]
                    ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
                    BxMM    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
                    BxMF    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxm.))
                    
                    ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
                    BxFF    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
                    BxFM    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxf.))
                    # sex-sex-ex- specific rates:
                    FxMM    <- Minf0(Mna0(BxMM / ExM))
                    FxFF    <- Minf0(Mna0(BxFF / ExF))
                    FxMF    <- Minf0(Mna0(BxMF / ExM))
                    FxFM    <- Minf0(Mna0(BxFM / ExF))
                    
                    r = exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
                            FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                            sigma = .sigma, maxit = .maxit, tol = .tol)
                    c(r = r, Ty = exTwoSexLinearTy(r, dxm = .dxm., dxf = .dxf., 
                                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                                    sigma = .sigma))
                    
                    
                }, .Bxymf = BxymfUS, .dxm = dxmUS, .dxf = dxfUS, .Ex = ExUS, .sigma = 1, .max = 5e2, .tol = 1e-15))

rES0 <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxymf, .dxm, .dxf, .Ex, .sigma, .maxit, .tol){
                    yri     <- as.integer(yr)
                    .dxm. <- .dxm[, yr]
                    .dxf. <- .dxf[, yr]
                    ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
                    BxMM    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
                    BxMF    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxm.))
                    
                    ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
                    BxFF    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
                    BxFM    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxf.))
                    # sex-sex-ex- specific rates:
                    FxMM    <- Minf0(Mna0(BxMM / ExM))
                    FxFF    <- Minf0(Mna0(BxFF / ExF))
                    FxMF    <- Minf0(Mna0(BxMF / ExM))
                    FxFM    <- Minf0(Mna0(BxFM / ExF))
                    
                    r = exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
                            FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                            sigma = .sigma, maxit = .maxit, tol = .tol)
                    c(r = r, Ty = exTwoSexLinearTy(r, dxm = .dxm., dxf = .dxf., 
                                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                                    sigma = .sigma))
                    
                    
                }, .Bxymf = BxymfES, .dxm = dxmES, .dxf = dxfES, .Ex = ExES, .sigma = 0, .max = 5e2, .tol = 1e-15))

rUS0 <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxymf, .dxm, .dxf, .Ex, .sigma, .maxit, .tol){
                    yri     <- as.integer(yr)
                    .dxm. <- .dxm[, yr]
                    .dxf. <- .dxf[, yr]
                    ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
                    BxMM    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
                    BxMF    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxm.))
                    
                    ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
                    BxFF    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
                    BxFM    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxf.))
                    # sex-sex-ex- specific rates:
                    FxMM    <- Minf0(Mna0(BxMM / ExM))
                    FxFF    <- Minf0(Mna0(BxFF / ExF))
                    FxMF    <- Minf0(Mna0(BxMF / ExM))
                    FxFM    <- Minf0(Mna0(BxFM / ExF))
                    
                    r = exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
                            FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                            sigma = .sigma, maxit = .maxit, tol = .tol)
                    c(r = r, Ty = exTwoSexLinearTy(r, dxm = .dxm., dxf = .dxf., 
                                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                                    sigma = .sigma))
                    
                    
                }, .Bxymf = BxymfUS, .dxm = dxmUS, .dxf = dxfUS, .Ex = ExUS, .sigma = 0, .max = 5e2, .tol = 1e-15))

ES <- data.frame(r0=rES0[, 1], r.5=rES.5[, 1], r1=rES1[, 1],T0=rES0[, 2], T.5=rES.5[, 2], T1=rES1[, 2],
        R0.0 = exp(rES0[,1]*rES0[,2]),  R0.5 = exp(rES.5[,1]*rES.5[,2]), R0.1 = exp(rES1[,1]*rES1[,2]))

colnames(ES) <- c("$r^{\\upsilon (\\sigma = 0)}$"  , "$r^{\\upsilon (\\sigma = .5)}$"  , "$r^{\\upsilon (\\sigma = 1)}$",
        "$T^{\\upsilon (\\sigma = 0)}$"  , "$T^{\\upsilon (\\sigma = .5)}$"  , "$T^{\\upsilon (\\sigma = 1)}$",
        "$R_0^{\\upsilon (\\sigma = 0)}$", "$R_0^{\\upsilon (\\sigma = .5)}$", "$R_0^{\\upsilon (\\sigma = 1)}$")
rownames(ES) <- yearsES
library(xtable)
print(xtable(ES, digits = c(0,4,4,4,2,2,2,3,3,3), align = c("c","c","c","c","c","c","c","c","c","c")),
        sanitize.colnames.function = identity, 
        file = "/home/triffe/git/DISS/latex/xtables/ex2sexlinearES.tex",floating=FALSE)


US <- data.frame(r0=rUS0[, 1], r.5=rUS.5[, 1], r1=rUS1[, 1],T0=rUS0[, 2], T.5=rUS.5[, 2], T1=rUS1[, 2],
        R0.0 = exp(rUS0[,1]*rUS0[,2]),  R0.5 = exp(rUS.5[,1]*rUS.5[,2]), R0.1 = exp(rUS1[,1]*rUS1[,2]))

colnames(US) <- c("$r^{\\upsilon (\\sigma = 0)}$"  , "$r^{\\upsilon (\\sigma = .5)}$"  , "$r^{\\upsilon (\\sigma = 1)}$",
        "$T^{\\upsilon (\\sigma = 0)}$"  , "$T^{\\upsilon (\\sigma = .5)}$"  , "$T^{\\upsilon (\\sigma = 1)}$",
        "$R_0^{\\upsilon (\\sigma = 0)}$", "$R_0^{\\upsilon (\\sigma = .5)}$", "$R_0^{\\upsilon (\\sigma = 1)}$")
rownames(US) <- yearsUS

print(xtable(US, digits = c(0,4,4,4,2,2,2,3,3,3), align = c("c","c","c","c","c","c","c","c","c","c")),
        sanitize.colnames.function = identity, 
        file = "/home/triffe/git/DISS/latex/xtables/ex2sexlinearUS.tex",floating=FALSE)

pdf("/home/triffe/git/DISS/latex/Figures/exLotka2sexlinear.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, US[, 2], type = 'n', ylim = c(-.016,.01),xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,-.016,2010,.01,col = gray(.95), border=NA),
                abline(h = seq(-.016,.01,by = .002), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.016,.01,by = .002),seq(-.016,.01,by = .002), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.016, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.0175, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,.0115, "r", cex = 1, xpd = TRUE)))
polygon(c(yearsUS, rev(yearsUS)), c(US[, 1], rev(US[, 3])), col = "#99999950")
lines(yearsUS, US[, 2],col = gray(.2), lwd = 2)

polygon(c(yearsES, rev(yearsES)), c(ES[, 1], rev(ES[, 3])), col = "#99999950", lty = 5)
lines(yearsES, ES[, 2],col = gray(.2), lwd = 2, lty = 5)

# sigma = 0
text(1978,-0.006,expression(sigma == 0))
segments(1978,-0.0055,1976,US["1976", 1])
segments(1978,-0.0055,1982,ES["1982", 1])


# sigma = 1
text(c(1996,1995),c(0.003,-.008),expression(sigma == 1))
segments(1996,0.0025,1995,US["1995", 3])
segments(1995,-0.0085,1994,ES["1994", 3])


legend(1970,-.01, lty = c(1,5), col = gray(.2), lwd = 2, bty = "n",
        legend = c("US","ES"), xpd = TRUE)
dev.off()



sum(rmUS[,1] < US[,3])
yearsUS[26]

rfUS[,1] > US[,1]
rfES[,1] > ES[,1]
plot(yearsES, rfES[,1] , type = 'l', col = "blue")
lines(yearsES, ES[,1], col = "red")

rmUS[,1] < US[,1]



r <- US[1,2]

N               <- length(FexFF)
dxM    <- dxF   <- matrix(0, ncol = N, nrow = N)
# remaining years go down rows. ages over columns
dxmi            <- dxm
dxfi            <- dxf
for (i in 1:N){
    dxM[i, 1:length(dxmi)  ] <- dxmi 
    dxmi                     <- dxmi[2:length(dxmi) ]
    
    dxF[i, 1:length(dxfi)  ] <- dxfi 
    dxfi                     <- dxfi[2:length(dxfi) ]
}     


b <-       1 / sum( sigma * rowSums(dxM %col% (1 / exp(-r * .a))) + 
                sigma * rowSums(dxF %col% (1 / exp(-r * .a))) +
                (1 - sigma) * rowSums(dxM %col% (1 / exp(-r * .a))) +
                (1 - sigma) * rowSums(dxF %col% (1 / exp(-r * .a)))  )
b <-       1 / sum(rowSums(dxM %col% (1 / exp(-r * .a))) + 
               rowSums(dxF %col% (1 / exp(-r * .a))) )

plot(0:110, b *  rowSums(dxF %col% (1 / exp(-r * .a))), type = 'l', col = "red")
lines(0:110,b *  rowSums(dxM %col% (1 / exp(-r * .a))),col = "blue")
sum(rowSums(dxF %col% (1 / exp(-r * .a)))) / sum(rowSums(dxM %col% (1 / exp(-r * .a))))
sigma <- .5
boys <- sum(sigma * rowSums(dxM %col% (1 / exp(-r * .a))) * FexMM) + 
        sum((1-sigma) * rowSums(dxF %col% (1 / exp(-r * .a))) * FexFM)
girls <- sum(sigma * rowSums(dxM %col% (1 / exp(-r * .a))) * FexMF) + 
        sum((1-sigma) * rowSums(dxF %col% (1 / exp(-r * .a))) * FexFF)
boys / girls

