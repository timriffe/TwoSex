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

sum(FxMM)
sum(FxMF)
sum(FxMM)
sum(FxMM)

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
        .a = .5:110.5, sigma = .5, maxit = 1e2, tol = 1e-15, r.start = 0.001){  
    
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
        #cat(r2,i,"\n")
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


#r <- exTwoSexLinearCoaleR(dxm = dxmUS[, "1975"], dxf = dxfUS[, "1975"],
#       FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF, sigma = .5)


rUS <- unlist(lapply(as.character(yearsUS), function(yr, .Bxymf, .dxm, .dxf, .Ex, .sigma){
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
            
            exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                    sigma = .sigma)
            
            
        }, .Bxymf = BxymfUS, .dxm = dxmUS, .dxf = dxfUS, .Ex = ExUS, .sigma = .5))

rES <- unlist(lapply(as.character(yearsES), function(yr, .Bxymf, .dxm, .dxf, .Ex, .sigma){
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
                    
                    exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
                            FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                            sigma = .sigma)
                    
                    
                }, .Bxymf = BxymfES, .dxm = dxmES, .dxf = dxfES, .Ex = ExES, .sigma = .5))

plot(yearsUS, rUS, type = 'l', ylim = c(-.015,.01))
lines(yearsES, rES, col = "red")

exTwoSexLinearTy <- function(r, dxm, dxf, FexFF, FexFM, FexMM, FexMF, .a = .5:110.5, sigma = .5){
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
}

exTwoSexLinearTy(r=r,dxm = dxmUS[, "1975"], dxf = dxfUS[, "1975"],
        FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF, sigma = .5) # sigma must be same as for r
(50.61+41.64) / 2