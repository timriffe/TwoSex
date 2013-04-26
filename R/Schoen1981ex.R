# mean function 2 sex ex
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

yearsUS <- 1969:2009
yearsES <- 1975:2009

# Bxy totals  (not used()
BxUS  <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
BxES  <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

# repeated below for better SRB assumptions
BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))
names(BxymfES) <- yearsES
names(BxymfUS) <- yearsUS
# exposures, as such, straight from HMD, all ages 0-110, long form
ExUS  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
dxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfES.Rdata"))) 
# make sum to 1
dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)

#------------------------------------------------------------
# compare with Lotka:
LxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmUS.Rdata"))) / 1e5
LxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfUS.Rdata"))) / 1e5
LxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmES.Rdata"))) / 1e5
LxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfES.Rdata"))) / 1e5


# harmonic mean of exposures:
HM <- compiler::cmpfun(function(x,y){
            Mna0(Minf0((2 * x * y) / (x + y)))
        })

# residual function
LotkaHminex <- compiler::cmpfun(function(r, dxm, dxf, FHf, FHm, SRB, M = HM, .a = .5:110.5){
          
            N               <- length(dxm)
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
            p.m <- SRB / (1 + SRB)
            p.f <- 1/ (1 + SRB)
            (1 - sum(outer(p.m * colSums(t(dxM) / exp(-r * .a)), 
                     p.f * colSums(t(dxF) / exp(-r * .a)),
                     M) * (FHf + FHm)))^2
                          
            
        })
yr  <- "1990"
dxm <- dxmUS[, yr]
dxf <- dxfUS[, yr]

BxyF <- ExpectedDxMxFmatrix(BxymfUS[[yr]][["Bxyf"]], dxm, dxf)
BxyM <- ExpectedDxMxFmatrix(BxymfUS[[yr]][["Bxym"]], dxm, dxf)

ExM <- rowSums(ExpectedDx(with(ExUS,Male[Year == as.integer(yr)]), dxm))
ExF <- rowSums(ExpectedDx(with(ExUS,Female[Year == as.integer(yr)]), dxf))
# get harmonic rates for boys, girls
FHf <- BxyF / outer(ExM, ExF, HM)
FHm <- BxyM / outer(ExM, ExF, HM)
SRB <- sum(BxyM) / sum(BxyF)
#optimize(LotkaHminex, interval = c(-.02,.02), dxm = dxm, dxf = dxf, FHf = FHf, FHm = FHm, SRB = SRB)


# using modified strategy of Coale
exMeanIt <- compiler::cmpfun(function(FHf, FHm, dxm, dxf, M = HM, .a = .5:110.5, maxit = 2e2, tol = 1e-15){
            # from Coale, Ansley J. (1957) A New Method for Calculating Lotka's r- the Intrinsic Rate of Growth in a Stable Population.
            # Population Studies, Vol. 11 no. 1, pp 92-94
            N               <- length(dxm)
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
            # first starting value for proportion male bzw female at birth
            p.m <- 1.05 / 2.05
            p.f <- 1 / 2.05
            # now see what SRB the data would produce given the above SRB and the given rates
            SRBi        <- sum(outer(p.m * rowSums(dxM), 
                                     p.f * rowSums(dxF), M) * FHm) / 
                           sum(outer(p.m * rowSums(dxM), 
                                     p.f * rowSums(dxF), M) * FHf)
            # reconvert to proportions
            p.m         <- SRBi / (1 + SRBi)
            p.f         <- 1 / (1 + SRBi)
            # R0, assumign r = 0
            R0          <- sum(outer(p.m * rowSums(dxM), 
                                     p.f * rowSums(dxF), M) * (FHf + FHm))
            # mean generation length (yers from death) assuming r = 0
            T.guess    <- sum(outer(.a * p.m * rowSums(dxM), 
                                    .a * p.f * rowSums(dxF), M) * (FHf + FHm)) / R0
            # starting value for r given above 
            r2 <- log(R0) /  T.guess
            for (i in 1:maxit){ # 10 is more than enough!
                r1      <- r2
                p.m     <- SRBi / (SRBi + 1) # at each iteration, regenerate proportions
                p.f     <- 1 / (SRBi + 1)
                # main unity equation, producing residual. Allow successive generations to grow according to r,
                # but compress [rowSums()]for total model exposure. No sex differentiation of births here, hence
                # (FHf + FHm)
                deltai  <- 1 - sum(outer(rowSums(dxM %col% (1 / exp(-r1 * .a))) * p.m, 
                                         rowSums(dxF %col% (1 / exp(-r1 * .a))) * p.f, M) * (FHf + FHm))

                # calibrate r according to the error produced by the Lotka equation
                r2      <- r1 - (deltai / (T.guess - (deltai / r1)))
                # use improved r (and old p.m, p.f) to update SRB (won't move much)
                SRBi    <-  sum(outer(rowSums(dxM %col% (1 / exp(-r2 * .a))) * p.m, 
                                        rowSums(dxF %col% (1 / exp(-r2 * .a)))* p.f, M) * FHm) / 
                            sum(outer(rowSums(dxM %col% (1 / exp(-r2 * .a))) * p.m, 
                                        rowSums(dxF %col% (1 / exp(-r2 * .a)))* p.f, M) * FHf)
                if (abs(deltai) < tol){ # if converged, break
                    break
                }
            }
            # spit back stable r and SRB.
            return(c(r=r2, SRB=SRBi, iterations = i))
        })
        

rUS <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxy, .dxm, .dxf, .Ex, .M){
                    
                    Exm <- rowSums(ExpectedDx(with(.Ex, Male[Year == as.integer(yr)]), .dxm[,yr]))
                    Exf <- rowSums(ExpectedDx(with(.Ex, Female[Year == as.integer(yr)]), .dxf[,yr]))
                    # harmonic mean of exposures
                    Hxy <- outer(Exm, Exf, .M)
                    # harmonic rates (divide by two since working with both sexes)
                    FxyHf <- Mna0(Minf0(ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxyf"]], .dxm[, yr], .dxf[, yr]) / Hxy))
                    FxyHm <- Mna0(Minf0(ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxym"]], .dxm[, yr], .dxf[, yr]) / Hxy))
                                      
                    LH <- exMeanIt(FHf = FxyHf, FHm = FxyHm, dxm =.dxm[,yr], dxf = .dxf[,yr], M = .M)
                    c(r = LH[1], SRB = LH[2], its = LH[3])
                }, .Bxy = BxymfUS, .dxm = dxmUS, .dxf = dxfUS, .Ex = ExUS, .M = HM))
rownames(rUS) <- yearsUS
rES <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxy, .dxm, .dxf, .Ex, .M){
                    
                    Exm <- rowSums(ExpectedDx(with(.Ex, Male[Year == as.integer(yr)]), .dxm[,yr]))
                    Exf <- rowSums(ExpectedDx(with(.Ex, Female[Year == as.integer(yr)]), .dxf[,yr]))
                    # harmonic mean of exposures
                    Hxy <- outer(Exm, Exf, .M)
                    # harmonic rates (divide by two since working with both sexes)
                    FxyHf <- Mna0(Minf0(ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxyf"]], .dxm[, yr], .dxf[, yr]) / Hxy))
                    FxyHm <- Mna0(Minf0(ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxym"]], .dxm[, yr], .dxf[, yr]) / Hxy))
                    
                    LH <- exMeanIt(FHf = FxyHf, FHm = FxyHm, dxm =.dxm[,yr], dxf = .dxf[,yr], M = .M)
                    c(r = LH[1], SRB = LH[2], its = LH[3])
                }, .Bxy = BxymfES, .dxm = dxmES, .dxf = dxfES, .Ex = ExES, .M = HM))
rownames(rES) <- yearsES


plot(yearsUS, rUS[,1], type = 'l',ylim = c(-.015,.01))
lines(yearsES, rES[,1], col = "red")

# lez get some stable ey structure:
r <- rUS[1,1]
SRB <- rUS[1,2]
dxm <- dxmUS[,"1969"]
dxf <- dxfUS[,"1969"]
exMstableex <- function(r, SRB, dxm, dxf, M = HM, .a = .5:110.5){
    p.m <- SRB / (1 + SRB)
    p.f <- 1 / (1 + SRB)
    
    N               <- length(dxm)
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
    
    b <-  1 / sum(colSums(t(dxM) * exp(-r * .a)) * p.m
                    + colSums(t(dxF) * exp(-r * .a)) * p.f)
    cbind(cym = b * p.m * colSums(t(dxM) * exp(-r * .a)),
          cyf = b * p.f * colSums(t(dxF) * exp(-r * .a)))
}


sum(cbind(cym = b * p.m * colSums(t(dxM) * exp(-r * .a)),
        cyf = b * p.f * colSums(t(dxF) * exp(-r * .a))))
