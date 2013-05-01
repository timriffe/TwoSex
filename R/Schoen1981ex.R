# mean function 2 sex ex
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

yearsUS <- 1969:2009
yearsES <- 1975:2009

# Bxy totals  (not used()
#BxUS  <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
#BxES  <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

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
#LotkaHminex <- compiler::cmpfun(function(r, dxm, dxf, FHf, FHm, SRB, M = HM, .a = .5:110.5){
#          
#            N               <- length(dxm)
#            dxM    <- dxF   <- matrix(0, ncol = N, nrow = N)
#            # remaining years go down rows. ages over columns
#            dxmi            <- dxm
#            dxfi            <- dxf
#            for (i in 1:N){
#                dxM[i, 1:length(dxmi)  ] <- dxmi 
#                dxmi                     <- dxmi[2:length(dxmi) ]
#                
#                dxF[i, 1:length(dxfi)  ] <- dxfi 
#                dxfi                     <- dxfi[2:length(dxfi) ]
#            }  
#            p.m <- SRB / (1 + SRB)
#            p.f <- 1/ (1 + SRB)
#            (1 - sum(outer(p.m * colSums(t(dxM) / exp(-r * .a)), 
#                     p.f * colSums(t(dxF) / exp(-r * .a)),
#                     M) * (FHf + FHm)))^2
#                          
#            
#        })
#yr  <- "1990"
#dxm <- dxmUS[, yr]
#dxf <- dxfUS[, yr]
#
#BxyF <- ExpectedDxMxFmatrix(BxymfUS[[yr]][["Bxyf"]], dxm, dxf)
#BxyM <- ExpectedDxMxFmatrix(BxymfUS[[yr]][["Bxym"]], dxm, dxf)
#
#ExM <- rowSums(ExpectedDx(with(ExUS,Male[Year == as.integer(yr)]), dxm))
#ExF <- rowSums(ExpectedDx(with(ExUS,Female[Year == as.integer(yr)]), dxf))
## get harmonic rates for boys, girls
#FHf <- BxyF / outer(ExM, ExF, HM)
#FHm <- BxyM / outer(ExM, ExF, HM)
#SRB <- sum(BxyM) / sum(BxyF)
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
        

rUShm <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxy, .dxm, .dxf, .Ex, .M){
                    
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
rownames(rUShm) <- yearsUS
rEShm <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxy, .dxm, .dxf, .Ex, .M){
                    
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
rownames(rEShm) <- yearsES



#plot(yearsUS, rUS[,1], type = 'l',ylim = c(-.015,.01))
#lines(yearsES, rES[,1], col = "red")

# lez get some stable ey structure:
do.cy <- FALSE
if (do.cy){
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


cyUShm <- do.call(rbind,lapply(as.character(yearsUS), function(yr, rSRB, .dxm, .dxf, .M){
                    
                    cbind(Year = as.integer(yr),exMstableex(r = rSRB[yr,1],
                            SRB = rSRB[yr,2],
                            dxm = .dxm[,yr],
                            dxf = .dxf[,yr],
                            M = .M))
                },rSRB = rUShm, .dxm = dxmUS, .dxf = dxfUS, .M = HM))
cyEShm <- do.call(rbind,lapply(as.character(yearsES), function(yr, rSRB, .dxm, .dxf, .M){
                    
                    cbind(Year = as.integer(yr),exMstableex(r = rSRB[yr,1],
                                    SRB = rSRB[yr,2],
                                    dxm = .dxm[,yr],
                                    dxf = .dxf[,yr],
                                    M = .M))
                },rSRB = rEShm, .dxm = dxmES, .dxf = dxfES, .M = HM))
#dim(cyEShm)
#head(cyUShm)
#library(Pyramid)
#for (yr in yearsUS){
#    Pyramid(males = cyUShm[cyUShm[,"Year"] == yr,"cym"], 
#            females = cyUShm[cyUShm[,"Year"] == yr,"cyf"], xlim = c(-1,1), main = yr, widths = rep(1,111))
#    Sys.sleep(1)
#}
#for (yr in yearsES){
#    Pyramid(males = cyEShm[cyEShm[,"Year"] == yr,"cym"], 
#            females = cyEShm[cyEShm[,"Year"] == yr,"cyf"], xlim = c(-1,1), main = yr, widths = rep(1,111))
#    Sys.sleep(1)
#}
}
# --------------------------------------------------
# judge stable ESFR vs original ESFR:
# 1975 and 2009, US and ES:

do.esfr.fig <- FALSE
if (do.esfr.fig){
USesfr <- do.call(rbind,lapply(as.character(yearsUS), function(yr, rSRB, .Bxy, .dxm, .dxf, .Ex, .M, .a = .5:110.5){

                    dxm <- .dxm[,yr]
                    dxf <- .dxf[,yr]
                    Exm <- rowSums(ExpectedDx(with(.Ex, Male[Year == as.integer(yr)]), dxm))
                    Exf <- rowSums(ExpectedDx(with(.Ex, Female[Year == as.integer(yr)]), dxf))
                    # harmonic mean of exposures
                    Hxy <- outer(Exm, Exf, .M)
                    Bxy <- ExpectedDxMxFmatrix(Mat = (.Bxy[[yr]][["Bxyf"]] + .Bxy[[yr]][["Bxym"]]), 
                            dxm = dxm, dxf = dxf)
                    Fxy <- Minf0(Mna0(Bxy / Hxy))
                    bxf <- rowSums(ExpectedDx(colSums(.Bxy[[yr]][["Bxyf"]] + .Bxy[[yr]][["Bxym"]]), dxf))
                    rowSums(Bxy) - bxf
                    # copied and pasted from structure function above
                    # get stable pop structure:
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
                    r       <- rSRB[yr,1]
                    SRB     <- rSRB[yr,2]
                    p.m     <- SRB / (1 + SRB)
                    p.f     <- 1 / (1 + SRB)
                    
                    b       <-  1 / sum(colSums(t(dxM) * exp(-r * .a)) * p.m
                                    + colSums(t(dxF) * exp(-r * .a)) * p.f)
                    cym     <- b * p.m * colSums(t(dxM) * exp(-r * .a))
                    cyf     <- b * p.f * colSums(t(dxF) * exp(-r * .a))
                    
                    bxy <- outer(cym, cyf, .M) * Fxy
                    # stable ESFR:
                    ESFRstm <- Minf0(Mna0(rowSums(bxy) / cym))
                    ESFRstf <- Minf0(Mna0(colSums(bxy) / cyf))
                    # original ESFR
                    ESFRm   <- Minf0(Mna0(rowSums(Bxy) / Exm))
                    ESFRf   <- Minf0(Mna0(colSums(Bxy) / Exf))
                    cbind(Year = as.integer(yr), ESFRstm = ESFRstm, ESFRm=ESFRm, ESFRstf = ESFRstf, ESFRf=ESFRf)
                }, .Bxy = BxymfUS, rSRB = rUShm, .dxm = dxmUS, .dxf = dxfUS, .Ex = ExUS, .M = HM))

ESesfr <- do.call(rbind,lapply(as.character(yearsES), function(yr, rSRB, .Bxy, .dxm, .dxf, .Ex, .M, .a = .5:110.5){
                    dxm <- .dxm[,yr]
                    dxf <- .dxf[,yr]
                    Exm <- rowSums(ExpectedDx(with(.Ex, Male[Year == as.integer(yr)]), dxm))
                    Exf <- rowSums(ExpectedDx(with(.Ex, Female[Year == as.integer(yr)]), dxf))
                    # harmonic mean of exposures
                    Hxy <- outer(Exm, Exf, .M)
                    Bxy <- ExpectedDxMxFmatrix((.Bxy[[yr]][["Bxyf"]]+.Bxy[[yr]][["Bxym"]]), dxm, dxf)
                    Fxy <- Minf0(Mna0(Bxy / Hxy))
                  
                    # copied and pasted from structure function above
                    # get stable pop structure:
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
                    r       <- rSRB[yr,1]
                    SRB     <- rSRB[yr,2]
                    p.m     <- SRB / (1 + SRB)
                    p.f     <- 1 / (1 + SRB)
                    
                    b       <-  1 / sum(colSums(t(dxM) * exp(-r * .a)) * p.m
                                    + colSums(t(dxF) * exp(-r * .a)) * p.f)
                    cym     <- b * p.m * colSums(t(dxM) * exp(-r * .a))
                    cyf     <- b * p.f * colSums(t(dxF) * exp(-r * .a))
                    
                    bxy <- outer(cym, cyf, .M) * Fxy
                    # stable ESFR:
                    ESFRstm <- Minf0(Mna0(rowSums(bxy) / cym))
                    ESFRstf <- Minf0(Mna0(colSums(bxy) / cyf))
                    # original ESFR
                    ESFRm   <- Minf0(Mna0(rowSums(Bxy) / Exm))
                    ESFRf   <- Minf0(Mna0(colSums(Bxy) / Exf))
                    cbind(Year = as.integer(yr), ESFRstm = ESFRstm, ESFRm=ESFRm, ESFRstf = ESFRstf, ESFRf=ESFRf)
                }, .Bxy = BxymfES, rSRB = rEShm, .dxm = dxmES, .dxf = dxfES, .Ex = ExES, .M = HM))

UScomphm <- do.call(rbind, lapply(yearsUS, function(yr, .esfr){
                    thisyr <- .esfr[.esfr[,1] == yr, 2:5]
                    TFRs <- colSums(thisyr)
                    mdiffcoef <- 1 - sum(pmin(thisyr[,1] / sum(thisyr[,1]), thisyr[,2] / sum(thisyr[,2])))
                    fdiffcoef <- 1 - sum(pmin(thisyr[,3] / sum(thisyr[,3]), thisyr[,4] / sum(thisyr[,4])))
                    c(TFRs,mdiffcoef=mdiffcoef,fdiffcoef=fdiffcoef)
                },.esfr = USesfr)) 
EScomphm <- do.call(rbind, lapply(yearsES, function(yr, .esfr){
                    thisyr <- .esfr[.esfr[,1] == yr, 2:5]
                    TFRs <- colSums(thisyr)
                    mdiffcoef <- 1 - sum(pmin(thisyr[,1] / sum(thisyr[,1]), thisyr[,2] / sum(thisyr[,2])))
                    fdiffcoef <- 1 - sum(pmin(thisyr[,3] / sum(thisyr[,3]), thisyr[,4] / sum(thisyr[,4])))
                    c(TFRs,mdiffcoef=mdiffcoef,fdiffcoef=fdiffcoef)
                },.esfr = ESesfr)) 
# 1)
# plot difference between stable and initial eTFR:
pdf("/home/triffe/git/DISS/latex/Figures/eTFRharmonic.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .3, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, UScomphm[,1] - UScomphm[,2], type = 'l', ylim = c(-.1, .1), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968, -.1, 2010, .1,col = gray(.95), border=NA),
                abline(h = seq(-.1, .1, by = .025), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.1, .1, by = .05),round(seq(-.1, .1, by = .05),2), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10), -.1, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.11, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1967, .11, "TFR diff", cex = 1, xpd = TRUE)))
lines(yearsUS,UScomphm[,3] - UScomphm[,4], col = gray(.5), lwd = 2)
lines(yearsES,EScomphm[,1] - EScomphm[,2],col = gray(.2), lty = 5)
lines(yearsES,EScomphm[,3] - EScomphm[,4], col = gray(.5), lwd = 2,lty = 5)

text(c(1990.918, 1991.334, 1990.849, 1991.057), 
        c(-0.02716403,  0.02811804, -0.08514278,  0.08137760),
        c("US males","US females","ES males","ES females"),pos = 4)
dev.off()




y <- 0:110
#plot(y, USesfr[USesfr[,1] == yr,2], type = 'l',col = "blue")
#lines(y, USesfr[USesfr[,1] == yr,3], col = "royalblue", lty =2)
#lines(y, USesfr[USesfr[,1] == yr,4], col = "red")
#lines(y, USesfr[USesfr[,1] == yr,5], col = "tomato", lty=2)

USFRus1975 <- USesfr[USesfr[,1] == 1975,]
ESFRus1975 <- ESesfr[ESesfr[,1] == 1975,]
pdf("/home/triffe/git/DISS/latex/Figures/eSFRharmonic.pdf", height = 5, width = 5)
par(mai = c(.3, .3, .3, .1), xaxs = "i", yaxs = "i", mfrow = c(2,2))
plot(y, USFRus1975[,"ESFRstm"], type = 'l', ylim = c(0, .1), xlim = c(0,111), axes = FALSE,
        col = gray(.2), lwd = 1.5, xlab = "", ylab = "", main = "US, 1975",
        panel.first = list(rect(0,0,111,0.1,col = gray(.95), border=NA),
                abline(h = seq(0,.1,by = .01), col = "white"),
                abline(v = seq(0, 110, by = 10), col = "white"),
                text(seq(0, 110, by = 10), 0,seq(0, 110, by = 10), pos = 1, cex = .7, xpd = TRUE),
                text(0,seq(0, .1, by = .01), seq(0, .1, by = .01), pos = 2, cex = .7, xpd = TRUE),
                text(55, -.007, expression(e[y]), cex = .8, pos = 1, xpd = TRUE),
                text(-19,.11, "Fertility Rate", cex = .8, xpd = TRUE, pos = 4)))
lines(y, USFRus1975[,"ESFRm"], lwd = 1, col = gray(.15), lty = 5)
lines(y, USFRus1975[,"ESFRstf"], lwd = 2, col = gray(.5), lty = 1)
lines(y, USFRus1975[,"ESFRf"], lwd = 1.5, col = gray(.3), lty = 5)
legend(58,.1, lty = c(1,5,1,5), col = gray(c(.2,.15,.5,.3)), lwd = c(1.5,1,2,1.5),bty = "n",
        legend = c("stable males", "initial males", "stable females", "initial females"), xpd = TRUE, cex = .7)
#dev.off()

par(mai = c(.3, .3, .3, .1), xaxs = "i", yaxs = "i")
plot(y, ESFRus1975[,"ESFRstm"], type = 'l', ylim = c(0, .1), xlim = c(0,111), axes = FALSE,
        col = gray(.2), lwd = 1.5, xlab = "", ylab = "", main = "Spain, 1975",
        panel.first = list(rect(0,0,111,0.1,col = gray(.95), border=NA),
                abline(h = seq(0,.1,by = .01), col = "white"),
                abline(v = seq(0, 110, by = 10), col = "white"),
                text(seq(0, 110, by = 10), 0,seq(0, 110, by = 10), pos = 1, cex = .7, xpd = TRUE),
                text(0,seq(0, .1, by = .01), seq(0, .1, by = .01), pos = 2, cex = .7, xpd = TRUE),
                text(55, -.007, expression(e[y]), cex = .8, pos = 1, xpd = TRUE),
                text(-19,.11, "Fertility Rate", cex = .8, xpd = TRUE, pos = 4)))
lines(y, ESFRus1975[,"ESFRm"], lwd = 1, col = gray(.15), lty = 5)
lines(y, ESFRus1975[,"ESFRstf"], lwd = 2, col = gray(.5), lty = 1)
lines(y, ESFRus1975[,"ESFRf"], lwd = 1.5, col = gray(.3), lty = 5)
legend(58,.1, lty = c(1,5,1,5), col = gray(c(.2,.15,.5,.3)), lwd = c(1.5,1,2,1.5),bty = "n",
        legend = c("stable males", "initial males", "stable females", "initial females"), xpd = TRUE, cex = .7)

# 2009
USFRus2009 <- USesfr[USesfr[,1] == 2009,]
ESFRus2009 <- ESesfr[ESesfr[,1] == 2009,]
par(mai = c(.3, .3, .3, .1), xaxs = "i", yaxs = "i")
plot(y, USFRus2009[,"ESFRstm"], type = 'l', ylim = c(0, .1), xlim = c(0,111), axes = FALSE,
        col = gray(.2), lwd = 1.5, xlab = "", ylab = "", main = "US, 2009",
        panel.first = list(rect(0,0,111,0.1,col = gray(.95), border=NA),
                abline(h = seq(0,.1,by = .01), col = "white"),
                abline(v = seq(0, 110, by = 10), col = "white"),
                text(seq(0, 110, by = 10), 0,seq(0, 110, by = 10), pos = 1, cex = .7, xpd = TRUE),
                text(0,seq(0, .1, by = .01), seq(0, .1, by = .01), pos = 2, cex = .7, xpd = TRUE),
                text(55, -.007, expression(e[y]), cex = .8, pos = 1, xpd = TRUE),
                text(-19,.11, "Fertility Rate", cex = .8, xpd = TRUE, pos = 4)))
lines(y, USFRus2009[,"ESFRm"], lwd = 1, col = gray(.15), lty = 5)
lines(y, USFRus2009[,"ESFRstf"], lwd = 2, col = gray(.5), lty = 1)
lines(y, USFRus2009[,"ESFRf"], lwd = 1.5, col = gray(.3), lty = 5)
legend(58,.1, lty = c(1,5,1,5), col = gray(c(.2,.15,.5,.3)), lwd = c(1.5,1,2,1.5),bty = "n",
        legend = c("stable males", "initial males", "stable females", "initial females"), xpd = TRUE, cex = .7)
#dev.off()

par(mai = c(.3, .3, .3, .1), xaxs = "i", yaxs = "i")
plot(y, ESFRus2009[,"ESFRstm"], type = 'l', ylim = c(0, .1), xlim = c(0,111), axes = FALSE,
        col = gray(.2), lwd = 1.5, xlab = "", ylab = "", main = "Spain, 2009",
        panel.first = list(rect(0,0,111,0.1,col = gray(.95), border=NA),
                abline(h = seq(0,.1,by = .01), col = "white"),
                abline(v = seq(0, 110, by = 10), col = "white"),
                text(seq(0, 110, by = 10), 0,seq(0, 110, by = 10), pos = 1, cex = .7, xpd = TRUE),
                text(0,seq(0, .1, by = .01), seq(0, .1, by = .01), pos = 2, cex = .7, xpd = TRUE),
                text(55, -.007, expression(e[y]), cex = .8, pos = 1, xpd = TRUE),
                text(-19,.11, "Fertility Rate", cex = .8, xpd = TRUE, pos = 4)))
lines(y, ESFRus2009[,"ESFRm"], lwd = 1, col = gray(.15), lty = 5)
lines(y, ESFRus2009[,"ESFRstf"], lwd = 2, col = gray(.5), lty = 1)
lines(y, ESFRus2009[,"ESFRf"], lwd = 1.5, col = gray(.3), lty = 5)
legend(58,.1, lty = c(1,5,1,5), col = gray(c(.2,.15,.5,.3)), lwd = c(1.5,1,2,1.5),bty = "n",
        legend = c("stable males", "initial males", "stable females", "initial females"), xpd = TRUE, cex = .7)
dev.off()

# diff coefs



pdf("/home/triffe/git/DISS/latex/Figures/exhmESFRdiffcoef.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, UScomphm[,"mdiffcoef"], type = 'l', ylim = c(0, .05), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 1.5, xlab = "", ylab = "",
        panel.first = list(rect(1968,0,2010,.05,col = gray(.95), border=NA),
                abline(h = seq(0,.05,by = .005), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(seq(1970, 2010, by = 10), 0,seq(1970, 2010, by = 10), pos = 1, cex = .7, xpd = TRUE),
                text(1968,seq(0, .05, by = .01), seq(0, .05, by = .01), pos = 2, cex = .7, xpd = TRUE),
                text(1990, -.003, "Year", cex = .8, pos = 1, xpd = TRUE),
                text(1963,.053, expression(theta~eSFR), cex = .8, xpd = TRUE, pos = 4)))
lines(yearsUS, UScomphm[,"fdiffcoef"], lwd = 2, col = gray(.4), lty = 1)
lines(yearsES, EScomphm[,"mdiffcoef"], lwd = 1, col = gray(.2), lty = 5)
lines(yearsES, EScomphm[,"fdiffcoef"], lwd = 2, col = gray(.4), lty = 5)
legend("topleft", lty = c(1,1,5,5), col = gray(c(.2,.4,.2,.4)), lwd = c(1,2,1,2),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
#
dev.off()



# end block
}

# ----------------------------
# transient dynamics?

do.transient <- FALSE
if (do.transient){
    
    
    # think about it
    
    
}





# --------------------------------------------------
# competition test:
# 1975
run.comp.code <- FALSE
if (run.comp.code){
yr  <- "1990"
dxm <- dxmUS[, yr]
dxf <- dxfUS[, yr]
#
BxyF <- ExpectedDxMxFmatrix(BxymfUS[[yr]][["Bxyf"]], dxm, dxf)
BxyM <- ExpectedDxMxFmatrix(BxymfUS[[yr]][["Bxym"]], dxm, dxf)
#
exm1 <- with(ExUS,Male[Year == as.integer(yr)])
exf1 <- with(ExUS,Female[Year == as.integer(yr)])
ExM1 <- rowSums(ExpectedDx(exm1, dxm))
ExF1 <- rowSums(ExpectedDx(exf1, dxf))
## get harmonic rates for boys, girls
FH <- (BxyF + BxyM) / outer(ExM1, ExF1, HM)

exm2 <- exm1
exm2[31] <- exm2[31] * 5
ExM2 <- rowSums(ExpectedDx(exm2, dxm))
ExF2 <- rowSums(ExpectedDx(exf1, dxf))

Bxy1 <- FH * outer(ExM1, ExF1, HM)
Bxy2 <- FH * outer(ExM2, ExF2, HM)

Fxm1 <- rowSums(Bxy1) / ExF1
Fxm2 <- rowSums(Bxy2) / ExF2
Fxf1 <- colSums(Bxy1) / ExF1
Fxf2 <- colSums(Bxy2) / ExF2

plot(0:110, Fxm1, type = 'l', col = "blue",lty=2, ylim=c(0,.08))
lines(0:110, Fxm2, col = "blue")
lines(0:110, Fxf1, col = "red", lty=2)
lines(0:110, Fxf2, col = "red")

sum(Fxm1)
sum(Fxm2)

ExM1.2 <- ExpectedDx(exm1, dxm)
ExM2.2 <- ExpectedDx(exm2, dxm)
Fxplit1 <- Fxm1 * (ExM1.2 / rowSums(ExM1.2))
Fxplit2 <- Fxm2 * (ExM2.2 / rowSums(ExM2.2))

fields::image.plot(Fxplit1-Fxplit2, zlim = c(-.0005,.0005))
}