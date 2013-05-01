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


IPFpred <- compiler::cmpfun(function(Bxy, Exm1, Exm2, Exf1, Exf2, marM = mean, tol = 1e-15, maxit = 200){
            
            # starting rates
            FxyF        <- Minf0(Mna0(colSums(Bxy) / Exf1))
            FxyM        <- Minf0(Mna0(rowSums(Bxy) / Exm1))
            # predicted counts
            Mpred       <- FxyM * Exm2
            Fpred       <- FxyF * Exf2
            # starting vals    
            Bxy2        <- Bxy1 <- Bxy
            
            # sums will differ
            Msum        <- sum(Mpred)
            Fsum        <- sum(Fpred)
            # take a mean of the male and female marginal totals
            Nsum        <- marM(c(Msum, Fsum))
            # rescale so that marginals match
            Mpredrsc    <- Mpred * (Nsum / Msum)
            Fpredrsc    <- Fpred * (Nsum / Fsum)
            # get a startin value for comparisons (assoc-free)
            BxyB        <- outer(Mpredrsc,Fpredrsc,"*") / Nsum
            for (i in maxit){
                # rows then cols
                Bxy1 <- Bxy1 * Minf0(Mna0(Mpredrsc / rowSums(Bxy1)))
                Bxy1 <- t(t(Bxy1) * Minf0(Mna0(Fpredrsc / colSums(Bxy1))))
                
                # cols then rows
                Bxy2 <- t(t(Bxy2) * Minf0(Mna0(Fpredrsc / colSums(Bxy2))))
                Bxy2 <- Bxy2 * Minf0(Mna0(Mpredrsc / rowSums(Bxy2)))
                
                # mean for next it
                BxyA <- (Bxy1 + Bxy2) / 2
                if (sum(abs(BxyA-BxyB))< tol){
                    break
                }
                BxyB <-Bxy2 <- Bxy1 <- BxyA
            }
            FxmPred <- Minf0(Mna0(rowSums(BxyA) / Exm2))
            FxfPred <- Minf0(Mna0(colSums(BxyA) / Exf2))
            list(FxmPred = FxmPred, FxfPred = FxfPred, i=i)
        })
        
rexIPFit <- compiler::cmpfun(function(Bxym, Bxyf, Exm, Exf, dxm, dxf, M = mean,
                .a = .5:110.5, maxit = 2e2, tol = 1e-15){
            
            # get dx matrices for zero growth
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
            # get starting values:
            
            SRB <- sum(Bxym) / sum(Bxyf)
            p.m <- SRB / (1 + SRB)
            p.f <- 1 / (1 + SRB)
            Bxy    <- Bxym + Bxyf
            Fxpred.hat <- IPFpred(Bxy, 
                    Exm1 = Exm, 
                    Exm2 = rowSums(dxM), # i.e. assuming r = 0
                    Exf1 = Exf, 
                    Exf2 = rowSums(dxF),
                    marM = M)
            R0.hat <- sum(p.m * rowSums(dxM) * Fxpred.hat[[1]] + p.f * rowSums(dxF) * Fxpred.hat[[2]]) / 2
            T.hat <- (sum(.a * p.m * rowSums(dxM) * Fxpred.hat[[1]] + .a * p.f * rowSums(dxF) * Fxpred.hat[[2]]) / 2) / R0.hat
            
            r.i <- log(R0.hat) / T.hat
            
            for (i in 1:maxit){
                # get adjusted rates for the given r.i
                Fxpredm <- IPFpred(Bxym, 
                        Exm1 = Exm, 
                        Exm2 = p.m * colSums(exp(-r.i * .a) * t(dxM)),  # now we scale for generation size
                        Exf1 = Exf, 
                        Exf2 = p.f * colSums(exp(-r.i * .a) * t(dxF)),
                        marM = M)
                Fxpredf <- IPFpred(Bxyf, 
                        Exm1 = Exm, 
                        Exm2 = p.m * colSums(exp(-r.i * .a) * t(dxM)), 
                        Exf1 = Exf, 
                        Exf2 = p.f * colSums(exp(-r.i * .a) * t(dxF)),
                        marM = M)
                # get residual
                # need to divide by two
                delta.i <- (2 - sum(p.m * colSums(exp(-r.i * .a) * t(dxM)) * (Fxpredm[[1]]+Fxpredf[[1]]) + 
                              p.f * colSums(exp(-r.i * .a) * t(dxF)) * (Fxpredm[[2]]+Fxpredf[[2]]))) / 2
                # update r
                r.i <- r.i - (delta.i / ( T.hat  - (delta.i / r.i)))
                # update SRB estimate. Same as above but for boy and girl births separately. take ratio as such
                # no need to re-update fertility prior to this. wont' speed things any
                SRB.i <- sum(p.m * colSums(exp(-r.i * .a) * t(dxM)) * Fxpredm[[1]] + 
                                        p.f * exp(-r.i * .a) * colSums(exp(-r.i * .a) * t(dxF)) * Fxpredm[[2]]) / 
                        sum(p.m * colSums(exp(-r.i * .a) * t(dxM)) * Fxpredf[[1]] + 
                                        p.f * colSums(exp(-r.i * .a) * t(dxF)) * Fxpredf[[2]])
                # convertto proportions male and female
                p.m <- (SRB.i / (1 + SRB.i))
                p.f <- (1 / (1 + SRB.i))
                
                if (abs(delta.i) <= tol){
                    break
                }
                
            }
            
            c(r = r.i, SRB = SRB.i, iter = i)
        })    

rUS <- do.call(rbind, lapply(as.character(yearsUS), function(yr, .Bxymf, .Ex, .dxm, .dxf, .M){
                    
                    Bxym <- ExpectedDxMxFmatrix(.Bxymf[[yr]][["Bxym"]], .dxm[,yr], .dxf[,yr])
                    Bxyf <- ExpectedDxMxFmatrix(.Bxymf[[yr]][["Bxyf"]], .dxm[,yr], .dxf[,yr])
                    Exm  <- rowSums(ExpectedDx(with(.Ex,Male[Year == as.integer(yr)]), .dxm[,yr]))
                    Exf  <- rowSums(ExpectedDx(with(.Ex,Female[Year == as.integer(yr)]), .dxf[,yr]))
                    rexIPFit(Bxym = Bxym, Bxyf = Bxyf, 
                             Exm = Exm, Exf = Exf, 
                             dxm = .dxm[,yr], dxf =  .dxf[,yr], M = .M)
                },.Bxymf = BxymfUS, .Ex = ExUS, .dxm = dxmUS, .dxf = dxfUS, .M = mean))
rES <- do.call(rbind, lapply(as.character(yearsES), function(yr, .Bxymf, .Ex, .dxm, .dxf, .M){
                    
                    Bxym <- ExpectedDxMxFmatrix(.Bxymf[[yr]][["Bxym"]], .dxm[,yr], .dxf[,yr])
                    Bxyf <- ExpectedDxMxFmatrix(.Bxymf[[yr]][["Bxyf"]], .dxm[,yr], .dxf[,yr])
                    Exm  <- rowSums(ExpectedDx(with(.Ex,Male[Year == as.integer(yr)]), .dxm[,yr]))
                    Exf  <- rowSums(ExpectedDx(with(.Ex,Female[Year == as.integer(yr)]), .dxf[,yr]))
                    rexIPFit(Bxym = Bxym, Bxyf = Bxyf, 
                            Exm = Exm, Exf = Exf, 
                            dxm = .dxm[,yr], dxf =  .dxf[,yr], M = .M)
                },.Bxymf = BxymfES, .Ex = ExES, .dxm = dxmES, .dxf = dxfES, .M = mean))

rownames(rUS) <- yearsUS
rownames(rES) <- yearsES

rUSipfhm <- do.call(rbind, lapply(as.character(yearsUS), function(yr, .Bxymf, .Ex, .dxm, .dxf, .M){
                    
                    Bxym <- ExpectedDxMxFmatrix(.Bxymf[[yr]][["Bxym"]], .dxm[,yr], .dxf[,yr])
                    Bxyf <- ExpectedDxMxFmatrix(.Bxymf[[yr]][["Bxyf"]], .dxm[,yr], .dxf[,yr])
                    Exm  <- rowSums(ExpectedDx(with(.Ex,Male[Year == as.integer(yr)]), .dxm[,yr]))
                    Exf  <- rowSums(ExpectedDx(with(.Ex,Female[Year == as.integer(yr)]), .dxf[,yr]))
                    rexIPFit(Bxym = Bxym, Bxyf = Bxyf, 
                            Exm = Exm, Exf = Exf, 
                            dxm = .dxm[,yr], dxf =  .dxf[,yr], M = .M)
                },.Bxymf = BxymfUS, .Ex = ExUS, .dxm = dxmUS, .dxf = dxfUS, .M = harmonic.mean))
rESipfhm <- do.call(rbind, lapply(as.character(yearsES), function(yr, .Bxymf, .Ex, .dxm, .dxf, .M){
                    
                    Bxym <- ExpectedDxMxFmatrix(.Bxymf[[yr]][["Bxym"]], .dxm[,yr], .dxf[,yr])
                    Bxyf <- ExpectedDxMxFmatrix(.Bxymf[[yr]][["Bxyf"]], .dxm[,yr], .dxf[,yr])
                    Exm  <- rowSums(ExpectedDx(with(.Ex,Male[Year == as.integer(yr)]), .dxm[,yr]))
                    Exf  <- rowSums(ExpectedDx(with(.Ex,Female[Year == as.integer(yr)]), .dxf[,yr]))
                    rexIPFit(Bxym = Bxym, Bxyf = Bxyf, 
                            Exm = Exm, Exf = Exf, 
                            dxm = .dxm[,yr], dxf =  .dxf[,yr], M = .M)
                },.Bxymf = BxymfES, .Ex = ExES, .dxm = dxmES, .dxf = dxfES, .M = harmonic.mean))
rownames(rUSipfhm) <- yearsUS
rownames(rESipfhm) <- yearsES
#plot((rUSipfhm[,1] - rUS[,1]) / ((rUSipfhm[,1] + rUS[,1])/2),type = 'l')
#plot(rUSipfhm[,1])
#lines(rUS[,1])
# this will take a minute:
# source("/home/triffe/git/DISS/R/ExLotka2SexLinear.R") # upper and lower bounds are dominance weighted 0 and 1,
# which are identical to the single-sex r estimates.
# produce the dominance weighted r estimates for comparison.
# plot(rESipfhm[,1] / rES[,1], ylim = c(.8,1.3))
# in order, 0, .5, 1 for sigmas
#plot(rUS[, 1] / US[, 2] )
#plot(US[,3] - rmUS[,1]) # rmUS came from exLotka1Sex.R

make.fig <- FALSE
if (make.fig){
pdf("/home/triffe/git/DISS/latex/Figures/exIPFr.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rUSipfhm[, 1], type = 'n', ylim = c(-.016,.01),xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,-.016,2010,.01,col = gray(.95), border=NA),
                abline(h = seq(-.016,.01,by = .002), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.016,.01,by = .002),seq(-.016,.01,by = .002), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.016, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.0175, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,.0115, "r", cex = 1, xpd = TRUE)))
polygon(c(yearsUS, rev(yearsUS)), c(US[, 1], rev(US[, 3])), col = "#99999950", border = NA)
lines(yearsUS, US[, 3],col = gray(.2), lwd = 1, lty = 5)
lines(yearsUS, US[, 1],col = gray(.2), lwd = 1, lty = 5)
lines(yearsUS, rUSipfhm[, 1],col = gray(.2), lwd = 2)
lines(yearsUS, rUS[, 1],col = gray(.1), lwd = 1)
#lines(yearsUS,rUShm[,1],col="red")
polygon(c(yearsES, rev(yearsES)), c(ES[, 1], rev(ES[, 3])), col = "#99999950", border = NA)
lines(yearsES, ES[, 3],col = gray(.2), lwd = 1, lty = 5)
lines(yearsES, ES[, 1],col = gray(.2), lwd = 1, lty = 5)
lines(yearsES, rESipfhm[, 1],col = gray(.2), lwd = 2, lty = 1)
lines(yearsES, rES[, 1],col = gray(.1), lwd = 1)
#lines(yearsES,rEShm[,1],col="red")
text(c(1995, 1992.372, 1995, 1988, 1983.303, 1988, 1983.303),
        c(0.0038379231,  0.0001698981,  0.0018907000, -0.0039962537, -0.0069850148,-0.0059887611,-0.010019060),
        c(expression(r^m~US),expression(r^f~US),expression(r^IPF~US),expression(r^m~ES),
        expression(r^f~ES),expression(r^{IPF(ar)}~ES), expression(r^{IPF(hm)~ES})), pos = c(4,1,4,4,1,4,1))
segments(c(1995.5,1995.349,1989,1988.564,1993.064,1983.787,1984.549),
        c(0.0037473546,0.0017095629,-0.0041321065,-0.0061698981,-0.0001018074,-0.0073020046,-0.010154913),
        c(1992.926,1994.588,1985.934,1986.764,1994.103,1985.034,1987.249),
        c( 0.0026152481,0.0009397305,-0.0051736444,-0.0072114361,0.0006227407,-0.0067585935,-0.008479395))
dev.off()

}

# ESFR stable:
do.esfr.ipf <- FALSE
if (do.esfr.ipf){
rUSesfr <- do.call(rbind, lapply(as.character(yearsUS), function(yr, .Bxymf, .Ex, .dxm, .dxf, .rSRB, .M, .a = .5:110.5){
                    # stable pop parameters
                    r               <- .rSRB[yr, 1]
                    SRB             <- .rSRB[yr, 2]
                    p.m             <- SRB / (1 + SRB)
                    p.f             <- 1 / (1 + SRB)
                    
                    dxm             <- .dxm[,yr]
                    dxf             <- .dxf[,yr]
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
                    # stable pops (not scaled to sum to 1)
                    Exm2 <- p.m * colSums(t(dxM) * exp(-r*.a))
                    Exf2 <- p.f * colSums(t(dxF) * exp(-r*.a))
   
                    #
                    Bxy    <- ExpectedDxMxFmatrix(.Bxymf[[yr]][["Bxym"]]+.Bxymf[[yr]][["Bxyf"]], .dxm[,yr], .dxf[,yr])
                  
                    Exm     <- rowSums(ExpectedDx(with(.Ex,Male[Year == as.integer(yr)]), .dxm[,yr]))
                    Exf     <- rowSums(ExpectedDx(with(.Ex,Female[Year == as.integer(yr)]), .dxf[,yr]))
                    
                    Fxm    <- Minf0(Mna0(rowSums(Bxy) / Exm))
                    Fxf    <- Minf0(Mna0(colSums(Bxy) / Exf))
                    esfrs <- IPFpred(Bxy = Bxy, 
                            Exm1 = Exm, Exm2 = Exm2, Exf1 = Exf, Exf2 = Exf2, marM = .M)
                    cbind(Year = as.integer(yr), Mst = esfrs$FxmPred, Fst = esfrs$FxfPred, Minit = Fxm, Finit = Fxf)
                },.Bxymf = BxymfUS, .Ex = ExUS, .dxm = dxmUS, .dxf = dxfUS, .rSRB = rUSipfhm,.M = mean))
rESesfr <- do.call(rbind, lapply(as.character(yearsES), function(yr, .Bxymf, .Ex, .dxm, .dxf, .rSRB, .M, .a = .5:110.5){
                    # stable pop parameters
                    r               <- .rSRB[yr, 1]
                    SRB             <- .rSRB[yr, 2]
                    p.m             <- SRB / (1 + SRB)
                    p.f             <- 1 / (1 + SRB)
                    
                    dxm             <- .dxm[,yr]
                    dxf             <- .dxf[,yr]
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
                    # stable pops (not scaled to sum to 1)
                    Exm2 <- p.m * colSums(t(dxM) * exp(-r*.a))
                    Exf2 <- p.f * colSums(t(dxF) * exp(-r*.a))
                    
                    #
                    Bxy    <- ExpectedDxMxFmatrix(.Bxymf[[yr]][["Bxym"]]+.Bxymf[[yr]][["Bxyf"]], .dxm[,yr], .dxf[,yr])
           
                    Exm     <- rowSums(ExpectedDx(with(.Ex,Male[Year == as.integer(yr)]), .dxm[,yr]))
                    Exf     <- rowSums(ExpectedDx(with(.Ex,Female[Year == as.integer(yr)]), .dxf[,yr]))
                    
                    Fxm    <- Minf0(Mna0(rowSums(Bxy) / Exm))
                    Fxf    <- Minf0(Mna0(colSums(Bxy) / Exf))
                    
                    esfrs <- IPFpred(Bxy = Bxy, 
                            Exm1 = Exm, Exm2 = Exm2, Exf1 = Exf, Exf2 = Exf2, marM = .M)
                    cbind(Year = as.integer(yr), Mst = esfrs$FxmPred, Fst = esfrs$FxfPred, Minit = Fxm, Finit = Fxf)
                },.Bxymf = BxymfES, .Ex = ExES, .dxm = dxmES, .dxf = dxfES, .rSRB = rESipfhm,.M = mean))

              
UScomp <- do.call(rbind, lapply(yearsUS, function(yr, .esfr){
                    thisyr <- .esfr[.esfr[,1] == yr, 2:5]
                    TFRs <- colSums(thisyr)
                    mdiffcoef <- 1 - sum(pmin(thisyr[,1] / sum(thisyr[,1]), thisyr[,3] / sum(thisyr[,3])))
                    fdiffcoef <- 1 - sum(pmin(thisyr[,2] / sum(thisyr[,2]), thisyr[,4] / sum(thisyr[,4])))
                    c(TFRs,mdiffcoef=mdiffcoef,fdiffcoef=fdiffcoef)
                },.esfr = rUSesfr))  
EScomp <- do.call(rbind, lapply(yearsES, function(yr, .esfr){
                    thisyr <- .esfr[.esfr[,1] == yr, 2:5]
                    TFRs <- colSums(thisyr)
                    mdiffcoef <- 1 - sum(pmin(thisyr[,1] / sum(thisyr[,1]), thisyr[,3] / sum(thisyr[,3])))
                    fdiffcoef <- 1 - sum(pmin(thisyr[,2] / sum(thisyr[,2]), thisyr[,4] / sum(thisyr[,4])))
                    c(TFRs,mdiffcoef=mdiffcoef,fdiffcoef=fdiffcoef)
                },.esfr = rESesfr))  

# plot differences in TFR:
pdf("/home/triffe/git/DISS/latex/Figures/exIPFTFRdiff.pdf",height=5, width=5)
par(mai = c(.5, .5, .3, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, UScomp[,1] - UScomp[,3], type = 'l', ylim = c(-.3, .5), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968, -.3, 2010, .5,col = gray(.95), border=NA),
                abline(h = seq(-.3, .5, by = .1), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.3, .5, by = .1),round(seq(-.3, .5, by = .1),1), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10), -.3, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.35, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1967, .54, "TFR diff", cex = 1, xpd = TRUE)))
lines(yearsUS,UScomp[,2] - UScomp[,4], col = gray(.5), lwd = 2)
lines(yearsES,EScomp[,1] - EScomp[,3],col = gray(.2), lty = 5)
lines(yearsES,EScomp[,2] - EScomp[,4], col = gray(.5), lwd = 2,lty = 5)

text(rep(1990,4), c(0.1099981, -0.1030401,  0.3203396, -0.216),
        c("US males","US females","ES males","ES females"),pos = 4)
dev.off()




# -----------------------
lines(EScomp[,5])
lines(EScomp[,6])


#y <- 0:110
#yr <- 1975
#
#for (yr in yearsES){
#
#USyr <- rESesfr[rESesfr[,1]==yr, ]
#plot(y, USyr[, 2], type = 'l', ylim = c(0,.1), col = "blue", main = yr)
#lines(y, USyr[, 3], col = "red")
#lines(y, USyr[, 4], col = "blue", lty = 2)
#lines(y, USyr[, 5], col = "red", lty = 2)
#locator(1) #manual forward click
#}

USFRus1975 <- rUSesfr[rUSesfr[,1]==1975, ]
ESFRus1975 <- rESesfr[rESesfr[,1]==1975, ]
y <- 0:110
pdf("/home/triffe/git/DISS/latex/Figures/eSFRIPF.pdf", height = 5, width = 5)
par(mai = c(.3, .3, .3, .1), xaxs = "i", yaxs = "i", mfrow = c(2,2))
plot(y, USFRus1975[,"Mst"], type = 'l', ylim = c(0, .1), xlim = c(0,111), axes = FALSE,
        col = gray(.2), lwd = 1.5, xlab = "", ylab = "", main = "US, 1975",
        panel.first = list(rect(0,0,111,0.1,col = gray(.95), border=NA),
                abline(h = seq(0,.1,by = .01), col = "white"),
                abline(v = seq(0, 110, by = 10), col = "white"),
                text(seq(0, 110, by = 10), 0,seq(0, 110, by = 10), pos = 1, cex = .7, xpd = TRUE),
                text(0,seq(0, .1, by = .01), seq(0, .1, by = .01), pos = 2, cex = .7, xpd = TRUE),
                text(55, -.007, expression(e[y]), cex = .8, pos = 1, xpd = TRUE),
                text(-19,.11, "Fertility Rate", cex = .8, xpd = TRUE, pos = 4)))
lines(y, USFRus1975[,"Minit"], lwd = 1, col = gray(.15), lty = 5)
lines(y, USFRus1975[,"Fst"], lwd = 2, col = gray(.5), lty = 1)
lines(y, USFRus1975[,"Finit"], lwd = 1.5, col = gray(.3), lty = 5)
legend(58,.1, lty = c(1,5,1,5), col = gray(c(.2,.15,.5,.3)), lwd = c(1.5,1,2,1.5),bty = "n",
        legend = c("stable males", "initial males", "stable females", "initial females"), xpd = TRUE, cex = .7)
#dev.off()

par(mai = c(.3, .3, .3, .1), xaxs = "i", yaxs = "i")
plot(y, ESFRus1975[,"Mst"], type = 'l', ylim = c(0, .1), xlim = c(0,111), axes = FALSE,
        col = gray(.2), lwd = 1.5, xlab = "", ylab = "", main = "Spain, 1975",
        panel.first = list(rect(0,0,111,0.1,col = gray(.95), border=NA),
                abline(h = seq(0,.1,by = .01), col = "white"),
                abline(v = seq(0, 110, by = 10), col = "white"),
                text(seq(0, 110, by = 10), 0,seq(0, 110, by = 10), pos = 1, cex = .7, xpd = TRUE),
                text(0,seq(0, .1, by = .01), seq(0, .1, by = .01), pos = 2, cex = .7, xpd = TRUE),
                text(55, -.007, expression(e[y]), cex = .8, pos = 1, xpd = TRUE),
                text(-19,.11, "Fertility Rate", cex = .8, xpd = TRUE, pos = 4)))
lines(y, ESFRus1975[,"Minit"], lwd = 1, col = gray(.15), lty = 5)
lines(y, ESFRus1975[,"Fst"], lwd = 2, col = gray(.5), lty = 1)
lines(y, ESFRus1975[,"Finit"], lwd = 1.5, col = gray(.3), lty = 5)
legend(58,.1, lty = c(1,5,1,5), col = gray(c(.2,.15,.5,.3)), lwd = c(1.5,1,2,1.5),bty = "n",
        legend = c("stable males", "initial males", "stable females", "initial females"), xpd = TRUE, cex = .7)

# 2009
USFRus2009 <- rUSesfr[rUSesfr[,1] == 2009,]
ESFRus2009 <- rESesfr[rESesfr[,1] == 2009,]
par(mai = c(.3, .3, .3, .1), xaxs = "i", yaxs = "i")
plot(y, USFRus2009[,"Mst"], type = 'l', ylim = c(0, .1), xlim = c(0,111), axes = FALSE,
        col = gray(.2), lwd = 1.5, xlab = "", ylab = "", main = "US, 2009",
        panel.first = list(rect(0,0,111,0.1,col = gray(.95), border=NA),
                abline(h = seq(0,.1,by = .01), col = "white"),
                abline(v = seq(0, 110, by = 10), col = "white"),
                text(seq(0, 110, by = 10), 0,seq(0, 110, by = 10), pos = 1, cex = .7, xpd = TRUE),
                text(0,seq(0, .1, by = .01), seq(0, .1, by = .01), pos = 2, cex = .7, xpd = TRUE),
                text(55, -.007, expression(e[y]), cex = .8, pos = 1, xpd = TRUE),
                text(-19,.11, "Fertility Rate", cex = .8, xpd = TRUE, pos = 4)))
lines(y, USFRus2009[,"Minit"], lwd = 1, col = gray(.15), lty = 5)
lines(y, USFRus2009[,"Fst"], lwd = 2, col = gray(.5), lty = 1)
lines(y, USFRus2009[,"Finit"], lwd = 1.5, col = gray(.3), lty = 5)
legend(58,.1, lty = c(1,5,1,5), col = gray(c(.2,.15,.5,.3)), lwd = c(1.5,1,2,1.5),bty = "n",
        legend = c("stable males", "initial males", "stable females", "initial females"), xpd = TRUE, cex = .7)
#dev.off()

par(mai = c(.3, .3, .3, .1), xaxs = "i", yaxs = "i")
plot(y, ESFRus2009[,"Mst"], type = 'l', ylim = c(0, .1), xlim = c(0,111), axes = FALSE,
        col = gray(.2), lwd = 1.5, xlab = "", ylab = "", main = "Spain, 2009",
        panel.first = list(rect(0,0,111,0.1,col = gray(.95), border=NA),
                abline(h = seq(0,.1,by = .01), col = "white"),
                abline(v = seq(0, 110, by = 10), col = "white"),
                text(seq(0, 110, by = 10), 0,seq(0, 110, by = 10), pos = 1, cex = .7, xpd = TRUE),
                text(0,seq(0, .1, by = .01), seq(0, .1, by = .01), pos = 2, cex = .7, xpd = TRUE),
                text(55, -.007, expression(e[y]), cex = .8, pos = 1, xpd = TRUE),
                text(-19,.11, "Fertility Rate", cex = .8, xpd = TRUE, pos = 4)))
lines(y, ESFRus2009[,"Minit"], lwd = 1, col = gray(.15), lty = 5)
lines(y, ESFRus2009[,"Fst"], lwd = 2, col = gray(.5), lty = 1)
lines(y, ESFRus2009[,"Finit"], lwd = 1.5, col = gray(.3), lty = 5)
legend(58,.1, lty = c(1,5,1,5), col = gray(c(.2,.15,.5,.3)), lwd = c(1.5,1,2,1.5),bty = "n",
        legend = c("stable males", "initial males", "stable females", "initial females"), xpd = TRUE, cex = .7)
dev.off()


}


