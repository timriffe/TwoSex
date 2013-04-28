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
                # need to divide by two
                delta.i <- (2 - sum(p.m * colSums(exp(-r.i * .a) * t(dxM)) * (Fxpredm[[1]]+Fxpredf[[1]]) + 
                              p.f * colSums(exp(-r.i * .a) * t(dxF)) * (Fxpredm[[2]]+Fxpredf[[2]]))) / 2
                r.i <- r.i - (delta.i / ( T.hat  - (delta.i / r.i)))
                
                SRB.i <- sum(p.m * colSums(exp(-r.i * .a) * t(dxM)) * Fxpredm[[1]] + 
                                        p.f * exp(-r.i * .a) * colSums(exp(-r.i * .a) * t(dxF)) * Fxpredm[[2]]) / 
                        sum(p.m * colSums(exp(-r.i * .a) * t(dxM)) * Fxpredf[[1]] + 
                                        p.f * colSums(exp(-r.i * .a) * t(dxF)) * Fxpredf[[2]])
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


plot(yearsES, rES[,1], col = "red", type = 'l', xlim = range(yearsUS))
lines(yearsUS, rUS[,1], col = "blue")
        abline(h=0)
        