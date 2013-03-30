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

Expect <- copiler::cmpfun(function(m1,m2){
    Minf0(Mna0(outer(m1,m2, "*") / (sum(m1+m2)/2)))
})
RatioAdj <- compiler::cmpfun(function(Ratio, BxyExp){
    (Ratio * BxyExp) * Minf0(Mna0((sum(BxyExp) / sum(RatioM * BxyExp))))
})

LotkaCrossRatioCoaleex <- compiler::cmpfun(function(BxyM, BxyF, 
                Exm, Exf, 
                dxm, dxf, 
                RatioM, RatioF,
                .a = .5:110.5, maxit = 2e2, tol = 1e-15){
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
            
            FexMM <- Minf0(Mna0(BxyM / Exm))
            FexMF <- Minf0(Mna0(BxyF / Exm))
            FexFF <- Minf0(Mna0(t(t(BxyF) / Exf)))
            FexFM <- Minf0(Mna0(t(t(BxyM) / Exf)))
           
            SRBi  <- sum( RatioAdj(RatioM, 
                                   Expect(rowSums(FexMM) * rowSums(dxM),
                                          colSums(FexFM) * rowSums(dxF)))) / 
                       sum(RatioAdj(RatioF, 
                                   Expect(rowSums(FexMF) * rowSums(dxM),
                                          colSums(FexFF) * rowSums(dxF))))
            pmi   <- SRBi / (1 + SRBi)
            pfi   <- 1 / (1 + SRBi)
            R0    <-  sum( RatioAdj(RatioM, 
                                    Expect(rowSums(FexMM) * rowSums(dxM) * pmi,
                                            colSums(FexFM) * rowSums(dxF) * pfi))) +
                    sum(RatioAdj(RatioF, 
                                    Expect(rowSums(FexMF) * rowSums(dxM) * pmi,
                                            colSums(FexFF) * rowSums(dxF) * pfi)))
           r2     <- log(R0) / 60
           
           for(i in 1:maxit){
               r1 <- r2
               deltai <- 1 - (sum( RatioAdj(RatioM, 
                                       Expect(rowSums(FexMM) * rowSums(dxM * exp(-r1 * .a)) * pmi,
                                               colSums(FexFM) * rowSums(dxF * exp(-r1 * .a)) * pfi))) +
                             sum(RatioAdj(RatioF, 
                                       Expect(rowSums(FexMF) * rowSums(dxM * exp(-r1 * .a)) * pmi,
                                               colSums(FexFF) * rowSums(dxF * exp(-r1 * .a)) * pfi))))
               r2    <- r1 - (deltai / (60 - (deltai / r1)))
               SRBi  <- sum( RatioAdj(RatioM, 
                                           Expect(rowSums(FexMM) * rowSums(dxM * exp(-r2 * .a)) * pmi,
                                                   colSums(FexFM) * rowSums(dxF * exp(-r2 * .a)) * pfi))) /
                        sum(RatioAdj(RatioF, 
                                           Expect(rowSums(FexMF) * rowSums(dxM * exp(-r2 * .a)) * pmi,
                                                   colSums(FexFF) * rowSums(dxF * exp(-r2 * .a)) * pfi)))
               pmi <- SRBi / (1 + SRBi)
               pfi <- 1 / (1 + SRBi)
               if (abs(deltai)<tol){
                   break
               }
           }
           c(r = r2, SRB = SRBi)
        })

rUScp <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxy, .dxm, .dxf, .Ex){
                    
                    dxm <- .dxm[,yr]
                    dxf <- .dxf[,yr]
                    Exm <- rowSums(ExpectedDx(with(.Ex, Male[Year == as.integer(yr)]), dxm))
                    Exf <- rowSums(ExpectedDx(with(.Ex, Female[Year == as.integer(yr)]), dxf))
                   
                    BxyM    <- ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxym"]], dxm, dxf) 
                    BxyF    <- ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxyf"]], dxm, dxf) 
                    RatioM  <- Minf0(Mna0(BxyM / (outer(rowSums(BxyM), colSums(BxyM)) / sum(BxyM))))
                    RatioF  <- Minf0(Mna0(BxyF / (outer(rowSums(BxyF), colSums(BxyF)) / sum(BxyF))))
                    
                    LotkaCrossRatioCoaleex(BxyM, BxyF, Exm, Exf, dxm, dxf, RatioM, RatioF)
                   
                }, .Bxy = BxymfUS, .dxm = dxmUS, .dxf = dxfUS, .Ex = ExUS))
rownames(rUScp) <- yearsUS
rEScp <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxy, .dxm, .dxf, .Ex){
                    
                    dxm <- .dxm[,yr]
                    dxf <- .dxf[,yr]
                    Exm <- rowSums(ExpectedDx(with(.Ex, Male[Year == as.integer(yr)]), dxm))
                    Exf <- rowSums(ExpectedDx(with(.Ex, Female[Year == as.integer(yr)]), dxf))
                    
                    BxyM    <- ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxym"]], dxm, dxf) 
                    BxyF    <- ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxyf"]], dxm, dxf) 
                    RatioM  <- Minf0(Mna0(BxyM / (outer(rowSums(BxyM), colSums(BxyM)) / sum(BxyM))))
                    RatioF  <- Minf0(Mna0(BxyF / (outer(rowSums(BxyF), colSums(BxyF)) / sum(BxyF))))
                    
                    LotkaCrossRatioCoaleex(BxyM, BxyF, Exm, Exf, dxm, dxf, RatioM, RatioF)
                    
                }, .Bxy = BxymfES, .dxm = dxmES, .dxf = dxfES, .Ex = ExES))
rownames(rEScp) <- yearsES
plot(yearsUS, rUScp[,1], type = 'l', ylim = c(-.015,.01))
lines(yearsUS, rUS[,1], lty = 2)
lines(yearsES, rEScp[,1], col = "red")
lines(yearsES, rES[,1], lty = 2, col = "red")





















