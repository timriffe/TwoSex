
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
LotkaHminex <- compiler::cmpfun(function(r, dxm, dxf, FHf, FHm, .a = .5:110.5){
          
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
            
            (1 - sum(outer(rowSums(dxM %col% (1 / exp(-r1 * .a))), 
                     rowSums(dxF %col% (1 / exp(-r1 * .a))),
                     HM) * (FHf + FHm)))^2
                          
            
        })

#optimize(LotkaHmin, interval = c(-.15,.15), Lxm = LxmUS[,10], Lxf = LxfUS[,10], 
#        FH = BxUS[[10]] /outer(with(ExUS, Male[Year == yearsUS[10]]), with(ExUS, Female[Year == yearsUS[10]]), HM), 
#        tol = 1e-15)$minimum

# using strategy of Coale
LotkaHRCoaleex <- compiler::cmpfun(function(FHf, FHm, dxm, dxf, .a = .5:110.5, maxit = 2e2, tol = 1e-15){
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
               
            r1 <- 0
            SRBi        <- sum(outer(rowSums(dxM), rowSums(dxF), HM) * FHm) / 
                    sum(outer(rowSums(dxM), rowSums(dxF), HM) * FHf)
            R0          <- sum(outer(rowSums(dxM) * SRBi / (SRBi + 1), rowSums(dxF) * 1 / (SRBi + 1), HM) * (FHf + FHm))
            # first assuming a mean generation time of 35
            r2 <- log(R0) / 60
            for (i in 1:maxit){ # 10 is more than enough!
                r1      <- r2
                p.m     <- SRBi / (SRBi+1)
                p.f     <- 1 / (SRBi+1)
                deltai  <- 1 - sum(outer(rowSums(dxM %col% (1 / exp(-r1 * .a))) * p.m, 
                                         rowSums(dxF %col% (1 / exp(-r1 * .a)))* p.f, HM) * (FHf + FHm))
                # the mean generation time self-corrects 
                # according to the error produced by the Lotka equation
                r2      <- r1 - (deltai / (60 - (deltai / r1)))
                
                SRBi    <-  sum(outer(rowSums(dxM %col% (1 / exp(-r2 * .a))) * p.m, 
                                        rowSums(dxF %col% (1 / exp(-r2 * .a)))* p.f, HM) * FHm) / 
                            sum(outer(rowSums(dxM %col% (1 / exp(-r2 * .a))) * p.m, 
                                        rowSums(dxF %col% (1 / exp(-r2 * .a)))* p.f, HM) * FHf)
                if (abs(r2 - r1) <= tol | zapsmall(abs(deltai)) <= tol){
                    break
                }
            }
            return(c(r2, SRBi))
        })
        
        yr <- "1975"
rUS <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxy, .dxm, .dxf, .Ex){
                    
                    Exm <- rowSums(ExpectedDx(with(.Ex, Male[Year == as.integer(yr)]), .dxm[,yr]))
                    Exf <- rowSums(ExpectedDx(with(.Ex, Female[Year == as.integer(yr)]), .dxf[,yr]))
                    # harmonic mean of exposures
                    Hxy <- outer(Exm, Exf, HM)
                    # harmonic rates (divide by two since working with both sexes)
                    FxyHf <- Mna0(Minf0(ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxyf"]], .dxm[, yr], .dxf[, yr]) / Hxy))
                    FxyHm <- Mna0(Minf0(ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxym"]], .dxm[, yr], .dxf[, yr]) / Hxy))
                                      
                    LH <- LotkaHRCoaleex(FHf = FxyHf, FHm = FxyHm, dxm =.dxm[,yr], dxf = .dxf[,yr])
                    c(r.H = LH[1], SRBH = LH[2])
                }, .Bxy = BxymfUS, .dxm = dxmUS, .dxf = dxfUS, .Ex = ExUS))
rownames(rUS) <- yearsUS
rES <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxy, .dxm, .dxf, .Ex){
                    
                    Exm <- rowSums(ExpectedDx(with(.Ex, Male[Year == as.integer(yr)]), .dxm[,yr]))
                    Exf <- rowSums(ExpectedDx(with(.Ex, Female[Year == as.integer(yr)]), .dxf[,yr]))
                    # harmonic mean of exposures
                    Hxy <- outer(Exm, Exf, HM)
                    # harmonic rates (divide by two since working with both sexes)
                    FxyHf <- Mna0(Minf0(ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxyf"]], .dxm[, yr], .dxf[, yr]) / Hxy))
                    FxyHm <- Mna0(Minf0(ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxym"]], .dxm[, yr], .dxf[, yr]) / Hxy))
                    
                    LH <- LotkaHRCoaleex(FHf = FxyHf, FHm = FxyHm, dxm =.dxm[,yr], dxf = .dxf[,yr])
                    c(r.H = LH[1], SRBH = LH[2])
                }, .Bxy = BxymfES, .dxm = dxmES, .dxf = dxfES, .Ex = ExES))
rownames(rES) <- yearsES


plot(yearsUS, rUS[,1], type = 'l',ylim = c(-.015,.01))
lines(yearsES, rES[,1], col = "red")


