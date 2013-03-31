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

Expect <- compiler::cmpfun(function(m1,m2,p=2){
    Minf0(Mna0(outer(m1,m2, "*") / stolarsky.mean(sum(m1), sum(m2), p = p)))
})
RatioAdj <- compiler::cmpfun(function(Ratio, BxyExp){
    (Ratio * BxyExp) * Minf0(Mna0((sum(BxyExp) / sum(Ratio * BxyExp))))
})

LotkaCrossRatioCoaleex <- compiler::cmpfun(function(BxyM, BxyF, 
                Exm, Exf, 
                dxm, dxf, 
                RatioM, RatioF, p = 2,
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
            
            # male and female rate matrices (male by female age), sex-specific rates.
            FexMM <- Minf0(Mna0(BxyM / Exm))      # father-son
            FexMF <- Minf0(Mna0(BxyF / Exm))      # father-daughter
            FexFF <- Minf0(Mna0(t(t(BxyF) / Exf)))# mother-daughter
            FexFM <- Minf0(Mna0(t(t(BxyM) / Exf)))# mother-son
           
            # a starting guess of the SRB, assuming r = 0 and equal weight to males and females
            SRBi  <- sum( RatioAdj(RatioM, 
                                   Expect(rowSums(FexMM) * rowSums(dxM),
                                          colSums(FexFM) * rowSums(dxF), p))) / 
                       sum(RatioAdj(RatioF, 
                                   Expect(rowSums(FexMF) * rowSums(dxM),
                                          colSums(FexFF) * rowSums(dxF), p)))
            # proportions male and female
            pmi   <- SRBi / (1 + SRBi)
            pfi   <- 1 / (1 + SRBi)
            # starting guess at NRR, assuming r = 0 and above guess at SRB
            R0    <-  sum( RatioAdj(RatioM, 
                                    Expect(rowSums(FexMM) * rowSums(dxM) * pmi,
                                            colSums(FexFM) * rowSums(dxF) * pfi, p))) +
                    sum(RatioAdj(RatioF, 
                                    Expect(rowSums(FexMF) * rowSums(dxM) * pmi,
                                            colSums(FexFF) * rowSums(dxF) * pfi, p)))
           # starting val for r, assuming above NRR and mean remaining life expectancy at reproduction of 60 years
           r2     <- log(R0) / 60
           
           for(i in 1:maxit){
               r1 <- r2
               # produce residual 
               deltai <- 1 - (sum( RatioAdj(RatioM, 
                                       Expect(rowSums(FexMM) * rowSums(dxM * exp(-r1 * .a)) * pmi,
                                               colSums(FexFM) * rowSums(dxF * exp(-r1 * .a)) * pfi, p))) +
                             sum(RatioAdj(RatioF, 
                                       Expect(rowSums(FexMF) * rowSums(dxM * exp(-r1 * .a)) * pmi,
                                               colSums(FexFF) * rowSums(dxF * exp(-r1 * .a)) * pfi, p))))
               # improve r
               r2    <- r1 - (deltai / (60 - (deltai / r1)))
               # improve SRB using improved R
               SRBi  <- sum( RatioAdj(RatioM, 
                                           Expect(rowSums(FexMM) * rowSums(dxM * exp(-r2 * .a)) * pmi,
                                                   colSums(FexFM) * rowSums(dxF * exp(-r2 * .a)) * pfi, p))) /
                        sum(RatioAdj(RatioF, 
                                           Expect(rowSums(FexMF) * rowSums(dxM * exp(-r2 * .a)) * pmi,
                                                   colSums(FexFF) * rowSums(dxF * exp(-r2 * .a)) * pfi, p)))
               # get back proportions to simplify above formulas
               pmi <- SRBi / (1 + SRBi)
               pfi <- 1 / (1 + SRBi)
               if (abs(deltai)<tol){
                   break
               }
           }
           c(r = r2, SRB = SRBi)
        })
# takes very long time
#DifferentPvals <- lapply(seq(-5,5,by=.1), function(p){
#         do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxy, .dxm, .dxf, .Ex, .p){
#                    
#                    dxm <- .dxm[,yr]
#                    dxf <- .dxf[,yr]
#                    Exm <- rowSums(ExpectedDx(with(.Ex, Male[Year == as.integer(yr)]), dxm))
#                    Exf <- rowSums(ExpectedDx(with(.Ex, Female[Year == as.integer(yr)]), dxf))
#                   
#                    BxyM    <- ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxym"]], dxm, dxf) 
#                    BxyF    <- ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxyf"]], dxm, dxf) 
#                    RatioM  <- Minf0(Mna0(BxyM / (outer(rowSums(BxyM), colSums(BxyM)) / sum(BxyM))))
#                    RatioF  <- Minf0(Mna0(BxyF / (outer(rowSums(BxyF), colSums(BxyF)) / sum(BxyF))))
#                    
#                    LotkaCrossRatioCoaleex(BxyM, BxyF, Exm, Exf, dxm, dxf, RatioM, RatioF, p = .p)
#                   
#                }, .Bxy = BxymfUS, .dxm = dxmUS, .dxf = dxfUS, .Ex = ExUS, .p = p))   
#        })
#plot(yearsUS, rUScp[,1], type = 'n', ylim = c(-.0005,.005))
#lapply(DifferentPvals, function(y){
#            lines(yearsUS, y[,1], col = "#00000004")
#        })
#lines(yearsUS, DifferentPvals[[1]][,1],col = "red")
#lines(yearsUS, DifferentPvals[[101]][,1],col = "blue")
#
#plot(yearsUS, DifferentPvals[[1]][,1]-DifferentPvals[[101]][,1], type = 'l', log = "xy" )

# LEARNED: the mean function used in the divisor in the expected value formula does not influece results much.
# In early US years the distance between p = -5 and p = 5 was greater than the year to year change. 
# after mid 1980s there was virtually no difference between p=-5 and p=5 (this interval is wide)


rUScp <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxy, .dxm, .dxf, .Ex, .p){
                    
                    dxm <- .dxm[,yr]
                    dxf <- .dxf[,yr]
                    Exm <- rowSums(ExpectedDx(with(.Ex, Male[Year == as.integer(yr)]), dxm))
                    Exf <- rowSums(ExpectedDx(with(.Ex, Female[Year == as.integer(yr)]), dxf))
                   
                    BxyM    <- ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxym"]], dxm, dxf) 
                    BxyF    <- ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxyf"]], dxm, dxf) 
                    RatioM  <- Minf0(Mna0(BxyM / (outer(rowSums(BxyM), colSums(BxyM)) / sum(BxyM))))
                    RatioF  <- Minf0(Mna0(BxyF / (outer(rowSums(BxyF), colSums(BxyF)) / sum(BxyF))))
                    
                    LotkaCrossRatioCoaleex(BxyM, BxyF, Exm, Exf, dxm, dxf, RatioM, RatioF, p = .p)
                   
                }, .Bxy = BxymfUS, .dxm = dxmUS, .dxf = dxfUS, .Ex = ExUS, .p = 2))
rownames(rUScp) <- yearsUS




rEScp <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxy, .dxm, .dxf, .Ex, .p){
                    
                    dxm <- .dxm[,yr]
                    dxf <- .dxf[,yr]
                    Exm <- rowSums(ExpectedDx(with(.Ex, Male[Year == as.integer(yr)]), dxm))
                    Exf <- rowSums(ExpectedDx(with(.Ex, Female[Year == as.integer(yr)]), dxf))
                    
                    BxyM    <- ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxym"]], dxm, dxf) 
                    BxyF    <- ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxyf"]], dxm, dxf) 
                    RatioM  <- Minf0(Mna0(BxyM / (outer(rowSums(BxyM), colSums(BxyM)) / sum(BxyM))))
                    RatioF  <- Minf0(Mna0(BxyF / (outer(rowSums(BxyF), colSums(BxyF)) / sum(BxyF))))
                    
                    LotkaCrossRatioCoaleex(BxyM, BxyF, Exm, Exf, dxm, dxf, RatioM, RatioF, p = .p)
                    
                }, .Bxy = BxymfES, .dxm = dxmES, .dxf = dxfES, .Ex = ExES, .p = 2))
rownames(rEScp) <- yearsES

plot(yearsUS, rUScp[,1], type = 'l', ylim = c(-.015,.01))
lines(yearsUS, rUS[,1], lty = 2)
#polygon(c(yearsES,rev(yearsES)), c(rEScpp5[,1],rev(rEScpm5[,1])), border = NA, col = "#55555550")
lines(yearsES, rEScp[,1], col = "red")
lines(yearsES, rES[,1], lty = 2, col = "red")
abline(h=0)
plot(yearsES,abs(rEScp[,1] / rES[,1]), type = 'l')
plot(yearsUS,abs(rUScp[,1] / rUS[,1]), type = 'l', ylim = c(0,1))


sign(diff(rEScp[,1])) == sign(diff(rES[,1]))
sign(diff(rUScp[,1])) == sign(diff(rUS[,1]))



HM(4,8)
stolarsky.mean.v(4,8,seq(-5,5,by=.01))

SHM <- optimize(function(p){
            (HM(4,8) - stolarsky.mean.v(4,8,p)) ^ 2 
        }, interval = c(-5,-3), tol = 1e-22)$min

stolarsky.mean(4,8,SHM) - HM(4,8)

-4.144574

(SHM)


stolarsky.mean(4,8,SHM)- HM(4,8)

0 =  (x - y) * ((2 * x * y) / (x + y)) - ((x^p - y^p)/p) ^ (1/(p-1))
Solve[((x^p - y^p) / (p*(x - y)))^(1 / (p - 1)) == (2 * x * y) / (x + y),p,Reals]
# would like to solve for p, as function of x and y
(2 * x * y) / (x + y) - ((x^p - y^p) / (p*(x - y)))^(1 / (p - 1)) = 0



(((4^p) - (8^p)) / (p*(4 - 8)))^(1 / (p - 1)) - (2 * 4 * 8) / (4 + 8) 




