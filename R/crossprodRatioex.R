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

Expect <- compiler::cmpfun(function(m1,m2,p=-2){
    Minf0(Mna0(outer(m1,m2, "*") / stolarsky.mean(sum(m1), sum(m2), p = p)))
})
RatioAdj <- compiler::cmpfun(function(Ratio, BxyExp){
    (Ratio * BxyExp) * Minf0(Mna0((sum(BxyExp) / sum(Ratio * BxyExp))))
})

rCrossRatioIt <- compiler::cmpfun(function(BxyM, BxyF, 
                Exm, Exf, 
                dxm, dxf, 
                RatioM, RatioF, p = -2,
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
           T.guess <-  (sum( RatioAdj(RatioM, 
                                   Expect(.a * rowSums(FexMM) * rowSums(dxM) * pmi,
                                           .a * colSums(FexFM) * rowSums(dxF) * pfi, p))) +
                   sum(RatioAdj(RatioF, 
                                   Expect(.a * rowSums(FexMF) * rowSums(dxM) * pmi,
                                           .a * colSums(FexFF) * rowSums(dxF) * pfi, p))) ) / R0
           # starting val for r, assuming above NRR and mean remaining life expectancy at reproduction of 60 years
           ri     <- log(R0) /  T.guess
           
           for(i in 1:maxit){
               # produce residual 
               deltai <- 1 - (sum( RatioAdj(RatioM, 
                                       Expect(rowSums(FexMM) * rowSums(dxM * exp(-ri * .a)) * pmi,
                                               colSums(FexFM) * rowSums(dxF * exp(-ri * .a)) * pfi, p))) +
                             sum(RatioAdj(RatioF, 
                                       Expect(rowSums(FexMF) * rowSums(dxM * exp(-ri * .a)) * pmi,
                                               colSums(FexFF) * rowSums(dxF * exp(-ri * .a)) * pfi, p))))
               # improve r
               ri    <- ri - (deltai / (T.guess - (deltai / ri)))
               # improve SRB using improved R
               SRBi  <- sum( RatioAdj(RatioM, 
                                           Expect(rowSums(FexMM) * rowSums(dxM * exp(-ri * .a)) * pmi,
                                                   colSums(FexFM) * rowSums(dxF * exp(-ri * .a)) * pfi, p))) /
                        sum(RatioAdj(RatioF, 
                                           Expect(rowSums(FexMF) * rowSums(dxM * exp(-ri * .a)) * pmi,
                                                   colSums(FexFF) * rowSums(dxF * exp(-ri * .a)) * pfi, p)))
               # get back proportions to simplify above formulas
               pmi <- SRBi / (1 + SRBi)
               pfi <- 1 / (1 + SRBi)
               
               if (abs(deltai)<tol){
                   break
               }
           }
           c(r = ri, SRB = SRBi, iter = i)
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
#                    rCrossRatioIt(BxyM, BxyF, Exm, Exf, dxm, dxf, RatioM, RatioF, p = .p)
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


rUScp2 <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxy, .dxm, .dxf, .Ex, .p){
                    
                    dxm <- .dxm[,yr]
                    dxf <- .dxf[,yr]
                    Exm <- rowSums(ExpectedDx(with(.Ex, Male[Year == as.integer(yr)]), dxm))
                    Exf <- rowSums(ExpectedDx(with(.Ex, Female[Year == as.integer(yr)]), dxf))
                   
                    BxyM    <- ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxym"]], dxm, dxf) 
                    BxyF    <- ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxyf"]], dxm, dxf) 
                    RatioM  <- Minf0(Mna0(BxyM / (outer(rowSums(BxyM), colSums(BxyM),"*") / sum(BxyM))))
                    RatioF  <- Minf0(Mna0(BxyF / (outer(rowSums(BxyF), colSums(BxyF),"*") / sum(BxyF))))
                    rCrossRatioIt(BxyM, BxyF, Exm, Exf, dxm, dxf, RatioM, RatioF, p = .p)
                }, .Bxy = BxymfUS, .dxm = dxmUS, .dxf = dxfUS, .Ex = ExUS, .p = -2))
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
                    
                    rCrossRatioIt(BxyM, BxyF, Exm, Exf, dxm, dxf, RatioM, RatioF, p = .p)
                    
                }, .Bxy = BxymfES, .dxm = dxmES, .dxf = dxfES, .Ex = ExES, .p = -2))
rownames(rEScp) <- yearsES

plot(yearsUS, rUScp[,1], type = 'l', ylim = c(-.015,.01))
#lines(yearsUS, rUS[,1], lty = 2)
#polygon(c(yearsES,rev(yearsES)), c(rEScpp5[,1],rev(rEScpm5[,1])), border = NA, col = "#55555550")
lines(yearsES, rEScp[,1], col = "red")
#lines(yearsES, rES[,1], lty = 2, col = "red")
abline(h=0)
plot(yearsES,abs(rEScp[,1] / rES[,1]), type = 'l')
plot(yearsUS,abs(rUScp[,1] / rUS[,1]), type = 'l', ylim = c(0,1))

#------------------------------------------------------------------
# 1) show example of Ratio object
# test IPF
#ExBxy2       <- ExpectedDxMxFmatrix( Mat=BxUS[["1970"]], dxm=dxmUS[,"1970"], dxf=dxfUS[,"1970"])
#expected2    <- outer(rowSums(ExBxy2), colSums(ExBxy2), "*") / sum(ExBxy2)
##


yr     <- "1975"
dxm    <- dxmUS[,yr] ; dxf <- dxfUS[,yr]
Bxy    <- ExpectedDxMxFmatrix(BxymfUS[[yr]][["Bxym"]] + BxymfUS[[yr]][["Bxyf"]], dxm, dxf) 
Ratio  <- Minf0(Mna0(Bxy / (outer(rowSums(Bxy), colSums(Bxy),"*") / sum(Bxy))))


# let's do an image, same 1975
colramp     <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"), space = "Lab")

# plotting pars
brks        <- seq(-6, 6, length.out = 51)
#brks <- seq(min(lExBxy, na.rm = TRUE), max(lExBxy, na.rm = TRUE), length.out = 51)
g.xy        <- seq(0, 100, by = 5)
gb.xy       <- seq(0, 100, by = 10)
levs        <- c(-6, 6, seq(500, 3000, by = 500))     # for contour plot
levs <- c(.01,.1,.5,.75,1.25,2,10,100,1000)
ExBxy[ExBxy == 0]       <- NA
expected[expected == 0] <- NA
pdf("/home/triffe/git/DISS/latex/Figures/exCPRatioExample.pdf", height = 5, width = 5)
par(mai = c(.4,.3,.3,.2), xaxs = "i", yaxs = "i")
image(x = ages + .5, y = ages + .5, MinfNA(log(Ratio)), 
        xlim = c(0, 101), ylim = c(0, 101), zlim = c(.5,7),
        col = rev(colramp(50)), breaks = brks, axes = FALSE, asp = 1,
        xlab = "", ylab = "", 
        panel.first = list(rect(0, 0, 101, 101, col = "#EEEEEE", xpd = TRUE, border = NA), 
                abline(h = g.xy, col = "white", lwd = .5),
                abline(v = g.xy, col = "white", lwd = .5),
                text(0, gb.xy, gb.xy, pos = 2, cex = .5, xpd = TRUE),
                text(gb.xy, 0, gb.xy, pos = 1, cex = .5, xpd = TRUE),
                segments(0, gb.xy, -1, gb.xy, xpd = TRUE),
                segments(gb.xy, 0, gb.xy, -1, xpd = TRUE)))
# contours
contour(x = ages + .5, y = ages + .5,  MinfNA(log(Ratio)), 
        levels = log(levs), labels = levs, add = TRUE)
# axis labels
fath <- "Father"
moth <- "Mother"
text(50,-4, bquote(.(fath) ~ e[y]), xpd = TRUE, cex = .7, pos =1)
text(-8,105,bquote(.(moth) ~ e[y]), xpd = TRUE, pos = 4, cex = .7)
#
segments(0,0,101,101,col = "#50505050")
dev.off()

(BirthWeightedAvg <- sum(Bxy * exp(abs(log(Ratio))), na.rm=TRUE) / sum(Bxy)) # very small

# -----------------------------------------------------------------------
# 2) demonstrate IPF superiority in predicting year t+1 joint birth distribution
IPFpred2 <- compiler::cmpfun(function(Bxy, Exm1, Exm2, Exf1, Exf2, marM = mean, tol = 1e-15, maxit = 20){
            
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
            BxyA
        })

Expect <- compiler::cmpfun(function(m1,m2,p=-2){
            Minf0(Mna0(outer(m1,m2, "*") / stolarsky.mean(sum(m1), sum(m2), p = p)))
        })
RatioAdj <- compiler::cmpfun(function(Ratio, BxyExp){
            (Ratio * BxyExp) * Minf0(Mna0((sum(BxyExp) / sum(Ratio * BxyExp))))
        })
yr <- 1975
CompUS <- do.call(rbind, lapply(yearsUS[-length(yearsUS)], function(yr, .Bxy, .Ex, .dxm, .dxf){
                    yrc         <- as.character(yr)
                    yrcp1       <- as.character(yr+1)
                    Bt          <- ExpectedDxMxFmatrix(.Bxy[[yrc]], .dxm[,yrc], .dxf[,yrc])
                    Btp1test    <- ExpectedDxMxFmatrix(.Bxy[[yrcp1]], .dxm[,yrcp1], .dxf[,yrcp1])
                    Exm1        <- rowSums(ExpectedDx(with(.Ex,Male[Year == yr]), .dxm[,yrc]))
                    Exm2        <- rowSums(ExpectedDx(with(.Ex,Male[Year == (yr + 1)]), .dxm[,yrcp1]))
                    Exf1        <- rowSums(ExpectedDx(with(.Ex,Female[Year == yr]), .dxf[,yrc]))
                    Exf2        <- rowSums(ExpectedDx(with(.Ex,Female[Year == (yr + 1)]), .dxf[,yrcp1]))
                    
                    IPFpred     <- IPFpred2(Bt, Exm1 = Exm1, Exm2 = Exm2, Exf1 = Exf1, Exf2 = Exf2, marM = mean)
                    Ratio       <- Minf0(Mna0(Bt / (outer(rowSums(Bt), colSums(Bt)) / sum(Bt))))
                   
                    Bpm         <- Minf0(Mna0(rowSums(Bt) / Exm1)) * Exm2
                    Bpf         <- Minf0(Mna0(colSums(Bt) / Exf1)) * Exf2
                    
                    Expec       <- Minf0(Mna0(outer(Bpm,Bpf, "*") / mean(c(sum(Bpm), sum(Bpf)))))
                    Pred        <- RatioAdj(Ratio, Expec)
                    
                    IPFpdf      <- IPFpred / sum(IPFpred)   # IPF predicted pdf
                    CRpdf       <- Pred / sum(Pred)         # CR predicted pdf
                    Test        <- Btp1test / sum(Btp1test) # pdf observed from t plus 1
                    
                    c(IPF = 1 - sum(pmin(IPFpdf, Test)), CR = 1 - sum(pmin(CRpdf, Test)), 
                            IPFmales = 1 - sum(pmin(rowSums(IPFpdf), rowSums(Test))),
                            IPFfem = 1 - sum(pmin(colSums(IPFpdf), colSums(Test))),
                            CRmales = 1 - sum(pmin(rowSums(CRpdf), rowSums(Test))),
                            CRfem = 1 - sum(pmin(colSums(CRpdf), colSums(Test))))
                }, .Bxy = BxUS, .Ex = ExUS, .dxm = dxmUS, .dxf = dxfUS))

CompES <- do.call(rbind, lapply(yearsES[-length(yearsES)], function(yr, .Bxy, .Ex, .dxm, .dxf){
                    yrc         <- as.character(yr)
                    yrcp1       <- as.character(yr+1)
                    Bt          <- ExpectedDxMxFmatrix(.Bxy[[yrc]], .dxm[,yrc], .dxf[,yrc])
                    Btp1test    <- ExpectedDxMxFmatrix(.Bxy[[yrcp1]], .dxm[,yrcp1], .dxf[,yrcp1])
                    Exm1        <- rowSums(ExpectedDx(with(.Ex,Male[Year == yr]), .dxm[,yrc]))
                    Exm2        <- rowSums(ExpectedDx(with(.Ex,Male[Year == (yr + 1)]), .dxm[,yrcp1]))
                    Exf1        <- rowSums(ExpectedDx(with(.Ex,Female[Year == yr]), .dxf[,yrc]))
                    Exf2        <- rowSums(ExpectedDx(with(.Ex,Female[Year == (yr + 1)]), .dxf[,yrcp1]))
                    
                    IPFpred     <- IPFpred2(Bt, Exm1 = Exm1, Exm2 = Exm2, Exf1 = Exf1, Exf2 = Exf2, marM = mean)
                    Ratio       <- Minf0(Mna0(Bt / (outer(rowSums(Bt), colSums(Bt)) / sum(Bt))))
                    
                    Bpm         <- Minf0(Mna0(rowSums(Bt) / Exm1)) * Exm2
                    Bpf         <- Minf0(Mna0(colSums(Bt) / Exf1)) * Exf2
                    
                    Expec       <- Minf0(Mna0(outer(Bpm,Bpf, "*") / mean(c(sum(Bpm), sum(Bpf)))))
                    Pred        <- RatioAdj(Ratio, Expec)
                    
                    IPFpdf      <- IPFpred / sum(IPFpred)   # IPF predicted pdf
                    CRpdf       <- Pred / sum(Pred)         # CR predicted pdf
                    Test        <- Btp1test / sum(Btp1test) # pdf observed from t plus 1
                    
                    c(IPF = 1 - sum(pmin(IPFpdf, Test)), CR = 1 - sum(pmin(CRpdf, Test)))
                }, .Bxy = BxES, .Ex = ExES, .dxm = dxmES, .dxf = dxfES))

head(CompUS)
plot(1969:2008, CompUS[,2] / CompUS[,1], type = 'l')
plot(1969:2008, CompUS[,5] / CompUS[,3], type = 'l')
lines(1969:2008, CompUS[,6] / CompUS[,4], col = 'red')
abline(h=1)
sum(CompUS[,5] < CompUS[,3])
sum(CompUS[,6] < CompUS[,4])

plot(1975:2008, CompES[,2] / CompES[,1], type = 'l')

sum(CompES[,2] < CompES[,1])
sum(CompUS[,2] < CompUS[,1])
length(yearsUS)

# ----------------------------------------------------------
# 3) check competition

yrc         <- as.character(yr)
Bt          <- ExpectedDxMxFmatrix(.Bxy[[yrc]], .dxm[,yrc], .dxf[,yrc])
Exf1        <- rowSums(ExpectedDx(with(.Ex,Female[Year == yr]), .dxf[,yrc]))
#Exf2        <- rowSums(ExpectedDx(with(.Ex,Female[Year == yr]), .dxf[,yrcp1]))
Exm1        <- rowSums(ExpectedDx(with(.Ex,Male[Year == yr]), .dxm[,yrc]))
#Exm2        <- rowSums(ExpectedDx(with(.Ex,Male[Year == yr]), .dxm[,yrcp1]))
Exm2        <- Exm1
Exm2[31]    <- Exm2[31] * 1.5
Exf2        <- Exf1
Ratio       <- Minf0(Mna0(Bt / (outer(rowSums(Bt), colSums(Bt)) / sum(Bt))))

Bpm         <- Minf0(Mna0(rowSums(Bt) / Exm1)) * Exm2
Bpf         <- Minf0(Mna0(colSums(Bt) / Exf1)) * Exf2

Expec       <- Minf0(Mna0(outer(Bpm,Bpf, "*") / mean(c(sum(Bpm), sum(Bpf)))))
Pred        <- RatioAdj(Ratio, Expec)

plot(0:110, rowSums(Pred), type = 'l', col = "red")
lines(0:110, rowSums(Bt), col = "blue")

plot(0:110, rowSums(Pred) / rowSums(Bt), type = 'l', ylim = c(.95,1.05))

plot(0:110, colSums(Pred) / colSums(Bt), type = 'l', ylim = c(.99,1.02))
abline(h=1)
