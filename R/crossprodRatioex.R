setwd("/home/triffe/git/DISS/")
source("R/UtilityFunctions.R")
source("R/MeanFunctions.R")

yearsUS <- 1969:2009
yearsES <- 1975:2009

# Bxy totals  (not used()
BxUS  <- local(get(load("Data/USbirths/USBxy0_110.Rdata")))
BxES  <- local(get(load("Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

# repeated below for better SRB assumptions
BxymfES <- local(get(load("Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("Data/USbirths/USBxymf0_110.Rdata")))
names(BxymfES) <- yearsES
names(BxymfUS) <- yearsUS
# exposures, as such, straight from HMD, all ages 0-110, long form
ExUS  <- local(get(load("Data/Exposures/USexp.Rdata")))
ExES  <- local(get(load("Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
dxmUS <- local(get(load("Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("Data/HMD_dx/dxfES.Rdata"))) 
# make sum to 1
dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)

Expect <- compiler::cmpfun(function(m1,m2,meanFun, ...){
    Minf0(Mna0(outer(m1,m2, "*") / meanFun(sum(m1), sum(m2), ...)))
})
RatioAdj <- compiler::cmpfun(function(Ratio, BxyExp){
    (Ratio * BxyExp) * Minf0(Mna0((sum(BxyExp) / sum(Ratio * BxyExp))))
})

rCrossRatioIt <- compiler::cmpfun(function(BxyM, BxyF, 
                Exm, Exf, 
                dxm, dxf, 
                RatioM, RatioF, meanFun = HM,
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
                                          colSums(FexFM) * rowSums(dxF), meanFun))) / 
                       sum(RatioAdj(RatioF, 
                                   Expect(rowSums(FexMF) * rowSums(dxM),
                                          colSums(FexFF) * rowSums(dxF), meanFun)))
            # proportions male and female
            pmi   <- SRBi / (1 + SRBi)
            pfi   <- 1 / (1 + SRBi)
            # starting guess at NRR, assuming r = 0 and above guess at SRB
            R0    <-  sum( RatioAdj(RatioM, 
                                    Expect(rowSums(FexMM) * rowSums(dxM) * pmi,
                                            colSums(FexFM) * rowSums(dxF) * pfi, meanFun))) +
                    sum(RatioAdj(RatioF, 
                                    Expect(rowSums(FexMF) * rowSums(dxM) * pmi,
                                            colSums(FexFF) * rowSums(dxF) * pfi, meanFun)))
           T.guess <-  (sum( RatioAdj(RatioM, 
                                   Expect(.a * rowSums(FexMM) * rowSums(dxM) * pmi,
                                           .a * colSums(FexFM) * rowSums(dxF) * pfi, meanFun))) +
                   sum(RatioAdj(RatioF, 
                                   Expect(.a * rowSums(FexMF) * rowSums(dxM) * pmi,
                                           .a * colSums(FexFF) * rowSums(dxF) * pfi, meanFun))) ) / R0
           # starting val for r, assuming above NRR and mean remaining life expectancy at reproduction of 60 years
           ri     <- log(R0) /  T.guess
           
           for(i in 1:maxit){
               # produce residual 
               deltai <- 1 - (sum( RatioAdj(RatioM, 
                                       Expect(rowSums(FexMM) * rowSums(dxM * exp(-ri * .a)) * pmi,
                                               colSums(FexFM) * rowSums(dxF * exp(-ri * .a)) * pfi, meanFun))) +
                             sum(RatioAdj(RatioF, 
                                       Expect(rowSums(FexMF) * rowSums(dxM * exp(-ri * .a)) * pmi,
                                               colSums(FexFF) * rowSums(dxF * exp(-ri * .a)) * pfi, meanFun))))
               # improve r
               ri    <- ri - (deltai / (T.guess - (deltai / ri)))
               # improve SRB using improved R
               SRBi  <- sum( RatioAdj(RatioM, 
                                           Expect(rowSums(FexMM) * rowSums(dxM * exp(-ri * .a)) * pmi,
                                                   colSums(FexFM) * rowSums(dxF * exp(-ri * .a)) * pfi, meanFun))) /
                        sum(RatioAdj(RatioF, 
                                           Expect(rowSums(FexMF) * rowSums(dxM * exp(-ri * .a)) * pmi,
                                                   colSums(FexFF) * rowSums(dxF * exp(-ri * .a)) * pfi, meanFun)))
               # get back proportions to simplify above formulas
               pmi <- SRBi / (1 + SRBi)
               pfi <- 1 / (1 + SRBi)
               
               if (abs(deltai)<tol){
                   break
               }
           }
           c(r = ri, SRB = SRBi, iter = i)
        })
        
rCrossRatiocy <- compiler::cmpfun(function(r, SRB, dxm, dxf, .a = .5:110.5){
    pmi <- SRB / (1+SRB)
    pfi <- 1 / (1+SRB)
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
    b <- 1 / sum(rowSums(dxM * exp(-r * .a)) * pmi + rowSums(dxF * exp(-r * .a)) * pfi)
    
    cym <- b * pmi * rowSums(dxM * exp(-r * .a))
    cyf <- b * pfi * rowSums(dxF * exp(-r * .a))
    cbind(cym, cyf)
})
rCrossRatiocyWide <- compiler::cmpfun(function(r, SRB, dxm, dxf, .a = .5:110.5){
            pmi <- SRB / (1+SRB)
            pfi <- 1 / (1+SRB)
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
            b <- 1 / sum(rowSums(dxM * exp(-r * .a)) * pmi + rowSums(dxF * exp(-r * .a)) * pfi)
            
            cym <- b * pmi * dxM * exp(-r * .a)
            cyf <- b * pfi * dxF * exp(-r * .a)
            list(cym=cym, cyf=cyf)
        })
HM <- compiler::cmpfun(function(x,y){
            Mna0(Minf0((2 * x * y) / (x + y)))
        })
# ---------------------------------------------------------------------------
# test sensitivity to 'p'
# this was done with a different Expec() function that had stolarsky.mean() inside it, iterating over p
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


# ----------------------------------------------------------------
# get r and SRB for US and Spain
rUScp <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxy, .dxm, .dxf, .Ex, .meanFun){
                    
                    dxm <- .dxm[,yr]
                    dxf <- .dxf[,yr]
                    Exm <- rowSums(ExpectedDx(with(.Ex, Male[Year == as.integer(yr)]), dxm))
                    Exf <- rowSums(ExpectedDx(with(.Ex, Female[Year == as.integer(yr)]), dxf))
                   
                    BxyM    <- ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxym"]], dxm, dxf) 
                    BxyF    <- ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxyf"]], dxm, dxf) 
                    RatioM  <- Minf0(Mna0(BxyM / (outer(rowSums(BxyM), colSums(BxyM),"*") / sum(BxyM))))
                    RatioF  <- Minf0(Mna0(BxyF / (outer(rowSums(BxyF), colSums(BxyF),"*") / sum(BxyF))))
                    rCrossRatioIt(BxyM, BxyF, Exm, Exf, dxm, dxf, RatioM, RatioF, meanFun = .meanFun)
                }, .Bxy = BxymfUS, .dxm = dxmUS, .dxf = dxfUS, .Ex = ExUS, .meanFun = HM))
rownames(rUScp) <- yearsUS

rEScp <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxy, .dxm, .dxf, .Ex, .meanFun){
                    
                    dxm <- .dxm[,yr]
                    dxf <- .dxf[,yr]
                    Exm <- rowSums(ExpectedDx(with(.Ex, Male[Year == as.integer(yr)]), dxm))
                    Exf <- rowSums(ExpectedDx(with(.Ex, Female[Year == as.integer(yr)]), dxf))
                    
                    BxyM    <- ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxym"]], dxm, dxf) 
                    BxyF    <- ExpectedDxMxFmatrix(.Bxy[[yr]][["Bxyf"]], dxm, dxf) 
                    RatioM  <- Minf0(Mna0(BxyM / (outer(rowSums(BxyM), colSums(BxyM)) / sum(BxyM))))
                    RatioF  <- Minf0(Mna0(BxyF / (outer(rowSums(BxyF), colSums(BxyF)) / sum(BxyF))))
                    
                    rCrossRatioIt(BxyM, BxyF, Exm, Exf, dxm, dxf, RatioM, RatioF, meanFun = .meanFun)
                    
                }, .Bxy = BxymfES, .dxm = dxmES, .dxf = dxfES, .Ex = ExES, .meanFun = HM))
rownames(rEScp) <- yearsES

#save(rUScp, file = "Data/results/exCPr/rUScp.Rdata")
#save(rEScp, file = "Data/results/exCPr/rEScp.Rdata")
rmUS <- local(get(load("Data/results/exSingleSex/rmUS.Rdata")))[,1]
rfUS <- local(get(load("Data/results/exSingleSex/rfUS.Rdata")))[,1]
rmES <- local(get(load("Data/results/exSingleSex/rmES.Rdata")))[,1]
rfES <- local(get(load("Data/results/exSingleSex/rfES.Rdata")))[,1]

# take a look:

pdf("latex/Figures/exCRr.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rUScp[,1], type = 'n', ylim = c(-.016,.01),xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,-.016,2010,.01,col = gray(.95), border=NA),
                abline(h = seq(-.016,.01,by = .002), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.016,.01,by = .002),seq(-.016,.01,by = .002), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.016, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.0175, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,.0115, "r", cex = 1, xpd = TRUE)))
polygon(c(yearsUS, rev(yearsUS)), c(rmUS, rev(rfUS)), col = "#99999950")
polygon(c(yearsES, rev(yearsES)), c(rmES, rev(rfES)), col = "#99999950", lty = 5)
lines(yearsUS, rUScp[,1],col = gray(.2), lwd = 2)
lines(yearsES, rEScp[,1],col = gray(.2), lwd = 2, lty = 5)
text(c(1992.5, 1995.626, 1995.8, 1977.5, 1972, 1974.3), 
        c(-0.013596517, -0.0098, -0.0065, -0.0025, -0.004, 0.001981268),
        c(expression(r^f~ES),expression(r^m~ES),expression(r^{RAdj-HM}~ES),
          expression(r^f~US),expression(r^{RAdj-HM}~US), expression(r^m~US)),  cex = .8, xpd = TRUE)
segments(1971.187, -0.0030000000,1974.026,-0.0002376602, col = gray(.4),lwd=.7)
dev.off()


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
pdf("latex/Figures/exCPRatioExample.pdf", height = 5, width = 5)
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


CompUS <- do.call(rbind, lapply(yearsUS[-length(yearsUS)], function(yr, .Bxy, .Ex, .dxm, .dxf){
                    # for indexing
                    yrc         <- as.character(yr)
                    yrcp1       <- as.character(yr+1)
                    # present year birth mat
                    Bt          <- ExpectedDxMxFmatrix(.Bxy[[yrc]], .dxm[,yrc], .dxf[,yrc])
                    # next year's birth mat
                    Btp1test    <- ExpectedDxMxFmatrix(.Bxy[[yrcp1]], .dxm[,yrcp1], .dxf[,yrcp1])
                    # male and female exposures 1 and 2
                    Exm1        <- rowSums(ExpectedDx(with(.Ex,Male[Year == yr]), .dxm[,yrc]))
                    Exm2        <- rowSums(ExpectedDx(with(.Ex,Male[Year == (yr + 1)]), .dxm[,yrcp1]))
                    Exf1        <- rowSums(ExpectedDx(with(.Ex,Female[Year == yr]), .dxf[,yrc]))
                    Exf2        <- rowSums(ExpectedDx(with(.Ex,Female[Year == (yr + 1)]), .dxf[,yrcp1]))
                    # prediction year t+1 from IPF
                    IPFpred     <- IPFpred2(Bt, Exm1 = Exm1, Exm2 = Exm2, Exf1 = Exf1, Exf2 = Exf2, marM = mean)
                    
                    # components for ratio adjustment method
                    Ratio       <- Minf0(Mna0(Bt / (outer(rowSums(Bt), colSums(Bt)) / sum(Bt))))
                    Bpm         <- Minf0(Mna0(rowSums(Bt) / Exm1)) * Exm2
                    Bpf         <- Minf0(Mna0(colSums(Bt) / Exf1)) * Exf2
                    # association-free, mean
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
                    
                    c(IPF = 1 - sum(pmin(IPFpdf, Test)), CR = 1 - sum(pmin(CRpdf, Test)), 
                            IPFmales = 1 - sum(pmin(rowSums(IPFpdf), rowSums(Test))),
                            IPFfem = 1 - sum(pmin(colSums(IPFpdf), colSums(Test))),
                            CRmales = 1 - sum(pmin(rowSums(CRpdf), rowSums(Test))),
                            CRfem = 1 - sum(pmin(colSums(CRpdf), colSums(Test))))
                }, .Bxy = BxES, .Ex = ExES, .dxm = dxmES, .dxf = dxfES))

head(CompUS)
plot(1969:2008, CompUS[,2] / CompUS[,1], type = 'l')
plot(1969:2008, CompUS[,5] / CompUS[,3], type = 'l')
lines(1969:2008, CompUS[,6] / CompUS[,4], col = 'red')
abline(h=1)
sum(CompUS[,5] < CompUS[,3])
sum(CompUS[,6] < CompUS[,4])
sum((CompUS[,2] / CompUS[,1]) < 1)

1-colMeans(CompUS)
1-colMeans(CompES)
plot(1975:2008, CompES[,2] / CompES[,1], type = 'l')

sum(CompES[,2] < CompES[,1])
sum(CompUS[,2] < CompUS[,1])
length(yearsUS)

# ----------------------------------------------------------
# 3) check competition

yr <- 1975
yrc         <- as.character(yr)
Bt          <- ExpectedDxMxFmatrix(BxUS[[yrc]], dxmUS[, yrc], dxfUS[, yrc])
Ratio       <- Minf0(Mna0(Bt / (outer(rowSums(Bt), colSums(Bt)) / sum(Bt))))
# age-classified rates:
Exf     <- with(ExUS,Female[Year == yr])
Exm     <- with(ExUS,Male[Year == yr])
Bxm     <- rowSums(BxUS[[yrc]])
Bxf     <- colSums(BxUS[[yrc]])
Fxm     <- Bxm / Exm
Fxf     <- Bxf / Exf
.a <- .5:110.5
a <- 0:110
#plot(a ,Fxf, type = 's', ylim = c(0,.14))
#lines(a ,Fxm, type = 's', col = "blue")

# how do these distribute out over ex? (last M for matrix)

Fym <- Minf0(Mna0(rowSums(Bt) / rowSums(ExpectedDx(Exm, dxmUS[,yrc]))))
Fyf <- Minf0(Mna0(colSums(Bt) / rowSums(ExpectedDx(Exf, dxfUS[,yrc]))))

FxmM    <- ExpectedDx(Fxm, dxmUS[,yrc])
FxfM    <- ExpectedDx(Fxf, dxfUS[,yrc])
# now what are the row-wise pdfs?
FxmPDF  <- Minf0(Mna0(FxmM / rowSums(FxmM)))  # use these to redist?
FxfPDF  <- Minf0(Mna0(FxfM / rowSums(FxfM)))  #

#plot(a, rowSums(FxmM), type = 's',col="blue")
#lines(a, rowSums(FxfM), type = 's')
# OK, now get stable pop:
r               <- rUScp[yrc,"r"]
SRB             <- rUScp[yrc,"SRB"]
cy1             <- rCrossRatiocyWide(r, SRB, dxmUS[,yrc], dxfUS[,yrc])
cym1            <- cy1[[1]]
cyf1            <- cy1[[2]]
# inflate one cohort of males
cym2            <- cym1
cym2[, 36]      <- cym2[,36] * 1.25
cyf2            <- cyf1

bpm1         <- Fym * rowSums(cym1)
bpf1         <- Fyf * rowSums(cyf1)
Expec1       <- Expect(bpm1,bpf1,p=-2)
Pred1        <- RatioAdj(Ratio, Expec1)

plot(a,colSums(rowSums(Pred1 / rowSums(cym1)) * FxmPDF),type= 'l',col = "blue")
lines(a,Fxm,lty=2,col="blue")
lines(a,colSums(colSums(Pred1 / rowSums(cyf1)) * FxfPDF),col="red")
lines(a,Fxf,lty=2,col = "red")
# but compare with ageified ex rates
plot(a,colSums(rowSums(Pred1 / rowSums(cym1)) * FxmPDF),type= 'l',col = "blue")
lines(a,colSums(Fym * FxmPDF),lty=2,col="blue")
lines(a,colSums(colSums(Pred1 / rowSums(cyf1)) * FxfPDF),col="red")
lines(a,colSums(Fyf * FxfPDF),lty=2,col = "red")

# now make pred 2:
bpm2         <- Fym * rowSums(cym2)
bpf2         <- Fyf * rowSums(cyf2)
Expec2       <- Expect(bpm2,bpf2,p=-2)
Pred2        <- RatioAdj(Ratio, Expec2)
mp1 <- rowSums(Pred1 / rowSums(cym1))
mp2 <- rowSums(Pred2) / rowSums(cym2) # new male rates
plot(a,mp1/mp2,type= 'l',col = "blue")

rowSums(Pred1) / cym1


mp1a <- colSums(Minf0(Mna0(rowSums(Pred1) / cym1)))
mp2a <- colSums(Minf0(Mna0(rowSums(Pred2) / cym2)))
plot(a,mp2a/mp1a,type= 'l',col = "blue", ylim = c(.99,1.01))
abline(h=1)


plot(a,colSums(mp * FxmPDF),type= 'l',col = "blue")
lines(a,colSums(rowSums(Pred1 / rowSums(cym1)) * FxmPDF), col = "red")


FxmPDF


#fields::image.plot(Minf0(Mna0((Pred2/sum(Pred2)) - (Pred1/sum(Pred1)))))


sum(hmm2)
sum(hmm1)
plot(a,colSums(hmm2)/sum(hmm2),type= 'l',col = "blue")
abline(v=30)
lines(a,colSums(hmm1)/sum(hmm1),lty=2,col="blue")

plot(Minf0(Mna0(colSums(hmm2)/colSums(hmm1))),type = 'l', ylim = c(.98,1.03))

plot(a,(colSums(hmm2)/sum(hmm2)) / (colSums(hmm1)/sum(hmm1)), type = 'l',col = "blue")
abline(h=1.25)

plot(a,(colSums(hmm2)/sum(hmm2)) / (colSums(hmm1)/sum(hmm1)), type = 'l',col = "blue", ylim = c(.95,1.05))
abline(h=1.25)

lines(a,colSums(colSums(Pred2 / rowSums(cyf2)) * FxfPDF),col="red")
lines(a,colSums(Fyf * FxfPDF),lty=2,col = "red")

plot(a,colSums(rowSums(Pred2 / rowSums(cym2)) * FxmPDF) / colSums(rowSums(Pred1 / rowSums(cym1)) * FxmPDF), type = 'l')

