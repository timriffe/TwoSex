source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
# BxES is 0:110, years 1975:2009
BxUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
BxES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata")))

BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))
BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata")))
# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
LxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmUS.Rdata"))) / 1e5
LxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfUS.Rdata"))) / 1e5
LxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmES.Rdata"))) / 1e5
LxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfES.Rdata"))) / 1e5
yearsUS <- 1969:2009
yearsES <- 1975:2009
# -------------------------------------------------------------
#

IPFpred <- compiler::cmpfun(function(Bxy, Exm1, Exm2, Exf1, Exf2, marM = mean, tol = 1e-15, maxit = 20){
            
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
            list(FxmPred = FxmPred, FxfPred = FxfPred)
        })
       
# Demonstrate competition property:


Bxy <- BxymfUS[["1975"]][["Bxym"]] + BxymfUS[["1975"]][["Bxyf"]]
#
Fxf <- colSums(Bxy) / with(ExUS, Female[Year == 1975])
Fxm <- rowSums(Bxy) / with(ExUS, Male[Year == 1975])
        
FxPred <- IPFpred(Bxy, 
        Exm1 = with(ExUS, Male[Year == 1975]), 
        Exm2 = with(ExUS, Male[Year == 1980]), 
        Exf1 = with(ExUS, Female[Year == 1975]), 
        Exf2 = with(ExUS, Female[Year == 1980]), marM = harmonic.mean)

## comp test. increase age 25 by 1.5 times
Exm2 <- with(ExUS, Male[Year == 1980])
Exm2[26] <- Exm2[26] * 1.5
FxPred2 <- IPFpred(Bxy, 
        Exm1 = with(ExUS, Male[Year == 1975]), 
        Exm2 = Exm2, 
        Exf1 = with(ExUS, Female[Year == 1975]), 
        Exf2 = with(ExUS, Female[Year == 1980]), marM = harmonic.mean)


MaleRatio   <- FxPred2[[1]] / FxPred[[1]]
FemaleRatio <- FxPred2[[2]] / FxPred[[2]]

pdf("/home/triffe/git/DISS/latex/Figures/IPFagecompetitiontest.pdf", height = 5, width = 5)
par(mai = c(.5,.5,.5,.3),xaxs = "i", yaxs = "i")
plot(0:110, MaleRatio, type = 'l', ylim = c(.95, 1.05), xlim = c(15,40), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(15,.95,40,1.05,col = gray(.95), border=NA),
                abline(h = seq(.95,1.05,by = .01), col = "white"),
                abline(v = seq(15, 40, by = 5), col = "white"),
                text(15, seq(.96,1.05,by = .01),seq(.96,1.05,by = .01),pos = 2, cex = .8, xpd = TRUE),
                text(seq(15, 40, by = 5),.95, seq(15, 40, by = 5), pos = 1, cex = .8, xpd = TRUE),
                text(27, .945, "Age", cex = 1, pos = 1, xpd = TRUE),
                text(15,1.06, "Rate Ratio", cex = 1, xpd = TRUE)))
lines(0:110, FemaleRatio, lwd = 2.5, col = gray(.5))

text(27.85882,0.977,"Male rate penalty",pos=4)
text(27.85882,1.008,"Female rate increase",pos=4)
dev.off()

#plot(FxPred[[1]], type = 'l', lty = 2, col = "blue")
#lines(FxPred[[2]], lty = 2, col = "red")
#
#lines(FxPred2[[1]], col = "blue")
#lines(FxPred2[[2]], col = "red")

# ------------------------------------------------------------
# now into stable pop land
# suboptimal because SRB not optimized
#rIPFmin <- compiler::cmpfun(function(r, Bxy, Exm, Exf, Lxm, Lxf, SRB = 1.05, .a = .5:110.5){
#    p.m <- SRB / (1 + SRB)
#    p.f <- 1 / (1 + SRB)
#    
#    Fxpred <- IPFpred(Bxy, 
#            Exm1 = Exm, 
#            Exm2 = exp(-r * .a) * Lxm, 
#            Exf1 = Exf, 
#            Exf2 = exp(-r * .a) * Lxf)
#    (2 - sum(p.m * exp(-r * .a) * Lxm * Fxpred[[1]] + p.f * exp(-r * .a) * Lxf * Fxpred[[2]]))^2
#})

#optimize(rIPFmin, interval = c(-.02,.02), Bxy = Bxy, 
#        Exm = Exm, Exf = Exf, 
#        Lxm = LxmUS[,"1969"], Lxf = LxfUS[,"1969"],
#        SRB = sum(BxymfUS[[1]][[1]]) / sum(BxymfUS[[1]][[2]]), tol = 1e-15)
#
#sum(Bxy)
#sum(Fxpred.hat[[1]])

rIPFit <- compiler::cmpfun(function(Bxym, Bxyf, Exm, Exf, Lxm, Lxf, M = mean,
                .a = .5:110.5, maxit = 2e3, tol = 1e-10){
            SRB <- sum(Bxym) / sum(Bxyf)
            p.m <- (SRB / (1 + SRB))
            p.f <- (1 / (1 + SRB))
            Bxy    <- Bxym + Bxyf
            Fxpred.hat <- IPFpred(Bxy, 
                    Exm1 = Exm, 
                    Exm2 = Lxm, 
                    Exf1 = Exf, 
                    Exf2 = Lxf,
                    marM = M)
            R0.hat <- sum(p.m * Lxm * Fxpred.hat[[1]] + p.f * Lxf * Fxpred.hat[[2]]) / 2
            T.hat <- (sum(.a * p.m * Lxm * Fxpred.hat[[1]] + .a * p.f * Lxf * Fxpred.hat[[2]]) / 2) / R0.hat
            
            r.i <- log(R0.hat) / T.hat
            
            for (i in 1:maxit){
                Fxpredm <- IPFpred(Bxym, 
                        Exm1 = Exm, 
                        Exm2 = p.m * exp(-r.i * .a) * Lxm, 
                        Exf1 = Exf, 
                        Exf2 = p.f * exp(-r.i * .a) * Lxf,
                        marM = M)
                Fxpredf <- IPFpred(Bxyf, 
                        Exm1 = Exm, 
                        Exm2 = p.m * exp(-r.i * .a) * Lxm, 
                        Exf1 = Exf, 
                        Exf2 = p.f * exp(-r.i * .a) * Lxf,
                        marM = M)
                
                delta.i <- (2 - sum(p.m * exp(-r.i * .a) * Lxm * (Fxpredm[[1]]+Fxpredf[[1]]) + 
                                p.f * exp(-r.i * .a) * Lxf * (Fxpredm[[2]]+Fxpredf[[2]]))) / 2
                r.i <- r.i - (delta.i / ( T.hat  - (delta.i / r.i)))
                
                SRB.i <- sum(p.m * exp(-r.i * .a) * Lxm * Fxpredm[[1]] + 
                                        p.f * exp(-r.i * .a) * Lxf * Fxpredm[[2]]) / 
                        sum(p.m * exp(-r.i * .a) * Lxm * Fxpredf[[1]] + 
                                        p.f * exp(-r.i * .a) * Lxf * Fxpredf[[2]])
                p.m <- (SRB.i / (1 + SRB.i))
                p.f <- (1 / (1 + SRB.i))
                
                if (abs(delta.i) <= tol){
                    break
                }
                
            }
            
            if (i == maxit){
                cat("Warning! maxit reached. result may not be converged")
            }
            c(r = r.i, SRB = SRB.i)
        })
#a <- rIPFit(Bxym, Bxyf, Exm, Exf, Lxm, Lxf, tol = 1e-15)
#        plot(a$r.vec, type = 'l')

rIPFUSavg <- do.call(rbind, lapply(as.character(yearsUS), function(yr, .Bxymf, .Ex, .Lxm, .Lxf, .M){
                    rIPFit(Bxym = .Bxymf[[yr]][["Bxym"]],
                           Bxyf = .Bxymf[[yr]][["Bxyf"]],
                           Exm = with(.Ex, Male[Year == as.integer(yr)]),
                           Exf = with(.Ex, Female[Year == as.integer(yr)]),
                           Lxm = .Lxm[,yr],
                           Lxf = .Lxf[,yr], M = .M, tol = 1e-15)
                },.Bxymf = BxymfUS, .Ex = ExUS, .Lxm = LxmUS, .Lxf = LxfUS, .M = mean))



rIPFUShm <- do.call(rbind, lapply(as.character(yearsUS), function(yr, .Bxymf, .Ex, .Lxm, .Lxf, .M){
                    rIPFit(Bxym = .Bxymf[[yr]][["Bxym"]],
                            Bxyf = .Bxymf[[yr]][["Bxyf"]],
                            Exm = with(.Ex, Male[Year == as.integer(yr)]),
                            Exf = with(.Ex, Female[Year == as.integer(yr)]),
                            Lxm = .Lxm[,yr],
                            Lxf = .Lxf[,yr], M = .M, tol = 1e-15)
                },.Bxymf = BxymfUS, .Ex = ExUS, .Lxm = LxmUS, .Lxf = LxfUS, .M = harmonic.mean))


rIPFESavg <- do.call(rbind, lapply(as.character(yearsES), function(yr, .Bxymf, .Ex, .Lxm, .Lxf, .M){
                    rIPFit(Bxym = .Bxymf[[yr]][["Bxym"]],
                            Bxyf = .Bxymf[[yr]][["Bxyf"]],
                            Exm = with(.Ex, Male[Year == as.integer(yr)]),
                            Exf = with(.Ex, Female[Year == as.integer(yr)]),
                            Lxm = .Lxm[,yr],
                            Lxf = .Lxf[,yr], , M = .M, tol = 1e-15)
                },.Bxymf = BxymfES, .Ex = ExES, .Lxm = LxmES, .Lxf = LxfES, .M = mean))
rIPFEShm <- do.call(rbind, lapply(as.character(yearsES), function(yr, .Bxymf, .Ex, .Lxm, .Lxf, .M){
                    rIPFit(Bxym = .Bxymf[[yr]][["Bxym"]],
                            Bxyf = .Bxymf[[yr]][["Bxyf"]],
                            Exm = with(.Ex, Male[Year == as.integer(yr)]),
                            Exf = with(.Ex, Female[Year == as.integer(yr)]),
                            Lxm = .Lxm[,yr],
                            Lxf = .Lxf[,yr], , M = .M, tol = 1e-15)
                },.Bxymf = BxymfES, .Ex = ExES, .Lxm = LxmES, .Lxf = LxfES, .M = harmonic.mean))
#plot(yearsUS,rIPFUSavg[,1],type = 'l')
#lines(yearsUS,rIPFUShm[,1])
#
#plot(yearsUS,rIPFUSavg[,1]/rIPFUShm[,1], ylim = c(.9,1.1))

#plot(yearsUS,rIPFUS[,1],type = 'l')
#polygon(c(yearsUS,rev(yearsUS)),c(rLotkaUS[, "r.m"],rev(rLotkaUS[, "r.f"])), border = NA, col = "#55555550")

rLotkaUS <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Lxm, .Lxf, .Bxy, .Ex){
                    Fxm <- Minf0(Mna0(rowSums(.Bxy[[yr]][["Bxym"]]) / with(.Ex, Male[Year == as.integer(yr)])))
                    Fxf <- Minf0(Mna0(colSums(.Bxy[[yr]][["Bxyf"]]) / with(.Ex, Female[Year == as.integer(yr)])))
                    r.m <- LotkaRCoale(fx = Fxm, Lx = .Lxm[,yr], x = .5:110.5)
                    r.f <- LotkaRCoale(fx = Fxf, Lx = .Lxf[,yr], x = .5:110.5)
                    c(r.m = r.m, r.f = r.f)
                },.Lxm = LxmUS, .Lxf = LxfUS, .Bxy = BxymfUS, .Ex = ExUS))

rLotkaES <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Lxm, .Lxf, .Bxy, .Ex){
                    Fxm <- Minf0(Mna0(rowSums(.Bxy[[yr]][["Bxym"]]) / with(.Ex, Male[Year == as.integer(yr)])))
                    Fxf <- Minf0(Mna0(colSums(.Bxy[[yr]][["Bxyf"]]) / with(.Ex, Female[Year == as.integer(yr)])))
                    r.m <- LotkaRCoale(fx = Fxm, Lx = .Lxm[,yr], x = .5:110.5)
                    r.f <- LotkaRCoale(fx = Fxf, Lx = .Lxf[,yr], x = .5:110.5)
                    c(r.m = r.m, r.f = r.f)
                },.Lxm = LxmES, .Lxf = LxfES, .Bxy = BxymfES, .Ex = ExES))

# ----------------------------------------------------------
# plot results. note that gupta not bracketted:
pdf("/home/triffe/git/DISS/latex/Figures/IPFager.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rIPFUShm[,1], type = 'n', ylim = c(-.02, .015), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968,-.02,2010,.015,col = gray(.95), border=NA),
                abline(h = seq(-.02, .015, by = .005), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.02, .015, by = .005),seq(-.02, .02, by = .005), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.02, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.0225, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1965.5,.018, "r", cex = 1, xpd = TRUE)))
# rm-rf regions
polygon(c(yearsUS,rev(yearsUS)),c(rLotkaUS[, "r.m"],rev(rLotkaUS[, "r.f"])), border = NA, col = "#55555550")

polygon(c(yearsES,rev(yearsES)),c(rLotkaES[, "r.m"],rev(rLotkaES[, "r.f"])), border = NA, col = "#55555550")
# US results
lines(yearsUS,rIPFUShm[,1], lwd = 2, col = gray(.2))
#lines(yearsUS, USrPollard, col = "red")
lines(yearsUS, rLotkaUS[, "r.m"], lwd = 1, col = gray(.2), lty = 5)
lines(yearsUS, rLotkaUS[, "r.f"], lwd = 1, col = gray(.2), lty = 5)
# add points for years not bracketted:
USind <- rIPFUShm[,1] > rLotkaUS[, "r.m"] & rIPFUShm[,1] < rLotkaUS[, "r.f"]
points(yearsUS[USind], rIPFUShm[USind,1], pch = 9, col = "black", xpd = TRUE)

# Spain results
lines(yearsES,rIPFEShm[,1], lwd = 2, col = gray(.2))
#lines(yearsES, ESrPollard, col = "red")
lines(yearsES, rLotkaES[, "r.m"], lwd = 1, col = gray(.2), lty = 5)
lines(yearsES, rLotkaES[, "r.f"], lwd = 1, col = gray(.2), lty = 5)
#ESind <- rIPFEShm[,1] > rLotkaES[, "r.m"] & rIPFEShm[,1] < rLotkaES[, "r.f"]
#points(yearsES[ESind], rGuptaES[ESind], pch = 9, col = "black", xpd = TRUE)

# label rm rf
text(c(1990, 1986.5, 1973, 1971),
        c(-0.0103542581, -0.0141970362,  0.0003650703, -0.0040170451),
        c(expression(r^m~ES),expression(r^f~ES),expression(r^m~US),expression(r^f~US)),
        cex = .8, pos = c(4,1,4,1))
legend(1995, .015, pch = 9, bty = "n", legend = c(expression(r^IPF(hm) > r^(single~sex))), xpd = TRUE)

dev.off()

# ----------------------------------------------------------
# test monotonicity
Bxy <- BxymfUS[["1975"]][["Bxym"]] + BxymfUS[["1975"]][["Bxyf"]]
#
Fxf <- colSums(Bxy) / with(ExUS, Female[Year == 1975])
Fxm <- rowSums(Bxy) / with(ExUS, Male[Year == 1975])

FxPred <- IPFpred(Bxy, 
        Exm1 = with(ExUS, Male[Year == 1975]), 
        Exm2 = with(ExUS, 2*Male[Year == 1975]), 
        Exf1 = with(ExUS, Female[Year == 1975]), 
        Exf2 = with(ExUS, 2*Female[Year == 1975]), marM = harmonic.mean)

plot(0:110,FxPred[[1]],type='l')
lines(0:110,Fxm,lty=2,col="blue")

plot(FxPred[[1]]-Fxm)

FxPred <- IPFpred(Bxy, 
        Exm1 = with(ExUS, Male[Year == 1975]), 
        Exm2 = with(ExUS, .5*Male[Year == 1975]), 
        Exf1 = with(ExUS, Female[Year == 1975]), 
        Exf2 = with(ExUS, .5*Female[Year == 1975]), marM = harmonic.mean)
plot(0:110,FxPred[[1]],type='l')
lines(0:110,Fxm,lty=2,col="blue")

plot(FxPred[[1]]-Fxm)