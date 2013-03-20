# Figures here found in PopulationStructure.tex
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")


# exposures, as such, straight from HMD, all ages 0-110, long form
PxUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxUS.Rdata")))
PxES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxES.Rdata")))
dxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfES.Rdata"))) 
# make sum to 1
dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)

yearsUS <- 1969:2009
yearsES <- 1975:2009

StabUS <- do.call(rbind,lapply(1969:2008, function(yr, .PxUS, .dxmUS, .dxfUS){
            .yr1 <- as.character(yr)
            .yr2 <- as.character(yr+1)
            Pxm <- .PxUS$Male[.PxUS$Year == yr]
            Pxf <- .PxUS$Female[.PxUS$Year == yr]
            Mex <- rowSums(ExpectedDx(Pxm, .dxmUS[,.yr1]))
            Fex <- rowSums(ExpectedDx(Pxf, .dxfUS[,.yr1]))
            
            cx1 <- c(Pxm,Pxf) / sum( c(Pxm,Pxf))
            cEx1 <- c(Mex,Fex) / sum( c(Mex,Fex))
            
            Pxm <- .PxUS$Male[.PxUS$Year == (yr+1)]
            Pxf <- .PxUS$Female[.PxUS$Year == (yr+1)]
            Mex <- rowSums(ExpectedDx(Pxm, .dxmUS[,.yr2]))
            Fex <- rowSums(ExpectedDx(Pxf, .dxfUS[,.yr2]))
            
            cx2 <- c(Pxm,Pxf) / sum( c(Pxm,Pxf))
            cEx2 <- c(Mex,Fex) / sum( c(Mex,Fex))
            
           1-c(Age = sum(pmin(cx1, cx2)), Ex = sum(pmin(cEx1, cEx2)))
            
        }, .PxUS = PxUS, .dxmUS = dxmUS, .dxfUS = dxfUS))
StabES <- do.call(rbind,lapply(1975:2008, function(yr, .PxES, .dxmES, .dxfES){
                    .yr1 <- as.character(yr)
                    .yr2 <- as.character(yr+1)
                    Pxm <- .PxES$Male[.PxES$Year == yr]
                    Pxf <- .PxES$Female[.PxES$Year == yr]
                    Mex <- rowSums(ExpectedDx(Pxm, .dxmES[,.yr1]))
                    Fex <- rowSums(ExpectedDx(Pxf, .dxfES[,.yr1]))
                    
                    cx1 <- c(Pxm,Pxf) / sum( c(Pxm,Pxf))
                    cEx1 <- c(Mex,Fex) / sum( c(Mex,Fex))
                    
                    Pxm <- .PxES$Male[.PxES$Year == (yr+1)]
                    Pxf <- .PxES$Female[.PxES$Year == (yr+1)]
                    Mex <- rowSums(ExpectedDx(Pxm, .dxmES[,.yr2]))
                    Fex <- rowSums(ExpectedDx(Pxf, .dxfES[,.yr2]))
                    
                    cx2 <- c(Pxm,Pxf) / sum( c(Pxm,Pxf))
                    cEx2 <- c(Mex,Fex) / sum( c(Mex,Fex))
                    
                    1-c(Age = sum(pmin(cx1, cx2)), Ex = sum(pmin(cEx1, cEx2)))
                    
                }, .PxES = PxES, .dxmES = dxmES, .dxfES = dxfES))

#
#plot(1969:2008, StabUS[,1], type = 'l', ylim = c(0,.02), col = "red")
#lines(1969:2008, StabUS[,2], col = "blue")
#lines(1975:2008, StabES[,1], col = "red", lty = 2)
#lines(1975:2008, StabES[,2], col = "blue", lty = 2)

pdf("/home/triffe/git/DISS/latex/Figures/PyramidStabilityThetaRatio.pdf", height = 5, width = 5)
par(mai = c(.5, .3, .6, .2),xaxs = "i", yaxs = "i")
plot(1969:2008, StabUS[,2]/StabUS[,1], type = 'l', ylim = c(0,.8), xlim = c(1968,2009), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968,0,2009,1,col = gray(.95), border=NA),
                abline(h = seq(0,.8,by = .1), col = "white"),
                abline(v = seq(1970, 2009, by = 5), col = "white"),
                text(1968, seq(0, .8, by = .1),seq(0, .8, by = .1), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2009, by = 10),0, seq(1970, 2009, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.05, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1968.5,.87, expression(frac(theta~e[x],theta~age)), cex = 1, xpd = TRUE)))

lines(1975:2008, StabES[,2]/StabES[,1], lwd = 2.5, col = gray(.5), lty = 5)

legend(1968,.78, lty = c(1,5), col = gray(c(.2,.5)), lwd = c(2,2.5),bty = "n",
        legend = c("US","Spain"), xpd = TRUE)
dev.off()


# ---------------------------------------------------------
# uncertainty in pyramid. self-sufficient code

PxUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxUS.Rdata")))
PxES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxES.Rdata")))
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
mxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxmUS.Rdata"))) 
mxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxfUS.Rdata"))) 
mxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxmES.Rdata"))) 
mxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxfES.Rdata"))) 

yearsUS <- 1969:2009
yearsES <- 1975:2009

PyramidProbUS <- lapply(as.character(yearsUS), function(yr, .Px, .Ex, .mxm, .mxf, .dxm, .dxf){
            
            Pxm <- .Px$Male[.Px$Year == as.integer(yr)]
            Pxf <- .Px$Female[.Px$Year == as.integer(yr)] 
            Exm <- .Ex$Male[.Ex$Year == as.integer(yr)]
            Exf <- .Ex$Female[.Ex$Year == as.integer(yr)]
            mxm <- .mxm[, yr]
            mxf <- .mxf[, yr]
            
            dxmSim <- apply(Mna0(Minf0(matrix(rpois(n = length(mxm)*1000, lambda = mxm * Exm), ncol = 1000) 
                                    / Exm)), 2, mx2dxHMD)
            dxfSim <- apply(Mna0(Minf0(matrix(rpois(n = length(mxf)*1000, lambda = mxf * Exf), ncol = 1000) 
                                    / Exf)), 2, mx2dxHMD)
            
            # redistributed and scaled to sum to 1
            Scale <- sum(c(Pxm, Pxf))
            Pymsc <- apply(dxmSim, 2, function(dx, .Pxm){
                        rowSums(ExpectedDx(.Pxm, dx))
                    }, .Pxm = Pxm) / Scale
            Pyfsc <- apply(dxfSim, 2, function(dx, .Pxf){
                        rowSums(ExpectedDx(.Pxf, dx))
                    }, .Pxf = Pxf) / Scale
          
            # now take the .025 and .975 quantiles by remaining years
            Qm <- t(apply(Pymsc, 1, function(x){
                        quantile(x, probs = c(.025, .975))
                    }))
            Qf <- t(apply(Pyfsc, 1, function(x){
                        quantile(x, probs = c(.025, .975))
                    }))
            Obsm <- rowSums(ExpectedDx(Pxm, .dxm[,yr]) / Scale)
            Obsf <- rowSums(ExpectedDx(Pxf, .dxf[,yr]) / Scale)
            
            cbind(Ml = Qm[,1], Mu = Qm[,2], Fl = Qf[,1], Fu = Qf[,2], Om = Obsm, Of = Obsf)
        },  .Px = PxUS, .Ex = ExUS, .mxm = mxmUS, .mxf = mxfUS, .dxm = dxmUS, .dxf = dxfUS)
PyramidProbES <- lapply(as.character(yearsES), function(yr, .Px, .Ex, .mxm, .mxf, .dxm, .dxf){
                    
                    Pxm <- .Px$Male[.Px$Year == as.integer(yr)]
                    Pxf <- .Px$Female[.Px$Year == as.integer(yr)] 
                    Exm <- .Ex$Male[.Ex$Year == as.integer(yr)]
                    Exf <- .Ex$Female[.Ex$Year == as.integer(yr)]
                    mxm <- .mxm[, yr]
                    mxf <- .mxf[, yr]
                    
                    dxmSim <- apply(Mna0(Minf0(matrix(rpois(n = length(mxm)*1000, lambda = mxm * Exm), ncol = 1000) 
                                                    / Exm)), 2, mx2dxHMD)
                    dxfSim <- apply(Mna0(Minf0(matrix(rpois(n = length(mxf)*1000, lambda = mxf * Exf), ncol = 1000) 
                                                    / Exf)), 2, mx2dxHMD)
                    
                    # redistributed and scaled to sum to 1
                    Scale <- sum(c(Pxm, Pxf))
                    Pymsc <- apply(dxmSim, 2, function(dx, .Pxm){
                                rowSums(ExpectedDx(.Pxm, dx))
                            }, .Pxm = Pxm) / Scale
                    Pyfsc <- apply(dxfSim, 2, function(dx, .Pxf){
                                rowSums(ExpectedDx(.Pxf, dx))
                            }, .Pxf = Pxf) / Scale
                    
                    # now take the .025 and .975 quantiles by remaining years
                    Qm <- t(apply(Pymsc, 1, function(x){
                                        quantile(x, probs = c(.025, .975))
                                    }))
                    Qf <- t(apply(Pyfsc, 1, function(x){
                                        quantile(x, probs = c(.025, .975))
                                    }))
                    Obsm <- rowSums(ExpectedDx(Pxm, .dxm[,yr]) / Scale)
                    Obsf <- rowSums(ExpectedDx(Pxf, .dxf[,yr]) / Scale)
                    
                    cbind(Ml = Qm[,1], Mu = Qm[,2], Fl = Qf[,1], Fu = Qf[,2], Om = Obsm, Of = Obsf)
                },  .Px = PxES, .Ex = ExES, .mxm = mxmES, .mxf = mxfES, .dxm = dxmES, .dxf = dxfES)
        
        
TestYr <- PyramidProbES[[1]]

#pdf("/home/triffe/git/DISS/latex/Figures/exPyramidUS.pdf", height = 5, width = 5)
par(mai = c(.6,.6,.3,.3), xaxs = "i", yaxs = "i")
plot(NULL, type = "n",axes = FALSE, xlab = "",ylab = "", xlim = c(-1, 1), ylim = c(0,111),
     panel.first = list(
        rect(-1, 0, 1, 111, col = gray(.95), border = NA),
        abline(v = seq(-1, 1, by = .2), col = "white"),
        abline(h = seq(0, 110, by = 10), col = "white"),
        text(seq(-1, 1, by = .2),0, xlabs, xpd = TRUE, pos = 1, cex = .7),
        text(-1, seq(0, 110, by = 10), seq(0, 110, by = 10), pos = 2, xpd = TRUE, cex = .7),
        text(-1.17, 116, expression(e[x]), xpd = TRUE, cex = 1),
        text(0, -12, "Percentage", xpd = TRUE, cex = 1),
      ))

PyramidOutline(TestYr[,"Ml"], TestYr[,"Fl"], scale =100, border = gray(.2), xpd = TRUE, lwd = .5) 
PyramidOutline(TestYr[,"Mu"], TestYr[,"Fu"], scale =100, border = gray(.2), xpd = TRUE, lwd = .5) 
text(c(-.5,-.55), c(50, 91), c("1975", "2009"), col = c("white", "black"), cex = 1.2, pos = 4)
segments(-.4, 88, -.3, 83)
#dev.off()

names(PyramidProbES) <- yearsES
names(PyramidProbUS) <- yearsUS
plot(0:110,PyramidProbES[["1975"]][,"Fu"]/PyramidProbES[["1975"]][,"Fl"], type ='l', ylim = c(1,1.1),
        col = gray(.4),lty=4)
lines(0:110,PyramidProbES[["1975"]][,"Mu"]/PyramidProbES[["1975"]][,"Ml"],col = gray(.4),lty=4)
lines(0:110,PyramidProbUS[["1975"]][,"Fu"]/PyramidProbUS[["1975"]][,"Fl"], col = gray(.2))
lines(0:110,PyramidProbUS[["1975"]][,"Mu"]/PyramidProbUS[["1975"]][,"Ml"], col = gray(.2))


MUSu <- (PyramidProbUS[["1975"]][,"Mu"]-PyramidProbUS[["1975"]][,"Om"]) / PyramidProbUS[["1975"]][,"Om"]
MUSl <- (PyramidProbUS[["1975"]][,"Ml"] - PyramidProbUS[["1975"]][,"Om"]) / PyramidProbUS[["1975"]][,"Om"]
FUSu <- (PyramidProbUS[["1975"]][,"Fu"]-PyramidProbUS[["1975"]][,"Of"]) / PyramidProbUS[["1975"]][,"Of"]
FUSl <- (PyramidProbUS[["1975"]][,"Fl"] - PyramidProbUS[["1975"]][,"Of"]) / PyramidProbUS[["1975"]][,"Of"]
MESu <- (PyramidProbES[["1975"]][,"Mu"]-PyramidProbES[["1975"]][,"Om"]) / PyramidProbES[["1975"]][,"Om"]
MESl <- (PyramidProbES[["1975"]][,"Ml"] - PyramidProbES[["1975"]][,"Om"]) / PyramidProbES[["1975"]][,"Om"]
FESu <- (PyramidProbES[["1975"]][,"Fu"]-PyramidProbES[["1975"]][,"Of"]) / PyramidProbES[["1975"]][,"Of"]
FESl <- (PyramidProbES[["1975"]][,"Fl"] - PyramidProbES[["1975"]][,"Of"]) / PyramidProbES[["1975"]][,"Of"]


pdf("/home/triffe/git/DISS/latex/Figures/exPyramiduncertainty1975.pdf", height = 5, width = 5)
par(mai = c(.5,.5,.5,.3), xaxs = "i", yaxs = "i")
ylim <- c(-.1,.1)
plot(NULL, type = "n",axes = FALSE, xlab = "",ylab = "", xlim = c(0, 100), ylim = ylim,
        panel.first = list(
                rect(0, ylim[1], 100, ylim[2], col = gray(.95), border = NA),
                abline(v = seq(0, 100, by = 10), col = "white"),
                abline(h = seq(ylim[1], ylim[2], by = .02), col = "white"),
                text(seq(0, 100, by = 10),ylim[1], seq(0, 100, by = 10), xpd = TRUE, pos = 1, cex = .7),
                text(0, seq(ylim[1], ylim[2], by = .02), paste(100*seq(ylim[1], ylim[2], by = .02),"%"), pos = 2, xpd = TRUE, cex = .7),
                text(50, ylim[1]-.02, "Remaining Years", xpd = TRUE, cex = 1),
                text(-10, .115, "Interval %", xpd = TRUE, cex = 1, pos=4)
        ))

lines(0:110,MUSu, col = gray(.2), lwd = 2)
lines(0:110,MUSl, col = gray(.2), lwd = 2)
lines(0:110,FUSu, col = gray(.2), lwd = 2)
lines(0:110,FUSl, col = gray(.2), lwd = 2)
lines(0:110,MESu, col = gray(.4), lwd = 2, lty = 4)
lines(0:110,MESl, col = gray(.4), lwd = 2, lty = 4)
lines(0:110,FESu, col = gray(.4), lwd = 2, lty = 4)
lines(0:110,FESl, col = gray(.4), lwd = 2, lty = 4)

text(c(82,82,82,82), c(.05,.04,-.04,-.05), c("US males","US females","ES females","ES males"), pos =2)
segments(c(82,82,82,82),c(.05,.04,-.04,-.05), c(98,98,93,93), c(MUSu[99], FUSu[99],FESl[94],MESl[94]))
text(c(40,40),c(-.03,.03),c("Lower","Upper"), cex = 1.2)
dev.off()
#
#plot((MUSu/FUSu)[1:10])
#plot((MESu/FESu)[1:10])
#plot((MUSu/FUSu)[90:100])
#plot((MESu/FESu)[90:100])




# -----------------------------------------------------------
# I'm very confused about this...
# self-sufficient code.
# take improvement into account deterministically
.Px = PxUS
.Ex = ExUS
.mxm = mxmUS
.mxf = mxfUS
.dxm = dxmUS
.dxf = dxfUS

mx <- mxmUS[,"1975"]
improv.prop <- .991
Px <- PxUS$Male[PxUS$Year == 1975]
ExpectedDxImpr <- function(Px, mx, improv.prop){
    N        <- length(mx)
    mxi      <- mx
    dxiImprov <- matrix(nrow = N, ncol = N)
    lxiImprov <- matrix(nrow = N, ncol = N)
    for (i in 1:N){
        lxiImprov[,i] <- rev(cumsum(rev(mx2dxHMD(mxi))))
        dxiImprov[,i] <- mx2dxHMD(mxi)
        mxi           <- mxi * improv.prop
    }

}
dx[1,]
dxiImprov[1,]
Pym1p <- rowSums(ExpectedDxImpr(PxUS$Male[PxUS$Year == 1975],  mxmUS[,"1975"], improv.prop = .995))
Pyf1p <- rowSums(ExpectedDxImpr(PxUS$Female[PxUS$Year == 1975],  mxfUS[,"1975"], improv.prop = .995))
Pyf <- rowSums(ExpectedDx(PxUS$Female[PxUS$Year == 1975], dxfUS[,"1975"]))
Pym <- rowSums(ExpectedDx(PxUS$Male[PxUS$Year == 1975], dxmUS[,"1975"]))
par(mai = c(.6,.6,.3,.3), xaxs = "i", yaxs = "i")
plot(NULL, type = "n",axes = FALSE, xlab = "",ylab = "", xlim = c(-1, 1), ylim = c(0,111),
        panel.first = list(
                rect(-1, 0, 1, 111, col = gray(.95), border = NA),
                abline(v = seq(-1, 1, by = .2), col = "white"),
                abline(h = seq(0, 110, by = 10), col = "white"),
                text(seq(-1, 1, by = .2),0, xlabs, xpd = TRUE, pos = 1, cex = .7),
                text(-1, seq(0, 110, by = 10), seq(0, 110, by = 10), pos = 2, xpd = TRUE, cex = .7),
                text(-1.17, 116, expression(e[x]), xpd = TRUE, cex = 1),
                text(0, -12, "Percentage", xpd = TRUE, cex = 1),
        ))
PyramidOutline(Pym1p, Pyf1p, scale =100, border = "red", xpd = TRUE, lwd = .5) 
PyramidOutline(Pym, Pyf, scale =100, border = "blue", xpd = TRUE, lwd = .5)






