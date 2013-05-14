# preamble:
setwd("/home/triffe/git/DISS/")
source("R/UtilityFunctions.R")
source("R/MeanFunctions.R")

# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
# BxES is 0:110, years 1975:2009
BxUS <- local(get(load("Data/USbirths/USBxy0_110.Rdata")))
BxES <- local(get(load("Data/ESbirths/ESBxy.Rdata")))

BxymfUS <- local(get(load("Data/USbirths/USBxymf0_110.Rdata")))
BxymfES <- local(get(load("Data/ESbirths/ESBxymf.Rdata")))
# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS <- local(get(load("Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("Data/Exposures/ESexp.Rdata")))
LxmUS <- local(get(load("Data/HMD_Lx/LxmUS.Rdata"))) / 1e5
LxfUS <- local(get(load("Data/HMD_Lx/LxfUS.Rdata"))) / 1e5
LxmES <- local(get(load("Data/HMD_Lx/LxmES.Rdata"))) / 1e5
LxfES <- local(get(load("Data/HMD_Lx/LxfES.Rdata"))) / 1e5
yearsUS <- 1969:2009
yearsES <- 1975:2009
#--------------------------
# Das Gupta 1972:


#Gupta1978Min <- function(r, Uxy, Vxy, lxm, lxf, mxy,.a = .5:110.5){
#    (1 - sum(mxy * (thetam * exp(-r * .a) * pm. * Uxy + t(thetaf * exp(-r * .a) * pf. * t(Vxy)))))^2
#}

Gupta1978It <- function(Uxy, Vxy, Lxm, Lxf, SRB, mxy,.a = .5:110.5, maxit = 200, tol = 1e-15){
    thetam <- SRB / (1 + SRB)
    thetaf <- 1 / (1 + SRB)
    R0.guess <- sum(mxy * (thetam * Lxm * Uxy + t(thetaf *  Lxf * t(Vxy))))
    T.guess <- sum(mxy * (.a * thetam * Lxm * Uxy + t(.a * thetaf *  Lxf * t(Vxy)))) / 
                      R0.guess 
    r.i <- log(R0.guess) / T.guess
    for (i in 1:maxit){
        deltai <- 1 - sum(mxy * (thetam * exp(-r.i * .a) * Lxm * Uxy + t(thetaf * exp(-r.i * .a) *  Lxf * t(Vxy))))
        r.i <- r.i - (deltai / ( T.guess  - (deltai / r.i)))
        if (abs(deltai) < tol){
            break
        }    
    }
    r.i
}



rGuptaUS <- unlist(lapply(as.character(yearsUS),function(yr, .Bxy, .Ex, .Lxm, .Lxf){
            Bm <- .Bxy[[yr]][["Bxym"]]
            Bf <- .Bxy[[yr]][["Bxyf"]]
            Bxy <- Bm + Bf
            SRB <- sum(Bm) / sum(Bf)
            PM <- with(.Ex, Male[Year == as.integer(yr)])
            PF <- with(.Ex, Female[Year == as.integer(yr)])               
            Uxy <- Minf0(Mna0(Bxy / rowSums(Bxy)))
            Vxy <- Minf0(Mna0(t(t(Bxy) / colSums(Bxy))))
          
            mxy <- Minf0(Mna0(Bxy /(Uxy * PM + t(t(Vxy) * PF))))
            Gupta1978It(Uxy = Uxy, Vxy = Vxy, Lxm = .Lxm[,yr], Lxf = .Lxf[,yr], SRB = SRB, mxy = mxy)
        }, .Bxy = BxymfUS, .Ex = ExUS, .Lxm = LxmUS, .Lxf = LxfUS))
rGuptaES <- unlist(lapply(as.character(yearsES),function(yr, .Bxy, .Ex, .Lxm, .Lxf){
                    Bm <- .Bxy[[yr]][["Bxym"]]
                    Bf <- .Bxy[[yr]][["Bxyf"]]
                    Bxy <- Bm + Bf
                    SRB <- sum(Bm) / sum(Bf)
                    PM <- with(.Ex, Male[Year == as.integer(yr)])
                    PF <- with(.Ex, Female[Year == as.integer(yr)])               
                    Uxy <- Minf0(Mna0(Bxy / rowSums(Bxy)))
                    Vxy <- Minf0(Mna0(t(t(Bxy) / colSums(Bxy))))
                    
                    mxy <- Minf0(Mna0(Bxy /(Uxy * PM + t(t(Vxy) * PF))))
                    Gupta1978It(Uxy = Uxy, Vxy = Vxy, Lxm = .Lxm[,yr], Lxf = .Lxf[,yr], SRB = SRB, mxy = mxy)
                }, .Bxy = BxymfES, .Ex = ExES, .Lxm = LxmES, .Lxf = LxfES))
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

#save(rGuptaUS, file = "Data/results/agerSRB/rGuptaUS.Rdata")
#save(rGuptaES, file = "Data/results/agerSRB/rGuptaES.Rdata")
#save(rLotkaUS, file = "Data/results/agerSRB/rLotkaUS.Rdata")
#save(rLotkaES, file = "Data/results/agerSRB/rLotkaES.Rdata")

sign(diff(rLotkaUS[,1]))==sign(diff(rLotkaUS[,2]))
sign(diff(rLotkaUS[,1]))==sign(diff(rGuptaUS))
sign(diff(rLotkaUS[,2]))==sign(diff(rGuptaUS))
sign(diff(rLotkaES[,1]))==sign(diff(rLotkaES[,2]))
sign(diff(rLotkaES[,1]))==sign(diff(rGuptaES))
sign(diff(rLotkaES[,2]))==sign(diff(rGuptaES))

# ----------------------------------------------------------
# plot results. note that gupta not bracketted:
pdf("latex/Figures/Gupta1978r.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rGuptaUS, type = 'n', ylim = c(-.02, .015), xlim = c(1968,2010), axes = FALSE,
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
lines(yearsUS,rGuptaUS, lwd = 2, col = gray(.2))
#lines(yearsUS, USrPollard, col = "red")
lines(yearsUS, rLotkaUS[, "r.m"], lwd = 1, col = gray(.2), lty = 5)
lines(yearsUS, rLotkaUS[, "r.f"], lwd = 1, col = gray(.2), lty = 5)
# add points for years not bracketted:
USind <- rGuptaUS > rLotkaUS[, "r.m"] & rGuptaUS < rLotkaUS[, "r.f"]
points(yearsUS[USind], rGuptaUS[USind], pch = 9, col = "black", xpd = TRUE)

# Spain results
lines(yearsES,rGuptaES, lwd = 2, col = gray(.2))
#lines(yearsES, ESrPollard, col = "red")
lines(yearsES, rLotkaES[, "r.m"], lwd = 1, col = gray(.2), lty = 5)
lines(yearsES, rLotkaES[, "r.f"], lwd = 1, col = gray(.2), lty = 5)
ESind <- rGuptaES > rLotkaES[, "r.m"] & rGuptaES < rLotkaES[, "r.f"]
points(yearsES[ESind], rGuptaES[ESind], pch = 9, col = "black", xpd = TRUE)

# label rm rf
text(c(1990, 1986.5, 1973, 1971),
        c(-0.0103542581, -0.0141970362,  0.0003650703, -0.0040170451),
        c(expression(r^m~ES),expression(r^f~ES),expression(r^m~US),expression(r^f~US)),
        cex = .8, pos = c(4,1,4,1))
#legend(1995, .015, pch = 9, bty = "n", legend = c(expression(r^Gupta > r^m)), xpd = TRUE)
dev.off()

# ----------------------------------------------------------
HM <- function(x,y){
    Minf0(Mna0((2 * x * y) / (x + y)))
}
AHM <- function(x,y){
    Minf0(Mna0(((x^2 + y^2) / 2) / ((x+y)/2)))
}
# compare harmonic:

harmonicIt <- function(Lxm, Lxf, SRB, mHM, .a = .5:110.5, maxit = 200, tol = 1e-15){
    thetam <- SRB / (1 + SRB)
    thetaf <- 1 / (1 + SRB)
    R0.guess <- sum(mHM * outer(thetam * Lxm, thetaf *  Lxf, HM))
    T.guess <- sum(mHM * outer(thetam * Lxm * .a, thetaf *  Lxf * .a, HM)) /  R0.guess
    r.i <- log(R0.guess) / T.guess
    for (i in 1:maxit){   
        deltai <- 1 - sum(mHM *  outer(thetam * exp(-r.i * .a) * Lxm, thetaf * exp(-r.i * .a) *  Lxf, HM))
        r.i <- r.i - (deltai / ( T.guess  - (deltai / r.i)))
        if (abs(deltai) < tol){
            break
        }    
    }
    r.i
}

antiharmonicIt <- function(Lxm, Lxf, SRB, mAHM, .a = .5:110.5, maxit = 200, tol = 1e-15){
    thetam <- SRB / (1 + SRB)
    thetaf <- 1 / (1 + SRB)
    R0.guess <- sum(mAHM *  outer(thetam * Lxm, thetaf*  Lxf, AHM))
    T.guess <-  sum(mAHM *  outer(thetam * Lxm * .a, thetaf *  Lxf* .a, AHM)) / R0.guess
    r.i <- log(R0.guess) / T.guess
    for (i in 1:maxit){   
        deltai <- 1 - sum(mAHM *  outer(thetam * exp(-r.i * .a) * Lxm, thetaf * exp(-r.i * .a) *  Lxf, AHM))
        r.i <- r.i - (deltai / ( T.guess  - (deltai / r.i))) 
        if (abs(deltai) < tol){
            break
        }    
    }
    r.i
}

rHMUS <- unlist(lapply(as.character(yearsUS),function(yr, .Bxy, .Ex, .Lxm, .Lxf){
                    Bm  <- .Bxy[[yr]][["Bxym"]]
                    Bf  <- .Bxy[[yr]][["Bxyf"]]
                    Bxy <- Bm + Bf
                    SRB <- sum(Bm) / sum(Bf)
                    PM  <- with(.Ex, Male[Year == as.integer(yr)])
                    PF  <- with(.Ex, Female[Year == as.integer(yr)])               
                    mm  <- Bxy / PM
                    mf  <- t(t(Bxy) / PF)
                    mHM <- HM(mm, mf)
                    harmonicIt(Lxm = .Lxm[,yr], Lxf = .Lxf[,yr], SRB = SRB, mHM = mHM)
                }, .Bxy = BxymfUS, .Ex = ExUS, .Lxm = LxmUS, .Lxf = LxfUS))
rAHMUS <- unlist(lapply(as.character(yearsUS),function(yr, .Bxy, .Ex, .Lxm, .Lxf){
                    Bm <- .Bxy[[yr]][["Bxym"]]
                    Bf <- .Bxy[[yr]][["Bxyf"]]
                    Bxy <- Bm + Bf
                    SRB <- sum(Bm) / sum(Bf)
                    PM <- with(.Ex, Male[Year == as.integer(yr)])
                    PF <- with(.Ex, Female[Year == as.integer(yr)])               
                    mm <- Bxy / PM
                    mf <- t(t(Bxy) / PF)
                    mAHM <- Minf0(Mna0(((mm^2 + mf^2) / 2)  / ((mm + mf) / 2)))
                    antiharmonicIt(Lxm = .Lxm[,yr], Lxf = .Lxf[,yr], SRB = SRB, mAHM = mAHM)
                }, .Bxy = BxymfUS, .Ex = ExUS, .Lxm = LxmUS, .Lxf = LxfUS))
rHMES <- unlist(lapply(as.character(yearsES),function(yr, .Bxy, .Ex, .Lxm, .Lxf){
                    Bm  <- .Bxy[[yr]][["Bxym"]]
                    Bf  <- .Bxy[[yr]][["Bxyf"]]
                    Bxy <- Bm + Bf
                    SRB <- sum(Bm) / sum(Bf)
                    PM  <- with(.Ex, Male[Year == as.integer(yr)])
                    PF  <- with(.Ex, Female[Year == as.integer(yr)])               
                    mm  <- Bxy / PM
                    mf  <- t(t(Bxy) / PF)
                    mHM <- HM(mm, mf)
                    harmonicIt(Lxm = .Lxm[,yr], Lxf = .Lxf[,yr], SRB = SRB, mHM = mHM)
                }, .Bxy = BxymfES, .Ex = ExES, .Lxm = LxmES, .Lxf = LxfES))
rAHMES <- unlist(lapply(as.character(yearsES),function(yr, .Bxy, .Ex, .Lxm, .Lxf){
                    Bm <- .Bxy[[yr]][["Bxym"]]
                    Bf <- .Bxy[[yr]][["Bxyf"]]
                    Bxy <- Bm + Bf
                    SRB <- sum(Bm) / sum(Bf)
                    PM <- with(.Ex, Male[Year == as.integer(yr)])
                    PF <- with(.Ex, Female[Year == as.integer(yr)])               
                    mm <- Bxy / PM
                    mf <- t(t(Bxy) / PF)
                    mAHM <- Minf0(Mna0(((mm^2 + mf^2) / 2)  / ((mm + mf) / 2)))
                    antiharmonicIt(Lxm = .Lxm[,yr], Lxf = .Lxf[,yr], SRB = SRB, mAHM = mAHM)
                }, .Bxy = BxymfES, .Ex = ExES, .Lxm = LxmES, .Lxf = LxfES))




par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rGuptaUS, type = 'n', ylim = c(-.02, .015), xlim = c(1968,2010), axes = FALSE,
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
lines(yearsUS,rGuptaUS, lwd = 2, col = gray(.2))
#lines(yearsUS, USrPollard, col = "red")
lines(yearsUS, rLotkaUS[, "r.m"], lwd = 1, col = gray(.2), lty = 5)
lines(yearsUS, rLotkaUS[, "r.f"], lwd = 1, col = gray(.2), lty = 5)
# add points for years not bracketted:
USind <- rGuptaUS > rLotkaUS[, "r.m"] | rGuptaUS < rLotkaUS[, "r.f"]
points(yearsUS[USind], rGuptaUS[USind], pch = 9, col = "black", xpd = TRUE)
lines(yearsUS, rHMUS, col = "red")
lines(yearsUS, rAHMUS, col = "blue")
# Spain results
lines(yearsES,rGuptaES, lwd = 2, col = gray(.2))
#lines(yearsES, ESrPollard, col = "red")
lines(yearsES, rLotkaES[, "r.m"], lwd = 1, col = gray(.2), lty = 5)
lines(yearsES, rLotkaES[, "r.f"], lwd = 1, col = gray(.2), lty = 5)
ESind <- rGuptaES > rLotkaES[, "r.m"] | rGuptaES < rLotkaES[, "r.f"]
points(yearsES[ESind], rGuptaES[ESind], pch = 9, col = "black", xpd = TRUE)
lines(yearsES, rHMES, col = "red")
lines(yearsES, rAHMES, col = "blue")
# label rm rf
text(c(1990, 1986.5, 1973, 1971),
        c(-0.0103542581, -0.0141970362,  0.0003650703, -0.0040170451),
        c(expression(r^m~ES),expression(r^f~ES),expression(r^m~US),expression(r^f~US)),
        cex = .8, pos = c(4,1,4,1))
legend(1995, .015, pch = 9, bty = "n", legend = c(expression(r^Gupta > r^m)), xpd = TRUE)

