
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

yearsUS <- 1969:2009
yearsES <- 1975:2009

# Bxy totals  (not used()
#BxUS  <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
#BxES  <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

# repeated below for better SRB assumptions
BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))
names(BxymfES) <- yearsES
names(BxymfUS) <- yearsUS
# exposures, as such, straight from HMD, all ages 0-110, long form
ExUS  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
#------------------------------------------------------------
# compare with Lotka:
LxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmUS.Rdata"))) / 1e5
LxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfUS.Rdata"))) / 1e5
LxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmES.Rdata"))) / 1e5
LxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfES.Rdata"))) / 1e5

#------------------------------------
yr <- "1969"
FxMM <- rowSums(BxymfUS[[yr]][["Bxym"]]) / with(ExUS,Male[Year == as.integer(yr)])
FxMF <- rowSums(BxymfUS[[yr]][["Bxyf"]]) / with(ExUS,Male[Year == as.integer(yr)])
FxFF <- colSums(BxymfUS[[yr]][["Bxyf"]]) / with(ExUS,Female[Year == as.integer(yr)])
FxFM <- colSums(BxymfUS[[yr]][["Bxym"]]) / with(ExUS,Female[Year == as.integer(yr)])



GoodmanMin <- function(r, Lxm, Lxf, FxMM, FxMF, FxFF, FxFM, SRB, sigma = .5, .a = .5:110.5){
    p.m <- SRB / (1 + SRB)
    p.f <- 1 / (1 + SRB)
    (1- (sum(sigma * p.m * Lxm * exp(-r * .a) * FxMM + (1-sigma) * p.f * Lxf * exp(-r * .a) * FxFM) +    # boy births
         sum(sigma * p.m * Lxm * exp(-r * .a) * FxMF + (1-sigma) * p.f * Lxf * exp(-r * .a) * FxFF) )) ^ 2 # girl births
}
SRB <- sum(BxymfUS[[yr]][["Bxym"]]) / sum(BxymfUS[[yr]][["Bxyf"]])
optimize(GoodmanMin, interval = c(-.02,.02),
        Lxm = LxmUS[,yr],
        Lxf = LxfUS[,yr],
        FxMM = FxMM, FxMF = FxMF, FxFF = FxFF, FxFM = FxFM,
        SRB = SRB, sigma = .5)

GoodmanIt <- compiler::cmpfun(function(Lxm, Lxf, FxMM, FxMF, FxFF, FxFM, 
                sigma = .5, .a = .5:110.5, tol = 1e-15, maxit = 2e2){
            
            SRBi        <- sum((1.05 / 2.05) * Lxm * sigma * FxMM + 
                            (1 / 2.05) * Lxf * (1 - sigma) * FxFM) / 
                           sum((1.05 / 2.05) * Lxm * sigma * FxMF + 
                            (1 / 2.05) * Lxf * (1 - sigma) * FxFF)
            # image(outer(Lxm * SRBi / (SRBi+1), Lxf * 1 / (SRBi+1), M) )
            R0          <- sum((SRBi / (SRBi + 1)) * Lxm * sigma * FxMM + 
                                    (1 / (SRBi + 1))  * Lxf * (1 - sigma) * FxFM) + 
                           sum((SRBi / (SRBi + 1))  * Lxm * sigma * FxMF + 
                                    (1 / (SRBi + 1))  * Lxf * (1 - sigma) * FxFF)
            T.guess <-    (sum(.a * (SRBi / (SRBi + 1)) * Lxm * sigma * FxMM + 
                                    .a * (1 / (SRBi + 1))  * Lxf * (1 - sigma) * FxFM) + 
                           sum(.a *(SRBi / (SRBi + 1))  * Lxm * sigma * FxMF + 
                                    .a * (1 / (SRBi + 1))  * Lxf * (1 - sigma) * FxFF)) / R0
            # first assuming a mean generation time of 35
            r.i <- log(R0) /  T.guess
            
            for (i in 1:maxit){
                
                p.m     <- (SRBi / (SRBi+1))
                p.f     <- (1 / (SRBi+1)) 
                deltai  <- 1 - sum(p.m * Lxm * exp(-r.i * .a) * sigma * (FxMM + FxMF) + 
                                p.f * Lxf * exp(-r.i * .a) * (1 - sigma) * (FxFM + FxFF))
                r.i     <- r.i - (deltai / (T.guess - (deltai / r.i)))
                
                SRBi    <-  sum(p.m * Lxm * exp(-r.i * .a) * sigma * FxMM + 
                                        p.f * Lxf * exp(-r.i * .a) * (1 - sigma) * FxFM) / 
                            sum(p.m * Lxm * exp(-r.i * .a) * sigma * FxMF + 
                                        p.f * Lxf * exp(-r.i * .a) * (1 - sigma) * FxFF)
                if (abs(deltai) < tol){
                    break
                }
            }
            c(r = r.i, SRB = SRBi)
        })

rUS <- do.call(rbind, lapply(as.character(yearsUS), function(yr, .Bxymf, .Ex, .Lxm, .Lxf){
                    
                    FxMM <- Minf0(Mna0(rowSums(.Bxymf[[yr]][["Bxym"]]) / with(.Ex,Male[Year == as.integer(yr)])))
                    FxMF <- Minf0(Mna0(rowSums(.Bxymf[[yr]][["Bxyf"]]) / with(.Ex,Male[Year == as.integer(yr)])))
                    FxFF <- Minf0(Mna0(colSums(.Bxymf[[yr]][["Bxyf"]]) / with(.Ex,Female[Year == as.integer(yr)])))
                    FxFM <- Minf0(Mna0(colSums(.Bxymf[[yr]][["Bxym"]]) / with(.Ex,Female[Year == as.integer(yr)])))
                    Lxm = .Lxm[,yr]
                    Lxf = .Lxf[,yr]
                    r.f <- LotkaRCoale(FxFF, Lxf, x = .5:110.5)
                    r.m <- LotkaRCoale(FxMM, Lxm, x = .5:110.5)
                    g0 <- GoodmanIt(Lxm = Lxm, Lxf = Lxf, 
                              FxMM = FxMM, FxMF = FxMF, FxFF = FxFF, FxFM = FxFM,
                              sigma = 0)
                    g1 <- GoodmanIt(Lxm = Lxm, Lxf = Lxf, 
                              FxMM = FxMM, FxMF = FxMF, FxFF = FxFF, FxFM = FxFM,
                              sigma = 1)
                    g0.5 <- GoodmanIt(Lxm = Lxm, Lxf = Lxf, 
                              FxMM = FxMM, FxMF = FxMF, FxFF = FxFF, FxFM = FxFM,
                              sigma = .5)
                    c(r.f = r.f, r.m = r.m, 
                            r.0 = g0["r"],
                            r.1 = g1["r"],
                            r0.5 = g0.5["r"],
                            SRB0 = g0["SRB"],
                            SRB1 = g1["SRB"],
                            SRB0.5 = g0.5["SRB"])
                }, .Bxymf = BxymfUS, .Ex = ExUS, .Lxm = LxmUS, .Lxf = LxfUS))

yr <- "1988"
rES <- do.call(rbind, lapply(as.character(1975:2009), function(yr, .Bxymf, .Ex, .Lxm, .Lxf){
                   
                    FxMM <- Minf0(Mna0(rowSums(.Bxymf[[yr]][["Bxym"]]) / with(.Ex,Male[Year == as.integer(yr)])))
                    FxMF <- Minf0(Mna0(rowSums(.Bxymf[[yr]][["Bxyf"]]) / with(.Ex,Male[Year == as.integer(yr)])))
                    FxFF <- Minf0(Mna0(colSums(.Bxymf[[yr]][["Bxyf"]]) / with(.Ex,Female[Year == as.integer(yr)])))
                    FxFM <- Minf0(Mna0(colSums(.Bxymf[[yr]][["Bxym"]]) / with(.Ex,Female[Year == as.integer(yr)])))
                    Lxm = .Lxm[,yr]
                    Lxf = .Lxf[,yr]
                    r.f <- LotkaRCoale(FxFF, Lxf, x = .5:110.5)
                    r.m <- LotkaRCoale(FxMM, Lxm, x = .5:110.5)
                    g0 <- GoodmanIt(Lxm = Lxm, Lxf = Lxf, 
                            FxMM = FxMM, FxMF = FxMF, FxFF = FxFF, FxFM = FxFM,
                            sigma = 0)
                    g1 <- GoodmanIt(Lxm = Lxm, Lxf = Lxf, 
                            FxMM = FxMM, FxMF = FxMF, FxFF = FxFF, FxFM = FxFM,
                            sigma = 1)
                    g0.5 <- GoodmanIt(Lxm = Lxm, Lxf = Lxf, 
                            FxMM = FxMM, FxMF = FxMF, FxFF = FxFF, FxFM = FxFM,
                            sigma = .5)
                    c(r.f = r.f, r.m = r.m, 
                            r.0 = g0["r"],
                            r.1 = g1["r"],
                            r0.5 = g0.5["r"],
                            SRB0 = g0["SRB"],
                            SRB1 = g1["SRB"],
                            SRB0.5 = g0.5["SRB"])
                
                }, .Bxymf = BxymfES, .Ex = ExES, .Lxm = LxmES, .Lxf = LxfES))

rownames(rUS) <- yearsUS
rownames(rES) <- yearsES

# plot it
pdf("/home/triffe/git/DISS/latex/Figures/Goodmanager.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rUS[,"r.f"], type = 'n', ylim = c(-.02, .015), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968,-.02,2010,.015,col = gray(.95), border=NA),
                abline(h = seq(-.02, .015, by = .005), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.02, .015, by = .005),seq(-.02, .02, by = .005), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.02, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.0225, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1965.5,.018, "r", cex = 1, xpd = TRUE)))
# rm-rf regions
polygon(c(yearsUS,rev(yearsUS)),c(rUS[, "r.m"],rev(rUS[, "r.f"])), border = NA, col = "#55555550")
polygon(c(yearsES,rev(yearsES)),c(rES[, "r.m"],rev(rES[, "r.f"])), border = NA, col = "#55555550")
# US results
lines(yearsUS, rUS[, "r.m"], lwd = 1, col = gray(.2), lty = 5)
lines(yearsUS, rUS[, "r.f"], lwd = 1, col = gray(.2), lty = 5)
lines(yearsUS, rUS[, "r0.5.r"], lwd = 2, col = gray(.2))
# Spain results
lines(yearsES, rES[, "r.m"], lwd = 1, col = gray(.2), lty = 5)
lines(yearsES, rES[, "r.f"], lwd = 1, col = gray(.2), lty = 5)
lines(yearsES, rES[, "r0.5.r"], lwd = 2, col = gray(.2))
# label rm rf
segments(1969.941, 0.011398784, 1969.526, 0.009996714)
segments(1988.910, -0.00835212, 1986.903, -0.00975419)
segments(1974.787, 0.0004870,1972.710,-0.0012198488)

text(c(1990, 1986.5, 1987, 1973, 1970,1976.5),
        c(-0.0103542581, -0.0141970362, -0.007803483, 0.0135, -0.008,0.002437726),
        c(expression(ES~sigma == 1),expression(ES~sigma == 0),expression(ES~sigma == .5),
                expression(US~sigma == 1),expression(US~sigma == 0),expression(US~sigma == .5)),
        cex = .8, pos = c(4,1,4,1),xpd=TRUE)
dev.off()

locator(2)