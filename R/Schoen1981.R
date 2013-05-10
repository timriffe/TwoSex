# Schoen cites 1977 rectangular model. Will try that later. First get basic Schoen model up and running


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

LM <- compiler::cmpfun(function(x,y){
           M <- (y - x) / (log(y) - log(x))
           M[is.nan(M)] <- x[is.nan(M)]
           Mna0(Minf0(M))
        })
GM <- compiler::cmpfun(function(x,y){
            Mna0(Minf0((x * y) ^ (1/2)))
        })

# the optimizer solution only will work with a fixed SRB
# residual function
#LotkaHmin <- compiler::cmpfun(function(r, Lxm, Lxf, FH, SRB,.a = .5:110.5){
#   theta.m <- SRB / (1 + SRB)
#   theta.f <- 1 / (1 + SRB)
#  (1 - sum(outer(theta.m * exp(-r * .a) * Lxm, theta.f * exp(-r * .a) * Lxf, HM) * FH)) ^ 2
#})
#optimize(LotkaHmin, interval = c(-.02,.01), Lxm = Lxm, Lxf = Lxf, FH = FH, SRB = SRB, tol = 1e-15)


# takes all iterations: doesn't coverge fast enough: better use optimizer for fixed SRB
#SchoenIt <- compiler::cmpfun(function(Lxm, Lxf, FH, SRB, .a = .5:110.5, maxit = 2e2, tol = 1e-15){
#            theta.m <- SRB / (1 + SRB)
#            theta.f <- 1 / (1 + SRB)
#            
#            R0.guess <- sum(outer(theta.m * Lxm, theta.f * Lxf, HM) * FH)
#            T.guess <- sum(outer(theta.m * Lxm * .a, theta.f * Lxf * .a, HM) * FH) / R0.guess
#            r.i <- log(R0.guess) / T.guess
#            delta.vec<-r.vec <- vector(length=maxit)
#            for (i in 1:maxit){
#                delta.i <- (1 - sum(outer(theta.m * exp(-r.i * .a) * Lxm, theta.f * exp(-r.i * .a) * Lxf, HM) * FH)) ^ 2
#                r.i <- r.i + (delta.i / ( T.guess  - (delta.i / r.i)))
#                if (abs(delta.i) < tol){
#                    break
#                }  
#                r.vec[i] <- r.i
#                delta.vec[i] <-  delta.i
#            }
#            cat(i,"\n")
#           r.i
#        })

# takes all iterations: doesn't coverge fast enough: better use optimizer for fixed SRB
#rUSshoen <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxy, .Lxm, .Lxf, .Ex){
#                    
#                    Exm  <- with(.Ex, Male[Year == as.integer(yr)])
#                    Exf  <- with(.Ex, Female[Year == as.integer(yr)])
#                    # harmonic mean of exposures
#                    Hxy  <- outer(Exm, Exf, HM)
#                    # harmonic rates (divide by two since working with both sexes)
#                    Bxyf <- .Bxy[[yr]][["Bxyf"]]
#                    Bxym <- .Bxy[[yr]][["Bxym"]]
#                    SRB  <- sum(Bxym) / sum(Bxyf)
#                    FxyH <- Mna0(Minf0((Bxyf + Bxym) / Hxy))
#
#                    Lxf <- .Lxf[, yr]
#                    Lxm <- .Lxm[, yr]
#                    
#                    Fxm <- Mna0(Minf0(rowSums(.Bxy[[yr]][["Bxym"]]) / Exm))
#                    Fxf <- Mna0(Minf0(colSums(.Bxy[[yr]][["Bxyf"]]) / Exf))
#                    
#                    LH <- SchoenIt(Lxm = Lxm, Lxf = Lxf, FH = FxyH, SRB = SRB)
#                    c(r.f = LotkaRCoale(Fxf, Lxf, x = .5:110.5),
#                            r.m = LotkaRCoale(Fxm, Lxm, x = .5:110.5),
#                            r.H = LH[1])
#                }, .Bxy = BxymfUS, .Lxm = LxmUS, .Lxf = LxfUS, .Ex = ExUS))
#rownames(rUSshoen) <- yearsUS

# using strategy of Coale
LotkaMCoale <- compiler::cmpfun(function(FHf, FHm, Lxm, Lxf, M = HM, x = .5:110.5, maxit = 2e2, tol = 1e-15){
    # from Coale, Ansley J. (1957) A New Method for Calculating Lotka's r- the Intrinsic Rate of Growth in a Stable Population.
    # Population Studies, Vol. 11 no. 1, pp 92-94
    
    # SRB that would result from r = 0
    SRBi        <- sum(outer(Lxm, Lxf, M) * FHm) / sum(outer(Lxm, Lxf, M) * FHf)
    # image(outer(Lxm * SRBi / (SRBi+1), Lxf * 1 / (SRBi+1), M) )
    R0          <- sum(outer(Lxm * SRBi / (SRBi+1), Lxf * 1 / (SRBi+1), M) * (FHf + FHm))
    # first assuming a mean generation time of 35
    r2 <- log(R0) / 35
    for (i in 1:maxit){ # 10 is more than enough!
        r1      <- r2
        p.m     <- SRBi / (SRBi+1)
        p.f     <- 1 / (SRBi+1)
        deltai  <- 1 - sum(outer(exp(-r1 * x) * Lxm * p.m, exp(-r1 * x) * Lxf * p.f, M) * (FHf + FHm)) 
        # the mean generation time self-corrects 
        # according to the error produced by the Lotka equation
        r2      <- r1 - (deltai / (35 - (deltai / r1)))
      
        SRBi    <-  sum(outer(exp(-r2 * x) * Lxm * p.m, exp(-r2 * x) * Lxf * p.f, M) * FHm) / 
                    sum(outer(exp(-r2 * x) * Lxm * p.m, exp(-r2 * x) * Lxf * p.f, M) * FHf)
        if (abs(r2 - r1) <= tol | zapsmall(abs(deltai)) <= tol){
            break
        }
    }
    return(c(r2, SRBi))
})
yr <- "1975"
rUSHM <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxy, .Lxm, .Lxf, .Ex, .M){
            
            Exm <- with(.Ex, Male[Year == as.integer(yr)])
            Exf <- with(.Ex, Female[Year == as.integer(yr)])
            # harmonic mean of exposures
            Hxy <- outer(Exm, Exf, .M)
            # harmonic rates (divide by two since working with both sexes)
            FxyHf <- Mna0(Minf0((.Bxy[[yr]][["Bxyf"]]) / Hxy))
            FxyHm <- Mna0(Minf0((.Bxy[[yr]][["Bxym"]]) / Hxy))
            
           
            # single sex rates
            Fxm <- Mna0(Minf0(rowSums(.Bxy[[yr]][["Bxym"]]) / Exm))
            Fxf <- Mna0(Minf0(colSums(.Bxy[[yr]][["Bxyf"]]) / Exf))
            
            Lxf <- .Lxf[, yr]
            Lxm <- .Lxm[, yr]
            
            LH <- LotkaMCoale(FHf = FxyHf, FHm = FxyHm, Lxm =Lxm, Lxf = Lxf, M = .M)
            c(r.f = LotkaRCoale(Fxf, Lxf, x = .5:110.5),
              r.m = LotkaRCoale(Fxm, Lxm, x = .5:110.5),
              r.H = LH[1], SRBH = LH[2])
        }, .Bxy = BxymfUS, .Lxm = LxmUS, .Lxf = LxfUS, .Ex = ExUS, .M = HM))
rownames(rUSHM) <- yearsUS
rESHM <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxy, .Lxm, .Lxf, .Ex, .M){
                    
                    Exm <- with(.Ex, Male[Year == as.integer(yr)])
                    Exf <- with(.Ex, Female[Year == as.integer(yr)])
                    # harmonic mean of exposures
                    Hxy <- outer(Exm, Exf, .M)
                    # harmonic rates (divide by two since working with both sexes)
                    FxyHf <- Mna0(Minf0((.Bxy[[yr]][["Bxyf"]]) / Hxy))
                    FxyHm <- Mna0(Minf0((.Bxy[[yr]][["Bxym"]]) / Hxy))
                    
                    
                    # single sex rates
                    Fxm <- Mna0(Minf0(rowSums(.Bxy[[yr]][["Bxym"]]) / Exm))
                    Fxf <- Mna0(Minf0(colSums(.Bxy[[yr]][["Bxyf"]]) / Exf))
                    
                    Lxf <- .Lxf[, yr]
                    Lxm <- .Lxm[, yr]
                    
                    LH <- LotkaMCoale(FHf = FxyHf, FHm = FxyHm, Lxm =Lxm, Lxf = Lxf, M = .M)
                    c(r.f = LotkaRCoale(Fxf, Lxf, x = .5:110.5),
                            r.m = LotkaRCoale(Fxm, Lxm, x = .5:110.5),
                            r.H = LH[1], SRBH = LH[2])
                }, .Bxy = BxymfES, .Lxm = LxmES, .Lxf = LxfES, .Ex = ExES, .M = HM))
rownames(rESHM) <- yearsES

rESLM <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxy, .Lxm, .Lxf, .Ex, .M){
                    
                    Exm <- with(.Ex, Male[Year == as.integer(yr)])
                    Exf <- with(.Ex, Female[Year == as.integer(yr)])
                    # harmonic mean of exposures
                    Hxy <- outer(Exm, Exf, .M)
                    # harmonic rates (divide by two since working with both sexes)
                    FxyHf <- Mna0(Minf0((.Bxy[[yr]][["Bxyf"]]) / Hxy))
                    FxyHm <- Mna0(Minf0((.Bxy[[yr]][["Bxym"]]) / Hxy))
                    
                    
                    # single sex rates
                    Fxm <- rowSums(Mna0(Minf0(.Bxy[[yr]][["Bxym"]] / Exm)))
                    Fxf <- colSums(Mna0(Minf0(t(t(.Bxy[[yr]][["Bxyf"]]) / Exf))))
                    
                    Lxf <- .Lxf[, yr]
                    Lxm <- .Lxm[, yr]
                    
                    LH <- LotkaMCoale(FHf = FxyHf, FHm = FxyHm, Lxm =Lxm, Lxf = Lxf, M = .M)
                    c(r.f = LotkaRCoale(Fxf, Lxf, x = .5:110.5),
                            r.m = LotkaRCoale(Fxm, Lxm, x = .5:110.5),
                            r.H = LH[1], SRBH = LH[2])
                }, .Bxy = BxymfES, .Lxm = LxmES, .Lxf = LxfES, .Ex = ExES, .M = LM))
rownames(rESLM) <- yearsES

rUSLM <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxy, .Lxm, .Lxf, .Ex, .M){
                    
                    Exm <- with(.Ex, Male[Year == as.integer(yr)])
                    Exf <- with(.Ex, Female[Year == as.integer(yr)])
                    # harmonic mean of exposures
                    Hxy <- outer(Exm, Exf, .M)
                    # harmonic rates (divide by two since working with both sexes)
                    FxyHf <- Mna0(Minf0((.Bxy[[yr]][["Bxyf"]]) / Hxy))
                    FxyHm <- Mna0(Minf0((.Bxy[[yr]][["Bxym"]]) / Hxy))
                    
                    
                    # single sex rates
                    Fxm <- Mna0(Minf0(rowSums(.Bxy[[yr]][["Bxym"]]) / Exm))
                    Fxf <- Mna0(Minf0(colSums(.Bxy[[yr]][["Bxyf"]]) / Exf))
                    
                    Lxf <- .Lxf[, yr]
                    Lxm <- .Lxm[, yr]
                    
                    LH <- LotkaMCoale(FHf = FxyHf, FHm = FxyHm, Lxm =Lxm, Lxf = Lxf, M = .M)
                    c(r.f = LotkaRCoale(Fxf, Lxf, x = .5:110.5),
                            r.m = LotkaRCoale(Fxm, Lxm, x = .5:110.5),
                            r.H = LH[1], SRBH = LH[2])
                }, .Bxy = BxymfUS, .Lxm = LxmUS, .Lxf = LxfUS, .Ex = ExUS, .M = LM))
rownames(rUSLM) <- yearsUS

rESGM <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxy, .Lxm, .Lxf, .Ex, .M){
                    
                    Exm <- with(.Ex, Male[Year == as.integer(yr)])
                    Exf <- with(.Ex, Female[Year == as.integer(yr)])
                    # harmonic mean of exposures
                    Hxy <- outer(Exm, Exf, .M)
                    # harmonic rates (divide by two since working with both sexes)
                    FxyHf <- Mna0(Minf0((.Bxy[[yr]][["Bxyf"]]) / Hxy))
                    FxyHm <- Mna0(Minf0((.Bxy[[yr]][["Bxym"]]) / Hxy))
                    
                    
                    # single sex rates
                    Fxm <- rowSums(Mna0(Minf0(.Bxy[[yr]][["Bxym"]] / Exm)))
                    Fxf <- colSums(Mna0(Minf0(t(t(.Bxy[[yr]][["Bxyf"]]) / Exf))))
                    
                    Lxf <- .Lxf[, yr]
                    Lxm <- .Lxm[, yr]
                    
                    LH <- LotkaMCoale(FHf = FxyHf, FHm = FxyHm, Lxm =Lxm, Lxf = Lxf, M = .M)
                    c(r.f = LotkaRCoale(Fxf, Lxf, x = .5:110.5),
                            r.m = LotkaRCoale(Fxm, Lxm, x = .5:110.5),
                            r.H = LH[1], SRBH = LH[2])
                }, .Bxy = BxymfES, .Lxm = LxmES, .Lxf = LxfES, .Ex = ExES, .M = GM))
rownames(rESGM) <- yearsES

rUSGM <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxy, .Lxm, .Lxf, .Ex, .M){
                    
                    Exm <- with(.Ex, Male[Year == as.integer(yr)])
                    Exf <- with(.Ex, Female[Year == as.integer(yr)])
                    # harmonic mean of exposures
                    Hxy <- outer(Exm, Exf, .M)
                    # harmonic rates (divide by two since working with both sexes)
                    FxyHf <- Mna0(Minf0((.Bxy[[yr]][["Bxyf"]]) / Hxy))
                    FxyHm <- Mna0(Minf0((.Bxy[[yr]][["Bxym"]]) / Hxy))
                    
                    
                    # single sex rates
                    Fxm <- Mna0(Minf0(rowSums(.Bxy[[yr]][["Bxym"]]) / Exm))
                    Fxf <- Mna0(Minf0(colSums(.Bxy[[yr]][["Bxyf"]]) / Exf))
                    
                    Lxf <- .Lxf[, yr]
                    Lxm <- .Lxm[, yr]
                    
                    LH <- LotkaMCoale(FHf = FxyHf, FHm = FxyHm, Lxm =Lxm, Lxf = Lxf, M = .M)
                    c(r.f = LotkaRCoale(Fxf, Lxf, x = .5:110.5),
                            r.m = LotkaRCoale(Fxm, Lxm, x = .5:110.5),
                            r.H = LH[1], SRBH = LH[2])
                }, .Bxy = BxymfUS, .Lxm = LxmUS, .Lxf = LxfUS, .Ex = ExUS, .M = GM))
rownames(rUSGM) <- yearsUS

rESMin <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxy, .Lxm, .Lxf, .Ex, .M){
                    
                    Exm <- with(.Ex, Male[Year == as.integer(yr)])
                    Exf <- with(.Ex, Female[Year == as.integer(yr)])
                    # harmonic mean of exposures
                    Hxy <- outer(Exm, Exf, .M)
                    # harmonic rates (divide by two since working with both sexes)
                    FxyHf <- Mna0(Minf0((.Bxy[[yr]][["Bxyf"]]) / Hxy))
                    FxyHm <- Mna0(Minf0((.Bxy[[yr]][["Bxym"]]) / Hxy))
                    
                    
                    # single sex rates
                    Fxm <- rowSums(Mna0(Minf0(.Bxy[[yr]][["Bxym"]] / Exm)))
                    Fxf <- colSums(Mna0(Minf0(t(t(.Bxy[[yr]][["Bxyf"]]) / Exf))))
                    
                    Lxf <- .Lxf[, yr]
                    Lxm <- .Lxm[, yr]
                    
                    LH <- LotkaMCoale(FHf = FxyHf, FHm = FxyHm, Lxm =Lxm, Lxf = Lxf, M = .M)
                    c(r.f = LotkaRCoale(Fxf, Lxf, x = .5:110.5),
                            r.m = LotkaRCoale(Fxm, Lxm, x = .5:110.5),
                            r.H = LH[1], SRBH = LH[2])
                }, .Bxy = BxymfES, .Lxm = LxmES, .Lxf = LxfES, .Ex = ExES, .M = pmin))
rownames(rESMin) <- yearsES

rUSMin <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxy, .Lxm, .Lxf, .Ex, .M){
                    
                    Exm <- with(.Ex, Male[Year == as.integer(yr)])
                    Exf <- with(.Ex, Female[Year == as.integer(yr)])
                    # harmonic mean of exposures
                    Hxy <- outer(Exm, Exf, .M)
                    # harmonic rates (divide by two since working with both sexes)
                    FxyHf <- Mna0(Minf0((.Bxy[[yr]][["Bxyf"]]) / Hxy))
                    FxyHm <- Mna0(Minf0((.Bxy[[yr]][["Bxym"]]) / Hxy))
                    
                    
                    # single sex rates
                    Fxm <- Mna0(Minf0(rowSums(.Bxy[[yr]][["Bxym"]]) / Exm))
                    Fxf <- Mna0(Minf0(colSums(.Bxy[[yr]][["Bxyf"]]) / Exf))
                    
                    Lxf <- .Lxf[, yr]
                    Lxm <- .Lxm[, yr]
                    
                    LH <- LotkaMCoale(FHf = FxyHf, FHm = FxyHm, Lxm =Lxm, Lxf = Lxf, M = .M)
                    c(r.f = LotkaRCoale(Fxf, Lxf, x = .5:110.5),
                            r.m = LotkaRCoale(Fxm, Lxm, x = .5:110.5),
                            r.H = LH[1], SRBH = LH[2])
                }, .Bxy = BxymfUS, .Lxm = LxmUS, .Lxf = LxfUS, .Ex = ExUS, .M = pmin))
rownames(rUSMin) <- yearsUS

#save(rUSHM, file = "/home/triffe/git/DISS/Data/results/agerSRB/rHMUS")
#save(rESHM, file = "/home/triffe/git/DISS/Data/results/agerSRB/rHMES")
#save(rUSGM, file = "/home/triffe/git/DISS/Data/results/agerSRB/rGMUS")
#save(rESGM, file = "/home/triffe/git/DISS/Data/results/agerSRB/rGMES")
#save(rUSLM, file = "/home/triffe/git/DISS/Data/results/agerSRB/rLMUS")
#save(rESLM, file = "/home/triffe/git/DISS/Data/results/agerSRB/rLMES")
#save(rUSMin, file = "/home/triffe/git/DISS/Data/results/agerSRB/rMinUS")
#save(rESMin, file = "/home/triffe/git/DISS/Data/results/agerSRB/rMinES")



# plot Harmonic mean and min results
pdf("/home/triffe/git/DISS/latex/Figures/HMager.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rUSHM[,"r.f"], type = 'n', ylim = c(-.02, .015), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968,-.02,2010,.015,col = gray(.95), border=NA),
                abline(h = seq(-.02, .015, by = .005), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.02, .015, by = .005),seq(-.02, .02, by = .005), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.02, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.0225, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1965.5,.018, "r", cex = 1, xpd = TRUE)))
# rm-rf regions
polygon(c(yearsUS,rev(yearsUS)),c(rUSHM[, "r.m"],rev(rUSHM[, "r.f"])), border = NA, col = "#55555550")
polygon(c(yearsES,rev(yearsES)),c(rESHM[, "r.m"],rev(rESHM[, "r.f"])), border = NA, col = "#55555550")
# US results
lines(yearsUS,rUSHM[, "r.H"], lwd = 2, col = gray(.2))
lines(yearsUS, rUSHM[, "r.m"], lwd = 1, col = gray(.2), lty = 5)
lines(yearsUS, rUSHM[, "r.f"], lwd = 1, col = gray(.2), lty = 5)
lines(yearsUS,rUSMin[, "r.H"], lwd = 1, col = "red")

# Spain results
lines(yearsES,rESHM[, "r.H"], lwd = 2, col = gray(.2))
lines(yearsES, rESHM[, "r.m"], lwd = 1, col = gray(.2), lty = 5)
lines(yearsES, rESHM[, "r.f"], lwd = 1, col = gray(.2), lty = 5)
lines(yearsES,rESMin[, "r.H"], lwd = 1, col = "red")

# label rm rf
text(c(1990, 1986.5, 1973, 1971),
        c(-0.0103542581, -0.0141970362,  0.0003650703, -0.0040170451),
        c(expression(r^m~ES),expression(r^f~ES),expression(r^m~US),expression(r^f~US)),
        cex = .8, pos = c(4,1,4,1))
legend(1995, .015, lty = 1, lwd = c(2,1), col = c(gray(.2),"red"), bty = "n", 
        legend = c(expression(r^H),expression(r^min)), xpd = TRUE)
dev.off()

# -------------------------------------


StableAgeH <- function(r, SRB, FHf, FHm, Lxm, Lxf, x = .5:110.5){
    p.m <- SRB / (1 + SRB)
    p.f <- 1 / (1 + SRB)
    b   <- 1 / sum(exp(-r * x) * (Lxm * p.m + Lxf * p.f))
    cam <- b * p.m * exp(-r * x) * Lxm
    caf <- b * p.f * exp(-r * x) * Lxf
    cbind(ca.m = cam, ca.f = caf)
}


plot(0:110, p.m * exp(-r * x) * Lxm, type = 'l')
lines(0:110, p.f * exp(-r * x) * Lxf)
