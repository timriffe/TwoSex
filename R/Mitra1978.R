# -----------------------------------------------------------------------
# Mitra 1978: [linear, OLS u, v]
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns

yearsES <- 1975:2009
yearsUS <- 1969:2009
BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata")))
names(BxymfES) <- yearsES
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))

# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
LxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmUS.Rdata"))) / 1e5
LxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfUS.Rdata"))) / 1e5
LxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmES.Rdata"))) / 1e5
LxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfES.Rdata"))) / 1e5

lxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_lx/lxmUS.Rdata"))) / 1e5
lxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_lx/lxfUS.Rdata"))) / 1e5
lxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_lx/lxmES.Rdata"))) / 1e5
lxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_lx/lxfES.Rdata"))) / 1e5


ages    <- 0:110

library(compiler)
Mitra1978initv0 <- cmpfun(function(r, .v0., .Bma., .Bfa., .Pma., .Pfa., .Mat., .Fat., .a.){
            # find r, given a starting .v0.
            .Hma. <-  Minf0(Mna0((.Bma. * .v0.) / .Mat.)) # eq 4
            .Hfa. <-  Minf0(Mna0((.Bfa. * (1 - .v0.)) / .Fat.)) # eq 5
            # residual from eq 17:
            (1 - sum(exp(-r * .a.)*(.Pma. * .Hma. + .Pfa. * .Hfa.), na.rm = TRUE)) ^ 2
        })

v0OLSmin <- function(.v0, .Bma, .Bfa, .Pma, .Pfa, .Mat, .Fat, .a, .SRB0){
    # find the v0 that minimizes the difference between the 
    # initial and stable unadjusted male and female rates
    # (optimize over logit of vt)
    .Hma <-  Minf0(Mna0((.Bma * .v0) / .Mat))  # eq 4
    .Hfa <-  Minf0(Mna0((.Bfa * (1 - .v0)) / .Fat))  # eq 5
    ri <- optimize(f = Mitra1978initv0, 
            interval = range(c(LotkaRCoale( Minf0(Mna0(.Bma / .Mat)), .Pma, .a), 
                            LotkaRCoale( Minf0(Mna0(.Bfa / .Fat)), .Pfa, .a))),
            .v0. = .v0, .Bma. = .Bma, .Bfa. = .Bfa, 
            .Pma. = .Pma, .Pfa. = .Pfa, 
            .Mat. = .Mat, .Fat. = .Fat,
            .a. = .a, tol = 1e-15)$minimum
    # eq 33
    #v.lim <- sum(exp(-r * a) * (.Pma * .Hma), na.rm = TRUE) / 
    #        (sum(exp(-r * a) * (.Pma * .Hma), na.rm = TRUE) + .SRB0 * sum(exp(-ri * .a) * (.Pfa * .Hfa), na.rm = TRUE) )
    #
    v.lim <- sum(exp(-ri * .a) * (.Pma * .Hma), na.rm = TRUE) # eq 34
    
    # initial rates:
    mma0 <-  Minf0(Mna0(.Bma / .Mat))
    mfa0 <-  Minf0(Mna0(.Bfa / .Fat))
    # stable rates:
    mma.lim <- .Hma / v.lim # eq 38
    mfa.lim <- .Hfa / (1 - v.lim)
    
    sum((mma.lim - mma0) ^2, na.rm = TRUE) + sum((mfa.lim - mfa0) ^2, na.rm = TRUE) # eq 45
}

Mitra1978OLS <- function(Bma, Bfa, Pma, Pfa, Mat, Fat){
    SRB0 <- sum(Bma, na.rm = TRUE) / sum(Bfa, na.rm = TRUE)
    a    <- 0:110 + .5
    v0 <- optimize(v0OLSmin, interval = c(0,1),
            .Bma = Bma, .Bfa = Bfa, 
            .Pma = Pma, .Pfa = Pfa, 
            .Mat = Mat, .Fat = Fat, .SRB0 = SRB0,
            .a = a, tol = 1e-15)$minimum
    
    Hma <-  Minf0(Mna0((Bma * v0) / Mat)) # from eq 4
    Hfa <-  Minf0(Mna0((Bfa * (1 - v0)) / Fat)) # from eq 5
    
    # use Coale's iterative solution for Lotka's r, single sex
    r.f <- LotkaRCoale( Minf0(Mna0(Bfa / Fat)), Pfa, a)
    r.m <- LotkaRCoale( Minf0(Mna0(Bma / Mat)), Pma, a)
    
    # (re) find r given by optimized v0 
    # [it was already used, but not returned in v0OLSmin()]
    r  <- optimize(f = Mitra1978initv0, 
            interval = range(c(r.m, r.f)),
            .v0. = v0, .Bma. = Bma, .Bfa. = Bfa, 
            .Pma. = Pma, .Pfa. = Pfa, 
            .Mat. = Mat, .Fat. = Fat, .a. = a,
            tol = 1e-15)$minimum
    # eq 33 (same as 34, simpler)
    #v.lim <- sum(exp(-r * a) * (Pma * Hma), na.rm = TRUE) / 
    #       (sum(exp(-r * a) * (Pma * Hma), na.rm = TRUE) + SRB0 * sum(exp(-r * a) * (Pfa * Hfa), na.rm = TRUE) )
    v.lim <- sum(exp(-r * a) * (Pma * Hma), na.rm = TRUE) # eq 34
    
    # return estimates
    c(v = v.lim, v0 = v0, r = r, r.m = r.m, r.f = r.f)
}

MitraOLSUSresults <- matrix(ncol = 5, nrow = length(yearsUS), dimnames = list(yearsUS,c("v","v0","r","r.m","r.f")))
for (i in 1:length(yearsUS)){
    yr <- as.character(yearsUS[i])
    MitraOLSUSresults[i, ] <- 
            Mitra1978OLS(   Bma = rowSums(BxymfUS[[yr]][["Bxym"]]), 
                    Bfa = colSums(BxymfUS[[yr]][["Bxyf"]]), 
                    Pma = LxmUS[, yr], 
                    Pfa = LxmUS[, yr], 
                    Mat = with(ExUS, Male[Year == yearsUS[i]]), 
                    Fat = with(ExUS, Female[Year == yearsUS[i]]))
}


# * makes no difference whether we optimize over v0 or expit(v0).
# no sensitivity issues, it would appear
plot(yearsUS, MitraOLSUSresults[, "v"], type = 'l', col = "blue", ylim = c(.4, .6))
lines(yearsUS, 1 - MitraOLSUSresults[, "v"], col = "pink")

plot(yearsUS, logit(MitraOLSUSresults[, "v0"]), type = 'l', col = "blue", ylim = c(-.1, .1))
lines(yearsUS, logit(1 - MitraOLSUSresults[, "v0"]), col = "pink")

plot(yearsUS, logit(MitraOLSUSresults[, "v0"]), type = 'l', col = "magenta", ylim = c(-.05,.2))
lines(yearsUS, logit(MitraOLSUSresults[, "v"]), col = "blue")

plot(yearsUS, MitraOLSUSresults[, "r.m"], type = 'l', col = "blue", ylim = c(-.01, .01))
lines(yearsUS, MitraOLSUSresults[, "r.f"], col = "pink")
lines(yearsUS,  MitraOLSUSresults[, "r"], col = "green")


# change strategy to slow optimization over a trajectory?


plot(v.vec, type = 'l')
plot(r.vec[2:length(r.vec)], type = 'l')


# now functionalize the above! 

v.r.stable <- compiler::cmpfun(function(.v0, .Lxm, .Lxf, .lxm, .lxf, .Bma, .Bfa, .Mat, 
                .Fat, tol = 1e-11, maxit = 1e3, 
                SRB  = sum(.Bma, na.rm = TRUE) / sum(.Bfa, na.rm = TRUE)
){
    N           <- length(.Lxm)
    Sxm         <- Minf0(Mna0(.Lxm[2:N] / .Lxm[1:(N - 1)]))
    Sxf         <- Minf0(Mna0(.Lxf[2:N] / .Lxf[1:(N - 1)]))
# Caswell eq 2.37 - infant mort to account for birth not making it to t+1
    infm        <- .lxm[1] * sqrt(.lxm[2])
    inff        <- .lxf[1] * sqrt(.lxf[2])
    # fertility-free Leslie matrices
    Lm          <- Leslie(rep(0, length(Sxm)), Mna0(Sxm))
    Lf          <- Leslie(rep(0, length(Sxf)),Mna0(Sxf))
    
    v.t         <- .v0 # v0 starting value - change
    Mat.1       <- .Mat
    Fat.1       <- .Fat
    Hma         <- Minf0(Mna0((.Bma * v.t) / .Mat))  # eq 4
    Hfa         <- Minf0(Mna0((.Bfa * (1 - v.t)) / .Fat))
    theta.m     <- Hma * .Lxm
    theta.f     <- Hfa * .Lxf
    
    v.vec       <- rep(NA, maxit)
    r.vec       <- rep(NA, maxit)
    p.vec       <- rep(NA, maxit)
    
    p.vec[1]    <- sum(Mat.1, na.rm = TRUE) + sum(Fat.1, na.rm = TRUE)
    v.vec[1]    <- v.t
    
    stable      <- FALSE
    i           <- 1
    while (i < maxit & !stable){
        i           <- i + 1
        Mat.2       <- c(Lm %*% Mat.1)
        Fat.2       <- c(Lf %*% Fat.1)
        
        #eq 8
        vR          <- (sum(Hma * Mat.2, na.rm = TRUE) / sum(Hfa * Fat.2, na.rm = TRUE)) / SRB
        v.t         <- vR / (1 + vR)
        
        # eq 9 , identical
        #v.i <- sum(Hma * Mat.i, na.rm = TRUE) / 
        #        (sum(Hma * Mat.i, na.rm = TRUE) + SRB * sum(Hfa * Fat.i, na.rm = TRUE)) 
        
        Bmi         <- sum(Hma * Mat.2, na.rm = TRUE) / v.t
        Bfi         <- sum(Hfa * Fat.2, na.rm = TRUE) / (1 - v.t)
        # Bmi / Bfi == SRB #TRUE
        # discount for infant mort (i.e. births not surviving to t+1
        Mat.2[1]    <- Bmi * infm
        Fat.2[1]    <- Bfi * inff
        
        # save parameters to return
        v.vec[i]    <- v.t 
        p.vec[i]    <- sum(Mat.2, na.rm = TRUE) + sum(Fat.2, na.rm = TRUE)
        r.vec[i]    <- log(p.vec[i] / p.vec[i - 1])
        stable      <- ifelse(i > 10, {abs(r.vec[i] - r.vec[i-1]) < tol}, FALSE)

        # reassign
        Mat.1       <- Mat.2
        Fat.1       <- Fat.2
    }
    
    if (i >= maxit){
        stop("maxit reached- increase maxit")
    }
    # really should just minimize variation in vt?
    cbind(v.vec, r.vec)
})

v.min <- function(.v0, .Lxm, .Lxf, .lxm, .lxf, .Bma, .Bfa, .Mat, .Fat, tol = 1e-11, 
        maxit = 1e3, what.min = "all", SRB){
    # minimize variance in rate trajectories .v0 <- .5
    v.r.traj <- v.r.stable(.v0 = .v0, .Lxm = .Lxm, .Lxf = .Lxf, 
            .lxm = .lxm, .lxf = .lxf, .Bma = .Bma, .Bfa = .Bfa, 
            .Mat = .Mat, .Fat = .Fat, tol = tol, maxit = maxit, SRB = SRB)[,"v.vec"]
    v.r.traj <- v.r.traj[!is.na(v.r.traj)]
    mm0 <- Minf0(Mna0(.Bma / .Mat))
    mf0 <- Minf0(Mna0(.Bfa / .Fat))
    
    .Hma         <- Minf0(Mna0((.Bma * .v0) / .Mat))  # eq 4
    .Hfa         <- Minf0(Mna0((.Bfa * (1 - .v0)) / .Fat))
    
    if (what.min == "all"){
        return(
          mean(colSums((outer(.Hma,  v.r.traj, "/") - mm0) ^ 2)) + 
                  mean(colSums((outer(.Hfa,  (1 - v.r.traj), "/") - mf0) ^ 2)))
    } else {
        return(
          sum((mm0 - .Hma / v.r.traj[length(v.r.traj)]) ^2) + 
                  sum((mf0 - .Hfa / (1 - v.r.traj[length(v.r.traj)])) ^2))
    }
    
}

v.r.OLS <- function(Lxm, Lxf, lxm, lxf, Bma, Bfa, Mat, Fat, 
        tol = 1e-11, maxit = 1e3, what.min = "firstfinal"){
    
    a <- 0:110 + .5
    v0 <- optimize(v.min, interval = c(0.2,.8), 
            .Lxm = Lxm, .Lxf = Lxf, 
            .lxm = lxm, .lxf = lxf, 
            .Bma = Bma, .Bfa = Bfa, 
            .Mat = Mat, .Fat = Fat, 
            maxit = maxit,  tol = tol,
            what.min = what.min, 
            SRB = sum(Bma, na.rm = TRUE) / sum(Bfa, na.rm = TRUE))$minimum
    trajectory <- v.r.stable(.v0 = v0, 
            .Lxm = Lxm, .Lxf = Lxf, 
            .lxm = lxm, .lxf = lxf, 
            .Bma = Bma, .Bfa = Bfa, 
            .Mat = Mat, .Fat = Fat, 
            tol = tol, maxit = 1e3,
            SRB = sum(Bma, na.rm = TRUE) / sum(Bfa, na.rm = TRUE))
    trajectory <- trajectory[!is.na(trajectory[, 1]), ]
    N          <- nrow(trajectory)
    invisible(list(v0 = v0, v.lim = trajectory[N, "v.vec"], 
                    r = trajectory[N, "r.vec"], 
                    r.f = LotkaRCoale(Minf0(Mna0(Bfa / Fat)), Lxf, a),
                    r.m = LotkaRCoale(Minf0(Mna0(Bma / Mat)), Lxm, a),
                    trajectory = trajectory))
}
# keep for testing:
#.Lxm <- Lxm     <- LxmUS[, yr]
#.Lxf <- Lxf     <- LxfUS[, yr]
#.lxm <- lxm     <- lxmUS[, yr]
#.lxf <- lxf     <- lxfUS[, yr]
#.Bma <- Bma     <- rowSums(BxymfUS[[yr]][["Bxym"]])
#.Bfa <- Bfa     <- colSums(BxymfUS[[yr]][["Bxyf"]])
#.Mat <- Mat     <- with(ExUS, Male[Year == yearsUS[i]])
#.Fat <- Fat     <- with(ExUS, Female[Year == yearsUS[i]])
#tol     <- 1e-11
#maxit   <- 1e5
#what.min <- "firstlast"



resultsUS <- list()
for (i in 1:length(yearsUS)){
    yr <- as.character(yearsUS[i])
    resultsUS[[yr]] <- v.r.OLS(Lxm = LxmUS[, yr], Lxf = LxfUS[, yr], 
            lxm = lxmUS[, yr], lxf = lxfUS[, yr], 
            Bma = rowSums(BxymfUS[[yr]][["Bxym"]]), 
            Bfa = colSums(BxymfUS[[yr]][["Bxyf"]]), 
            Mat = with(ExUS, Male[Year == yearsUS[i]]), 
            Fat = with(ExUS, Female[Year == yearsUS[i]]), 
            tol = 1e-08, maxit = 1e5, what.min = "firstlast")
    
}

resultsES <- list()
for (i in 1:length(yearsES)){
    yr <- as.character(yearsES[i])
    resultsES[[yr]] <- v.r.OLS(Lxm =  Mna0(LxmES[, yr]), Lxf =  Mna0(LxfES[, yr]), 
            lxm =  Mna0(lxmES[, yr]), lxf =  Mna0(lxfES[, yr]), 
            Bma = Mna0(rowSums(BxymfES[[yr]][["Bxym"]])), 
            Bfa = Mna0(colSums(BxymfES[[yr]][["Bxyf"]])), 
            Mat = Mna0(with(ExES, Male[Year == yearsES[i]])), 
            Fat = Mna0(with(ExES, Female[Year == yearsES[i]])), 
            tol = 1e-11, maxit = 1e4, what.min = "firstlast")
    
}

names(resultsUS[[1]])
r.stableUS <- unlist(lapply(resultsUS, "[[", "r"))
plot(yearsUS, MitraOLSUSresults[, "r"], type = 'l')
lines(yearsUS, r.stableUS, col = "blue")



v0.stableUS <- unlist(lapply(resultsUS, "[[", "v0"))
v.stableUS <- unlist(lapply(resultsUS, "[[", "v.lim"))
rf.stableUS <- unlist(lapply(resultsUS, "[[","r.f"))
rm.stableUS <- unlist(lapply(resultsUS, "[[","r.m"))


plot(yearsUS, v0.stableUS, type = 'l', ylim = range(c(v0.stableUS, v.stableUS)))
lines(yearsUS, v.stableUS, type = 'l', col = "blue")


plot(yearsUS, r.stableUS, type = 'l', ylim = range(c(rf.stableUS, rm.stableUS)))
lines(yearsUS, rf.stableUS, col = "red")
lines(yearsUS, rm.stableUS, col = "blue")

plot(unlist(lapply(resultsUS, function(x){nrow(x$trajectory)})))


v.stableES <- unlist(lapply(resultsES, "[[", "v.lim"))
v0.stableES <- unlist(lapply(resultsES, "[[", "v0"))
plot(yearsES,v.stableES, ylim = range(c(v.stableES, v0.stableES)), type = 'l', col = "red")
lines(yearsES,v0.stableES, col = "blue")
abline(h = .5)
legend("topright", lty = 1, col = c("red","blue"))

plot(yearsES, unlist(lapply(resultsES, "[[", 1)), type = 'l')
r.stableES <- unlist(lapply(resultsES, "[[", 3))
rf.stableES <- unlist(lapply(resultsES, "[[",4))
rm.stableES <- unlist(lapply(resultsES, "[[",5))
plot(yearsES, r.stableES, type = 'l', ylim = range(c(rf.stableES, rm.stableES)))
lines(yearsES, rf.stableES, col = "red")
lines(yearsES, rm.stableES, col = "blue")


plot(unlist(lapply(resultsES, function(x){nrow(x$trajectory)})))


for (i in 1:length(yearsES)){
    yr <- as.character(yearsES[i])
    plot(resultsES[[yr]]$trajectory[,"r.vec"][1:300],type='l',main=yr, ylim = c(-.02,.02))
    Sys.sleep(.6)
}


