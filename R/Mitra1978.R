# -----------------------------------------------------------------------
# Mitra 1978: [linear, OLS u, v]
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
# BxES is 0:110, years 1975:2009
BxUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy10_65.Rdata")))
BxES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))

# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
LxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmUS.Rdata"))) / 1e5
LxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfUS.Rdata"))) / 1e5
LxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmES.Rdata"))) / 1e5
LxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfES.Rdata"))) / 1e5


yearsES <- 1975:2009
yearsUS <- 1969:2009
ages    <- 0:110

library(compiler)
Mitra1978initv0 <- cmpfun(function(r, .v0., .Bma., .Bfa., .Pma., .Pfa., .Mat., .Fat., .a.){
            # find r, given a starting .v0.
            .Hma. <- (.Bma. * .v0.) / .Mat. # eq 4
            .Hfa. <- (.Bfa. * (1 - .v0.)) / .Fat. # eq 5
            # residual from eq 17:
            (1 - sum(exp(-r * .a.)*(.Pma. * .Hma. + .Pfa. * .Hfa.), na.rm = TRUE)) ^ 2
        })

v0OLSmin <- function(.v0, .Bma, .Bfa, .Pma, .Pfa, .Mat, .Fat, .a){
    # find the v0 that minimizes the difference between the 
    # initial and stable unadjusted male and female rates
    # (optimize over logit of vt)
    .Hma <- (.Bma * .v0) / .Mat  # eq 4
    .Hfa <- (.Bfa * (1 - .v0)) / .Fat  # eq 5
    ri <- optimize(f = Mitra1978initv0, 
            interval = range(c(LotkaRCoale(.Bma / .Mat, .Pma, .a), 
                            LotkaRCoale(.Bfa / .Fat, .Pfa, .a))),
            .v0. = .v0, .Bma. = .Bma, .Bfa. = .Bfa, 
            .Pma. = .Pma, .Pfa. = .Pfa, 
            .Mat. = .Mat, .Fat. = .Fat,
            .a. = .a, tol = 1e-15)$minimum
    #
    v.lim <- sum(exp(-ri * .a) * (.Pma * .Hma), na.rm = TRUE) # eq 34
    
    # initial rates:
    mma0 <- .Bma / .Mat
    mfa0 <- .Bfa / .Fat
    # stable rates:
    mma.lim <- .Hma / v.lim # eq 38
    mfa.lim <- .Hfa / (1 - v.lim)
    
    sum((mma.lim - mma0) ^2, na.rm = TRUE) + sum((mfa.lim - mfa0) ^2, na.rm = TRUE) # eq 45
}

Mitra1978OLS <- function(Bma, Bfa, Pma, Pfa, Mat, Fat){
    a  <- 0:110 + .5
    v0 <- optimize(v0OLSmin, interval = c(0,1),
            .Bma = Bma, .Bfa = Bfa, 
            .Pma = Pma, .Pfa = Pfa, 
            .Mat = Mat, .Fat = Fat, 
            .a = a, tol = 1e-15)$minimum
    
    Hma <- (Bma * v0) / Mat # from eq 4
    Hfa <- (Bfa * (1 - v0)) / Fat # from eq 5
    
    # use Coale's iterative solution for Lotka's r, single sex
    r.f <- LotkaRCoale(Bfa / Fat, Pfa, a)
    r.m <- LotkaRCoale(Bma / Mat, Pma, a)
    
    # (re) find r given by optimized v0 
    # [it was already used, but not returned in v0OLSmin()]
    r  <- optimize(f = Mitra1978initv0, 
            interval = range(c(r.m, r.f)),
            .v0. = v0, .Bma. = Bma, .Bfa. = Bfa, 
            .Pma. = Pma, .Pfa. = Pfa, 
            .Mat. = Mat, .Fat. = Fat, .a. = a,
            tol = 1e-15)$minimum
    
    v.lim <- sum(exp(-r * a)*(Pma * Hma), na.rm = TRUE) # eq 34
    
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


