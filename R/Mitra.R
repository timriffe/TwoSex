

# from Mitra (1976)
# Mitra used 1966 US data to compare with Das Gupta. Can't reproduce his results
# but I can apply the method.
# ----------------------------------------------
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
# BxES is 0:110, years 1975:2009
BxUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy10_65.Rdata")))
BxES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf10_65.Rdata")))
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf10_65.Rdata")))

# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
LxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmUS.Rdata")))[11:66, ] / 1e5
LxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfUS.Rdata")))[11:50, ] / 1e5
LxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmES.Rdata")))[11:66, ] / 1e5
LxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfES.Rdata")))[11:50, ] / 1e5


yearsES <- 1975:2009
yearsUS <- 1969:2009
ages    <- 10:65

agem <- 10:65
agef <- 10:49
# example year to get started:

# US 1975:

Bxym <- BxymfUS[["1975"]][["Bxym"]][, 1:40]
Bxyf <- BxymfUS[["1975"]][["Bxyf"]][, 1:40]
# just take ages 10-49 for females here



Exm <- with(ExUS, Male[Year == 1975 & Age >= 10 & Age <= 65])
Eyf <- with(ExUS, Female[Year == 1975 & Age >= 10 & Age <= 49])
Lxm <- LxmUS[,"1975"]
Lyf <- LxfUS[,"1975"]
# only interested in truly reproductive ages here, not ages 10:65
Ftot <- sum(Eyf)
Mtot <- sum(Exm)
Ptot <- sum(Ftot,Mtot)

mm <- rowSums(Bxym) / Exm
mf <- colSums(Bxyf) / Eyf


Km <- mm / (Ftot / Ptot)
all(round(mm - Km * (Ftot / Ptot), digits = 12) == 0) # eq 2
#plot(10:65, Km, type = 'l')
#lines(10:65, mm)

Kf <- mf / (Mtot / Ptot)
all(round(mf - Kf * (Mtot / Ptot), digits = 12) == 0) # eq 3

ut <- Ftot / Ptot # eq 5
vt <- Mtot / Ptot # eq 7
(ut + vt == 1)    # eq 8
gm <- Lxm * Km    # eq 6
gf <- Lyf * Kf

R0m <- sum(mm * Lxm)
R0f <- sum(mf * Lxf)


eq10optim <- function(r, gf, gm, agef = 10.5:49.5, agem = 10.5:65.5){
    (1 - (1 / sum(exp(-r * agem) * gm) + 1 / sum(exp(-r * agef) * gf))) ^ 2
}

# produces residual to minimize:
# eq10optim(pars, gf, gm) 

# eq23, first approximation of two-sex r, to start optimization:
rstart <- (Ptot * (R0m / Ftot) * (R0f / Mtot) - R0m / Ftot - R0f / Mtot) / 
            (wmean(agem + .5, gm) / Mtot + wmean(agef + .5, gf) / Ftot)

  
rstar <- optimize(f = eq10optim, 
                  # unique r is bounded by male and female r, thus we specify for faster convergence
                  interval = range(c(LotkaRCoale(mm, Lxm, agem + .5),
                                     LotkaRCoale(mf, Lxf, agef + .5)
                                )), 
                  gf = gf, 
                  gm = gm, 
                  agef = agef + .5, 
                  agem = agem + .5,
                  tol = 1e-11)$minimum
 
rMitra <- compiler::cmpfun(function(Bxym, Bxyf, Exm, Eyf, Lxm, Lyf, agem, agef){
    Ftot    <- sum(Eyf)
    Mtot    <- sum(Exm)
    Ptot    <- sum(Ftot, Mtot)
    mm      <- rowSums(Bxym) / Exm
    mf      <- colSums(Bxyf) / Eyf
    Km      <- mm / (Ftot / Ptot)
    Kf      <- mf / (Mtot / Ptot)
    gm      <- Lxm * Km    # eq 6
    gf      <- Lyf * Kf
    
    optimize(f = eq10optim, 
            # unique r is bounded by male and female r, thus we specify for faster convergence
            interval = range(c(LotkaRCoale(mm, Lxm, agem),
                            LotkaRCoale(mf, Lxf, agef)
                    )), 
            gf = gf, 
            gm = gm, 
            agef = agef, 
            agem = agem,
            tol = 1e-11)$minimum
})

rMitra(Bxym,Bxyf,Exm,Eyf,Lxm,Lyf, agem + .5, agef + .5)

# calculate ES and US single sex and Mitra r estimates:

rmfUS <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .BxymfUS, .ExUS, .LxmUS, .LxfUS, agem = 10:65, agef = 10:49){
            Bxym <- .BxymfUS[[yr]][["Bxym"]][, 1:40]
            Bxyf <- .BxymfUS[[yr]][["Bxyf"]][, 1:40]
            Exm  <- with(.ExUS, Male[Year == as.integer(yr) & Age >= 10 & Age <= 65])
            Eyf  <- with(.ExUS, Female[Year == as.integer(yr) & Age >= 10 & Age <= 49])
            Lxm  <- .LxmUS[, yr]
            Lyf  <- .LxfUS[, yr]
            mm      <- rowSums(Bxym) / Exm
            mf      <- colSums(Bxyf) / Eyf
            c( r.mf = rMitra(Bxym, Bxyf, Exm, Eyf, Lxm, Lyf, agem + .5, agef + .5),
               r.m = LotkaRCoale(mm, Lxm, agem + .5),
               r.f = LotkaRCoale(mf, Lxf, agef))
        }, .BxymfUS = BxymfUS, .ExUS = ExUS, .LxmUS = LxmUS, .LxfUS = LxfUS))
# TODO: ES throws error, fix
rmfES <- do.call(rbind,lapply(as.character(yearsES), function(yr, .BxymfES, .ExES, .LxmES, .LxfES, agem = 10:65, agef = 10:49){
                    Bxym <- .BxymfES[[yr]][["Bxym"]][, 1:40]
                    Bxyf <- .BxymfES[[yr]][["Bxyf"]][, 1:40]
                    Exm  <- with(.ExES, Male[Year == as.integer(yr) & Age >= 10 & Age <= 65])
                    Eyf  <- with(.ExES, Female[Year == as.integer(yr) & Age >= 10 & Age <= 49])
                    Lxm  <- .LxmES[, yr]
                    Lyf  <- .LxfES[, yr]
                    mm   <- rowSums(Bxym) / Exm
                    mf   <- colSums(Bxyf) / Eyf
                    c( r.mf = rMitra(Bxym, Bxyf, Exm, Eyf, Lxm, Lyf, agem + .5, agef + .5),
                            r.m = LotkaRCoale(mm, Lxm, agem + .5),
                            r.f = LotkaRCoale(mf, Lxf, agef))
                }, .BxymfES = BxymfES, .ExES = ExES, .LxmES = LxmES, .LxfES = LxfES))

plot(yearsUS, rmfUS[,"r.mf"], type = 'l')
abline(h=0)

