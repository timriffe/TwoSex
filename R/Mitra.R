

# from Mitra (1976) [non-linear]
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
# ------------------------------------------------
# run on while series:          
eq10optim <- function(r, gf, gm, agef = 10.5:49.5, agem = 10.5:65.5){
    (1 - (1 / sum(exp(-r * agem) * gm) + 1 / sum(exp(-r * agef) * gf))) ^ 2
}

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
                            LotkaRCoale(mf, Lyf, agef)
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
               r.f = LotkaRCoale(mf, Lyf, agef + .5))
        }, .BxymfUS = BxymfUS, .ExUS = ExUS, .LxmUS = LxmUS, .LxfUS = LxfUS))
# TODO: ES throws error, fix
names(BxymfES) <- yearsES
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
                            r.f = LotkaRCoale(mf, Lyf, agef + .5))
                }, .BxymfES = BxymfES, .ExES = ExES, .LxmES = LxmES, .LxfES = LxfES))

plot(yearsUS, rmfUS[,"r.mf"], type = 'l', ylim = c(-.02,.02))
polygon(c(yearsUS,rev(yearsUS)), c(rmfUS[,"r.m"],rev(rmfUS[,"r.f"])), col = "#CCCCCC50")
lines(yearsUS,rmfUS[,"r.m"], col = gray(.2), lwd = 2)
lines(yearsUS,rmfUS[,"r.f"], col = gray(.4), lwd = 2)



abline(h=0)
lines(yearsES, rmfES[,"r.mf"], col = "red")
lines(yearsUS,rmfUS[,"r.m"], col = gray(.2), lwd = 2)
lines(yearsUS,rmfUS[,"r.f"], col = gray(.4), lwd = 2)




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
ages    <- 10:65

agem <- 10:65
agef <- 10:49

# test year USA 1970:

Bma <- rowSums(BxymfUS[["1970"]][["Bxym"]])
Bfa <- colSums(BxymfUS[["1970"]][["Bxyf"]])

Mat <- with(ExUS, Male[Year == 1970])
Fat <- with(ExUS, Female[Year == 1970])

# eq 4
# Bm = [Hma*Mat] / vt
# where vt is a single number
# Bf = [Hfa*Fat] / ut
# Hm
vt <- .5
ut <- 1 - vt

# vt male prop
# ut female
# ut + vt = 1, and vary so keep SRB constant over time.
# Hma, Hfa constant proportions (like rates)

Hma_equal <- (Bma * vt) / Mat
Hfa_equal <- (Bfa * ut) / Fat

Bm <- sum(Mat * Hma_equal) / vt 
Bf <- sum(Fat * Hfa_equal) / ut
Bt <- Bm + Bf # eq 4,5
Bt == sum(c(Bma, Bfa)) # see, works

(SRB <- sum(Bm) / sum(Bf)) # eq 7
# a.k.a. eq. 8
SRB * vt / ut == sum(Hma_equal * Mat) / sum(Hfa_equal * Fat)
# eq 9
vt == sum(Hma_equal * Mat) / (sum(Hma_equal * Mat)+ SRB *sum(Hfa_equal * Fat))

# thus vt, ut change over time given an initial value choice
# and constant Hma, Hfa vectors

# initial values of vt, ut can be derived from their stable values
Pma <- LxmUS[,"1970"]
Pfa <- LxfUS[,"1970"]
theta_m <- Pma * Hma_equal # eq 13
theta_f <- Pfa * Hma_equal

# root eq, 17
a <- 0:110 + .5
r <- .001
sum(exp(-r * a)*(theta_m + theta_f))
# true r is the one that sets above to 1
# where theta already contains a fixed u,v pair
Mitra1978min <- function(r, vt, Bma, Bfa, Hma, Hfa, Pma, Pfa, a =  0:110 + .5){
    Hma <- (Bma * vt) / Mat
    Hfa <- (Bfa * (1 - vt)) / Fat
    (1 - sum(exp(-r * a)*(Pma * Hma + Pfa * Hfa))) ^ 2
}
re.5 <- optimize(f = Mitra1978min, interval = range(c(LotkaRCoale(Bma / Mat, Pma, a), LotkaRCoale(Bfa / Fat, Pfa, a))),
        vt = .5, Bma = Bma, Bfa = Bfa, Hma = Hma, Hfa = Hfa, Pma = Pma, Pfa = Pfa, tol = 1e-11)
re.6 <- optimize(f = Mitra1978min, interval = range(c(LotkaRCoale(Bma / Mat, Pma, a), LotkaRCoale(Bfa / Fat, Pfa, a))),
                vt = .6, Bma = Bma, Bfa = Bfa, Hma = Hma, Hfa = Hfa, Pma = Pma, Pfa = Pfa, tol = 1e-11)

sum(Hma * Mat) / (sum(Hma * Mat)+ SRB *sum(Hfa * Fat))

# needed in iteration to get vt
getvt <- function(Hma, Hfa, Mat, Fat, SRB){
    # returns new vt:
    sum(Hma * Mat) / (sum(Hma * Mat) + SRB * sum(Hfa * Fat))
}


# eq 32:
getMat <- function(r, Bmt, Pma, a = 0:110 + .5){
    Bmt * exp(-r*a) * Pma
}
Mat.r <- getMat(re.5$minimum, sum(Bma), Pma)
Fat.r <- getMat(re.5$minimum, sum(Bfa), Pfa)
sum(Mat.r)
sum(Fat.r)
BirthRate <- function(r, P, a = 0:110 + .5){
    1 / sum(exp(-r * a) * P)
}
BirthRate(re.5$minimum, Pma) -
BirthRate(re.5$minimum, Pfa)

# now choosing vt that minimizes final
        
pars <- c(r = 0.008, vt = .5)
Mitra1978min2 <- function(pars, Bma, Bfa, Ma0, Fa0, Pma, Pfa, a =  0:110 + .5){
    Hma <- (Bma * pars["vt"]) / Mat
    Hfa <- (Bfa * (1 - pars["vt"])) / Fat
    (1 - sum(exp(-pars["r"] * a)*(Pma * Hma + Pfa * Hfa))) ^ 2
}

pars <- optim(pars, Mitra1978min2, Bma = Bma, Bfa = Bfa, Ma0 = Mat, Fa0 = Fat, Pma = Pma, Pfa = Pfa, a =  0:110 + .5)$par
optim(pars, Mitra1978min2, Bma = Bma, Bfa = Bfa, Ma0 = Mat, Fa0 = Fat, Pma = Pma, Pfa = Pfa, a =  0:110 + .5)$par

# ideal: find stable mma, mfa, such that difference from original
#        mma, mfa is minimized,
# eq 39:
mm0 <- Bma / Mat 


net_parenthood <- theta_m + theta_f
net_maternity  <- Pfa * Bfa / Fat
net_paternity  <- Pma * Bma / Mat

plot(0:110, net_maternity, type = 'l', col = "pink")
lines(0:110, net_paternity,col = "blue")
lines(0:110, net_parenthood,col = "green")

# not intermediate?
sum(net_parenthood)
sum(net_maternity)
sum(net_paternity)

# example rates from paper:
#Male <- c(10.1, 76.4, 89.6, 56.9, 30.0, 12.7, 4.9, 2.7)
#Female <- c(32.7, 83.0, 67.5, 34.4, 16.0, 4.4, 0.3, 0)
#a5 <- seq(17.5,52.5,by=5)
#wmean(a5, Male)

# where we let vt vary over logit space to force between 0 and 1
getBmat <- function(Hma, Mat, vt){
    (Hma * Mat) / expit(vt)
}
getBfat <- function(Hfa, Fat, vt){
    (Hfa * Fat) / (1 - expit(vt))
}
# so vt can take any value between -Inf and Inf, this forces (6) to hold true
getvtopt <- function(Hma, Hfa, Mat, Fat, SRB = 1.054){
    vtoptim <- function(.vt, .Hma, .Hfa, .Mat, .Fat, .SRB){
        ((sum(getBmat(.Hma, .Mat, .vt)) / sum(getBfat(.Hfa, .Fat, .vt))) - .SRB) ^ 2
    }
    expit(optimize(vtoptim, interval = c(-10,10), .Hma = Hma, .Hfa = Hfa, .Mat = Mat, .Fat = Fat, .SRB = SRB)$minimum)
} 

getvt <- function(Hma, Hfa, Mat, Fat, SRB){
    # returns new vt:
    sum(Hma * Mat) / (sum(Hma * Mat) + SRB * sum(Hfa * Fat))
}
getvt(Hma, Hfa, Mat, Fat, SRB = 1.054)
getvtopt(Hma, Hfa, Mat, Fat, SRB = 1.054)

Hma <- Hma_equal
Hfa <- Hfa_equal

# repeated, using r and a given vt
Mitra1978minrv <- function(pars = c(r = 0, v0 = .5), Bma, Bfa, Pma, Pfa, a =  0:110 + .5){
    Hma <- (Bma * pars["v0"]) / Mat
    Hfa <- (Bfa * (1 - pars["v0"])) / Fat
    (1 - sum(exp(-pars["r"] * a) * (Pma * Hma + Pfa * Hfa))) ^ 2
}
optim(pars, Mitra1978minrv, Bma = Bma, Bfa = Bfa, Pma = Pma, Pfa = Pfa, a =  0:110 + .5 )
Mitra1978min2 <- function()
