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
# get Lx estimates for R0, r
dxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfES.Rdata"))) 
# make sum to 1
dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)

#------------------------------------------------------------
# compare with Lotka:
LxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmUS.Rdata"))) / 1e5
LxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfUS.Rdata"))) / 1e5
LxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmES.Rdata"))) / 1e5
LxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfES.Rdata"))) / 1e5



Bxy <- BxUS[[1]]

Exm <- with(ExUS, Male[Year == yearsUS[1]])
Exf <- with(ExUS, Female[Year == yearsUS[1]])

# harmonic mean of exposures:
HM <- function(x,y){Mna0(Minf0((2*x*y)/(x+y)))}
Hxy <- outer(Exm, Exf, HM)

# corresponding fertility rates:
FHxy <- Bxy / Hxy
ages <- 0:110

Fxf <- t(t(Bxy) / Exf)
Fxm <- Bxy / Exm


sum(FHxy * Exm)
sum(t(FHxy) * Exf)
# marginals are NOT constrained to intermediate values.
#plot(ages, colSums(FHxy), type = 'l', col = "red", lty =2)
#lines(ages, colSums(Fxf), col = "red")
#
#lines(ages, rowSums(Fxm), col = "blue")
#lines(ages, rowSums(FHxy), col = "blue", lty = 2)

Lxm <- LxmUS[,1]
Lxf <- LxfUS[,1]

#Lotka-ish equation:
.a <- .5:110.5
r <- .01
sum(exp(-r * .a) * Lxm * FHxy)
sum(exp(-r * .a) * Lxf * t(FHxy))

# residual function
LotkaHmin <- function(r, Lxm, Lxf, FH, .a = .5:110.5){
  (1 - sum(outer(exp(-r * .a) * Lxm, exp(-r * .a) * Lxf, HM) * FH)) ^ 2
}

optimize(LotkaHmin, interval = c(-.15,.15), Lxm = LxmUS[,10], Lxf = LxfUS[,10], 
        FH = BxUS[[10]] /outer(with(ExUS, Male[Year == yearsUS[10]]), with(ExUS, Female[Year == yearsUS[10]]), HM), 
        tol = 1e-15)$minimum
Fxm <- rowSums(BxUS[[10]] / with(ExUS, Male[Year == yearsUS[10]]))
Fxf <- colSums(t(t(BxUS[[10]]) / with(ExUS, Female[Year == yearsUS[10]])))
LotkaRCoale(Fxm,  LxmUS[,10], x = .5:110.5)
LotkaRCoale(Fxf,  LxfUS[,10], x = .5:110.5) - LotkaRCoale(Fxm,  LxmUS[,10], x = .5:110.5)

LotkaHRCoale <- function(FH, Lxm, Lxf, x = .5:110.5, maxit = 2e2, tol = 1e-15){
    # from Coale, Ansley J. (1957) A New Method for Calculating Lotka's r- the Intrinsic Rate of Growth in a Stable Population.
    # Population Studies, Vol. 11 no. 1, pp 92-94
    R0 <- sum(outer(Lxm, Lxf, HM) * FH)
    # first assuming a mean generation time of 29
    r2 <- log(R0) / 35
    for (i in 1:maxit){ # 10 is more than enough!
        r1 <- r2
        deltai <- 1 - sum(outer(exp(-r1 * x) * Lxm, exp(-r1 * x) * Lxf, HM) * FH) 
        # the mean generation time self-corrects 
        # according to the error produced by the Lotka equation
        r2 <- r1 - (deltai / (35 - (deltai / r1)))
        if (abs(r2 - r1) <= tol | zapsmall(abs(deltai)) <= tol){
            break
        }
    }
    return(r2)  
}

rUS <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxy, .Lxm, .Lxf, .Ex){
            
            Exm <- with(.Ex, Male[Year == as.integer(yr)])
            Exf <- with(.Ex, Female[Year == as.integer(yr)])
            # harmonic mean of exposures
            Hxy <- outer(Exm, Exf, HM)
            # harmonic rates (divide by two since working with both sexes)
            FxyH <- Mna0(Minf0((.Bxy[[yr]][["Bxym"]] + .Bxy[[yr]][["Bxyf"]]) / Hxy)) / 2
            # single sex rates
            Fxm <- rowSums(Mna0(Minf0(.Bxy[[yr]][["Bxym"]] / Exm)))
            Fxf <- colSums(Mna0(Minf0(t(t(.Bxy[[yr]][["Bxyf"]]) / Exf))))
            
            Lxf <- .Lxf[, yr]
            Lxm <- .Lxm[, yr]
            
            c(r.f = LotkaRCoale(Fxf, Lxf, x = .5:110.5),
              r.m = LotkaRCoale(Fxm, Lxm, x = .5:110.5),
              r.H = LotkaHRCoale(FH = FxyH, Lxm =Lxm, Lxf = Lxf))
        }, .Bxy = BxymfUS, .Lxm = LxmUS, .Lxf = LxfUS, .Ex = ExUS))
rownames(rUS) <- yearsUS


rES <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxy, .Lxm, .Lxf, .Ex){
                    
                    Exm <- with(.Ex, Male[Year == as.integer(yr)])
                    Exf <- with(.Ex, Female[Year == as.integer(yr)])
                    # harmonic mean of exposures
                    Hxy <- outer(Exm, Exf, HM)
                    # harmonic rates (divide by two since working with both sexes)
                    FxyH <- Mna0(Minf0((.Bxy[[yr]][["Bxym"]] + .Bxy[[yr]][["Bxyf"]]) / Hxy)) / 2
                    # single sex rates
                    Fxm <- rowSums(Mna0(Minf0(.Bxy[[yr]][["Bxym"]] / Exm)))
                    Fxf <- colSums(Mna0(Minf0(t(t(.Bxy[[yr]][["Bxyf"]]) / Exf))))
                    
                    Lxf <- .Lxf[, yr]
                    Lxm <- .Lxm[, yr]
                    
                    c(r.f = LotkaRCoale(Fxf, Lxf, x = .5:110.5),
                            r.m = LotkaRCoale(Fxm, Lxm, x = .5:110.5),
                            r.H = LotkaHRCoale(FH = FxyH, Lxm =Lxm, Lxf = Lxf))
                }, .Bxy = BxymfES, .Lxm = LxmES, .Lxf = LxfES, .Ex = ExES))
rownames(rES) <- yearsES

plot(yearsUS, rUS[,"r.f"], type = 'l', col = "red", ylim = c(-.02,.013))
lines(yearsUS, rUS[,"r.m"], col = "blue")
lines(yearsUS, rUS[,"r.H"], col = "green")

lines(yearsES, rES[,"r.f"], col = "red", lty=2)
lines(yearsES, rES[,"r.m"], col = "blue", lty=2)
lines(yearsES, rES[,"r.H"], col = "green", lty=2)





