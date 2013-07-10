setwd("/home/triffe/git/DISS/")
source("R/UtilityFunctions.R")
source("R/MeanFunctions.R")

yearsUS <- 1969:2009
yearsES <- 1975:2009

# Bxy totals  (not used()
BxUS  <- local(get(load("Data/USbirths/USBxy0_110.Rdata")))
BxES  <- local(get(load("Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

# repeated below for better SRB assumptions
BxymfES <- local(get(load("Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("Data/USbirths/USBxymf0_110.Rdata")))
names(BxymfES) <- yearsES
names(BxymfUS) <- yearsUS
# exposures, as such, straight from HMD, all ages 0-110, long form
ExUS  <- local(get(load("Data/Exposures/USexp.Rdata")))
ExES  <- local(get(load("Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
dxmUS <- local(get(load("Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("Data/HMD_dx/dxfES.Rdata"))) 
# make sum to 1
dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)

#------------------------------------------------------------
# compare with Lotka:
LxmUS <- local(get(load("Data/HMD_Lx/LxmUS.Rdata"))) / 1e5
LxfUS <- local(get(load("Data/HMD_Lx/LxfUS.Rdata"))) / 1e5
LxmES <- local(get(load("Data/HMD_Lx/LxmES.Rdata"))) / 1e5
LxfES <- local(get(load("Data/HMD_Lx/LxfES.Rdata"))) / 1e5
#------------------------------------------------------------

# Function Definitions:
#------------------------------------------------------------
# minimizer function for 1 sex ex-perspective renewal function:
# use with: optimize()
exOneSexMin <- compiler::cmpfun(function(r, dx, Fex, .a = .5:110.5){
    # get the overlapped / staggered dx structure
    N    <- length(Fex)
    dxM  <- matrix(0, ncol = N, nrow = N)
    dxi  <- dx
    for (i in 1:N){
        dxM[i, 1:length(dxi)  ] <- dxi 
        dxi <- dxi[2:length(dxi) ]
    }   
    1 - sum(rowSums(t(t(dxM) * exp(-r * .a))) * Fex)
})

exOneSexCoaleR <- function(Fex, dx, .a = .5:110.5, maxit = 1e2, tol = 1e-15, r.start = 0.001){  
    
    N    <- length(Fex)
    dxM  <- matrix(0, ncol = N, nrow = N)
    dxi  <- dx
    for (i in 1:N){
        dxM[i, 1:length(dxi)  ] <- dxi 
        dxi <- dxi[2:length(dxi) ]
    }     
    R0      <- sum(dxM * Fex)
    T.guess <- sum(.a * dxM * Fex) / R0 # assuming r = 0
    r.i      <- log(R0) / T.guess
   
    # be careful to discount Fex by SRB appropriately for males / females
    # prior to specification
    # Based on Coale (1957)
    for (i in 1:maxit){ # 15 is more than enough!
       #cat(r2,i,"\n")
        deltai <- 1 - sum(rowSums(t(t(dxM) * exp(-r.i * .a))) * Fex)
        # the mean generation time self-corrects 
        # according to the error produced by the Lotka equation
        r.i <- r.i - (deltai / (T.guess - (deltai / r.i))) 
        if (abs(deltai) <= tol){
            break
        }
    }
    return(r.i)  
}

ex1SexStableAge <- function(r, Fex, dx, .a = .5:110.5){
    N    <- length(Fex)
    dxM  <- matrix(0, ncol = N, nrow = N)
    dxi  <- dx
    for (i in 1:N){
        dxM[i, 1:length(dxi)  ] <- dxi 
        dxi <- dxi[2:length(dxi) ]
    }  
    # birth rate
    b <- 1 / sum(t(dxM) * exp(-r * .a))
    b * rowSums(t(t(dxM) * exp(-r * .a)))
}

exOneSexTy <- function(r, Fex, dx, .a = .5:110.5){
    N    <- length(Fex)
    dxM  <- matrix(0, ncol = N, nrow = N)
    dxi  <- dx
    for (i in 1:N){
        dxM[i, 1:length(dxi)  ] <- dxi 
        dxi <- dxi[2:length(dxi) ]
    }     
    sum(rowSums(.a*t(t(dxM) * exp(-r * .a))) * Fex) /
            sum(rowSums(t(t(dxM) * exp(-r * .a))) * Fex)   
}


rmUS <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bx, .Ex, .dx, MF, rc,Bxymf){       
                    Ex <- .Ex[.Ex$Year == as.integer(yr),MF]
                    .ey <- rowSums(ExpectedDx(Ex, .dx[, yr]))
                    .by <-  rowSums(ExpectedDx(rc(.Bx[[yr]][[Bxymf]], na.rm = TRUE), .dx[, yr]))
                    .fy <- .by / .ey
                    r <- exOneSexCoaleR(Fex = .fy, dx =.dx[, yr], tol =1e-12)
                    c(r = r, Ty = exOneSexTy(r, Fex = .fy, dx =.dx[, yr]))
                }, .Bx = BxymfUS, .Ex = ExUS, .dx = dxmUS, Bxymf = "Bxym", MF = "Male", rc = rowSums))

rfUS <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bx, .Ex, .dx, MF, rc,Bxymf){       
                    Ex <- .Ex[.Ex$Year == as.integer(yr),MF]
                    .ey <- rowSums(ExpectedDx(Ex, .dx[, yr]))
                    .by <-  rowSums(ExpectedDx(rc(.Bx[[yr]][[Bxymf]], na.rm = TRUE), .dx[, yr]))
                    .fy <- Minf0(Mna0(.by / .ey))
                    r <- exOneSexCoaleR(Fex = .fy  , dx = Mna0(.dx[, yr]), tol =1e-12)
                    c(r = r, Ty = exOneSexTy(r, Fex = .fy, dx = Mna0(.dx[, yr])))
                }, .Bx = BxymfUS, .Ex = ExUS, .dx = dxfUS, Bxymf = "Bxyf", MF = "Female", rc = colSums))

rmES <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bx, .Ex, .dx, MF, rc,Bxymf){       
                    Ex <- .Ex[.Ex$Year == as.integer(yr),MF]
                    .ey <- rowSums(ExpectedDx(Ex, Mna0(.dx[, yr])))
                    .by <-  rowSums(ExpectedDx(rc(.Bx[[yr]][[Bxymf]], na.rm = TRUE), Mna0(.dx[, yr])))
                    .fy <-  Minf0(Mna0(.by / .ey))
                    r <- exOneSexCoaleR(Fex = .fy, dx =.dx[, yr], tol =1e-12)
                    c(r = r, Ty = exOneSexTy(r, Fex = .fy, dx = Mna0(.dx[, yr])))
                }, .Bx = BxymfES, .Ex = ExES, .dx = dxmES,  Bxymf = "Bxym",MF = "Male", rc = rowSums))

rfES <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bx, .Ex, .dx, MF, rc,Bxymf){       
                    Ex <- .Ex[.Ex$Year == as.integer(yr),MF]
                    .ey <- rowSums(ExpectedDx(Ex,  Mna0(.dx[, yr])))
                    .by <-  rowSums(ExpectedDx(rc(.Bx[[yr]][[Bxymf]], na.rm = TRUE),  Mna0(.dx[, yr])))
                    .fy <-  Minf0(Mna0(.by / .ey))
                    r <- exOneSexCoaleR(Fex = .fy  , dx =.dx[, yr], tol =1e-12)
                    c(r = r, Ty = exOneSexTy(r, Fex = .fy, dx = Mna0(.dx[, yr])))
                }, .Bx = BxymfES, .Ex = ExES, .dx = dxfES,Bxymf = "Bxyf", MF = "Female", rc = colSums))

#save(rmUS, file = "Data/results/exSingleSex/rmUS.Rdata")
#save(rfUS, file = "Data/results/exSingleSex/rfUS.Rdata")
#save(rmES, file = "Data/results/exSingleSex/rmES.Rdata")
#save(rfES, file = "Data/results/exSingleSex/rfES.Rdata")


#--------------------------------------------------------------
# Tables for appendix 'zAppendix.exSingleSexLotka.tex'
rmUS <- local(get(load("Data/results/exSingleSex/rmUS.Rdata")))
rfUS <- local(get(load(file = "Data/results/exSingleSex/rfUS.Rdata")))
rmES <- local(get(load(file = "Data/results/exSingleSex/rmES.Rdata")))
rfES <- local(get(load(file = "Data/results/exSingleSex/rfES.Rdata")))
exRepUSm <- cbind(rmUS, R0 = exp(rmUS[,1]*rmUS[,2]) )
exRepUSf <- cbind(rfUS, R0 = exp(rfUS[,1]*rfUS[,2]) )
exRepESm <- cbind(rmES, R0 = exp(rmES[,1]*rmES[,2]) )
exRepESf <- cbind(rfES, R0 = exp(rfES[,1]*rfES[,2]) )
colnames(exRepESf) <-colnames(exRepESm) <-colnames(exRepUSf) <-colnames(exRepUSm) <- c("$r$","$T$","$R_0$")
rownames(exRepUSm) <- rownames(exRepUSf) <- yearsUS
rownames(exRepESm) <- rownames(exRepESf) <- yearsES
library(xtable)
print(xtable(exRepUSm, digits = c(0,4,2,3), align = c("c","c","c","c")),
        sanitize.colnames.function = identity, 
        file = "latex/xtables/exRepUSm.tex",floating=FALSE)
print(xtable(exRepUSf, digits = c(0,4,2,3), align = c("c","c","c","c")),
        sanitize.colnames.function = identity, 
        file = "latex/xtables/exRepUSf.tex",floating=FALSE)
print(xtable(exRepESm, digits = c(0,4,2,3), align = c("c","c","c","c")),
        sanitize.colnames.function = identity, 
        file = "latex/xtables/exRepESm.tex",floating=FALSE)
print(xtable(exRepESf, digits = c(0,4,2,3), align = c("c","c","c","c")),
        sanitize.colnames.function = identity, 
        file = "latex/xtables/exRepESf.tex",floating=FALSE)
#--------------------------------------------------------------
# calculate age-classified Lotka r:

rmLUS <- unlist(lapply(as.character(yearsUS), function(yr, .Bx, .Ex, .Lx, Bxymf, MF, rc){       
                    Fx  <- Minf0(Mna0(rc(.Bx[[yr]][[Bxymf]], na.rm = TRUE) / .Ex[.Ex$Year == as.integer(yr),MF]))
                    LotkaRCoale(Fx, Lx = .Lx[,yr], x = .5:110.5)
                }, .Bx = BxymfUS, .Ex = ExUS, .Lx = LxmUS, Bxymf = "Bxym", MF = "Male", rc = rowSums))

rfLUS <- unlist(lapply(as.character(yearsUS), function(yr, .Bx, .Ex, .Lx, Bxymf, MF, rc){       
                    Fx  <- Minf0(Mna0(rc(.Bx[[yr]][[Bxymf]], na.rm = TRUE) / .Ex[.Ex$Year == as.integer(yr),MF]))
                    LotkaRCoale(Fx, Lx = .Lx[,yr], x = .5:110.5)
                }, .Bx = BxymfUS, .Ex = ExUS, .Lx = LxfUS, Bxymf = "Bxyf", MF = "Female", rc = colSums))

rmLES <- unlist(lapply(as.character(yearsES), function(yr, .Bx, .Ex, .Lx, Bxymf, MF, rc){       
                    Fx  <- Minf0(Mna0(rc(.Bx[[yr]][[Bxymf]], na.rm = TRUE) / .Ex[.Ex$Year == as.integer(yr),MF]))
                    LotkaRCoale(Fx, Lx = .Lx[,yr], x = .5:110.5)
                }, .Bx = BxymfES, .Ex = ExES, .Lx = LxmES, Bxymf = "Bxym", MF = "Male", rc = rowSums))

rfLES <- unlist(lapply(as.character(yearsES), function(yr, .Bx, .Ex, .Lx, Bxymf, MF, rc){       
                    Fx  <- Minf0(Mna0(rc(.Bx[[yr]][[Bxymf]], na.rm = TRUE) / .Ex[.Ex$Year == as.integer(yr),MF]))
                    LotkaRCoale(Fx, Lx = .Lx[,yr], x = .5:110.5)
                }, .Bx = BxymfES, .Ex = ExES, .Lx = LxfES, Bxymf = "Bxyf", MF = "Female", rc = colSums))


#plot(yearsUS, exRepUSm[,2], type = 'l', col = "blue", ylim = c(40,60))
#lines(yearsUS, exRepUSm[,2], col ="red")
#lines(yearsES, exRepUSm[,2], col ="red", lty = 2)
#lines(yearsES, exRepUSm[,2], col ="blue", lty = 2)
#
#plot(yearsUS,R0mUS, type = 'l', ylim = c(.8,1.4))
#lines(yearsUS, exp(rfUS[,1]*rfUS[,2])  , type = 'l')
#
#plot(yearsUS, exRepUSm[,1], type = 'l', col = "blue", ylim = c(-.02,.015))
#lines(yearsUS, exRepUSf[,1], col ="red")
#lines(yearsES, exRepESf[,1], col ="red", lty = 2)
#lines(yearsES, exRepESm[,1], col ="blue", lty = 2)
#abline(h=0)
#lines(yearsUS, rmLUS, col = "royalblue", lwd = 2)
#lines(yearsUS, rfLUS, col ="pink", lwd = 2)
#lines(yearsES, rfLES, col ="pink", lty = 2, lwd = 2)
#lines(yearsES, rmLES, col ="royalblue", lty = 2, lwd = 2)
#abline(h=0)

# -----------------------------------------------------
# Figure
#pdf("latex/Figures/exLotka1sex.pdf", height = 5, width = 5)
#par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
#plot(yearsUS, rmUS[, 1], type = 'n', ylim = c(-.016,.01),xlim = c(1968,2010), axes = FALSE,
#        xlab = "", ylab = "",
#        panel.first = list(rect(1968,-.016,2010,.01,col = gray(.95), border=NA),
#                abline(h = seq(-.016,.01,by = .002), col = "white"),
#                abline(v = seq(1970, 2010, by = 5), col = "white"),
#                text(1968, seq(-.016,.01,by = .002),seq(-.016,.01,by = .002), pos = 2, cex = .8, xpd = TRUE),
#                text(seq(1970, 2010, by = 10),-.016, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
#                text(1990, -.0175, "Year", cex = 1, pos = 1, xpd = TRUE),
#                text(1966,.0115, "r", cex = 1, xpd = TRUE)))
#
#lines(yearsUS, rmUS[, 1],col = gray(.2), lwd = 2)
#lines(yearsUS, rfUS[, 1],col = gray(.5), lwd = 2.5)
#lines(yearsES, rmES[, 1],col = gray(.2), lwd = 2, lty = 5)
#lines(yearsES, rfES[, 1],col = gray(.5), lwd = 2.5, lty = 5)
#
#legend(1970,-.01, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
#        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
#dev.off()

#-----------------------------------------------
# add normal Lotka r to series
# it's pretty noisy, but will pass
pdf("latex/Figures/exLotka1sex2.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rmUS[, 1], type = 'n', ylim = c(-.02,.011),xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,-.02,2010,.011,col = gray(.95), border=NA),
                abline(h = seq(-.02,.011,by = .0025), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.02,.011,by = .005),seq(-.02,.011,by = .005), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.02, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.022, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,.0125, "r", cex = 1, xpd = TRUE)))

lines(yearsUS, rmUS[, 1],col = gray(.2), lwd = 2.5)
lines(yearsUS, rfUS[, 1],col = gray(.5), lwd = 3)
lines(yearsES, rmES[, 1],col = gray(.2), lwd = 2.5, lty = 5)
lines(yearsES, rfES[, 1],col = gray(.5), lwd = 3, lty = 5)

lines(yearsUS, rmLUS,col = gray(.2), lwd = 1)
lines(yearsUS, rfLUS,col = grey(.5), lwd = 1)
lines(yearsES, rmLES,col = gray(.2), lwd = 1, lty = 5)
lines(yearsES, rfLES,col = gray(.5), lwd = 1, lty = 5)

legend(1969,-.0085, lty = c(1,1,5,5,1,1,5,5), 
        col = gray(c(.2,.5,.2,.5,.2,.5,.2,.5)), 
        lwd = c(2.5,3,2.5,3,1,1,1,1),
        bty = "n",
        legend = c(expression(US~males~e[y]), expression(US~females~e[y]), 
                expression(ES~males~e[y]), expression(ES~females~e[y]),
                "US males ageexOneSexMin", "US females age", "ES males age", "ES females age"), 
        xpd = TRUE, cex = .8)
dev.off()




# Lotka stable age:

LotkaStableAge <- compiler::cmpfun(function(r, Lx, Fx,.a = .5:110.5){
    (1 / sum(exp(-r * .a) * Lx)) * exp(-r * .a) * Lx
})

# --------------------------------------------------------------
# year t ex-structure with its own stable structure. Stability:
# verbose code!! separate calcs for 'true' and simulated thetas
rownames(rmUS) <- yearsUS
rownames(rfUS) <- yearsUS
rownames(rmES) <- yearsES
rownames(rfES) <- yearsES
names(rmLUS) <- yearsUS
names(rfLUS) <- yearsUS
names(rmLES) <- yearsES
names(rfLES) <- yearsESexOneSexMin

PxUS <- local(get(load("Data/HMD_Px/PxUS.Rdata")))
PxES <- local(get(load("Data/HMD_Px/PxES.Rdata")))
#-----------------------------
yr <- "1980"
DiffCoefryUSm <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bx, .Ex, .Px, .dx, .Lx, rmat, rL){  

                    Bx      <- rowSums(.Bx[[yr]][["Bxym"]], na.rm = TRUE)
                    Px      <- .Px$Male[.Px$Year == as.integer(yr)]
                    Ex      <- .Ex$Male[.Ex$Year == as.integer(yr)]
                    
                    Py      <- rowSums(ExpectedDx(Px, .dx[, yr]))
                    Ey      <- rowSums(ExpectedDx(Ex, .dx[, yr]))
                    By      <- rowSums(ExpectedDx(Bx, .dx[, yr]))
                    Fex     <- Minf0(Mna0(By / Ey))
                    # stable structure
                    cyst    <- ex1SexStableAge(r = rmat[yr,"r"], Fex = Fex, dx = .dx[, yr])
                    
                    # preesnt structure
                    cy      <- Py / sum(Py)
                    
                    cast    <- LotkaStableAge(r = rL[yr], Lx = .Lx[,yr], Fx = Bx / Ex)
                    ca.now  <- Px / sum(Px)
                    # difference coef
                    c(extheta = 1-sum(pmin(cyst, cy)), Lotkatheta =  1-sum(pmin(cast, ca.now)))
                }, .Bx = BxymfUS, .Ex = ExUS, .Px = PxUS, .Lx = LxmUS, .dx = dxmUS, rmat = rmUS, rL = rmLUS))

DiffCoefryUSf <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bx, .Ex, .Px, .dx, .Lx, rmat, rL){  
                   
               
                    Bx      <- colSums(.Bx[[yr]][["Bxyf"]], na.rm = TRUE)
                    Px      <- .Px$Female[.Px$Year == as.integer(yr)]
                    Ex      <- .Ex$Female[.Ex$Year == as.integer(yr)]
                    
                    Py      <- rowSums(ExpectedDx(Px, .dx[, yr]))
                    Ey      <- rowSums(ExpectedDx(Ex, .dx[, yr]))
                    By      <- rowSums(ExpectedDx(Bx, .dx[, yr]))
                    Fex     <- Minf0(Mna0(By / Ey))
                    # stable structure
                    cyst    <- ex1SexStableAge(r = rmat[yr,"r"], Fex = Fex, dx = .dx[, yr])
                    
                    # preesnt structure
                    cy      <-exOneSexMin Py / sum(Py)
                    
                    cast    <- LotkaStableAge(r = rL[yr], Lx = .Lx[,yr], Fx = Bx / Ex)
                    ca.now  <- Px / sum(Px)
                    # difference coef
                    c(extheta = 1-sum(pmin(cyst, cy)), Lotkatheta =  1-sum(pmin(cast, ca.now)))
                }, .Bx = BxymfUS, .Ex = ExUS, .Px = PxUS, .Lx = LxfUS, .dx = dxfUS, rmat = rfUS, rL = rfLUS))

DiffCoefryESm <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bx, .Ex, .Px, .dx, .Lx, rmat, rL){  
                   
               
                    Bx      <- rowSums(.Bx[[yr]][["Bxym"]], na.rm = TRUE)
                    Px      <- .Px$Male[.Px$Year == as.integer(yr)]
                    Ex      <- .Ex$Male[.Ex$Year == as.integer(yr)]
                    
                    Py      <- rowSums(ExpectedDx(Px, .dx[, yr]))
                    Ey      <- rowSums(ExpectedDx(Ex, .dx[, yr]))
                    By      <- rowSums(ExpectedDx(Bx, .dx[, yr]))
                    Fex     <- Minf0(Mna0(By / Ey))
                    # stable structure
                    cyst    <- ex1SexStableAge(r = rmat[yr,"r"], Fex = Fex, dx = .dx[, yr])
                    
                    # preesnt structure
                    cy      <- Py / sum(Py)
                    
                    cast    <- LotkaStableAge(r = rL[yr], Lx = .Lx[,yr], Fx = Bx / Ex)
                    ca.now  <- Px / sum(Px)
                    # difference coef
                    c(extheta exOneSexMin= 1-sum(pmin(cyst, cy)), Lotkatheta =  1-sum(pmin(cast, ca.now)))
                }, .Bx = BxymfES, .Ex = ExES, .Px = PxES, .Lx = LxmES, .dx = dxmES, rmat = rmES, rL = rmLES))

DiffCoefryESf <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bx, .Ex, .Px, .dx, .Lx, rmat, rL){  
                   
               
                    Bx      <- colSums(.Bx[[yr]][["Bxyf"]], na.rm = TRUE)
                    Px      <- .Px$Female[.Px$Year == as.integer(yr)]
                    Ex      <- .Ex$Female[.Ex$Year == as.integer(yr)]
                    
                    Py      <- rowSums(ExpectedDx(Px, .dx[, yr]))
                    Ey      <- rowSums(ExpectedDx(Ex, .dx[, yr]))
                    By      <- rowSums(ExpectedDx(Bx, .dx[, yr]))
                    Fex     <- Minf0(Mna0(By / Ey))
                    # stable structure
                    cyst    <- ex1SexStableAge(r = rmat[yr,"r"], Fex = Fex, dx = .dx[, yr])
                    
                    # preesnt exOneSexMinstructure
                    cy      <- Py / sum(Py)
                    
                    cast    <- LotkaStableAge(r = rL[yr], Lx = .Lx[,yr], Fx = Bx / Ex)
                    ca.now  <- Px / sum(Px)
                    # difference coef
                    c(extheta = 1-sum(pmin(cyst, cy)), Lotkatheta =  1-sum(pmin(cast, ca.now)))
                }, .Bx = BxymfES, .Ex = ExES, .Px = PxES, .Lx = LxfES, .dx = dxfES, rmat = rfES, rL = rfLES))


mxmUS <- local(get(load("Data/HMD_mux/muxmUS.Rdata"))) 
mxfUS <- local(get(load("Data/HMD_mux/muxfUS.Rdata"))) 
mxmES <- local(get(load("Data/HMD_mux/muxmES.Rdata"))) 
mxfES <- local(get(load("Data/HMD_mux/muxfES.Rdata"))) 

library(parallel)
# SAVED:
# US
# Males
#DiffCoefryUSthetam <- do.call(rbind,mclapply(as.character(yearsUS), function(yr, .Bx, .Ex, .Px, .mx, rmat, rL){  
#                 
#                    Px      <- .Px$Male[.Px$Year == as.integer(yr)]
#                    Ex      <- .Ex$Male[.Ex$Year == as.integer(yr)]
#                    Bx      <- rowSums(.Bx[[yr]][["Bxym"]])
#                    mx      <- .mx[, yr]
#                                   
#                    mxsim <- matrix(rpois(n=length(mx)*1000, lambda = mx * Ex), ncol = 1000) / Ex
#                    
#                    Lxsim <-  apply(mxsim, 2, function(x){
#                                mx2LxHMD(x)}
#                    )
#                    
#                    dxSim <- apply(mxsim, 2, function(x){
#                                        mx2dxHMD(x)}
#                                   )
#                    Pysim <- apply(dxSim, 2, function(dx, .Px){
#                                rowSums(ExpectedDx(.Px, dx))
#                            }, .Px = Px)
#                    Cysim <- t(t(Pysim) / colSums(Pysim))
#                    
#                    Eysim <- apply(dxSim, 2, function(dx, .Ex){
#                                rowSums(ExpectedDx(.Ex, dx))
#                            }, .Ex = Ex)
#                    BirthsSim <- matrix(rpois(n=length(Bx)*1000, lambda = Bx), ncol = 1000)
#                    
#                    BySim <- exOneSexMinBirthsSim * 0 
#               
#                    for (i in 1:1000){
#                        BySim[,i] <- rowSums(ExpectedDx(BirthsSim[,i], dxSim[,i]))
#                    }
#                    FexSim    <- Minf0(Mna0(BySim / Eysim))
#                  
#                    FxSim     <- Minf0(Mna0(BirthsSim / Ex))
#                     # stable structure
#                    castsim <- cystsim <- FexSim * 0
#                    for (i in 1:1000){
#                        cystsim[,i]    <- ex1SexStableAge(r = rmat[yr,"r"], Fex = FexSim[,i], dx = dxSim[,i])
#                        castsim[,i]    <- LotkaStableAge(rL[yr], Lx = Lxsim[,i], Fx = FxSim[,i])
#                    }
#                    
#                    # preesnt structure
#                      
#                    # difference coef
#                    thetas  <- 1-colSums(pmin(cystsim, Cysim))
#                    thetasL <- 1-colSums(pmin(castsim, Px / sum(Px)))
#                    ratio   <- thetasL / thetas
#                    c(ExQ = quantile(thetas, probs = c(.025,.975)), LQ = quantile(thetasL, probs = c(.025,.975)), ratio = quantile(ratio, probs = c(.025,.975)))
#                }, .Bx = BxymfUS, .Ex = ExUS, .Px = PxUS, .mx = mxmUS, rL = rmLUS, rmat = rmUS, mc.cores = 2))
#save(DiffCoefryUSthetam,file="Data/rDecompResults/DiffCoefryUSthetam.Rdata")
## Females
#DiffCoefryUSthetaf <- do.call(rbind,mclapply(as.character(yearsUS), function(yr, .Bx, .Ex, .Px, .mx, rmat, rL){  
#                 
#                    Px      <- .Px$Female[.Px$Year == as.integer(yr)]
#                    Ex      <- .Ex$Female[.Ex$Year == as.integer(yr)]
#                    Bx      <- colSums(.Bx[[yr]][["Bxyf"]])
#                    mx      <- .mx[, yr]
#                                   
#                    mxsim <- matrix(rpois(n=length(mx)*1000, lambda = mx * Ex), ncol = 1000) / Ex
#                    
#                    Lxsim <-  apply(mxsim, 2, function(x){
#                                mx2LxHMD(x)}
#                    )
#                    
#                    dxSim <- apply(mxsim, 2, function(x){
#                                        mx2dxHMD(x)}
#                                   )
#                    Pysim <- apply(dxSim, 2, function(dx, .Px){
#                                rowSums(ExpectedDx(.Px, dx))
#                            }, .Px = Px)
#                    Cysim <- t(t(Pysim) / colSums(Pysim))
#                    
#                    Eysim <- apply(dxSim, 2, function(dx, .Ex){
#                                rowSums(ExpectedDx(.Ex, dx))
#                            }, .Ex = Ex)
#                    BirthsSim <- matrix(rpois(n=length(Bx)*1000, lambda = Bx), ncol = 1000)
#                    
#                    BySim <- BirthsSim * 0 
#               
#                    for (i in 1:1000){
#                        BySim[,i] <- rowSums(ExpectedDx(BirthsSim[,i], dxSim[,i]))
#                    }
#                    FexSim    <- Minf0(Mna0(BySim / Eysim))
#                  
#                    FxSim     <- Minf0(Mna0(BirthsSim / Ex))
#                     # stable structure
#                    castsim <- cystsim <- FexSim * 0
#                    for (i in 1:1000){
#                        cystsim[,i]    <- ex1SexStableAge(r = rmat[yr,"r"], Fex = FexSim[,i], dx = dxSim[,i])
#                        castsim[,i]    <- LotkaStableAge(rL[yr], Lx = Lxsim[,i], Fx = FxSim[,i])
#                    }
#                    
#                    # preesnt structure
#                      
#                    # difference coef
#                    thetas  <- 1-colSums(pmin(cystsim, Cysim))
#                    thetasL <- 1-colSums(pmin(castsim, Px / sum(Px)))
#                    ratio   <- thetasL / thetas
#                    c(ExQ = quantile(thetas, probs = c(.025,.975)), LQ = quantile(thetasL, probs = c(.025,.975)), ratio = quantile(ratio, probs = c(.025,.975)))
#                    
#                }, .Bx = BxymfUS, .Ex = ExUS, .Px = PxUS, .mx = mxfUS, rL = rfLUS, rmat = rfUS, mc.cores = 2))
#save(DiffCoefryUSthetaf,file="Data/rDecompResults/DiffCoefryUSthetaf.Rdata")
## Spain
## Males
#DiffCoefryESthetam <- do.call(rbind,mclapply(as.character(yearsES), function(yr, .Bx, .Ex, .Px, .mx, rmat, rL){  
#                    Px      <- .Px$Male[.Px$Year == as.integer(yr)]
#                    Ex      <- .Ex$Male[.Ex$Year == as.integer(yr)]
#                    Bx      <- rowSums(.Bx[[yr]][["Bxym"]])
#                    mx      <- .mx[, yr]
#                                   
#                    mxsim <- Minf0(Mna0(matrix(rpois(n=length(mx)*1000, lambda = mx * Ex), ncol = 1000) / Ex))
#                    
#                    Lxsim <-  apply(mxsim, 2, function(x){
#                                mx2LxHMD(x)}
#                    )
#                
#                    dxSim <- apply(mxsim, 2, function(x){
#                                        mx2dxHMD(x)}
#                                   )
#                    Pysim <- apply(dxSim, 2, function(dx, .Px){
#                                rowSums(ExpectedDx(.Px, dx))
#                            }, .Px = Px)
#                    Cysim <- t(t(Pysim) / colSums(Pysim))
#                    
#                    Eysim <- apply(dxSim, 2, function(dx, .Ex){
#                                rowSums(ExpectedDx(.Ex, dx))
#                            }, .Ex = Ex)
#                    BirthsSim <- matrix(rpois(n=length(Bx)*1000, lambda = Bx), ncol = 1000)
#                    
#                    BySim <- BirthsSim * 0 
#               
#                    for (i in 1:1000){
#                        BySim[,i] <- rowSums(ExpectedDx(BirthsSim[,i], dxSim[,i]))
#                    }
#                    FexSim    <- Minf0(Mna0(BySim / Eysim))
#                  
#                    FxSim     <- Minf0(Mna0(BirthsSim / Ex))
#                     # stable structure
#                    castsim <- cystsim <- FexSim * 0
#                    for (i in 1:1000){
#                        cystsim[,i]    <- ex1SexStableAge(r = rmat[yr,"r"], Fex = FexSim[,i], dx = dxSim[,i])
#                        castsim[,i]    <- LotkaStableAge(rL[yr], Lx = Lxsim[,i], Fx = FxSim[,i])
#                    }
#                    # preesnt structure
#                    
#                    # difference coef
#                    thetas  <- 1-colSums(pmin(cystsim, Cysim))
#                    thetasL <- 1-colSums(pmin(castsim, Px / sum(Px, na.rm = TRUE)))
#                    ratio   <- thetasL / thetas
#                    c(ExQ = quantile(thetas, probs = c(.025,.975)), 
#                            LQ = quantile(thetasL, probs = c(.025,.975)), 
#                            ratio = quantile(ratio, probs = c(.025,.975)))
#                    
#                }, .Bx = BxymfES, .Ex = ExES, .Px = PxES, .mx = mxmES, rL = rmLES, rmat = rmES, mc.cores = 2))
#save(DiffCoefryESthetam,file="Data/rDecompResults/DiffCoefryESthetam.Rdata")
## Females
#DiffCoefryESthetaf <- do.call(rbind,mclapply(as.character(yearsES), function(yr, .Bx, .Ex, .Px, .mx, rmat, rL){  
#                 
#                    Px      <- .Px$Female[.Px$Year == as.integer(yr)]
#                    Ex      <- .Ex$Female[.Ex$Year == as.integer(yr)]
#                    Bx      <- colSums(.Bx[[yr]][["Bxyf"]])
#                    mx      <- .mx[, yr]
#                                   
#                    mxsim <- matrix(rpois(n=length(mx)*1000, lambda = mx * Ex), ncol = 1000) / Ex
#                    
#                    Lxsim <-  apply(mxsim, 2, function(x){
#                                mx2LxHMD(x)}
#                    )
#                    
#                    dxSim <- apply(mxsim, 2, function(x){
#                                        mx2dxHMD(x)}
#                                   )
#                    Pysim <- apply(dxSim, 2, function(dx, .Px){
#                                rowSums(ExpectedDx(.Px, dx))
#                            }, .Px = Px)
#                    Cysim <- t(t(Pysim) / colSums(Pysim))
#                    
#                    Eysim <- apply(dxSim, 2, function(dx, .Ex){
#                                rowSums(ExpectedDx(.Ex, dx))
#                            }, .Ex = Ex)
#                    BirthsSim <- matrix(rpois(n=length(Bx)*1000, lambda = Bx), ncol = 1000)
#                    
#                    BySim <- BirthsSim * 0 
#               
#                    for (i in 1:1000){
#                        BySim[,i] <- rowSums(ExpectedDx(BirthsSim[,i], dxSim[,i]))
#                    }
#                    FexSim    <- Minf0(Mna0(BySim / Eysim))
#                  
#                    FxSim     <- Minf0(Mna0(BirthsSim / Ex))
#                     # stable structure
#                    castsim <- cystsim <- FexSim * 0
#                    for (i in 1:1000){
#                        cystsim[,i]    <- ex1SexStableAge(r = rmat[yr,"r"], Fex = FexSim[,i], dx = dxSim[,i])
#                        castsim[,i]    <- LotkaStableAge(rL[yr], Lx = Lxsim[,i], Fx = FxSim[,i])
#                    }
#                    
#                    # preesnt structure
#                      
#                    # difference coef
#                    thetas  <- 1-colSums(pmin(cystsim, Cysim))
#                    thetasL <- 1-colSums(pmin(castsim, Px / sum(Px)))
#                    ratio   <- thetasL / thetas
#                    c(ExQ = quantile(thetas, probs = c(.025,.975)), LQ = quantile(thetasL, probs = c(.025,.975)), ratio = quantile(ratio, probs = c(.025,.975)))
#                    
#                }, .Bx = BxymfES, .Ex = ExES, .Px = PxES, .mx = mxfES, rL = rfLES, rmat = rfES, mc.cores = 2))
#save(DiffCoefryESthetaf,file="Data/rDecompResults/DiffCoefryESthetaf.Rdata")

DiffCoefryUSthetam <- local(get(load("Data/rDecompResults/DiffCoefryUSthetam.Rdata")))
DiffCoefryUSthetaf <- local(get(load("Data/rDecompResults/DiffCoefryUSthetaf.Rdata")))
DiffCoefryESthetam <- local(get(load("Data/rDecompResults/DiffCoefryESthetam.Rdata")))
DiffCoefryESthetaf <- local(get(load("Data/rDecompResults/DiffCoefryESthetaf.Rdata")))

#plot(NULL, type = "n", xlim = c(1968,2010) , ylim = c(0,.2))
#lines(yearsUS, DiffCoefryUSm, col = gray(.2), lwd = 1, lty = 1)                
#lines(yearsUS, DiffCoefryUSf, col = gray(.2), lwd = 1, lty = 1)
#lines(yearsES, DiffCoefryESm, col = gray(.4), lwd = 1.2, lty = 4)                
#lines(yearsES, DiffCoefryESf, col = gray(.4), lwd = 1.2, lty = 4)

# plot it!
pdf("latex/Figures/exPyramidPresentvsStableDivergence.pdf", height = 4.5, width = 4.5)
par(mai = c(.5,.5,.5,.3), xaxs = "i", yaxs = "i")
plot(yearsUS, DiffCoefryUSm[,1], type = 'n', ylim = c(0,.3), xlim = c(1968, 2010), 
        axes = FALSE, xlab = "", ylab = "",
        panel.first = list(rect(1968, 0, 2010, .3, col = gray(.95), border = NA),
                abline(v = seq(1970, 2010, by = 5),col = "white"),
                abline(h = seq(0, .3, by = .025),col = "white"),
                text(1990,-.02, "Year", xpd = TRUE),
                text(1965,.22, expression(theta), xpd = TRUE),
                text(seq(1970, 2010, by = 5), 0, seq(1970, 2010, by = 5), pos = 1 ,cex = .7, xpd = TRUE),
                text(1968, seq(0, .2, by = .05), seq(0, .2, by = .05), pos = 2, cex = .7, xpd = TRUE)
        ))
polygon(c(yearsUS,rev(yearsUS)), c(DiffCoefryUSthetam[,1],rev(DiffCoefryUSthetam[,2])), 
        border = gray(.2), lwd = .5, col = "#BBBBBB40")
polygon(c(yearsUS,rev(yearsUS)), c(DiffCoefryUSthetaf[,1],rev(DiffCoefryUSthetaf[,2])), 
        border = gray(.2), lwd = .5, col = "#BBBBBB40")

polygon(c(yearsES,rev(yearsES)), c(DiffCoefryESthetam[,1],rev(DiffCoefryESthetam[,2])), 
        border = gray(.2), lwd = .5, col = "#BBBBBB40")
polygon(c(yearsES,rev(yearsES)), c(DiffCoefryESthetaf[,1],rev(DiffCoefryESthetaf[,2])), 
        border = gray(.2), lwd = .5, col = "#BBBBBB40")
# Lotka same thing
#polygon(c(yearsUS,rev(yearsUS)), c(DiffCoefryUSthetam[,3],rev(DiffCoefryUSthetam[,4])), 
#        border = "blue", lwd = .5, col = "#BBBBBB40")
#polygon(c(yearsUS,rev(yearsUS)), c(DiffCoefryUSthetaf[,3],rev(DiffCoefryUSthetaf[,4])), 
#        border = "blue", lwd = .5, col = "#BBBBBB40")
#
#polygon(c(yearsES,rev(yearsES)), c(DiffCoefryESthetam[,3],rev(DiffCoefryESthetam[,4])), 
#        border = "red", lwd = .5, col = "#BBBBBB40")
#polygon(c(yearsES,rev(yearsES)), c(DiffCoefryESthetaf[,3],rev(DiffCoefryESthetaf[,4])), 
#        border = "red", lwd = .5, col = "#BBBBBB40")

# change text
text(c(1993, 1993,1985,1986), c(0.09198817, 0.04217548, 0.17349984, 0.11114689),
        c(expression(paste(theta," US females")), 
          expression(paste(theta," US males")),
          expression(paste(theta," ES females")), 
          expression(paste(theta," ES males"))), pos = c(4,4,2,4))

segments(c(1993, 1993,1985,1986), c(0.09198817, 0.04217548, 0.17349984, 0.11114689),
        c(1994, 1992, 1986, 1985), c(DiffCoefryUSf[yearsUS == 1994,1],
                DiffCoefryUSm[yearsUS == 1992,1],
                DiffCoefryESf[yearsES == 1986,1],
                DiffCoefryESm[yearsES == 1985,1]))
dev.off()

pdf("latex/Figures/exPyramidthetaratio.pdf", height = 4.5, width = 4.5)
par(mai = c(.5,.5,.5,.3), xaxs = "i", yaxs = "i")
plot(yearsUS, DiffCoefryUSm[,1], type = 'n', ylim = c(1,2.5), xlim = c(1968, 2010), 
        axes = FALSE, xlab = "", ylab = "",
        panel.first = list(rect(1968, 1, 2010, 2.5, col = gray(.95), border = NA),
                abline(v = seq(1970, 2010, by = 5),col = "white"),
                abline(h = seq(1, 2.5, by = .2),col = "white"),
                text(1990,.85, "Year", xpd = TRUE),
                text(1966,2.6, expression(frac(theta~Lotka,theta~e[y])), xpd = TRUE),
                text(seq(1970, 2010, by = 5), 1, seq(1970, 2010, by = 5), pos = 1 ,cex = .7, xpd = TRUE),
                text(1968, seq(1, 2.5, by = .2), seq(1, 2.5, by = .2), pos = 2, cex = .7, xpd = TRUE)
        ))
polygon(c(yearsUS,rev(yearsUS)), c(DiffCoefryUSthetam[,5],rev(DiffCoefryUSthetam[,6])), 
        border = gray(.2), lwd = .5, col = "#BBBBBB40")
polygon(c(yearsUS,rev(yearsUS)), c(DiffCoefryUSthetaf[,5],rev(DiffCoefryUSthetaf[,6])), 
        border = gray(.2), lwd = .5, col = "#BBBBBB40")

polygon(c(yearsES,rev(yearsES)), c(DiffCoefryESthetam[,5],rev(DiffCoefryESthetam[,6])), 
        border = gray(.2), lwd = .5, col = "#BBBBBB40", lty = 4)
polygon(c(yearsES,rev(yearsES)), c(DiffCoefryESthetaf[,5],rev(DiffCoefryESthetaf[,6])), 
        border = gray(.2), lwd = .5, col = "#BBBBBB40", lty = 4)

# indicate
text(c(1985, 1980, 1990.5, 1990), c( 2, 2.267286, 1.2, 1.486132),
        c("US females", "US males", "ES females", "ES males"), pos = c(1,2,4,2))

segments(c(1982, 1980, 1991, 1990), c(2, 2.267286, 1.245777, 1.486132),
        c(1981, 1981, 1990, 1991), c(DiffCoefryUSf[yearsUS == 1981,2]/DiffCoefryUSf[yearsUS == 1981,1],
                DiffCoefryUSm[yearsUS == 1981,2]/DiffCoefryUSm[yearsUS == 1981,1],
                DiffCoefryESf[yearsES == 1990,2]/DiffCoefryESf[yearsES == 1990,1],
                DiffCoefryESm[yearsES == 1991,2]/DiffCoefryESm[yearsES == 1991,1]))

dev.off()
#
#plot(NULL, type = "n", xlim = c(1968,2010) , ylim = c(0,.3))
#lines(yearsUS, DiffCoefryUSm[,1], col = gray(.2), lwd = 1, lty = 1)                
#lines(yearsUS, DiffCoefryUSf[,1], col = gray(.2), lwd = 1, lty = 1)
#lines(yearsES, DiffCoefryESm[,1], col = gray(.4), lwd = 1.2, lty = 4)                
#lines(yearsES, DiffCoefryESf[,1], col = gray(.4), lwd = 1.2, lty = 4)
#
#
#lines(yearsUS, DiffCoefryUSm[,2], col = "blue", lwd = 1, lty = 1)                
#lines(yearsUS, DiffCoefryUSf[,2], col = "blue", lwd = 1, lty = 1)
#lines(yearsES, DiffCoefryESm[,2], col = "red", lwd = 1.2, lty = 4)                
#lines(yearsES, DiffCoefryESf[,2], col = "red", lwd = 1.2, lty = 4)
#
#cor(diff(DiffCoefryUSm[,1]),diff(DiffCoefryUSm[,2]))
#cor(diff(DiffCoefryUSf[,1]),diff(DiffCoefryUSf[,2]))
#cor(diff(DiffCoefryESm[,1]),diff(DiffCoefryESm[,2]))
#cor(diff(DiffCoefryESf[,1]),diff(DiffCoefryESf[,2]))

# --------------------------------------------------------------------------------------------------
lapply(as.character(yearsUS), function(yr, .Lxm ,.rmLUS))
rmUS
names(rmLUS) <- yearsUS
yr <- "1975"
TmLUS <- unlist(lapply(as.character(yearsUS), function(yr, .Bx, .Ex, .Lx, Bxymf, MF, rc,.rmLUS){       
                    Fx  <- Minf0(Mna0(rc(.Bx[[yr]][[Bxymf]], na.rm = TRUE) / .Ex[.Ex$Year == as.integer(yr),MF]))
                    sum((.5:110.5)^2*Fx*.Lx[,yr]*exp(-.rmLUS[yr]*.5:110.5))/
                    sum((.5:110.5)*Fx*.Lx[,yr]*exp(-.rmLUS[yr]*.5:110.5))
                }, .Bx = BxymfUS, .Ex = ExUS, .Lx = LxmUS, Bxymf = "Bxym", 
                MF = "Male", rc = rowSums, .rmLUS = rmLUS))
names(TmLUS) <- yearsUS

plot(exp(rmLUS*TmLUS), ylim = c(.7,1.5))
lines(exp(rmUS[,1]*rmUS[,2]))
rmLUS - rmUS[,1]


# ----------------------------------------------------
# check solution unique and concave:
yr <- "1975"
Eym <- rowSums(ExpectedDx(with(ExUS,Male[ExUS$Year == as.integer(yr)]), dxmUS[, yr]))

Bym <-  rowSums(ExpectedDx(rowSums(BxymfUS[[yr]][["Bxym"]], na.rm = TRUE), dxmUS[, yr]))
Fym <- Minf0(Mna0(Bym / Eym))
rvals <- seq(-.05,.1,by=.001)
resids <- sapply(rvals,exOneSexMin, dx = dxmUS[, yr], Fex = Fym)
plot(rvals,resids,type = 'l', panel.first=list(abline(h=0,col="red")))




