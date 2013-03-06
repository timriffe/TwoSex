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
#------------------------------------------------------------

# Function Definitions:
#------------------------------------------------------------
# minimizer function for 1 sex ex-perspective renewal function:
# use with: optimize()
exOneSexMin <- function(r, dx, Fex, .a = .5:110.5){
    # get the overlapped / staggered dx structure
    N               <- length(Fex)
    dxM  <- matrix(0, ncol = N, nrow = N)
    # remaining years go down rows. ages over columns
    dxi  <- dx
    for (i in 1:N){
        dxM[i, 1:length(dxi)  ] <- dxi 
        dxi <- dxi[2:length(dxi) ]
    }     
    (1 - sum(rowSums(dxM %col% (1 / exp(-r * .a))) * Fex)) ^ 2
}

exOneSexCoaleR <- function(Fex, dx, .a = .5:110.5, maxit = 1e2, tol = 1e-15, r.start = 0.001){  
    
    N    <- length(Fex)
    dxM  <- matrix(0, ncol = N, nrow = N)
    dxi  <- dx
    for (i in 1:N){
        dxM[i, 1:length(dxi)  ] <- dxi 
        dxi <- dxi[2:length(dxi) ]
    }     
    R0      <- sum(dxM*Fex)
    T.guess <- wmean(.a,rowSums(dxM)*Fex) # assuming r = 0
    r2      <- log(R0) / T.guess
   
    # be careful to discount Fex by SRB appropriately for males / females
    # prior to specification
    # Based on Coale (1957)
    for (i in 1:maxit){ # 15 is more than enough!
       #cat(r2,i,"\n")
        r1 <- r2
        deltai <- 1 - sum(rowSums(dxM %col% (1 / exp(-r1 * .a))) * Fex)
        # the mean generation time self-corrects 
        # according to the error produced by the Lotka equation
        r2 <- r1 - (deltai / (T.guess - (deltai / r1))) 
        if (abs(r2 - r1) <= tol | zapsmall(abs(deltai)) <= tol){
            break
        }
    }
    return(r2)  
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
    b <- 1 / sum(rowSums(dxM %col% (1 / exp(-r * .a))))
    b * rowSums(dxM %col% (1 / exp(-r * .a)))
}

exOneSexTy <- function(r, Fex, dx, .a = .5:110.5){
    N    <- length(Fex)
    dxM  <- matrix(0, ncol = N, nrow = N)
    dxi  <- dx
    for (i in 1:N){
        dxM[i, 1:length(dxi)  ] <- dxi 
        dxi <- dxi[2:length(dxi) ]
    }     
    wmean(.a, rowSums(dxM %col% (1 / exp(-r * .a))) * Fex)
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

#--------------------------------------------------------------
# Tables for appendix 'zAppendix.exSingleSexLotka.tex'
exRepUSm <- cbind(rmUS, R0 = exp(rmUS[,1]*rmUS[,2]) )
exRepUSf <- cbind(rfUS, R0 = exp(rfUS[,1]*rfUS[,2]) )
exRepESm <- cbind(rmES, R0 = exp(rmES[,1]*rmES[,2]) )
exRepESf <- cbind(rfES, R0 = exp(rfES[,1]*rfES[,2]) )
colnames(exRepESf) <-colnames(exRepESm) <-colnames(exRepUSf) <-colnames(exRepUSm) <- c("$r$","$T^y$","$R_0$")
rownames(exRepUSm) <- rownames(exRepUSf) <- yearsUS
rownames(exRepESm) <- rownames(exRepESf) <- yearsES

print(xtable(exRepUSm, digits = c(0,4,2,3), align = c("c","c","c","c")),
        sanitize.colnames.function = identity, 
        file = "/home/triffe/git/DISS/latex/xtables/exRepUSm.tex",floating=FALSE)
print(xtable(exRepUSf, digits = c(0,4,2,3), align = c("c","c","c","c")),
        sanitize.colnames.function = identity, 
        file = "/home/triffe/git/DISS/latex/xtables/exRepUSf.tex",floating=FALSE)
print(xtable(exRepESm, digits = c(0,4,2,3), align = c("c","c","c","c")),
        sanitize.colnames.function = identity, 
        file = "/home/triffe/git/DISS/latex/xtables/exRepESm.tex",floating=FALSE)
print(xtable(exRepESf, digits = c(0,4,2,3), align = c("c","c","c","c")),
        sanitize.colnames.function = identity, 
        file = "/home/triffe/git/DISS/latex/xtables/exRepESf.tex",floating=FALSE)
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


plot(yearsUS, exRepUSm[,2], type = 'l', col = "blue", ylim = c(40,60))
lines(yearsUS, exRepUSm[,2], col ="red")
lines(yearsES, exRepUSm[,2], col ="red", lty = 2)
lines(yearsES, exRepUSm[,2], col ="blue", lty = 2)

plot(yearsUS,R0mUS, type = 'l', ylim = c(.8,1.4))
lines(yearsUS, exp(rfUS[,1]*rfUS[,2])  , type = 'l')

plot(yearsUS, exRepUSm[,1], type = 'l', col = "blue", ylim = c(-.02,.015))
lines(yearsUS, exRepUSf[,1], col ="red")
lines(yearsES, exRepESf[,1], col ="red", lty = 2)
lines(yearsES, exRepESm[,1], col ="blue", lty = 2)
abline(h=0)
lines(yearsUS, rmLUS, col = "royalblue", lwd = 2)
lines(yearsUS, rfLUS, col ="pink", lwd = 2)
lines(yearsES, rfLES, col ="pink", lty = 2, lwd = 2)
lines(yearsES, rmLES, col ="royalblue", lty = 2, lwd = 2)
abline(h=0)

# -----------------------------------------------------
# Figure
pdf("/home/triffe/git/DISS/latex/Figures/exLotka1sex.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rmUS[, 1], type = 'n', ylim = c(-.016,.01),xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,-.016,2010,.01,col = gray(.95), border=NA),
                abline(h = seq(-.016,.01,by = .002), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.016,.01,by = .002),seq(-.016,.01,by = .002), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.016, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.0175, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,.0115, "r", cex = 1, xpd = TRUE)))

lines(yearsUS, rmUS[, 1],col = gray(.2), lwd = 2)
lines(yearsUS, rfUS[, 1],col = gray(.5), lwd = 2.5)
lines(yearsES, rmES[, 1],col = gray(.2), lwd = 2, lty = 5)
lines(yearsES, rfES[, 1],col = gray(.5), lwd = 2.5, lty = 5)

legend(1970,-.01, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
dev.off()
