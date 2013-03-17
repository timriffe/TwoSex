# 1) make function to estimate lotka's R using basic inputs of dx, SRB and Fx
# 2) decompose difference between male and female r using DecompHoriuchi()
# 3) plot in 2 ways: 1) proportion due to each factor over time 2) total from each factor over time.
Mna0 <- function(M){
    M[is.na(M)]  <- 0
    M[is.nan(M)] <- 0
    M
}

Minf0 <- function(M){
    M[is.infinite(M)]  <- 0
    M
}
#source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")
yearsUS <- 1969:2009
yearsES <- 1975:2009

BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))
names(BxymfES) <- yearsES
# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
dxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfES.Rdata"))) 

dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)

mxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxmUS.Rdata"))) 
mxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxfUS.Rdata"))) 
mxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxmES.Rdata"))) 
mxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxfES.Rdata"))) 


mx2LxHMD <- compiler::cmpfun(function(mx){
    mx                  <- Mna0(as.numeric(mx))
    
    # mean proportion of interval passed at death
    ax                  <- mx * 0 + .5                      # ax = .5, pg 38 MPv5
    
    ax[1]   <- ((0.045 + 2.684 * mx[1]) + (0.053 + 2.800 * mx[1])) / 2 # hack, hard to pass in sex variable

    qx                  <- mx / (1 + (1 - ax) * mx)          # Eq 60 MPv5 (identity)
# ---------------------------------------------------------------------------------
# set open age qx to 1
    i.openage           <- 111 # removed argument OPENAGE
    qx[i.openage]       <- ifelse(is.na(qx[i.openage]), NA, 1)
    ax[i.openage]       <- 1 / mx[i.openage]                   
# ---------------------------------------------------------------------------------
# define remaining lifetable columns:
    px                  <- 1 - qx                                                                                 # Eq 64 MPv5
    px[is.nan(px)]      <- 0 # skips BEL NAs, as these are distinct from NaNs
# lx needs to be done columnwise over px, argument 2 refers to the margin.
    lx                  <- c(1, cumprod(px[1:(i.openage-1)]))
    # NA should only be possible if there was a death with no Exp below age 80- impossible, but just to be sure
    # lx[is.na(lx)]   <- 0 # removed for BEL testing        
    dx                  <- lx * qx                                                                                # Eq 66 MPv5
    Lx                  <- lx - (1 - ax) * dx                                                         # Eq 67 MPv5
    Lx[i.openage]       <- lx[i.openage] * ax[i.openage]
    Lx
})

LotkardxFxSRB1 <- compiler::cmpfun(function(rates, .a = .5:110.5, T.guess = 30){
            N       <- length(.a)
           
            fx      <- rates[1:N]
            sig     <- rates[(N+1):(2*N)]
            dx      <- rates[(2*N+1):(3*N)]
            
            lx      <- rev(cumsum(rev(dx)))         # identity
            Lx      <- (lx[1:(N-1)] + lx[2:N]) / 2  # hack 1
            Lx[1]   <- Lx[1] - dx[1] / 5            # hack 2
            Lx[N]   <- lx[N] / 2                    # hack 3
            # from Coale, Ansley J. (1957) A New Method for Calculating Lotka's r- the Intrinsic Rate of Growth in a Stable Population.
            # Population Studies, Vol. 11 no. 1, pp 92-94
            R0 <- sum(sig * fx * Lx)
            # first assuming a mean generation time of 29
            ri <- log(R0)/30
            
            for (i in 1:15){ # 10 is more than enough!
                deltai <- sum(exp(-ri * .a) * sig * fx * Lx) - 1
                # the mean generation time self-corrects 
                # according to the error produced by the Lotka equation
                ri <- ri + (deltai / (T.guess - (deltai / ri)))
            }
            return(ri)  
        })

LotkarmxFxSRB1 <- compiler::cmpfun(function(rates, .a = .5:110.5, T.guess = 30){
            N       <- length(.a)
            
            fx      <- rates[1:N]
            sig     <- rates[(N+1):(2*N)]
            mx      <- rates[(2*N+1):(3*N)]
            
            Lx      <- mx2LxHMD(mx = mx)
            # from Coale, Ansley J. (1957) A New Method for Calculating Lotka's r- the Intrinsic Rate of Growth in a Stable Population.
            # Population Studies, Vol. 11 no. 1, pp 92-94
            R0 <- sum(sig * fx * Lx)
            # first assuming a mean generation time of 29
            ri <- log(R0)/30
            
            for (i in 1:15){ # 10 is more than enough!
                deltai <- sum(exp(-ri * .a) * sig * fx * Lx) - 1
                # the mean generation time self-corrects 
                # according to the error produced by the Lotka equation
                ri <- ri + (deltai / (T.guess - (deltai / ri)))
            }
            return(ri)  
        })
library(DecompHoriuchi)
library(parallel)
#USdecompR <- do.call(rbind, lapply(as.character(yearsUS), function(yr, .Bxymf, .Ex, .dxm, .dxf){
#            
#            BTm  <- (rowSums(.Bxymf[[yr]][["Bxym"]]) + rowSums(.Bxymf[[yr]][["Bxyf"]]))
#            .Fxm.  <- BTm / .Ex$Male[.Ex$Year == as.integer(yr)]
#            .sigm. <- Mna0(rowSums(.Bxymf[[yr]][["Bxym"]]) / BTm)
#            .dxm.  <- .dxm[,yr]
#          
#            BTf  <- (colSums(.Bxymf[[yr]][["Bxym"]]) + colSums(.Bxymf[[yr]][["Bxyf"]]))
#            .Fxf.  <- BTf / .Ex$Female[.Ex$Year == as.integer(yr)]
#            .sigf. <- Mna0(colSums(.Bxymf[[yr]][["Bxyf"]]) / BTf)
#            .dxf.  <- .dxf[,yr]
#            
#            Dec <- DecompContinuousOrig(func = LotkardxFxSRB1, 
#                    rates2 = c(.Fxm., .sigm., .dxm.), 
#                    rates1 = c(.Fxf., .sigf., .dxf.), N = 100)
#            c(Fert = sum(Dec[1:N]), SRB = sum(Dec[(N+1):(2*N)]), Mort = sum(Dec[(2*N+1):(3*N)]))
#        }, .Bxymf = BxymfUS, .Ex = ExUS, .dxm = dxmUS, .dxf = dxfUS))
#ESdecompR <- do.call(rbind, lapply(as.character(yearsES), function(yr, .Bxymf, .Ex, .dxm, .dxf){
#                    BTm  <- (rowSums(.Bxymf[[yr]][["Bxym"]]) + rowSums(.Bxymf[[yr]][["Bxyf"]]))
#                    .Fxm.  <-  Minf0(Mna0(BTm / .Ex$Male[.Ex$Year == as.integer(yr)]))
#                    .sigm. <- Mna0(rowSums(.Bxymf[[yr]][["Bxym"]]) / BTm)
#                    .dxm.  <- .dxm[,yr]
#                    
#                    BTf  <- (colSums(.Bxymf[[yr]][["Bxym"]]) + colSums(.Bxymf[[yr]][["Bxyf"]]))
#                    .Fxf.  <-  Minf0(Mna0(BTf / .Ex$Female[.Ex$Year == as.integer(yr)]))
#                    .sigf. <- Mna0(colSums(.Bxymf[[yr]][["Bxyf"]]) / BTf)
#                    .dxf.  <- .dxf[,yr]
#                    
#                    Dec <- DecompContinuousOrig(func = LotkardxFxSRB1, 
#                            rates2 = c(.Fxm., .sigm., .dxm.), 
#                            rates1 = c(.Fxf., .sigf., .dxf.), N = 100)
#                    c(Fert = sum(Dec[1:N]), SRB = sum(Dec[(N+1):(2*N)]), Mort = sum(Dec[(2*N+1):(3*N)]))
#                }, .Bxymf = BxymfES, .Ex = ExES, .dxm = dxmES, .dxf = dxfES))
#save(USdecompR, file = "/home/triffe/git/DISS/Data/rDecompResults/USdecompR.Rdata")
#save(ESdecompR, file = "/home/triffe/git/DISS/Data/rDecompResults/ESdecompR.Rdata")

# these two decompositions were done on the WORLDFAM server, using N = 500, since for some reason the 
# residual error was high in a few years. It should now be negligible.

#USdecompmxR <- do.call(rbind, parallel::mclapply(as.character(yearsUS), function(yr, .Bxymf, .Ex, .mxm, .mxf){
#            N <- 111
#            BTm  <- (rowSums(.Bxymf[[yr]][["Bxym"]]) + rowSums(.Bxymf[[yr]][["Bxyf"]]))
#            .Fxm.  <- BTm / .Ex$Male[.Ex$Year == as.integer(yr)]
#            .sigm. <- Mna0(rowSums(.Bxymf[[yr]][["Bxym"]]) / BTm)
#            .mxm.  <- .mxm[,yr]
#          
#            BTf  <- (colSums(.Bxymf[[yr]][["Bxym"]]) + colSums(.Bxymf[[yr]][["Bxyf"]]))
#            .Fxf.  <- BTf / .Ex$Female[.Ex$Year == as.integer(yr)]
#            .sigf. <- Mna0(colSums(.Bxymf[[yr]][["Bxyf"]]) / BTf)
#            .mxf.  <- .mxf[,yr]
#            
#            Dec <- DecompContinuousOrig(func = LotkarmxFxSRB1, 
#                    rates2 = c(.Fxm., .sigm., .mxm.), 
#                    rates1 = c(.Fxf., .sigf., .mxf.), N = 200)
#            c(Fert = sum(Dec[1:N]), SRB = sum(Dec[(N+1):(2*N)]), Mort = sum(Dec[(2*N+1):(3*N)]))
#        }, .Bxymf = BxymfUS, .Ex = ExUS, .mxm = mxmUS, .mxf = mxfUS))
#ESdecompmxR <- do.call(rbind, parallel::mclapply(as.character(yearsES), function(yr, .Bxymf, .Ex, .mxm, .mxf){
#                    N <- 111
#                    BTm  <- (rowSums(.Bxymf[[yr]][["Bxym"]]) + rowSums(.Bxymf[[yr]][["Bxyf"]]))
#                    .Fxm.  <- BTm / .Ex$Male[.Ex$Year == as.integer(yr)]
#                    .sigm. <- Mna0(rowSums(.Bxymf[[yr]][["Bxym"]]) / BTm)
#                    .mxm.  <- .mxm[,yr]
#                    
#                    BTf  <- (colSums(.Bxymf[[yr]][["Bxym"]]) + colSums(.Bxymf[[yr]][["Bxyf"]]))
#                    .Fxf.  <- BTf / .Ex$Female[.Ex$Year == as.integer(yr)]
#                    .sigf. <- Mna0(colSums(.Bxymf[[yr]][["Bxyf"]]) / BTf)
#                    .mxf.  <- .mxf[,yr]
#                    
#                    Dec <- DecompContinuousOrig(func = LotkarmxFxSRB1, 
#                            rates2 = c(.Fxm., .sigm., .mxm.), 
#                            rates1 = c(.Fxf., .sigf., .mxf.), N = 200)
#                    c(Fert = sum(Dec[1:N]), SRB = sum(Dec[(N+1):(2*N)]), Mort = sum(Dec[(2*N+1):(3*N)]))
#                }, .Bxymf = BxymfES, .Ex = ExES, .mxm = mxmES, .mxf = mxfES))
#USdecompR <- local(get(load("/home/triffe/git/DISS/Data/rDecompResults/USdecompR.Rdata")))
#ESdecompR <- local(get(load("/home/triffe/git/DISS/Data/rDecompResults/ESdecompR.Rdata")))
USdecompR <- local(get(load("/home/triffe/git/DISS/Data/rDecompResults/USdecompmxR.Rdata")))
ESdecompR <- local(get(load("/home/triffe/git/DISS/Data/rDecompResults/ESdecompmxR.Rdata")))
# determine axes compatible with output from both countries
Neg <- USdecompmxR
Neg[Neg > 0] <- 0
Pos <- USdecompmxR
Pos[Pos < 0] <- 0
2010-1968

pdf("/home/triffe/git/DISS/latex/Figures/DecomprUS.pdf", height = 5, width = 5)
par(mai = c(.5,.5,.5,.2), xaxs = "i", yaxs = "i")
plot(NULL, type = "n", xlab = "", ylab = "", xlim = c(-1, 42), ylim = c(-.0025,.006), 
        axes = FALSE)
rect(-1,-.0025,42,.006,col = gray(.95),border = NA)
abline(h = seq(-.002,.006,by=.001), col = "white")
abline(v = seq(1,41,by=5), col = "white")
barplot(t(Neg), add = TRUE, space = 0, border = NA, col = paste0(gray(c(.6,.4,.2)),"BB"),width = 1,axes = FALSE)
barplot(t(Pos), add = TRUE, space = 0, border = NA, col = paste0(gray(c(.6,.4,.2)),"BB"),width = 1,axes = FALSE)

text(seq(1,41,by=5),-.0025,seq(1970,2010,by=5),pos=1,cex=.8,xpd=TRUE)
text(-1,seq(-.002,.006,by=.001),seq(-.002,.006,by=.001),cex=.8,pos=2,xpd=TRUE)
text(.5,c(-.0004,.0015,.0042), c("Mortality","Fertility","SRB"), cex = 1.5,pos = 4, col = "white")
text(20,-.0032,"Year",xpd=TRUE)
text(-2,.0065,"Contribution\nto difference in r", pos = 4, xpd = TRUE)
dev.off()

# Spain
Neg <- ESdecompR
Neg[Neg > 0] <- 0
Pos <- ESdecompR
Pos[Pos < 0] <- 0

pdf("/home/triffe/git/DISS/latex/Figures/DecomprES.pdf", height = 5, width = 5)
par(mai = c(.5,.5,.5,.2), xaxs = "i", yaxs = "i")
plot(NULL, type = "n", xlab = "", ylab = "", xlim = c(-7, 36), ylim = c(-.0025,.006), 
        axes = FALSE)
rect(-7,-.0025,36,.006,col = gray(.95),border = NA)
abline(h = seq(-.002,.006,by=.001), col = "white")
abline(v = seq(-5,35,by=5), col = "white")
barplot(t(Neg), add = TRUE, space = 0, border = NA, col = paste0(gray(c(.6,.4,.2)),"BB"),width = 1,axes = FALSE)
barplot(t(Pos), add = TRUE, space = 0, border = NA, col = paste0(gray(c(.6,.4,.2)),"BB"),width = 1,axes = FALSE)

text(seq(-5,35,by=5),-.0025,seq(1970,2010,by=5),pos=1,cex=.8,xpd=TRUE)
text(-7,seq(-.002,.006,by=.001),seq(-.002,.006,by=.001),cex=.8,pos=2,xpd=TRUE)
text(10,c(-0.000275,0.001,0.003), c("Mortality","Fertility","SRB"), cex = 1.5,pos = 4, col = "white")

text(15,-.0032,"Year",xpd=TRUE)
text(-9,.0065,"Contribution\nto difference in r", pos = 4, xpd = TRUE)
dev.off()


plot(rowSums(abs(ESdecompR)))
plot(rowSums(abs(USdecompR)))

