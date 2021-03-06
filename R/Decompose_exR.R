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

mxmUS <- local(get(load("Data/HMD_mux/muxmUS.Rdata"))) 
mxfUS <- local(get(load("Data/HMD_mux/muxfUS.Rdata"))) 
mxmES <- local(get(load("Data/HMD_mux/muxmES.Rdata"))) 
mxfES <- local(get(load("Data/HMD_mux/muxfES.Rdata"))) 


# update function to determine sex gap
# should run on mx and not dx:

exOneSexCoaleRdec1 <- compiler::cmpfun(function(rates, .a = .5:110.5, maxit = 2e2, tol = 1e-11, 
                mx2dxHMD, Mna0, wmean){  
    N        <- length(.a)
    # extract rates
    Fex      <- rates[1:N]
    sig      <- rates[(N+1):(2*N)]
    mx       <- rates[(2*N+1):(3*N)]
    
    dx       <- mx2dxHMD(mx, Mna0)
    
    dxM      <- matrix(0, ncol = N, nrow = N)
    dxi      <- dx
    for (i in 1:N){
        dxM[i, 1:length(dxi)  ] <- dxi 
        dxi  <- dxi[2:length(dxi) ]
    }     
    R0       <- sum(dxM*(sig*Fex))
    T.guess  <- wmean(.a,rowSums(dxM)*(sig*Fex)) # assuming r = 0
    r2       <- log(R0) / T.guess
    
    # be careful to discount Fex by SRB appropriately for males / females
    # prior to specification
    # Based on Coale (1957)
    for (i in 1:maxit){ # 15 is more than enough!
        #cat(r2,i,"\n")
        r1     <- r2
        deltai <- 1 - sum(rowSums(t(t(dxM) / (1 / exp(-r1 * .a)))) * (sig*Fex))
        # the mean generation time self-corrects 
        # according to the error produced by the Lotka equation
        r2     <- r1 - (deltai / (T.guess - (deltai / r1))) 
        if (abs(r2 - r1) <= tol | zapsmall(abs(deltai)) <= tol){
            break
        }
    }
    r2
})

# this one takes shape into account separately:
exOneSexCoaleRdec2 <- compiler::cmpfun(function(rates, .a = .5:110.5, maxit = 2e2, tol = 1e-11){  
            N        <- length(.a)
            # extract rates
            Fexex      <- rates[1:N]
            sig      <- rates[(N+1):(2*N)]
            mx       <- rates[(2*N+1):(3*N)]
            eTFR     <- rates[length(rates)]
            
            dx       <- mx2dxHMD(mx)
            
            dxM      <- matrix(0, ncol = N, nrow = N)
            dxi      <- dx
            for (i in 1:N){
                dxM[i, 1:length(dxi)  ] <- dxi 
                dxi  <- dxi[2:length(dxi) ]
            }     
            R0       <- sum(dxM*(sig*Fexex*eTFR))
            T.guess  <- wmean(.a,rowSums(dxM)*(sig*Fexex*eTFR)) # assuming r = 0
            r2       <- log(R0) / T.guess
            
            # be careful to discount Fex by SRB appropriately for males / females
            # prior to specification
            # Based on Coale (1957)
            for (i in 1:maxit){ # 15 is more than enough!
                #cat(r2,i,"\n")
                r1     <- r2
                deltai <- 1 - sum(rowSums(t(t(dxM) / (1 / exp(-r1 * .a)))) * (sig*Fexex*eTFR))
                # the mean generation time self-corrects 
                # according to the error produced by the Lotka equation
                r2     <- r1 - (deltai / (T.guess - (deltai / r1))) 
                if (abs(r2 - r1) <= tol | zapsmall(abs(deltai)) <= tol){
                    break
                }
            }
            r2
        })




# all functions needed to be passed explicitly inside because parLapply uses snow framework,
# which fires up r sessions that need to be fed all relevant materials.
# executed on WORLDFAM server. Takes several hours to run, but produces small error. 
#cl <- makeCluster(4)
#USdecompExR <- do.call(rbind, parLapply(cl, as.character(yearsUS), 
#				function(yr, .Bxymf, .Ex, .mxm, .mxf,  mx2dxHMD, Mna0, Minf0, wmean, ExpectedDx, exOneSexCoaleRdec1, DecompContinuousOrig){
#					N <- 111
#					
#					# 1) get mx for year
#					.mxf.     <- .mxf[,yr]
#					.mxm.     <- .mxm[,yr]
#					
#					dxf1      <- mx2dxHMD(.mxf.,Mna0)
#					dxm1      <- mx2dxHMD(.mxm., Mna0)
#					
#                   # these differ because dxm used to distribute girl births: compared to male stable exposure: not redistributed
#                   # by female dx. Vice versa below
#					BexMM     <- rowSums(ExpectedDx(rowSums(.Bxymf[[yr]][["Bxym"]]), dxm1)) 
#					BexMF     <- rowSums(ExpectedDx(rowSums(.Bxymf[[yr]][["Bxyf"]]), dxm1))
#					.Fexm.    <- Mna0(Minf0((BexMM+BexMF) /  rowSums(ExpectedDx(.Ex$Male[.Ex$Year == as.integer(yr)], dxm1))))
#					.sigexm.  <- Mna0(Minf0(BexMM / (BexMM + BexMF)))
#					
#					
#					BexFM     <- rowSums(ExpectedDx(colSums(.Bxymf[[yr]][["Bxym"]]), dxf1)) 
#					BexFF     <- rowSums(ExpectedDx(colSums(.Bxymf[[yr]][["Bxyf"]]), dxf1))
#					.Fexf.    <- Mna0(Minf0((BexFM+BexFF) /  rowSums(ExpectedDx(.Ex$Female[.Ex$Year == as.integer(yr)], dxf1))))
#					.sigexf.  <- Mna0(Minf0(BexFF / (BexFF + BexFM)))
#					rates2    <- c(.Fexm., .sigexm., .mxm.)
#					rates1    <- c(.Fexf., .sigexf., .mxf.)
#					
#					Dec <- DecompContinuousOrig(func = exOneSexCoaleRdec1, 
#							rates2 = rates2, 
#							rates1 = rates1, N = 300, 
#							mx2dxHMD = mx2dxHMD, Mna0 = Mna0, wmean = wmean)
#					c(Fert = sum(Dec[1:N]), SRB = sum(Dec[(N+1):(2*N)]), Mort = sum(Dec[(2*N+1):(3*N)]))
#				}, .Bxymf = BxymfUS, .Ex = ExUS, .mxm = mxmUS, .mxf = mxfUS, 
#				mx2dxHMD = mx2dxHMD, Mna0 = Mna0, Minf0 = Minf0, wmean = wmean, 
#				ExpectedDx = ExpectedDx, DecompContinuousOrig = DecompContinuousOrig, 
#				exOneSexCoaleRdec1 = exOneSexCoaleRdec1))
#stopCluster(cl)
##
#cl <- makeCluster(4)
#ESdecompExR <- do.call(rbind, parLapply(cl, as.character(yearsES), 
#				function(yr, .Bxymf, .Ex, .mxm, .mxf,  mx2dxHMD, Mna0, Minf0, wmean, ExpectedDx, exOneSexCoaleRdec1, DecompContinuousOrig){
#					N <- 111
#					
#					# 1) get mx for year
#					.mxf.     <- .mxf[,yr]
#					.mxm.     <- .mxm[,yr]
#					
#					dxf1      <- mx2dxHMD(.mxf.,Mna0)
#					dxm1      <- mx2dxHMD(.mxm., Mna0)
#					
#					BexMM     <- rowSums(ExpectedDx(rowSums(.Bxymf[[yr]][["Bxym"]]), dxm1)) 
#					BexMF     <- rowSums(ExpectedDx(rowSums(.Bxymf[[yr]][["Bxyf"]]), dxm1))
#					.Fexm.    <- Mna0(Minf0((BexMM+BexMF) /  rowSums(ExpectedDx(.Ex$Male[.Ex$Year == as.integer(yr)], dxm1))))
#					.sigexm.  <- Mna0(Minf0(BexMM / (BexMM + BexMF)))
##					
#					
#					BexFM     <- rowSums(ExpectedDx(colSums(.Bxymf[[yr]][["Bxym"]]), dxf1)) 
#					BexFF     <- rowSums(ExpectedDx(colSums(.Bxymf[[yr]][["Bxyf"]]), dxf1))
#					.Fexf.    <- Mna0(Minf0((BexFM+BexFF) /  rowSums(ExpectedDx(.Ex$Female[.Ex$Year == as.integer(yr)], dxf1))))
#					.sigexf.  <- Mna0(Minf0(BexFF / (BexFF + BexFM)))
#					rates2    <- c(.Fexm., .sigexm., .mxm.)
#					rates1    <- c(.Fexf., .sigexf., .mxf.)
#					
#					Dec <- DecompContinuousOrig(func = exOneSexCoaleRdec1, 
#							rates2 = rates2, 
#							rates1 = rates1, N = 300, 
#							mx2dxHMD = mx2dxHMD, Mna0 = Mna0, wmean = wmean)
#					c(Fert = sum(Dec[1:N]), SRB = sum(Dec[(N+1):(2*N)]), Mort = sum(Dec[(2*N+1):(3*N)]))
#				}, .Bxymf = BxymfES, .Ex = ExES, .mxm = mxmES, .mxf = mxfES, 
#				mx2dxHMD = mx2dxHMD, Mna0 = Mna0, Minf0 = Minf0, wmean = wmean, 
#				ExpectedDx = ExpectedDx, DecompContinuousOrig = DecompContinuousOrig, 
#				exOneSexCoaleRdec1 = exOneSexCoaleRdec1))
#stopCluster(cl)

USdecompR <- local(get(load("Data/rDecompResults/USdecompExR.Rdata")))
ESdecompR <- local(get(load("Data/rDecompResults/ESdecompExR.Rdata")))
# determine axes compatible with output from both countries
Neg <- USdecompR
Neg[Neg > 0] <- 0
Pos <- USdecompR
Pos[Pos < 0] <- 0

pdf("latex/Figures/DecomprExUS.pdf", height = 5, width = 5)
par(mai = c(.5,.5,.5,.2), xaxs = "i", yaxs = "i")
plot(NULL, type = "n", xlab = "", ylab = "", xlim = c(-1, 42), ylim = c(-.004,.007), 
        axes = FALSE)
rect(-1,-.004,42,.007,col = gray(.95),border = NA)
abline(h = seq(-.004,.007,by=.001), col = "white")
abline(v = seq(1,41,by=5), col = "white")
barplot(t(Neg), add = TRUE, space = 0, border = NA, col = paste0(gray(c(.6,.4,.2)),"BB"),width = 1,axes = FALSE)
barplot(t(Pos), add = TRUE, space = 0, border = NA, col = paste0(gray(c(.6,.4,.2)),"BB"),width = 1,axes = FALSE)

text(seq(1,41,by=5),-.004,seq(1970,2010,by=5),pos=1,cex=.8,xpd=TRUE)
text(-1,seq(-.004,.007,by=.001),seq(-.004,.007,by=.001),cex=.8,pos=2,xpd=TRUE)
text(10,c(-0.0008636214, 0.0004774893, 0.0033583963), c("Mortality","Fertility","SRB"), cex = 1.5,pos = 4, col = "white")

text(20,-.005,"Year",xpd=TRUE)
text(-2,.0075,"Contribution\nto difference in r", pos = 4, xpd = TRUE)
dev.off()

# Spain
Neg <- ESdecompR
Neg[Neg > 0] <- 0
Pos <- ESdecompR
Pos[Pos < 0] <- 0

pdf("latex/Figures/DecomprExES.pdf", height = 5, width = 5)
par(mai = c(.5,.5,.5,.2), xaxs = "i", yaxs = "i")
plot(NULL, type = "n", xlab = "", ylab = "", xlim = c(-7, 36), ylim = c(-.004,.007), 
        axes = FALSE)
rect(-7,-.004,36,.007,col = gray(.95),border = NA)
abline(h = seq(-.004,.007,by=.001), col = "white")
abline(v = seq(-5,35,by=5), col = "white")
barplot(t(Neg), add = TRUE, space = 0, border = NA, col = paste0(gray(c(.6,.4,.2)),"BB"),width = 1,axes = FALSE)
barplot(t(Pos), add = TRUE, space = 0, border = NA, col = paste0(gray(c(.6,.4,.2)),"BB"),width = 1,axes = FALSE)

text(seq(-5,35,by=5),-.004,seq(1970,2010,by=5),pos=1,cex=.8,xpd=TRUE)
text(-7,seq(-.004,.007,by=.001),seq(-.004,.007,by=.001),cex=.8,pos=2,xpd=TRUE)
text(16,c(-0.001009793,  0.001825698,  0.005216793), c("Mortality","Fertility","SRB"), cex = 1.5,pos = 4, col = "white")
text(15,-.005,"Year",xpd=TRUE)
text(-8,.0075,"Contribution\nto difference in r", pos = 4, xpd = TRUE)
dev.off()

plot(rowSums(abs(ESdecompR)))
USdecompRL <- local(get(load("Data/rDecompResults/USdecompmxR.Rdata")))
ESdecompRL <- local(get(load("Data/rDecompResults/ESdecompmxR.Rdata")))

plot(yearsUS, rowSums(abs(USdecompRL)), type = 'l', col = "blue", ylim = c(0,.01))
lines(yearsES, rowSums(abs(ESdecompRL)), col = "red")
lines(yearsUS, rowSums(abs(USdecompR)), col = "blue", lty = 2)
lines(yearsES, rowSums(abs(ESdecompR)), col = "red", lty = 2)


plot(yearsUS, USdecompRL[,1], type = 'l', col = "blue", ylim = c(0,.01))


# confirm that decomposition adds properly: 
#(rmUS[,1] - rfUS[,1]) - rowSums(USdecompR)
#(rmES[,1] - rfES[,1]) - rowSums(ESdecompR)

plot(yearsUS, USdecompR[,1], type = 'l', col = "#11FF33", lwd =2, ylim = c(-.005,.005))
lines(yearsUS, USdecompR[,2], col = "#AABBFF", lwd = 2)
lines(yearsUS, USdecompR[,3], col = "#FF1111", lwd = 2)
#
lines(yearsES, ESdecompR[,1],  col = "#11FF33", lwd =2, lty = 4)
lines(yearsES, ESdecompR[,2], col = "#AABBFF", lwd = 2, lty = 4)
lines(yearsES, ESdecompR[,3], col = "#FF1111", lwd = 2, lty = 4)
#

plot(yearsUS, USdecompR[,"Fert"] + USdecompR[,"Mort"], type = 'l', ylim = c(-.003,.004))
lines(yearsUS, USdecompRL[,"Fert"] + USdecompRL[,"Mort"])
cor(diff(USdecompR[,"Fert"] + USdecompR[,"Mort"]), diff(USdecompRL[,"Fert"]))
cor(USdecompR[,"Fert"], USdecompR[,"Mort"])

###############################################
library(parallel)
library(DecompHoriuchi)
# this code transferred to Coale or Galton, as it takes too long to run...
# a single year takes 20-30 min on my old laptop...
#ESdecompExR2 <- do.call(rbind, mclapply(as.character(yearsES), function(yr, .Bxymf, .Ex, .mxm, .mxf){
#                    N <- 111
#                    
#                    # 1) get mx for year
#                    .mxf.     <- .mxf[,yr]
#                    .mxm.     <- .mxm[,yr]
#                    
#                    dxf1      <- mx2dxHMD(.mxf.)
#                    dxm1      <- mx2dxHMD(.mxm.)
#                    
#                    BexMM     <- rowSums(ExpectedDx(rowSums(.Bxymf[[yr]][["Bxym"]]), dxm1)) 
#                    BexMF     <- rowSums(ExpectedDx(rowSums(.Bxymf[[yr]][["Bxyf"]]), dxm1))
#                    .Fexm.    <- Mna0(Minf0((BexMM+BexMF) /  rowSums(ExpectedDx(.Ex$Male[.Ex$Year == as.integer(yr)], dxm1))))
#                    .sigexm.  <- Mna0(Minf0(BexMM / (BexMM + BexMF)))
#                    eTFRm     <- sum(.Fexm.)
#                   .Fexexm.   <- .Fexm. / eTFRm
#                    
#                    BexFM     <- rowSums(ExpectedDx(colSums(.Bxymf[[yr]][["Bxym"]]), dxf1)) 
#                    BexFF     <- rowSums(ExpectedDx(colSums(.Bxymf[[yr]][["Bxyf"]]), dxf1))
#                    .Fexf.    <- Mna0(Minf0((BexFM+BexFF) /  rowSums(ExpectedDx(.Ex$Female[.Ex$Year == as.integer(yr)], dxf1))))
#                    .sigexf.  <- Mna0(Minf0(BexFF / (BexFF + BexFM)))
#                    eTFRf     <- sum(.Fexf.)
#                    .Fexexf.   <- .Fexf. / eTFRf
#                    
#                    rates2    <- c(.Fexexm., .sigexm., .mxm., eTFRm)
#                    rates1    <- c(.Fexexf., .sigexf., .mxf., eTFRf)
#                    
#                    Dec <- DecompContinuousOrig(func = exOneSexCoaleRdec2, 
#                            rates2 = rates2, 
#                            rates1 = rates1, N = 300)
#            c(Fertpdf = sum(Dec[1:N]), SRB = sum(Dec[(N+1):(2*N)]), Mort = sum(Dec[(2*N+1):(3*N)]), TFR = Dec[length(Dec)])
#                }, .Bxymf = BxymfES, .Ex = ExES, .mxm = mxmES, .mxf = mxfES))

# NOTE: US was computed analagously to the above (swap all instances of 'ES' to 'US' in named variables)
# computed on UCD demog Coale old server. results emailed back to laptop.


USdecompR <- local(get(load("Data/rDecompResults/USdecompExR2.Rdata")))
ESdecompR <- local(get(load("Data/rDecompResults/ESdecompExR2.Rdata")))


head(USdecompR)
# determine axes compatible with output from both countries
Neg <- USdecompR
Neg[Neg > 0] <- 0
Pos <- USdecompR
Pos[Pos < 0] <- 0

pdf("latex/Figures/DecomprExUS2.pdf", height = 5, width = 5)
par(mai = c(.5,.5,.5,.2), xaxs = "i", yaxs = "i")
plot(NULL, type = "n", xlab = "", ylab = "", xlim = c(-1, 42), ylim = c(-.004,.007), 
        axes = FALSE)
rect(-1,-.004,42,.007,col = gray(.95),border = NA)
abline(h = seq(-.004,.007,by=.001), col = "white")
abline(v = seq(1,41,by=5), col = "white")
barplot(t(Neg), add = TRUE, space = 0, border = NA, col = paste0(gray(c(.8,.6,.4,.2)),"BB"),width = 1,axes = FALSE)
barplot(t(Pos), add = TRUE, space = 0, border = NA, col = paste0(gray(c(.8,.6,.4,.2)),"BB"),width = 1,axes = FALSE)

text(seq(1,41,by=5),-.004,seq(1970,2010,by=5),pos=1,cex=.8,xpd=TRUE)
text(-1,seq(-.004,.007,by=.001),seq(-.004,.007,by=.001),cex=.8,pos=2,xpd=TRUE)
text(c(11, 11, 11, 12,32),c(0.000522905,  0.002898587, -0.001163063,  0.005350904,-.0027), 
        c("Fertility shape", "SRB", "Mortality", "eTFR","eTFR"), cex = 1.3,pos = 4, col = c(gray(.2),"white","white",gray(.2),gray(.2)))

segments(12.604584, 0.005159316, 9.953926, 0.004431285)
segments(32.97280, -0.002465856, 31.50796, -0.001833618)

text(20,-.005,"Year",xpd=TRUE)
text(-2,.0075,"Contribution\nto difference in r", pos = 4, xpd = TRUE)
dev.off()

# Spain
Neg <- ESdecompR
Neg[Neg > 0] <- 0
Pos <- ESdecompR
Pos[Pos < 0] <- 0

pdf("latex/Figures/DecomprExES2.pdf", height = 5, width = 5)
par(mai = c(.5,.5,.5,.2), xaxs = "i", yaxs = "i")
plot(NULL, type = "n", xlab = "", ylab = "", xlim = c(-7, 36), ylim = c(-.004,.007), 
        axes = FALSE)
rect(-7,-.004,36,.007,col = gray(.95),border = NA)
abline(h = seq(-.004,.007,by=.001), col = "white")
abline(v = seq(-5,35,by=5), col = "white")
barplot(t(Neg), add = TRUE, space = 0, border = NA, col = paste0(gray(c(.8,.6,.4,.2)),"BB"),width = 1,axes = FALSE)
barplot(t(Pos), add = TRUE, space = 0, border = NA, col = paste0(gray(c(.8,.6,.4,.2)),"BB"),width = 1,axes = FALSE)

text(seq(-5,35,by=5),-.004,seq(1970,2010,by=5),pos=1,cex=.8,xpd=TRUE)
text(-7,seq(-.004,.007,by=.001),seq(-.004,.007,by=.001),cex=.8,pos=2,xpd=TRUE)
text(c(12,12,12,  6, 29),c(0.0009635557,  0.0040289517, -0.0015079198,  0.0058107131, -0.0033), 
        c("Fertility shape","SRB","Mortality","eTFR","eTFR"), cex = 1.3,pos = 4, col = c(gray(.2),"white","white",gray(.2),gray(.2)))
segments(9.185489,0.005580808 ,10.441064,0.005044364)
segments(32.20436,-0.003002300,32.90191,-0.002389221)
text(15,-.005,"Year",xpd=TRUE)
text(-8,.0075,"Contribution\nto difference in r", pos = 4, xpd = TRUE)
dev.off()

# -----------------------------------------
# no point to lapply, but was easy to copy and paste
ES1990 <- lapply("1990", function(yr, .Bxymf, .Ex, .mxm, .mxf){
                    N <- 111
                    
                    # 1) get mx for year
                    .mxf.     <- .mxf[,yr]
                    .mxm.     <- .mxm[,yr]
                    
                    dxf1      <- mx2dxHMD(.mxf.)
                    dxm1      <- mx2dxHMD(.mxm.)
                    
                    BexMM     <- rowSums(ExpectedDx(rowSums(.Bxymf[[yr]][["Bxym"]]), dxm1)) 
                    BexMF     <- rowSums(ExpectedDx(rowSums(.Bxymf[[yr]][["Bxyf"]]), dxm1))
                    .Fexm.    <- Mna0(Minf0((BexMM+BexMF) /  rowSums(ExpectedDx(.Ex$Male[.Ex$Year == as.integer(yr)], dxm1))))
                    .sigexm.  <- Mna0(Minf0(BexMM / (BexMM + BexMF)))
                    eTFRm     <- sum(.Fexm.)
                   .Fexexm.   <- .Fexm. / eTFRm
                    
                    BexFM     <- rowSums(ExpectedDx(colSums(.Bxymf[[yr]][["Bxym"]]), dxf1)) 
                    BexFF     <- rowSums(ExpectedDx(colSums(.Bxymf[[yr]][["Bxyf"]]), dxf1))
                    .Fexf.    <- Mna0(Minf0((BexFM+BexFF) /  rowSums(ExpectedDx(.Ex$Female[.Ex$Year == as.integer(yr)], dxf1))))
                    .sigexf.  <- Mna0(Minf0(BexFF / (BexFF + BexFM)))
                    eTFRf     <- sum(.Fexf.)
                    .Fexexf.   <- .Fexf. / eTFRf
                    
                    rates2    <- c(.Fexexm., .sigexm., .mxm., eTFRm)
                    rates1    <- c(.Fexexf., .sigexf., .mxf., eTFRf)
                    
                    DecompContinuousOrig(func = exOneSexCoaleRdec2, 
                            rates2 = rates2, 
                            rates1 = rates1, N = 300)
          
                }, .Bxymf = BxymfES, .Ex = ExES, .mxm = mxmES, .mxf = mxfES)
ES1990 <- unlist(ES1990)
N <- 111
Fertpdf <- ES1990[1:N]
SRB     <- ES1990[(N+1):(2*N)]
Mort    <- ES1990[(2*N+1):(3*N)]
ES1990[length(ES1990)]
All <- rbind(SRB,Mort,Fertpdf)
Neg <- All
Neg[Neg > 0] <- 0
Pos <- All
Pos[Pos < 0] <- 0

pdf("latex/Figures/DecomprExES1990.pdf", height = 5, width = 5)
par(mai = c(.5,.5,.5,.2), xaxs = "i", yaxs = "i")
plot(NULL, type = "n", xlab = "", ylab = "", xlim = c(-1, 111), ylim = c(-.0005,.0004), 
        axes = FALSE)
rect(-1,-.0005,111,.0004,col = gray(.95),border = NA)
abline(h = seq(-.0005,.0004,by=.0001), col = "white")
abline(v = seq(0,110,by=10), col = "white")
barplot(Neg, add = TRUE, space = 0, border = NA, col = paste0(gray(c(.6,.4,.8)),"BB"),width = 1,axes = FALSE)
barplot(Pos, add = TRUE, space = 0, border = NA, col = paste0(gray(c(.6,.4,.8)),"BB"),width = 1,axes = FALSE)

text(seq(0,110,by=10),-.0005,seq(0,110,by=10),pos=1,cex=.8,xpd=TRUE)
text(-1,seq(-.0005,.0004,by=.0001),seq(-.0005,.0004,by=.0001),cex=.8,pos=2,xpd=TRUE)
text(c(5, 64, 26),c( 5.815971e-05,  8.010516e-05, -1.111337e-04), 
        c("Fertility shape","SRB","Mortality"), cex = 1.3,pos = 4, col = c(gray(.2),gray(.2),gray(.2)))
segments(43.88283,-8.605324e-05 ,48.24328,-2.648702e-05)
segments(67.86529,6.442984e-05,62.95979,2.367401e-05)
text(50,-.0006,"remaining years",xpd=TRUE)
text(-1,.00045,"Contribution\nto difference in r", pos = 4, xpd = TRUE)
dev.off()

