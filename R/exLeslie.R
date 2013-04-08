source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

yearsUS <- 1969:2009
yearsES <- 1975:2009

BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))
names(BxymfES) <- yearsES
names(BxymfUS) <- yearsUS

# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
dxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfES.Rdata"))) 

dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)

PxUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxUS.Rdata"))) 
PxES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxUS.Rdata"))) 

# get male and female lambdas
DTLTUUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_TLTU/DTLTUUS.Rdata"))) 
DTLTUES <- local(get(load("/home/triffe/git/DISS/Data/HMD_TLTU/DTLTUES.Rdata"))) 


getLambda <- function(DTLTU, sex = "Female"){
    D0 <- reshape2::acast(DTLTU[DTLTU$Age == 0,],Year~Cohort,sum,value.var=sex)
    D0[D0==0] <- NA
    D0 <- do.call(cbind,apply(D0,2,function(x){
                        if (sum(!is.na(x))==2){
                            x[!is.na(x)]
                        }
                    }))
    lambda        <- D0[1,] / colSums(D0)
    years         <- as.integer(names(lambda))
    years         <-c(years,years[length(years)]+1)
    lambda        <- c(lambda,lambda[length(lambda)])
    names(lambda) <- years
    lambda
}

lambdamUS <- getLambda(DTLTUUS, "Male")
lambdafUS <- getLambda(DTLTUUS, "Female")
lambdamES <- getLambda(DTLTUES, "Male")
lambdafES <- getLambda(DTLTUES, "Female")

# take a look at lambda:
# pretty stable, but clear pattern.
# don't bother including in thesis, just mention
# range observed.
#plot(yearsUS, lambdamUS,  col = "blue", ylim = c(.84,.91), pch = 19)
#points(yearsUS, lambdafUS, col = "red", pch = 19)
#points(yearsES, lambdamES, col = "blue", pch = 17)
#points(yearsES, lambdafES, col = "red", pch = 17)
#lines(lowess(lambdafES~yearsES), col = "red", lty = 2)
#lines(lowess(lambdamES~yearsES), col = "blue", lty = 2)
#lines(lowess(lambdafUS~yearsUS), col = "red")
#lines(lowess(lambdamUS~yearsUS), col = "blue")
#




#Lex2 <- MakeLExTwoSexProjMatrix()

# ------------------------------------------------------------------
# make 1-sex Leslie matrices for males and females
Y_1sexUS <- lapply(as.character(yearsUS), function(yr, .dxm, .dxf, .Ex, .Bxymf, .lambdaf, .lambdam){
            yri     <- as.integer(yr)
            .dxm.   <- .dxm[, yr]
            .dxf.   <- .dxf[, yr]
            ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
            BxMM    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
            
            ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
            BxFF    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
         
            # sex-sex-ex- specific rates:
            FxMM    <- BxMM / ExM
            FxFF    <- BxFF / ExF
  
            
            list(Yf = MakeLExOneSexProjMatrix(FxFF, .dxf., .lambdaf[yr]),
            Ym = MakeLExOneSexProjMatrix(FxMM, .dxm., .lambdam[yr]))
        }, .dxm = dxmUS, .dxf = dxfUS, 
        .Ex = ExUS, .Bxymf = BxymfUS, 
        .lambdaf = lambdafUS, .lambdam = lambdafES)

names(Y_1sexUS) <- yearsUS
# repeat for Spain
Y_1sexES <- lapply(as.character(yearsES), function(yr, .dxm, .dxf, .Ex, .Bxymf, .lambdaf, .lambdam){
            yri     <- as.integer(yr)
            .dxm.   <- .dxm[, yr]
            .dxf.   <- .dxf[, yr]
            ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
            BxMM    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
            
            ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
            BxFF    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
            
            # sex-sex-ex- specific rates:
            FxMM    <- BxMM / ExM
            FxFF    <- BxFF / ExF
            
            
            list(Yf = MakeLExOneSexProjMatrix(FxFF, .dxf., .lambdaf[yr]),
                    Ym = MakeLExOneSexProjMatrix(FxMM, .dxm., .lambdam[yr]))
        }, .dxm = dxmES, .dxf = dxfES, 
        .Ex = ExES, .Bxymf = BxymfES, 
        .lambdaf = lambdafES, .lambdam = lambdafES)
names(Y_1sexES) <- yearsES

library(popbio)



rES1sex <- do.call(rbind,lapply(Y_1sexES, function(x){
            c(r.f = log(max(Re(eigen(x[["Yf"]])$values))),
            r.m = log(max(Re(eigen(x[["Ym"]])$values))))
        }))
rUS1sex <- do.call(rbind,lapply(Y_1sexUS, function(x){
                    c(r.f = log(max(Re(eigen(x[["Yf"]])$values))),
                            r.m = log(max(Re(eigen(x[["Ym"]])$values))))
                }))

#plot(yearsUS, rUS1sex[,"r.f"], type = 'l', col = "red", ylim = c(-.02,.02))
#lines(yearsUS, rUS1sex[,"r.m"], col = "blue")
#lines(yearsES, rES1sex[,"r.f"], col = "red",lty=2)
#lines(yearsES, rES1sex[,"r.m"], col = "blue",lty=2)

# make 2-sex matrices
Y_2sexUS <- lapply(as.character(yearsUS), function(yr, .dxm, .dxf, .Ex, .Bxymf, .lambdaf, .lambdam, .sigma){
            yri     <- as.integer(yr)
            .dxm.   <- .dxm[, yr]
            .dxf.   <- .dxf[, yr]
            ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
            BxMM    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
            BxMF    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxm.))
            
            ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
            BxFF    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
            BxFM    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxf.))
            
            # sex-sex-ex- specific rates:
            FxMM    <- Mna0(Minf0(BxMM / ExM))
            FxFF    <- Mna0(Minf0(BxFF / ExF))
            FxMF    <- Mna0(Minf0(BxMF / ExM))
            FxFM    <- Mna0(Minf0(BxFM / ExF))
            
            MakeLExTwoSexProjMatrix(dxm = .dxm., dxf = .dxf., 
                    FexFF = FxFF, FexFM = FxFM, 
                    FexMM = FxMM, FexMF = FxMF, 
                    lambdaM = .lambdam[yr], lambdaF = .lambdaf[yr], sigma = .sigma)
        },
        .dxm = dxmUS, .dxf = dxfUS,
        .Ex = ExUS, .Bxymf = BxymfUS, 
        .lambdaf = lambdafUS, .lambdam = lambdafUS, .sigma= .5)

names(Y_2sexUS) <- yearsUS
# repeat for Spain
Y_2sexES <- lapply(as.character(yearsES), function(yr, .dxm, .dxf, .Ex, .Bxymf, .lambdaf, .lambdam, .sigma){
            yri     <- as.integer(yr)
            .dxm.   <- .dxm[, yr]
            .dxf.   <- .dxf[, yr]
            ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
            BxMM    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
            BxMF    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxm.))
            
            ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
            BxFF    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
            BxFM    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxf.))
            
            # sex-sex-ex- specific rates:
            FxMM    <- Mna0(Minf0(BxMM / ExM))
            FxFF    <- Mna0(Minf0(BxFF / ExF))
            FxMF    <- Mna0(Minf0(BxMF / ExM))
            FxFM    <- Mna0(Minf0(BxFM / ExF))
            
            MakeLExTwoSexProjMatrix(dxm = .dxm., dxf = .dxf., 
                    FexFF = FxFF, FexFM = FxFM, 
                    FexMM = FxMM, FexMF = FxMF, 
                    lambdaM = .lambdam[yr], lambdaF = .lambdaf[yr], sigma = .sigma)
        },
        .dxm = dxmES, .dxf = dxfES,
        .Ex = ExES, .Bxymf = BxymfES, 
        .lambdaf = lambdafES, .lambdam = lambdafES, .sigma= .5)
names(Y_2sexES) <- yearsES


# exercise: measure time to one sex is double size of other:

doubleUS <- unlist(lapply(as.character(yearsUS), function(yr, .Y_1sex, .Px, .dxm, .dxf, maxit = 1e4){
            
           Pyfi <-  ExpectedDx( with(.Px, Female[Year == as.integer(yr)]), .dxf[,yr])
           Pymi <-  ExpectedDx( with(.Px, Male[Year == as.integer(yr)]), .dxm[,yr])
            
           for (i in 1:maxit){
               Pyfi <- .Y_1sex[[yr]][["Yf"]] %*% Pyfi
               Pymi <- .Y_1sex[[yr]][["Ym"]] %*% Pymi
               SR <- sum(Pymi) / sum(Pyfi)
               if (SR > 2 | SR < .5){
                   break
               }
           }
           i
        }, .Y_1sex = Y_1sexUS, .dxm = dxmUS, .dxf = dxfUS, .Px = PxUS))
doubleES <- unlist(lapply(as.character(yearsES), function(yr, .Y_1sex, .Px, .dxm, .dxf, maxit = 1e4){
                    
                    Pyfi <-  ExpectedDx( with(.Px, Female[Year == as.integer(yr)]), .dxf[,yr])
                    Pymi <-  ExpectedDx( with(.Px, Male[Year == as.integer(yr)]), .dxm[,yr])
                    
                    for (i in 1:maxit){
                        Pyfi <- .Y_1sex[[yr]][["Yf"]] %*% Pyfi
                        Pymi <- .Y_1sex[[yr]][["Ym"]] %*% Pymi
                        SR <- sum(Pymi) / sum(Pyfi)
                        if (SR > 2 | SR < .5){
                            break
                        }
                    }
                    i
                }, .Y_1sex = Y_1sexES, .dxm = dxmES, .dxf = dxfES, .Px = PxES))

# repeating for Lotka divergence:
pxmUS <- apply(dxmUS, 2, dx2pxLes)
pxfUS <- apply(dxfUS, 2, dx2pxLes)
pxmES <- apply(dxmES, 2, dx2pxLes)
pxfES <- apply(dxfES, 2, dx2pxLes)
yrs2double <- compiler::cmpfun(function(Pxm, Pyf, Lf, Lm, maxit = 5000){
            Pxmi  <- Pxm
            Pyfi  <- Pyf
            SR    <- sum(Pxmi) / sum(Pyfi) 
            yrt   <- 0
            while (SR < 2 & SR > .5 & yrt < maxit){
                Pxmi <- c(Lm %*% Pxmi)
                Pyfi <- c(Lf %*% Pyfi)
                SR   <- sum(Pxmi) / sum(Pyfi) 
                yrt  <- yrt + 1
            }
            yrt
        })

doubleESL <- unlist(lapply(as.character(yearsES), function(yr, .BxymfES, .PxES, .ExES = ExES, .pxmES, .pxfES, age = 0:110){        
                    Bxym <- Mna0(.BxymfES[[yr]][["Bxym"]])
                    Bxyf <- Mna0(.BxymfES[[yr]][["Bxyf"]])
                    Exm  <- with(.ExES, Male[Year == as.integer(yr)])
                    Eyf  <- with(.ExES, Female[Year == as.integer(yr)])
                    Pxm  <- with(.PxES, Male[Year == as.integer(yr)])
                    Pyf  <- with(.PxES, Female[Year == as.integer(yr)])
                    
                    mm   <- Mna0(rowSums(Bxym) / Exm)
                    mf   <- Mna0(colSums(Bxyf) / Eyf)
                    
                    yrs2double(Pxm, Pyf, Lf = Leslie(mf, pxfUS[,yr]), Lm = Leslie(mm, pxmUS[,yr]), maxit = 500000)
                }, .BxymfES = BxymfES, .PxES = PxES, .ExES = ExES, .pxmES = pxmES, .pxfES = pxfES))
doubleUSL <- unlist(lapply(as.character(yearsUS), function(yr, .BxymfUS, .PxUS, .ExUS = ExUS, .pxmUS, .pxfUS, age = 0:110){        
                    Bxym <- Mna0(.BxymfUS[[yr]][["Bxym"]])
                    Bxyf <- Mna0(.BxymfUS[[yr]][["Bxyf"]])
                    Exm  <- with(.ExUS, Male[Year == as.integer(yr)])
                    Eyf  <- with(.ExUS, Female[Year == as.integer(yr)])
                    Pxm  <- with(.PxUS, Male[Year == as.integer(yr)])
                    Pyf  <- with(.PxUS, Female[Year == as.integer(yr)])
                    
                    mm   <- Mna0(rowSums(Bxym) / Exm)
                    mf   <- Mna0(colSums(Bxyf) / Eyf)
                    
                    yrs2double(Pxm, Pyf, Lf = Leslie(mf, pxfUS[,yr]), Lm = Leslie(mm, pxmUS[,yr]), maxit = 500000)
                }, .BxymfUS = BxymfUS, .PxUS = PxUS, .ExUS = ExUS, .pxmUS = pxmUS, .pxfUS = pxfUS))



TicksMaj <- 10 ^ (2:5)
TicksMagLab <- c("100","1000","10000","10000")
TickMin     <- log(c(t(outer(TicksMaj, 2:9))))

pdf("/home/triffe/git/DISS/latex/Figures/ExrSRdoubling.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, log(doubleUS), type = 'l', ylim = log(c(100,150000)), xlim = c(1968, 2010),
        axes = FALSE, xlab = "", ylab = "", col = gray(.2), lwd = 2, 
        panel.first = list(rect(1968, log(100), 2010, log(150000), col = gray(.95), border=NA),
                abline(h = log(TicksMaj), col = "white"),
                segments(1968, TickMin, 1969, TickMin,col = "white"),
                segments(2009, TickMin, 2010, TickMin,  col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, log(TicksMaj), TicksMagLab, pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10), log(100), seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, 4.1, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1962, 12.4, "years to\nSR > 2 or < .5", cex = 1, xpd = TRUE, pos = 4)))
lines(yearsES, log(doubleES), lty = 1, col = gray(.5), lwd = 2.5)
lines(yearsUS, log(doubleUSL), lty = 5, col = gray(.2), lwd = 2)
lines(yearsES, log(doubleESL), lty = 5, col = gray(.5), lwd = 2.5)
legend(1970,11.7, lty = c(1,5,1,5), col = gray(c(.2,.2,.5,.5)), lwd = c(2,2,2.5,2.5),bty = "n",
        legend = c(expression(US~e[x]), "US age",expression(ES~e[x]), "ES age"), xpd = TRUE)
dev.off()

# ----------------------------------------------------------------
# compare damping ratios for convergence
# indicator of speed to convergence:
DampES <- do.call(rbind,lapply(Y_1sexES, function(x){
                    c(FD = eigen.analysis(x[["Yf"]])$damping.ratio, MD = eigen.analysis(x[["Ym"]])$damping.ratio)
                }))
DampUS <- do.call(rbind,lapply(Y_1sexUS, function(x){
                    c(FD = eigen.analysis(x[["Yf"]])$damping.ratio, MD = eigen.analysis(x[["Ym"]])$damping.ratio)
                }))
plot(yearsUS, DampUS[,"FD"], type = 'l', ylim = c(1,1.1))
lines(yearsUS, DampUS[,"MD"])
lines(yearsES, DampES[,"FD"])
lines(yearsES, DampES[,"MD"])



ESL <- lapply(as.character(yearsES), 
                function(yr, .BxymfES, .PxES, .ExES = ExES, .pxmES, .pxfES, age = 0:110){        
                    Bxym <- Mna0(.BxymfES[[yr]][["Bxym"]])
                    Bxyf <- Mna0(.BxymfES[[yr]][["Bxyf"]])
                    Exm  <- with(.ExES, Male[Year == as.integer(yr)])
                    Eyf  <- with(.ExES, Female[Year == as.integer(yr)])
                    Pxm  <- with(.PxES, Male[Year == as.integer(yr)])
                    Pyf  <- with(.PxES, Female[Year == as.integer(yr)])
                    
                    mm   <- Mna0(rowSums(Bxym) / Exm)
                    mf   <- Mna0(colSums(Bxyf) / Eyf)
                    
                    list(Lf = Leslie(mf, pxfUS[,yr]),
                    Lm = Leslie(mm, pxmUS[,yr]))
                    
                }, .BxymfES = BxymfES, .PxES = PxES, .ExES = ExES, .pxmES = pxmES, .pxfES = pxfES)
USL <- lapply(as.character(yearsUS), 
                function(yr, .BxymfUS, .PxUS, .ExUS = ExUS, .pxmUS, .pxfUS, age = 0:110){        
                    Bxym <- Mna0(.BxymfUS[[yr]][["Bxym"]])
                    Bxyf <- Mna0(.BxymfUS[[yr]][["Bxyf"]])
                    Exm  <- with(.ExUS, Male[Year == as.integer(yr)])
                    Eyf  <- with(.ExUS, Female[Year == as.integer(yr)])
                    Pxm  <- with(.PxUS, Male[Year == as.integer(yr)])
                    Pyf  <- with(.PxUS, Female[Year == as.integer(yr)])
                    
                    mm   <- Mna0(rowSums(Bxym) / Exm)
                    mf   <- Mna0(colSums(Bxyf) / Eyf)
                    
                    list(Lf = Leslie(mf, pxfUS[,yr]),
                            Lm = Leslie(mm, pxmUS[,yr]))
                    
                }, .BxymfUS = BxymfUS, .PxUS = PxUS, .ExUS = ExUS, .pxmUS = pxmUS, .pxfUS = pxfUS)
names(USL) <- yearsUS
names(ESL) <- yearsES
DampESL <- do.call(rbind,lapply(ESL, function(x){
                    c(FD = eigen.analysis(x[["Lf"]])$damping.ratio, MD = eigen.analysis(x[["Lm"]])$damping.ratio)
                }))
DampUSL <- do.call(rbind,lapply(USL, function(x){
                    c(FD = eigen.analysis(x[["Lf"]])$damping.ratio, MD = eigen.analysis(x[["Lm"]])$damping.ratio)
                }))
        
# Plot Damping ratios:
pdf("/home/triffe/git/DISS/latex/Figures/Damping.pdf", height = 5, width = 5)
par(mai = c(.5,.5,.5,.2), xaxs = "i", yaxs = "i")
plot(yearsES, DampESL[,1], type = 'l', ylim = c(1.01,1.08), lty = 4, col = gray(.4), lwd = 1.2,
        axes = FALSE, xlim = c(1968,2010), xlab = "", ylab = "",
        panel.first = list(rect(1968, 1.01, 2010, 1.08, col = gray(.95), border=NA),
                abline(h = seq(1.01,1.08,by=.005), col = "white"),   
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(1.01,1.08,by=.01), seq(1.01,1.08,by=.01), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10), 1.01, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, 1, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1963, 1.085, "Damping ratio", cex = 1, xpd = TRUE, pos = 4)))
lines(yearsES, DampESL[,2], lty = 4, lwd = 1.2)
lines(yearsES, DampES[,1], lty = 4, lwd = 3, col = gray(.4))
lines(yearsES, DampES[,2], lty = 4, lwd = 3)

lines(yearsUS, DampUSL[,1], col = gray(.4), lwd = 1.2)
lines(yearsUS, DampUSL[,2], lwd = 1.2)
lines(yearsUS, DampUS[,1], lwd = 3, col = gray(.4))
lines(yearsUS, DampUS[,2], lwd = 3)

text(rep(1994,4),c(1.060118, 1.071, 1.032564, 1.041464),
        c(expression(US^M~e[y]),
                expression(US^F~e[y]),expression(US^M~age),expression(US^F~age)), cex = .8)
text(c(1985, 1977, 1985, 1990),c(1.055363, 1.072, 1.019397, 1.025493),c(expression(ES^M~e[y]),
                expression(ES^F~e[y]),expression(ES^M~age),expression(ES^F~age)), cex = .8)
dev.off()

# Total Oscillation along the path to stability:
EScohenL <- do.call(rbind,lapply(as.character(yearsES), function(yr, .L, .Px){
            Pxm  <- with(.Px, Male[Year == as.integer(yr)])
            Pxf  <- with(.Px, Female[Year == as.integer(yr)])
            
            FProj <- pop.projection(.L[[yr]][["Lf"]],  Pxf, 500)
            FCat <- t(t(FProj$stage.vectors) / colSums(FProj$stage.vectors))
            
            MProj <- pop.projection(.L[[yr]][["Lm"]],  Pxm, 500)
            MCat <- t(t(MProj$stage.vectors) / colSums(MProj$stage.vectors))
            
            
            c(FD1 = sum(abs(colSums(FCat - FProj$stable.stage))), FD2 = sum(abs(FCat - FProj$stable.stage)), 
              MD1 = sum(abs(colSums(MCat - MProj$stable.stage))), MD2 = sum(abs(MCat - MProj$stable.stage)))
        }, .L = ESL, .Px = PxES))
#
UScohenL <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .L, .Px){
                    Pxm  <- with(.Px, Male[Year == as.integer(yr)])
                    Pxf  <- with(.Px, Female[Year == as.integer(yr)])
                    
                    FProj <- pop.projection(.L[[yr]][["Lf"]],  Pxf, 500)
                    FCat <- t(t(FProj$stage.vectors) / colSums(FProj$stage.vectors))
                    
                    MProj <- pop.projection(.L[[yr]][["Lm"]],  Pxm, 500)
                    MCat <- t(t(MProj$stage.vectors) / colSums(MProj$stage.vectors))
                    
                    
                    c(FD1 = sum(abs(colSums(FCat - FProj$stable.stage))), FD2 = sum(abs(FCat - FProj$stable.stage)), 
                            MD1 = sum(abs(colSums(MCat - MProj$stable.stage))), MD2 = sum(abs(MCat - MProj$stable.stage)))
                }, .L = USL, .Px = PxUS))

EScohenY <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Y, .Px, .dxm, .dxf){
                   
                    Pxm  <- with(.Px, Male[Year == as.integer(yr)])
                    Pxf  <- with(.Px, Female[Year == as.integer(yr)])
                    
                    Pym  <- rowSums(ExpectedDx(Pxm, .dxm[,yr]))
                    Pyf  <- rowSums(ExpectedDx(Pxf, .dxf[,yr]))
                    
                    FProj <- pop.projection(.Y[[yr]][["Yf"]],  Pyf, 500)
                    FCat <- t(t(FProj$stage.vectors) / colSums(FProj$stage.vectors))
                    
                    MProj <- pop.projection(.Y[[yr]][["Ym"]],  Pym, 500)
                    MCat <- t(t(MProj$stage.vectors) / colSums(MProj$stage.vectors))
                    
                    c(FD1 = sum(abs(colSums(FCat - FProj$stable.stage))), FD2 = sum(abs(FCat - FProj$stable.stage)), 
                            MD1 = sum(abs(colSums(MCat - MProj$stable.stage))), MD2 = sum(abs(MCat - MProj$stable.stage)))
                }, .Y = Y_1sexES, .Px = PxES, .dxm = dxmES, .dxf = dxfES))

UScohenY <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Y, .Px, .dxm, .dxf){
                    Pxm  <- with(.Px, Male[Year == as.integer(yr)])
                    Pxf  <- with(.Px, Female[Year == as.integer(yr)])
                    
                    Pym  <- rowSums(ExpectedDx(Pxm, .dxm[,yr]))
                    Pyf  <- rowSums(ExpectedDx(Pxf, .dxf[,yr]))
                    
                    FProj <- pop.projection(.Y[[yr]][["Yf"]],  Pyf, 500)
                    FCat <- t(t(FProj$stage.vectors) / colSums(FProj$stage.vectors))
                    
                    MProj <- pop.projection(.Y[[yr]][["Ym"]],  Pym, 500)
                    MCat <- t(t(MProj$stage.vectors) / colSums(MProj$stage.vectors))
                    
                    
                    c(FD1 = sum(abs(colSums(FCat - FProj$stable.stage))), FD2 = sum(abs(FCat - FProj$stable.stage)), 
                            MD1 = sum(abs(colSums(MCat - MProj$stable.stage))), MD2 = sum(abs(MCat - MProj$stable.stage)))
                }, .Y = Y_1sexUS, .Px = PxUS, .dxm = dxmUS, .dxf = dxfUS))



# D2 plot
pdf("/home/triffe/git/DISS/latex/Figures/CohenD2.pdf", height = 5, width = 5)
par(mai = c(.5,.5,.5,.2), xaxs = "i", yaxs = "i")
plot(yearsES, EScohenL[,2], type = 'l', ylim = c(0,28), lty = 4, col = gray(.4), lwd = 1.2,
        axes = FALSE, xlim = c(1968,2010), xlab = "", ylab = "",
        panel.first = list(rect(1968, 0, 2010, 28, col = gray(.95), border=NA),
                abline(h = seq(0,28,by=2.5), col = "white"),   
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(0,28,by=5), seq(0,28,by=5), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10), 0, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -2, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1963, 30, "Total Oscillation", cex = 1, xpd = TRUE, pos = 4)))
lines(yearsES, EScohenL[,4], lty = 4, lwd = 1.2)
lines(yearsES, EScohenY[,2], lty = 4, lwd = 3, col = gray(.4))
lines(yearsES, EScohenY[,4], lty = 4, lwd = 3)

lines(yearsUS, UScohenL[,2], col = gray(.4), lwd = 1.2)
lines(yearsUS, UScohenL[,4], lwd = 1.2)
lines(yearsUS, UScohenY[,2], lwd = 3, col = gray(.4))
lines(yearsUS, UScohenY[,4], lwd = 3)


text(c(1988, 1988, 1976, 1976),c(2.003155,4.295235,10.635031,14.195071),c(expression(US^M~e[y]),
        expression(US^F~e[y]),expression(US^M~age),expression(US^F~age)), cex = .8)
text(c(1995, 1995, 1996, 1996),c(9.562143, 13.561091, 21.7, 26.289451),c(expression(ES^M~e[y]),
                expression(ES^F~e[y]),expression(ES^M~age),expression(ES^F~age)), cex = .8)
dev.off()

#pdf("/home/triffe/git/DISS/latex/Figures/CohenD1.pdf", height = 5, width = 5)
#par(mai = c(.5,.5,.5,.2), xaxs = "i", yaxs = "i")
#plot(yearsES, EScohenL[,1], type = 'l', ylim = c(0,28), lty = 4, col = gray(.4), lwd = 1.2,
#        axes = FALSE, xlim = c(1968,2010), xlab = "", ylab = "",
#        panel.first = list(rect(1968, 0, 2010, 28, col = gray(.95), border=NA),
#                abline(h = seq(0,28,by=2.5), col = "white"),   
#                abline(v = seq(1970, 2010, by = 5), col = "white"),
#                text(1968, seq(0,28,by=5), seq(0,28,by=5), pos = 2, cex = .8, xpd = TRUE),
#                text(seq(1970, 2010, by = 10), 0, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
#                text(1990, -2, "Year", cex = 1, pos = 1, xpd = TRUE),
#                text(1963, 30, "Total Oscillation", cex = 1, xpd = TRUE, pos = 4)))
#lines(yearsES, EScohenL[,3], lty = 4, lwd = 1.2)
#lines(yearsES, EScohenY[,1], lty = 4, lwd = 3, col = gray(.4))
#lines(yearsES, EScohenY[,3], lty = 4, lwd = 3)
#
#lines(yearsUS, UScohenL[,1], col = gray(.4), lwd = 1.2)
#lines(yearsUS, UScohenL[,3], lwd = 1.2)
#lines(yearsUS, UScohenY[,1], lwd = 3, col = gray(.4))
#lines(yearsUS, UScohenY[,3], lwd = 3)


#text(c(1988, 1988, 1976, 1976),c(2.003155,4.295235,10.635031,14.195071),c(expression(US^M~e[y]),
#                expression(US^F~e[y]),expression(US^M~age),expression(US^F~age)), cex = .8)
#text(c(1995, 1995, 1996, 1996),c(9.562143, 13.561091, 21.7, 26.289451),c(expression(ES^M~e[y]),
#                expression(ES^F~e[y]),expression(ES^M~age),expression(ES^F~age)), cex = .8)
#dev.off()
#library(popbio)
#citation(popbio)


DampES2 <- unlist(lapply(Y_2sexES, function(x){
                    dimnames(x) <- list(c(paste0("m",0:110),paste0("f",0:110)),
                            c(paste0("m",0:110),paste0("f",0:110)))
                    eigen.analysis(x)$damping.ratio
                }))
DampUS2 <- unlist(lapply(Y_2sexUS, function(x){
                    dimnames(x) <- list(c(paste0("m",0:110),paste0("f",0:110)),
                            c(paste0("m",0:110),paste0("f",0:110)))
                    eigen.analysis(x)$damping.ratio
                }))


# comparing damping ratios
pdf("/home/triffe/git/DISS/latex/Figures/Damping2.pdf", height = 5, width = 5)
par(mai = c(.5,.5,.5,.2), xaxs = "i", yaxs = "i")
plot(yearsES, DampES[,1], type = 'l', ylim = c(1.05,1.075), lty = 4, col = gray(.4), lwd = 1.2,
        axes = FALSE, xlim = c(1968,2010), xlab = "", ylab = "",
        panel.first = list(rect(1968, 1.05, 2010, 1.075, col = gray(.95), border=NA),
                abline(h = seq(1.05,1.075,by=.0025), col = "white"),   
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(1.05,1.075,by=.005), seq(1.05,1.075,by=.005), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10), 1.05, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, 1.0483, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1963, 1.077, "Damping ratio", cex = 1, xpd = TRUE, pos = 4)))
lines(yearsES, DampES[,2], lty = 4, lwd = 1.2, col = gray(.2))

lines(yearsUS, DampUS[,1], lwd = 1.2, col = gray(.4))
lines(yearsUS, DampUS[,2], lwd = 1.2, col = gray(.2))

lines(yearsUS, DampUS2,lwd=3)
lines(yearsES, DampES2,lwd=3,lty=4)
text(c(2002, 2003,2002),c(1.058886,1.067116,1.064),
        c(expression(US^M),
                expression(US^F),expression(US^2~sex~sigma==.5)), cex = .8)
text(c(1984, 1974, 1977),c(1.057536, 1.070381, 1.074),c(expression(ES^M),
                expression(ES^F),expression(ES^2~sex~sigma==.5)), cex = .8)
dev.off()

EScohenY2 <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Y, .Px, .dxm, .dxf){
                    Pxm  <- with(.Px, Male[Year == as.integer(yr)])
                    Pxf  <- with(.Px, Female[Year == as.integer(yr)])
                    
                    Pym  <- rowSums(ExpectedDx(Pxm, .dxm[,yr]))
                    Pyf  <- rowSums(ExpectedDx(Pxf, .dxf[,yr]))
                    
                    Proj <- pop.projection(.Y[[yr]],  c(Pym,Pyf), 500)
                    Cat <- t(t(Proj$stage.vectors) / colSums(Proj$stage.vectors))
                    
             c(D1 = sum(abs(colSums(Cat - Proj$stable.stage))), D2 = sum(abs(Cat - Proj$stable.stage))) 
                          
                }, .Y = Y_2sexES, .Px = PxES, .dxm = dxmES, .dxf = dxfES))
UScohenY2 <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Y, .Px, .dxm, .dxf){
                    Pxm  <- with(.Px, Male[Year == as.integer(yr)])
                    Pxf  <- with(.Px, Female[Year == as.integer(yr)])
                    
                    Pym  <- rowSums(ExpectedDx(Pxm, .dxm[,yr]))
                    Pyf  <- rowSums(ExpectedDx(Pxf, .dxf[,yr]))
                    
                    Proj <- pop.projection(.Y[[yr]],  c(Pym,Pyf), 500)
                    Cat <- t(t(Proj$stage.vectors) / colSums(Proj$stage.vectors))
                    
                    c(D1 = sum(abs(colSums(Cat - Proj$stable.stage))), D2 = sum(abs(Cat - Proj$stable.stage))) 
                    
                }, .Y = Y_2sexUS, .Px = PxUS, .dxm = dxmUS, .dxf = dxfUS))

pdf("/home/triffe/git/DISS/latex/Figures/CohenD22sex.pdf", height = 5, width = 5)
par(mai = c(.5,.5,.5,.2), xaxs = "i", yaxs = "i")
plot(yearsES, EScohenY[,2], type = 'l', ylim = c(0,12), lty = 4, col = gray(.4), lwd = 1.2,
        axes = FALSE, xlim = c(1968,2010), xlab = "", ylab = "",
        panel.first = list(rect(1968, 0, 2010, 12, col = gray(.95), border=NA),
                abline(h = seq(0,12,by=2), col = "white"),   
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(0,12,by=2), seq(0,12,by=2), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10), 0, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.5, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1963, 12.8, "Total Oscillation", cex = 1, xpd = TRUE, pos = 4)))
lines(yearsES, EScohenY[,4], lty = 4, lwd = 1.2, col = gray(.2))
lines(yearsES, EScohenY2[,2],  lty = 4, lwd = 3)

lines(yearsUS, UScohenY[,2], lwd = 1.2, col = gray(.4))
lines(yearsUS, UScohenY[,4], lwd = 1.2, col = gray(.2))
lines(yearsUS, UScohenY2[,2], lwd = 3)

text(c(1988, 1988,1998),c(2,4.3,4.3),c(expression(US^M~e[y]),
                expression(US^F~e[y]),expression(US~2~sex~sigma==.5)), cex = .8)
text(c(1988, 1988,1978),c(6.5, 11.225107,9.344068),c(expression(ES^M~e[y]),
                expression(ES^F~e[y]),expression(ES~2~sex~sigma==.5)), cex = .8)
segments(1998,4.1, 1998, UScohenY2[yearsUS==1998,2])
segments(1978,9.1, 1985, EScohenY2[yearsES==1985,2])
dev.off()