setwd("/home/triffe/git/DISS/")
source("R/UtilityFunctions.R")
source("R/MeanFunctions.R")

# BxES is 0:110, years 1975:2009
BxUS  <- local(get(load("Data/USbirths/USBxy0_110.Rdata")))
BxES  <- local(get(load("Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65
PxUS <- local(get(load("Data/HMD_Px/PxUS.Rdata")))
PxES <- local(get(load("Data/HMD_Px/PxES.Rdata")))
yearsUS <- 1969:2009
yearsES <- 1975:2009

BxymfES <- local(get(load("Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("Data/USbirths/USBxymf0_110.Rdata")))
names(BxymfES) <- yearsES
names(BxymfUS) <- yearsUS

# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS  <- local(get(load("Data/Exposures/USexp.Rdata")))
ExES  <- local(get(load("Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
dxmUS <- local(get(load("Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("Data/HMD_dx/dxfES.Rdata"))) 

dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)

# some sensitivity tests
exTwoSexLinearCoaleR <- compiler::cmpfun(function(dxm, dxf, FexFF, FexFM, FexMM, FexMF, 
                .a = .5:110.5, sigma = .5, maxit = 2e3, tol = 1e-15){  
            N               <- length(dxm)
            dxM    <- dxF   <- matrix(0, ncol = N, nrow = N)
            # remaining years go down rows. ages over columns
            dxmi            <- dxm
            dxfi            <- dxf
            for (i in 1:N){
                dxM[i, 1:length(dxmi)  ] <- dxmi 
                dxmi                     <- dxmi[2:length(dxmi) ]
                
                dxF[i, 1:length(dxfi)  ] <- dxfi 
                dxfi                     <- dxfi[2:length(dxfi) ]
            }  
            
             # assuming r = 0
            # first guess at SRB
            p.m  <- 1.05 / 2.05
            p.f  <- 1 / 2.05
            SRBi <- sum(sigma * p.m * rowSums(dxM) * FexMM + (1 - sigma) * p.f * rowSums(dxF) * FexFM) / 
                    sum(sigma * p.m * rowSums(dxM) * FexMF + (1 - sigma) * p.f * rowSums(dxF) * FexFF)
            # derive proportions male and female
            p.m <- SRBi / (1 + SRBi)
            p.f <- 1 / (1 + SRBi)
            # guess at R0
            R0      <-   (sigma * sum(p.m * rowSums(dxM) * (FexMF + FexMM)) + 
                        (1 - sigma) * sum(p.f * rowSums(dxF) * (FexFF + FexFM))) 
            T.guess <- sum(.a * sigma * p.m * rowSums(dxM) * (FexMF + FexMM) + 
                        .a * (1 - sigma) * p.f * rowSums(dxF) * (FexFF + FexFM)) / R0
            r.i      <- log(R0) / T.guess
            
            # be careful to discount Fex by SRB appropriately for males / females
            # prior to specification sigma
            # Based on modification of Coale (1957)
            for (i in 1:maxit){ # 
                deltai <- 1 - sum(
                        sigma * rowSums(p.m * t(t(dxM) * exp(-r.i * .a))) * (FexMM + FexMF) + # male - female   
                        (1 - sigma) * rowSums(p.f * t(t(dxF) * exp(-r.i * .a))) * (FexFF + FexFM)  # female - female
                ) 
                # the mean generation time self-corrects 
                # according to the error produced by the Lotka equation
                r.i <- r.i - (deltai / (T.guess - (deltai / r.i))) 
          
                # update SRB
                SRBi <- sum(rowSums(sigma * p.m * t(t(dxM) * exp(-r.i * .a))) * FexMM + 
                                 rowSums( (1 - sigma) * p.f * t(t(dxF) * exp(-r.i * .a))) * FexFM) / 
                        sum(rowSums(sigma * p.m * t(t(dxM) * exp(-r.i * .a))) * FexMF + 
                                 rowSums( (1 - sigma) * p.f * t(t(dxF) * exp(-r.i * .a))) * FexFF)
                p.m <- SRBi / (1 + SRBi)
                p.f <- 1 / (1 + SRBi)
                if (abs(deltai) < tol){
                    break
                }
                
            }
            return(c(r = r.i, SRB = SRBi))  

            })
exTwoSexLinearTy <- compiler::cmpfun(function(r, SRB, dxm, dxf, 
                FexFF, FexFM, FexMM, FexMF, .a = .5:110.5, sigma = .5){
            N               <- length(FexFF)
            dxM    <- dxF   <- matrix(0, ncol = N, nrow = N)
            # remaining years go down rows. ages over columns
            dxmi            <- dxm
            dxfi            <- dxf
            for (i in 1:N){
                dxM[i, 1:length(dxmi)  ] <- dxmi 
                dxmi                     <- dxmi[2:length(dxmi) ]
                
                dxF[i, 1:length(dxfi)  ] <- dxfi 
                dxfi                     <- dxfi[2:length(dxfi) ]
            }     
            p.m <- SRB / (1 + SRB)
            p.f <- 1 / (1 + SRB)
            
            sum(sigma * rowSums(.a * p.m * t(t(dxM) * exp(-r * .a))) * (FexMM + FexMF) +
                (1-sigma) * rowSums(.a * p.f * t(t(dxF) * exp(-r * .a))) * (FexFF + FexFM)) / 
            sum(sigma * rowSums(p.m * t(t(dxM) * exp(-r * .a))) * (FexMM + FexMF) +
                (1-sigma) * rowSums(p.f * t(t(dxF) * exp(-r * .a))) * (FexFF + FexFM))
        })

exTwoSexStableAge <- compiler::cmpfun(function(r, SRB, dxm, dxf, .a = .5:110.5){
    N    <- length(dxm)
    dxM  <- matrix(0, ncol = N, nrow = N)
    dxF  <- matrix(0, ncol = N, nrow = N)
    dxmi  <- dxm
    dxfi  <- dxf
    for (i in 1:N){
        dxM[i, 1:length(dxmi)  ] <- dxmi 
        dxmi <- dxmi[2:length(dxmi) ]
        dxF[i, 1:length(dxfi)  ] <- dxfi 
        dxfi <- dxfi[2:length(dxfi) ]
    }  
    # birth rate
    b <- 1 / (sum(rowSums((SRB / (1 + SRB)) * t(t(dxM)* exp(-r * .a)))) +
                sum(rowSums((1 / (1 + SRB)) * t(t(dxF) * exp(-r * .a))))
                )
    cym <- b * (SRB / (1 + SRB)) * rowSums(t(t(dxM)* exp(-r * .a)))
    cyf <- b * (1 / (1 + SRB)) * rowSums(t(t(dxF)* exp(-r * .a)))
    cbind(cym = cym, cyf = cyf)
        })

US <-do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxymf, .dxm, .dxf, .Ex){
                            yri     <- as.integer(yr)
                            .dxm. <- .dxm[, yr]
                            .dxf. <- .dxf[, yr]
                            ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
                            BxMM    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
                            BxMF    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxm.))
                            
                            ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
                            BxFF    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
                            BxFM    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxf.))
                            # sex-sex-ex- specific rates:
                            FxMM    <- BxMM / ExM
                            FxFF    <- BxFF / ExF
                            FxMF    <- BxMF / ExM
                            FxFM    <- BxFM / ExF
                            
                            r.0 <- exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
                                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                                    sigma = 0) 
                            r.5 <- exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
                                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                                    sigma = .5) 
                            r.1 <- exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
                                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                                    sigma = 1) 
                            Ty.0 <- exTwoSexLinearTy(r=r.0[1], SRB = r.0[2], dxm = .dxm., dxf = .dxf., 
                                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                                    sigma = 0)
                            Ty.5 <- exTwoSexLinearTy(r=r.5[1], SRB = r.5[2], dxm = .dxm., dxf = .dxf., 
                                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                                    sigma = .5)
                            Ty.1 <- exTwoSexLinearTy(r=r.1[1], SRB = r.1[2], dxm = .dxm., dxf = .dxf., 
                                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                                    sigma = 1)
                            R0.0 <- exp(r.0[1] * Ty.0)
                            R0.5 <- exp(r.5[1] * Ty.5)
                            R0.1 <- exp(r.1[1] * Ty.1)
                            c(r.0 = r.0[1], r.5 = r.5[1], r.1 = r.1[1], 
                              Ty.0 = Ty.0, Ty.5 = Ty.5, Ty.1 = Ty.1, 
                              R0.0 = R0.0, R0.5 = R0.5, R0.1 = R0.1,S.0 = r.0[2],S.5 = r.5[2],S.1 = r.1[2])
                        }, .Bxymf = BxymfUS, .dxm = dxmUS, .dxf = dxfUS, .Ex = ExUS))
rownames(US) <- yearsUS
        
ES <-do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxymf, .dxm, .dxf, .Ex){
                            yri     <- as.integer(yr)
                            .dxm. <- .dxm[, yr]
                            .dxf. <- .dxf[, yr]
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
                            
                            r.0 <- exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
                                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                                    sigma = 0) 
                            r.5 <- exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
                                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                                    sigma = .5) 
                            r.1 <- exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
                                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                                    sigma = 1) 
                            Ty.0 <- exTwoSexLinearTy(r=r.0[1], SRB = r.0[2], dxm = .dxm., dxf = .dxf., 
                                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                                    sigma = 0)
                            Ty.5 <- exTwoSexLinearTy(r=r.5[1], SRB = r.5[2], dxm = .dxm., dxf = .dxf., 
                                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                                    sigma = .5)
                            Ty.1 <- exTwoSexLinearTy(r=r.1[1], SRB = r.1[2], dxm = .dxm., dxf = .dxf., 
                                    FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                                    sigma = 1)
                            R0.0 <- exp(r.0[1] * Ty.0)
                            R0.5 <- exp(r.5[1] * Ty.5)
                            R0.1 <- exp(r.1[1] * Ty.1)
                            c(r.0 = r.0[1], r.5 = r.5[1], r.1 = r.1[1], 
                              Ty.0 = Ty.0, Ty.5 = Ty.5, Ty.1 = Ty.1, 
                              R0.0 = R0.0, R0.5 = R0.5, R0.1 = R0.1,S.0 = r.0[2],S.5 = r.5[2],S.1 = r.1[2])
                        }, .Bxymf = BxymfES, .dxm = dxmES, .dxf = dxfES, .Ex = ExES))
rownames(ES) <- yearsES

#rSigmaUS <- US
#rSigmaES <- ES
#save(rSigmaUS, file = "Data/results/exGoodmanr/rSigmaUS.Rdata")
#save(rSigmaES, file = "Data/results/exGoodmanr/rSigmaES.Rdata")

# manual block on table creation
#save.tables <- FALSE
#if (save.tables){                    
#dimnames(ES) <- list(yearsES, c("$r^{\\upsilon (\\sigma = 0)}$"  , "$r^{\\upsilon (\\sigma = .5)}$"  , "$r^{\\upsilon (\\sigma = 1)}$",
#                "$T^{\\upsilon (\\sigma = 0)}$"  , "$T^{\\upsilon (\\sigma = .5)}$"  , "$T^{\\upsilon (\\sigma = 1)}$",
#                "$R_0^{\\upsilon (\\sigma = 0)}$", "$R_0^{\\upsilon (\\sigma = .5)}$", "$R_0^{\\upsilon (\\sigma = 1)}$"))
#dimnames(US) <- list(yearsUS, c("$r^{\\upsilon (\\sigma = 0)}$"  , "$r^{\\upsilon (\\sigma = .5)}$"  , "$r^{\\upsilon (\\sigma = 1)}$",
#                "$T^{\\upsilon (\\sigma = 0)}$"  , "$T^{\\upsilon (\\sigma = .5)}$"  , "$T^{\\upsilon (\\sigma = 1)}$",
#                "$R_0^{\\upsilon (\\sigma = 0)}$", "$R_0^{\\upsilon (\\sigma = .5)}$", "$R_0^{\\upsilon (\\sigma = 1)}$"))
#
#library(xtable)
#print(xtable(ES, digits = c(0,4,4,4,2,2,2,3,3,3), align = c("c","c","c","c","c","c","c","c","c","c")),
#        sanitize.colnames.function = identity, 
#        file = "latex/xtables/ex2sexlinearES.tex",floating=FALSE)
#
#print(xtable(US, digits = c(0,4,4,4,2,2,2,3,3,3), align = c("c","c","c","c","c","c","c","c","c","c")),
#        sanitize.colnames.function = identity, 
#        file = "latex/xtables/ex2sexlinearUS.tex",floating=FALSE)
#}
#

# manual block on figure creation
make.fig <- FALSE
if (male.fig){
pdf("latex/Figures/exLotka2sexlinear.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, US[, 2], type = 'n', ylim = c(-.016,.01),xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,-.016,2010,.01,col = gray(.95), border=NA),
                abline(h = seq(-.016,.01,by = .002), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.016,.01,by = .002),seq(-.016,.01,by = .002), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.016, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.0175, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,.0115, "r", cex = 1, xpd = TRUE)))
polygon(c(yearsUS, rev(yearsUS)), c(US[, 1], rev(US[, 3])), col = "#99999950")
lines(yearsUS, US[, 2],col = gray(.2), lwd = 2)
lines(yearsUS, US[, 3],col = "blue", lwd = 1)
lines(yearsUS, US[, 1],col = "red", lwd = 1)

polygon(c(yearsES, rev(yearsES)), c(ES[, 1], rev(ES[, 3])), col = "#99999950", lty = 5)
lines(yearsES, ES[, 2],col = gray(.2), lwd = 2, lty = 5)
lines(yearsES, ES[, 3],col = "blue", lwd = 1, lty = 5)
lines(yearsES, ES[, 1],col = "red", lwd = 1, lty = 5)

legend(1970,-.007, 
        lty = c(1,1,1,5,5,5), 
        col = c("red",gray(.2),"blue","red",gray(.2),"blue"), 
        lwd = c(1,2,1,1,2,1), bty = "n",
        legend = c(expression(US~sigma == 0),
                expression(US~sigma == .5),
                expression(US~sigma == 1),
                expression(ES~sigma == 0),
                expression(ES~sigma == .5),
                expression(ES~sigma == 1)), xpd = TRUE)
dev.off()
}
# TODO: change text?
#US[,1] < US[,3]
#ES[,1] < ES[,3]
#
#
#rfUS[,1] < US[,1]
#rmUS[,1] > US[,3]
#
#rfES[,1] < ES[,1]
#rmES[,1] > ES[,3]

#US <-do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxymf, .dxm, .dxf, .Ex){
#                    yri     <- as.integer(yr)
#                    .dxm. <- .dxm[, yr]
#                    .dxf. <- .dxf[, yr]
#                    ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
#                    BxMM    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
#                    BxMF    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxm.))
#                    
#                    ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
#                    BxFF    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
#                    BxFM    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxf.))
#                    # sex-sex-ex- specific rates:
#                    FxMM    <- BxMM / ExM
#                    FxFF    <- BxFF / ExF
#                    FxMF    <- BxMF / ExM
#                    FxFM    <- BxFM / ExF
#                    
#                    r.5 <- exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
#                            FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
#                            sigma = .5) 
#             
#                }, .Bxymf = BxymfUS, .dxm = dxmUS, .dxf = dxfUS, .Ex = ExUS))

# ---------------------------------------------

# This is a manual block on the remainder of the code here, so that the above can be source()ed
# the below does some stable pyramids, transient analysis, etc. several Figures produced
do.rest <- FALSE
if (do.rest){

StableStructUS.5 <- lapply(as.character(yearsUS), function(yr, .Bx, .Ex, .Px, .dxm, .dxf, .sigma){  
                    yri     <- as.integer(yr)
                    .dxm.   <- .dxm[, yr]
                    .dxf.   <- .dxf[, yr]
                    ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
                    BxMM    <- rowSums(ExpectedDx( rowSums(.Bx[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
                    BxMF    <- rowSums(ExpectedDx( rowSums(.Bx[[yr]][["Bxyf"]], na.rm = TRUE), .dxm.))
                    
                    ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
                    BxFF    <- rowSums(ExpectedDx( colSums(.Bx[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
                    BxFM    <- rowSums(ExpectedDx( colSums(.Bx[[yr]][["Bxym"]], na.rm = TRUE), .dxf.))
                    # sex-sex-ex- specific rates:
                    FxMM    <- Minf0(Mna0(BxMM / ExM))
                    FxFF    <- Minf0(Mna0(BxFF / ExF))
                    FxMF    <- Minf0(Mna0(BxMF / ExM))
                    FxFM    <- Minf0(Mna0(BxFM / ExF))
                    # 2000: 0.0005731307 1.0480159308 
                    r.5     <- exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
                            FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                            sigma = .sigma) 
            # stable structure
                    cyst    <- exTwoSexStableAge(r = r.5[1], SRB = r.5[2], dxm = .dxm., dxf = .dxf.)
                    PyM     <- rowSums(ExpectedDx( with(.Px, Male[Year == yri]), .dxm.))
                    PyF     <- rowSums(ExpectedDx( with(.Px, Female[Year == yri]), .dxf.))
                    # preesnt structure
                    cyinit  <- cbind(PyM,PyF) / sum(PyM+PyF)

            # difference coef
            list(extheta = 1-sum(pmin(cyst, cyinit)), cyst = cyst, cyinit = cyinit, r.5 = r.5)
                }, .Bx = BxymfUS, .Ex = ExUS, .Px = PxUS,  .dxm = dxmUS, .dxf = dxfUS, .sigma = .5)
        
StableStructES.5 <- lapply(as.character(yearsES), function(yr, .Bx, .Ex, .Px, .dxm, .dxf, .sigma){  
                    yri     <- as.integer(yr)
                    .dxm.   <- .dxm[, yr]
                    .dxf.   <- .dxf[, yr]
                    ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
                    BxMM    <- rowSums(ExpectedDx( rowSums(.Bx[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
                    BxMF    <- rowSums(ExpectedDx( rowSums(.Bx[[yr]][["Bxyf"]], na.rm = TRUE), .dxm.))
                    
                    ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
                    BxFF    <- rowSums(ExpectedDx( colSums(.Bx[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
                    BxFM    <- rowSums(ExpectedDx( colSums(.Bx[[yr]][["Bxym"]], na.rm = TRUE), .dxf.))
                    # sex-sex-ex- specific rates:
                    FxMM    <- Minf0(Mna0(BxMM / ExM))
                    FxFF    <- Minf0(Mna0(BxFF / ExF))
                    FxMF    <- Minf0(Mna0(BxMF / ExM))
                    FxFM    <- Minf0(Mna0(BxFM / ExF))
                    
                    r.5     <- exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
                            FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
                            sigma = .sigma) 
                    # stable structure
                    cyst    <- exTwoSexStableAge(r = r.5[1], SRB = r.5[2], dxm = .dxm., dxf = .dxf.)
                    PyM     <- rowSums(ExpectedDx( with(.Px, Male[Year == yri]), .dxm.))
                    PyF     <- rowSums(ExpectedDx( with(.Px, Female[Year == yri]), .dxf.))
                    # preesnt structure
                    cyinit  <- cbind(PyM,PyF) / sum(PyM+PyF)
                    
                    # difference coef
                    list(extheta = 1-sum(pmin(cyst, cyinit)), cyst = cyst, cyinit = cyinit, r.5 = r.5)
                }, .Bx = BxymfES, .Ex = ExES, .Px = PxES,  .dxm = dxmES, .dxf = dxfES, .sigma = .5)
# interesting: male dominance implies closer structure to present...
#StableStructES.1 <- lapply(as.character(yearsES), function(yr, .Bx, .Ex, .Px, .dxm, .dxf, .sigma){  
#                    yri     <- as.integer(yr)
#                    .dxm.   <- .dxm[, yr]
#                    .dxf.   <- .dxf[, yr]
#                    ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
#                    BxMM    <- rowSums(ExpectedDx( rowSums(.Bx[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
#                    BxMF    <- rowSums(ExpectedDx( rowSums(.Bx[[yr]][["Bxyf"]], na.rm = TRUE), .dxm.))
#                    
#                    ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
#                    BxFF    <- rowSums(ExpectedDx( colSums(.Bx[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
#                    BxFM    <- rowSums(ExpectedDx( colSums(.Bx[[yr]][["Bxym"]], na.rm = TRUE), .dxf.))
#                    # sex-sex-ex- specific rates:
#                    FxMM    <- Minf0(Mna0(BxMM / ExM))
#                    FxFF    <- Minf0(Mna0(BxFF / ExF))
#                    FxMF    <- Minf0(Mna0(BxMF / ExM))
#                    FxFM    <- Minf0(Mna0(BxFM / ExF))
#                    
#                    r.5     <- exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
#                            FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
#                            sigma = .sigma) 
#                    # stable structure
#                    cyst    <- exTwoSexStableAge(r = r.5[1], SRB = r.5[2], dxm = .dxm., dxf = .dxf.)
#                    PyM     <- rowSums(ExpectedDx( with(.Px, Male[Year == yri]), .dxm.))
#                    PyF     <- rowSums(ExpectedDx( with(.Px, Female[Year == yri]), .dxf.))
#                    # preesnt structure
#                    cyinit  <- cbind(PyM,PyF) / sum(PyM+PyF)
#                    
#                    # difference coef
#                    list(extheta = 1-sum(pmin(cyst, cyinit)), cyst = cyst, cyinit = cyinit)
#                }, .Bx = BxymfES, .Ex = ExES, .Px = PxES,  .dxm = dxmES, .dxf = dxfES, .sigma = 1)
#StableStructES.0 <- lapply(as.character(yearsES), function(yr, .Bx, .Ex, .Px, .dxm, .dxf, .sigma){  
#                    yri     <- as.integer(yr)
#                    .dxm.   <- .dxm[, yr]
#                    .dxf.   <- .dxf[, yr]
#                    ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
#                    BxMM    <- rowSums(ExpectedDx( rowSums(.Bx[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
#                    BxMF    <- rowSums(ExpectedDx( rowSums(.Bx[[yr]][["Bxyf"]], na.rm = TRUE), .dxm.))
#                    
#                    ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
#                    BxFF    <- rowSums(ExpectedDx( colSums(.Bx[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
#                    BxFM    <- rowSums(ExpectedDx( colSums(.Bx[[yr]][["Bxym"]], na.rm = TRUE), .dxf.))
#                    # sex-sex-ex- specific rates:
#                    FxMM    <- Minf0(Mna0(BxMM / ExM))
#                    FxFF    <- Minf0(Mna0(BxFF / ExF))
#                    FxMF    <- Minf0(Mna0(BxMF / ExM))
#                    FxFM    <- Minf0(Mna0(BxFM / ExF))
#                    
#                    r.5     <- exTwoSexLinearCoaleR(dxm = .dxm., dxf = .dxf., 
#                            FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,
#                            sigma = .sigma) 
#                    # stable structure
#                    cyst    <- exTwoSexStableAge(r = r.5[1], SRB = r.5[2], dxm = .dxm., dxf = .dxf.)
#                    PyM     <- rowSums(ExpectedDx( with(.Px, Male[Year == yri]), .dxm.))
#                    PyF     <- rowSums(ExpectedDx( with(.Px, Female[Year == yri]), .dxf.))
#                    # preesnt structure
#                    cyinit  <- cbind(PyM,PyF) / sum(PyM+PyF)
#                    
#                    # difference coef
#                    list(extheta = 1-sum(pmin(cyst, cyinit)), cyst = cyst, cyinit = cyinit)
#                }, .Bx = BxymfES, .Ex = ExES, .Px = PxES,  .dxm = dxmES, .dxf = dxfES, .sigma = 0)
#        

# objects manually created from exLotka1Sex.R
#DiffCoefryUSm[,1]
#DiffCoefryUSf[,1]
#DiffCoefryESm[,1]
#DiffCoefryESf[,1]

thetaUS <- unlist(lapply(StableStructUS.5,"[[",1))
thetaES <- unlist(lapply(StableStructES.5,"[[",1))

#plot(yearsUS, thetaUS, type = 'l')
#lines(yearsUS, DiffCoefryUSm[,1], col = "blue")
#lines(yearsUS, DiffCoefryUSf[,1], col = "red")

mxmUS <- local(get(load("Data/HMD_mux/muxmUS.Rdata"))) 
mxfUS <- local(get(load("Data/HMD_mux/muxfUS.Rdata"))) 
mxmES <- local(get(load("Data/HMD_mux/muxmES.Rdata"))) 
mxfES <- local(get(load("Data/HMD_mux/muxfES.Rdata"))) 

library(parallel)

DiffCoefryUSm

UScy <- lapply(StableStructUS.5,"[[",2)
USin <- lapply(StableStructUS.5,"[[",3)

EScy <- lapply(StableStructES.5,"[[",2)
ESin <- lapply(StableStructES.5,"[[",3)

names(UScy) <- yearsUS
names(USin) <- yearsUS
names(EScy) <- yearsES
names(ESin) <- yearsES
xlabs <- c("1.0%","0.8%","0.6%","0.4%","0.2%","0%","0.2%","0.4%","0.5%","0.6%","0.8%","1.0%")

pdf("latex/Figures/exLotka2sexlinear1975USpyr.pdf", height = 5, width = 5)
par(mai = c(.6,.6,.3,.3), xaxs = "i", yaxs = "i")
plot(NULL, type = "n",axes = FALSE, xlab = "",ylab = "", xlim = c(-1, 1), ylim = c(0,111),
        panel.first = list(
                rect(-1, 0, 1, 111, col = gray(.95), border = NA),
                abline(v = seq(-1, 1, by = .2), col = "white"),
                abline(h = seq(0, 110, by = 10), col = "white"),
                text(seq(-1, 1, by = .2),0, xlabs, xpd = TRUE, pos = 1, cex = .7),
                text(-1, seq(0, 110, by = 10), seq(0, 110, by = 10), pos = 2, xpd = TRUE, cex = .7),
                text(-1.17, 116, expression(e[y]), xpd = TRUE, cex = 1),
                text(0, -12, "Percentage", xpd = TRUE, cex = 1)
        ))

PyramidOutline(UScy[["1975"]][,1], UScy[["1975"]][,2], scale =100, border = gray(.2), xpd = TRUE, lwd = .5, col = "#55555550") 
PyramidOutline(USin[["1975"]][,1], USin[["1975"]][,2], scale =100, border = gray(.2), xpd = TRUE, lwd = .5) 
text(c(-.63,-.55), c(50, 91), c(expression(Stable~sigma==0.5),"Initial"), col = c("white", "black"), cex = 1.2, pos = 4)
segments(-.4, 88, -0.2485961, 80.14833)
dev.off()

pdf("latex/Figures/exLotka2sexlinear2009USpyr.pdf", height = 5, width = 5)
par(mai = c(.6,.6,.3,.3), xaxs = "i", yaxs = "i")
plot(NULL, type = "n",axes = FALSE, xlab = "",ylab = "", xlim = c(-1, 1), ylim = c(0,111),
        panel.first = list(
                rect(-1, 0, 1, 111, col = gray(.95), border = NA),
                abline(v = seq(-1, 1, by = .2), col = "white"),
                abline(h = seq(0, 110, by = 10), col = "white"),
                text(seq(-1, 1, by = .2),0, xlabs, xpd = TRUE, pos = 1, cex = .7),
                text(-1, seq(0, 110, by = 10), seq(0, 110, by = 10), pos = 2, xpd = TRUE, cex = .7),
                text(-1.17, 116, expression(e[y]), xpd = TRUE, cex = 1),
                text(0, -12, "Percentage", xpd = TRUE, cex = 1)
        ))

PyramidOutline(UScy[["2009"]][,1], UScy[["2009"]][,2], scale =100, border = gray(.2), xpd = TRUE, lwd = .5, col = "#55555550") 
PyramidOutline(USin[["2009"]][,1], USin[["2009"]][,2], scale =100, border = gray(.2), xpd = TRUE, lwd = .5) 
text(c(-.63,-.55), c(50, 91), c(expression(Stable~sigma==0.5),"Initial"), col = c("white", "black"), cex = 1.2, pos = 4)
segments(-.4, 88,  -0.3022095, 83)
dev.off()
# Spain 1975 and 2009:

pdf("latex/Figures/exLotka2sexlinear1975ESpyr.pdf", height = 5, width = 5)
par(mai = c(.6,.6,.3,.3), xaxs = "i", yaxs = "i")
plot(NULL, type = "n",axes = FALSE, xlab = "",ylab = "", xlim = c(-1, 1), ylim = c(0,111),
        panel.first = list(
                rect(-1, 0, 1, 111, col = gray(.95), border = NA),
                abline(v = seq(-1, 1, by = .2), col = "white"),
                abline(h = seq(0, 110, by = 10), col = "white"),
                text(seq(-1, 1, by = .2),0, xlabs, xpd = TRUE, pos = 1, cex = .7),
                text(-1, seq(0, 110, by = 10), seq(0, 110, by = 10), pos = 2, xpd = TRUE, cex = .7),
                text(-1.17, 116, expression(e[y]), xpd = TRUE, cex = 1),
                text(0, -12, "Percentage", xpd = TRUE, cex = 1)
        ))

PyramidOutline(EScy[["1975"]][,1], EScy[["1975"]][,2], scale =100, border = gray(.2), xpd = TRUE, lwd = .5, col = "#55555550") 
PyramidOutline(ESin[["1975"]][,1], ESin[["1975"]][,2], scale =100, border = gray(.2), xpd = TRUE, lwd = .5) 
text(c(-.66,.55), c(50, 91), c(expression(Stable~sigma==0.5),"Initial"), col = c("white", "black"), cex = 1.2, pos = 4)
segments(0.6293243,  88.70528,0.5153957, 79.00740)
dev.off()

pdf("latex/Figures/exLotka2sexlinear2009ESpyr.pdf", height = 5, width = 5)
par(mai = c(.6,.6,.3,.3), xaxs = "i", yaxs = "i")
plot(NULL, type = "n",axes = FALSE, xlab = "",ylab = "", xlim = c(-1, 1), ylim = c(0,111),
        panel.first = list(
                rect(-1, 0, 1, 111, col = gray(.95), border = NA),
                abline(v = seq(-1, 1, by = .2), col = "white"),
                abline(h = seq(0, 110, by = 10), col = "white"),
                text(seq(-1, 1, by = .2),0, xlabs, xpd = TRUE, pos = 1, cex = .7),
                text(-1, seq(0, 110, by = 10), seq(0, 110, by = 10), pos = 2, xpd = TRUE, cex = .7),
                text(-1.17, 116, expression(e[y]), xpd = TRUE, cex = 1),
                text(0, -12, "Percentage", xpd = TRUE, cex = 1)
        ))

PyramidOutline(EScy[["2009"]][,1], EScy[["2009"]][,2], scale =100, border = gray(.2), xpd = TRUE, lwd = .5, col = "#55555550") 
PyramidOutline(ESin[["2009"]][,1], ESin[["2009"]][,2], scale =100, border = gray(.2), xpd = TRUE, lwd = .5) 
text(c(-.64,-.87), c(35, 64), c(expression(Stable~sigma==0.5),"Initial"), col = c("white", "black"), cex = 1.2, pos = 4)
dev.off()

#for (yr in as.character(yearsUS)){
#par(mai = c(.6,.6,.3,.3), xaxs = "i", yaxs = "i")
#plot(NULL, type = "n",axes = FALSE, xlab = "",ylab = "", xlim = c(-1, 1), ylim = c(0,111),main=yr,
#        panel.first = list(
#                rect(-1, 0, 1, 111, col = gray(.95), border = NA),
#                abline(v = seq(-1, 1, by = .2), col = "white"),
#                abline(h = seq(0, 110, by = 10), col = "white"),
#                text(seq(-1, 1, by = .2),0, xlabs, xpd = TRUE, pos = 1, cex = .7),
#                text(-1, seq(0, 110, by = 10), seq(0, 110, by = 10), pos = 2, xpd = TRUE, cex = .7),
#                text(-1.17, 116, expression(e[y]), xpd = TRUE, cex = 1),
#                text(0, -12, "Percentage", xpd = TRUE, cex = 1)
#        ))
#PyramidOutline(UScy[[yr]][,1], UScy[[yr]][,2], 
#        scale =100, border = gray(.2), xpd = TRUE, lwd = .5, col = "#55555550") 
#Sys.sleep(.3)
#}

# ----------------------------------------------------------
# stable distribution (cy) under different levels of r (-.01 - .01),
# using dx from 1975 USA.
r.5US <- do.call(rbind,lapply(StableStructUS.5,"[[",4))
r.5ES <- do.call(rbind,lapply(StableStructES.5,"[[",4))
rownames(r.5US) <- yearsUS
rownames(r.5ES) <- yearsES

yr <- "1975"
pdf("latex/Figures/exLotka2sexlinearPyrDiffr.pdf", height = 5, width = 5)
par(mai = c(.6,.6,.3,.3), xaxs = "i", yaxs = "i")
plot(NULL, type = "n",axes = FALSE, xlab = "",ylab = "", xlim = c(-1, 1), ylim = c(0,111),
        panel.first = list(
                rect(-1, 0, 1, 111, col = gray(.95), border = NA),
                abline(v = seq(-1, 1, by = .2), col = "white"),
                abline(h = seq(0, 110, by = 10), col = "white"),
                text(seq(-1, 1, by = .2),0, xlabs, xpd = TRUE, pos = 1, cex = .7),
                text(-1, seq(0, 110, by = 10), seq(0, 110, by = 10), pos = 2, xpd = TRUE, cex = .7),
                text(-1.17, 116, expression(e[y]), xpd = TRUE, cex = 1),
                text(0, -12, "Percentage", xpd = TRUE, cex = 1)
        ))
rvals <- seq(-.01,.01,by = .005)
for(rs in 1:length(rvals)){
    pop <- exTwoSexStableAge(r = rvals[rs], SRB = r.5US[yr,2], dxm = dxmUS[,yr], dxf = dxfUS[,yr])
    PyramidOutline(pop[,1], pop[,2], 
            scale =100, border = gray(.4), xpd = TRUE, lwd = .5) 
#        text(-100*pop[4+rs*2,1],4+rs*2-.5,paste0(rvals[rs],"%"),cex = .9)
#        text(-100*pop[52+rs*2,1],52+rs*2-.5,paste0(rvals[rs],"%"),cex = .9)
}
text(c(-0.82, -0.61, -0.58),
        c(8,  9.5, 11),rvals[c(1,3,5)],cex = .7,pos=c(2,2,4),xpd=TRUE)
text(c(-0.51, -0.66, -0.66),
        c(51.5,  56, 60),rvals[c(1,3,5)],cex = .7,pos=c(4,4,2),xpd=TRUE)
dev.off()

# what about under different mortalit assumptions?
mxmUS <- local(get(load("Data/HMD_mux/muxmUS.Rdata"))) 
mxfUS <- local(get(load("Data/HMD_mux/muxfUS.Rdata"))) 
yr <- "1975"
pdf("latex/Figures/exLotka2sexlinearPyrDiffdx.pdf", height = 5, width = 5)
par(mai = c(.6,.6,.3,.3), xaxs = "i", yaxs = "i")
plot(NULL, type = "n",axes = FALSE, xlab = "",ylab = "", xlim = c(-1, 1), ylim = c(0,111),
        panel.first = list(
                rect(-1, 0, 1, 111, col = gray(.95), border = NA),
                abline(v = seq(-1, 1, by = .2), col = "white"),
                abline(h = seq(0, 110, by = 10), col = "white"),
                text(seq(-1, 1, by = .2),0, xlabs, xpd = TRUE, pos = 1, cex = .7),
                text(-1, seq(0, 110, by = 10), seq(0, 110, by = 10), pos = 2, xpd = TRUE, cex = .7),
                text(-1.17, 116, expression(e[y]), xpd = TRUE, cex = 1),
                text(0, -12, "Percentage", xpd = TRUE, cex = 1)
        ))
imp <- seq(.8,1.2,by=.2)
for(rs in 1:length(imp)){
    pop <- exTwoSexStableAge(r = 0, SRB = r.5US[yr,2], 
            dxm = mx2dxHMD(mxmUS[,yr]*imp[rs]), 
            dxf = mx2dxHMD(mxfUS[,yr]*imp[rs]))
    PyramidOutline(pop[,1], pop[,2], 
            scale =100, border = gray(.4), xpd = TRUE, lwd = .5) 
#        text(-100*pop[4+rs*2,1],4+rs*2-.5,paste0(rvals[rs],"%"),cex = .9)
#        text(-100*pop[52+rs*2,1],52+rs*2-.5,paste0(rvals[rs],"%"),cex = .9)
}
text(c(-.65,-.7),c(15,15),c("80%","120%"),pos=c(4,2),cex = .9)
text(c(-0.38, -0.15),c(80,73),c("80%","120%"),cex=.9)
dev.off()



#plot(0:110,mx2dxHMD(mxmUS[,yr]*.8),type = 'l')
#lines(0:110,mx2dxHMD(mxmUS[,yr]*1.2),col = "red")
#
#wmean(.5:110.5,dxmUS[,yr])
#wmean(.5:110.5,mx2dxHMD(mxmUS[,yr]*.8))
#wmean(.5:110.5,mx2dxHMD(mxmUS[,yr]*1.2))

}