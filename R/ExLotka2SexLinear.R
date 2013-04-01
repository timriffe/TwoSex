source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxES is 0:110, years 1975:2009
BxUS  <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
BxES  <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

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


ExM     <- rowSums(ExpectedDx( with(ExUS, Male[Year == 1975]), dxmUS[, "1975"]))
BxMM    <- rowSums(ExpectedDx( rowSums(BxymfUS[["1975"]][["Bxym"]], na.rm = TRUE), dxmUS[, "1975"]))
BxMF    <- rowSums(ExpectedDx( rowSums(BxymfUS[["1975"]][["Bxyf"]], na.rm = TRUE), dxmUS[, "1975"]))

BxM    <- rowSums(ExpectedDx( rowSums(BxUS[["1975"]], na.rm = TRUE), dxmUS[, "1975"]))
BxF    <- rowSums(ExpectedDx( colSums(BxUS[["1975"]], na.rm = TRUE), dxfUS[, "1975"]))
ExF     <- rowSums(ExpectedDx( with(ExUS, Female[Year == 1975]), dxfUS[, "1975"]))
BxFF    <- rowSums(ExpectedDx( colSums(BxymfUS[["1975"]][["Bxyf"]], na.rm = TRUE), dxfUS[, "1975"]))
BxFM    <- rowSums(ExpectedDx( colSums(BxymfUS[["1975"]][["Bxym"]], na.rm = TRUE), dxfUS[, "1975"]))
# some rates:
FxMM    <- BxMM / ExM
FxFF    <- BxFF / ExF
FxMF    <- BxMF / ExM
FxFM    <- BxFM / ExF
# minimizer function for 1 sex ex-perspective renewal function:


#exTwoSexLinearMin <- compiler::cmpfun(function(r, dxm, dxf, FexFF, FexFM, FexMM, FexMF, .a = .5:110.5, sigma = .5){
#    # all input vectors (except r) must be same length!
#    # get the overlapped / staggered dx structure
#    N               <- length(FexFF)
#    dxM  <- dxF   <- matrix(0, ncol = N, nrow = N)
#    # remaining years go down rows. ages over columns
#    dxmi  <- dxm
#    dxfi  <- dxf
#    for (i in 1:N){
#        dxM[i, 1:length(dxmi)  ] <- dxmi 
#        dxmi <- dxmi[2:length(dxmi) ]
#        
#        dxF[i, 1:length(dxfi)  ] <- dxfi 
#        dxfi <- dxfi[2:length(dxfi) ]
#    }    
#    
#    
#    (2 - sum(sigma * rowSums(dxM %col% (1 / exp(-r * .a))) * (FexMF + FexMM)# male - female
#         (1-sigma) * rowSums(dxF %col% (1 / exp(-r * .a))) * (FexFF + FexFM) # female - female
#             )
#           ) ^ 2
#})

exTwoSexLinearCoaleR <- compiler::cmpfun(function(dxm, dxf, FexFF, FexFM, FexMM, FexMF, 
                .a = .5:110.5, sigma = .5, maxit = 2e2, tol = 1e-15){  
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
            
                      
             # assuming r = 0
            # first guess at SRB
            SRBi <- sum(rowSums(dxM) * FexMM + rowSums(dxF) * FexFM) / 
                    sum(rowSums(dxM) * FexMF + rowSums(dxF) * FexFF)
            # derive proportions male and female
            p.m <- SRBi / (1 + SRBi)
            p.f <- 1 / (1 + SRBi)
            # guess at R0
            R0      <-   (sigma * sum(p.m * dxM * (FexMF + FexMM)) + 
                        (1 - sigma) * sum(p.f * dxF * (FexFF + FexFM))) 
            T.guess <- wmean(.a,
                    sigma * p.m * rowSums(dxM) * (FexMM + FexMF) +
                            (1 - sigma) * p.f * rowSums(dxF) * (FexFF + FexFM)
            )
            r2      <- log(R0) / T.guess
            
            # be careful to discount Fex by SRB appropriately for males / females
            # prior to specification
            # Based on Coale (1957)
            for (i in 1:maxit){ # 15 is more than enough!
                #cat(r2,i,"\n")
                r1 <- r2
                deltai <- 1 - sum(
                        sigma * rowSums(p.m * dxM %col% (1 / exp(-r1 * .a))) * (FexMM + FexMF) + # male - female   
                                (1-sigma) * rowSums(p.f * dxF %col% (1 / exp(-r1 * .a))) * (FexFF + FexFM)  # female - female
                )
                # the mean generation time self-corrects 
                # according to the error produced by the Lotka equation
                r2 <- r1 - (deltai / (T.guess - (deltai / r1))) 
                # update SRB
                SRBi <- sum(rowSums(p.m * dxM %col% (1 / exp(-r2 * .a))) * FexMM + 
                                 rowSums(p.f * dxF %col% (1 / exp(-r2 * .a))) * FexFM) / 
                        sum(rowSums(p.m * dxM %col% (1 / exp(-r2 * .a))) * FexMF + 
                                 rowSums(p.f * dxF %col% (1 / exp(-r2 * .a))) * FexFF)
                p.m <- SRBi / (1 + SRBi)
                p.f <- 1 / (1 + SRBi)
                if (abs(r2 - r1) <= tol | zapsmall(abs(deltai)) <= tol){
                    break
                }
                
            }
            if (i == maxit){
                cat("WARNING: max iterations reached, r may not be solution")
            }
            return(c(r=r2, SRB=SRBi))  

            })
exTwoSexLinearTy <- compiler::cmpfun(function(r, SRB, dxm, dxf, FexFF, FexFM, FexMM, FexMF, .a = .5:110.5, sigma = .5){
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
            p.m <- SRB / (1+SRB)
            p.f <- 1 / (1+SRB)
            wmean(.a,
                    sigma * rowSums(p.m * dxM %col% (1 / exp(-r * .a))) * (FexMM + FexMF) +
                            (1 - sigma) * rowSums(p.f*dxF %col% (1 / exp(-r * .a))) * (FexFF + FexFM)
            )
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
                              R0.0 = R0.0, R0.5 = R0.5, R0.1 = R0.1)
                        }, .Bxymf = BxymfUS, .dxm = dxmUS, .dxf = dxfUS, .Ex = ExUS))
dimnames(US) <- list(yearsUS, c("$r^{\\upsilon (\\sigma = 0)}$"  , "$r^{\\upsilon (\\sigma = .5)}$"  , "$r^{\\upsilon (\\sigma = 1)}$",
        "$T^{\\upsilon (\\sigma = 0)}$"  , "$T^{\\upsilon (\\sigma = .5)}$"  , "$T^{\\upsilon (\\sigma = 1)}$",
        "$R_0^{\\upsilon (\\sigma = 0)}$", "$R_0^{\\upsilon (\\sigma = .5)}$", "$R_0^{\\upsilon (\\sigma = 1)}$"))

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
                              R0.0 = R0.0, R0.5 = R0.5, R0.1 = R0.1)
                        }, .Bxymf = BxymfES, .dxm = dxmES, .dxf = dxfES, .Ex = ExES))

dimnames(ES) <- list(yearsES, c("$r^{\\upsilon (\\sigma = 0)}$"  , "$r^{\\upsilon (\\sigma = .5)}$"  , "$r^{\\upsilon (\\sigma = 1)}$",
                "$T^{\\upsilon (\\sigma = 0)}$"  , "$T^{\\upsilon (\\sigma = .5)}$"  , "$T^{\\upsilon (\\sigma = 1)}$",
                "$R_0^{\\upsilon (\\sigma = 0)}$", "$R_0^{\\upsilon (\\sigma = .5)}$", "$R_0^{\\upsilon (\\sigma = 1)}$"))
        
library(xtable)
print(xtable(ES, digits = c(0,4,4,4,2,2,2,3,3,3), align = c("c","c","c","c","c","c","c","c","c","c")),
        sanitize.colnames.function = identity, 
        file = "/home/triffe/git/DISS/latex/xtables/ex2sexlinearES.tex",floating=FALSE)

print(xtable(US, digits = c(0,4,4,4,2,2,2,3,3,3), align = c("c","c","c","c","c","c","c","c","c","c")),
        sanitize.colnames.function = identity, 
        file = "/home/triffe/git/DISS/latex/xtables/ex2sexlinearUS.tex",floating=FALSE)

#

pdf("/home/triffe/git/DISS/latex/Figures/exLotka2sexlinear.pdf", height = 5, width = 5)
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

US[,1] < US[,3]
ES[,1] < ES[,3]


rfUS[,1] < US[,1]
rmUS[,1] > US[,3]

rfES[,1] < ES[,1]
rmES[,1] > ES[,3]
# ---------------------------------------------
# tests
r     <- US[1, 2]







yr <- "1980"
yri     <- as.integer(yr)
ExM     <- rowSums(ExpectedDx( with(ExUS, Male[Year == yri]), dxm))
BxMM    <- rowSums(ExpectedDx( rowSums(BxymfUS[[yr]][["Bxym"]], na.rm = TRUE), dxm))
BxMF    <- rowSums(ExpectedDx( rowSums(BxymfUS[[yr]][["Bxyf"]], na.rm = TRUE), dxm))

ExF     <- rowSums(ExpectedDx( with(ExUS, Female[Year == yri]), dxf))
BxFF    <- rowSums(ExpectedDx( colSums(BxymfUS[[yr]][["Bxyf"]], na.rm = TRUE), dxf))
BxFM    <- rowSums(ExpectedDx( colSums(BxymfUS[[yr]][["Bxym"]], na.rm = TRUE), dxf))
# sex-sex-ex- specific rates:
FxMM    <- BxMM / ExM
FxFF    <- BxFF / ExF
FxMF    <- BxMF / ExM
FxFM    <- BxFM / ExF
plot(BxFM[1:90]/BxFF[1:90], type = 'l', col = "red", ylim = c(1.05,1.06))
lines(BxMM[1:90]/BxMF[1:90], col = "blue")

dxm <- dxmUS[,yr]
dxf <- dxfUS[,yr]
N               <- 111
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
.a <- .5:110.5
sigma <- .4
r <- exTwoSexLinearCoaleR(dxm = dxm, dxf = dxf, FexFF = FxFF, FexFM = FxFM, FexMM = FxMM, FexMF = FxMF,sigma = .4) 
gF <- sum((1-sigma)*rowSums(dxF %col% (1 / exp(-r * .a))) * FxFF)
bF <- sum((1-sigma)*rowSums(dxF %col% (1 / exp(-r * .a))) * FxFM)
gM <- sum(sigma*rowSums(dxM %col% (1 / exp(-r * .a))) * FxMF)
bM <- sum(sigma*rowSums(dxM %col% (1 / exp(-r * .a))) * FxMM)
# try recursion to get SRB

boys  <-(sum((1-sigma)/2*rowSums(dxF %col% (1 / exp(-r * .a))) * FxFM) + sum(sigma/2*rowSums(dxM %col% (1 / exp(-r * .a))) * FxMM))
girls <- (sum((1-sigma)/2*rowSums(dxF %col% (1 / exp(-r * .a))) * FxFF) + sum(sigma/2*rowSums(dxM %col% (1 / exp(-r * .a))) * FxMF)) 
SRBi <- boys / girls
srbvec <- vector(length=10)
for (i in 1:10){
    srbvec[i] <- SRBi
    p.f <- 1 / (1 + SRBi)
    p.m <- SRBi / (1 + SRBi)
    
    gi <- sum((1-sigma)*p.f*rowSums(dxF %col% (1 / exp(-r * .a))) * FxFF) + sum(sigma*p.m*rowSums(dxM %col% (1 / exp(-r * .a))) * FxMF)
    bi <- sum((1-sigma)*p.f*rowSums(dxF %col% (1 / exp(-r * .a))) * FxFM) + sum(sigma*p.m*rowSums(dxM %col% (1 / exp(-r * .a))) * FxMM)
    
    SRBi <- bi/gi
}
boys
bi/gi == SRBi
bi + gi
plot(srbvec, type='l')

b <- 1 / sum((SRBi / (1+SRBi)) * rowSums(dxM %col% (1 / exp(-r * .a))) +
        (1 / (1+SRBi)) * rowSums(dxF %col% (1 / exp(-r * .a))) )

g <- sum((1-sigma)/2*rowSums(dxF %col% (1 / exp(-r * .a))) * FxFF) + sum(sigma/2*rowSums(dxM %col% (1 / exp(-r * .a))) * FxMF)
b <- sum((1-sigma)/2*rowSums(dxF %col% (1 / exp(-r * .a))) * FxFM) + sum(sigma/2*rowSums(dxM %col% (1 / exp(-r * .a))) * FxMM)
g+b

sum(b * (1 / (1+SRBi)) * rowSums(dxF %col% (1 / exp(-r * .a))))+
sum(b * (SRBi / (1+SRBi)) * rowSums(dxM %col% (1 / exp(-r * .a))))

sum(b * ((SRBi / (1+SRBi)) *rowSums(dxM %col% (1 / exp(-r * .a))) + (1 / (1+SRBi)) * rowSums(dxF %col% (1 / exp(-r * .a)))))

b == bi

b / g
SRB <- b/g

bmi <- (SRBi / (1+SRBi)) * (1 / sum((sigma) * rowSums(dxM %col% (1 / exp(-r * .a))) +
                ((1 - sigma)) * rowSums(dxF %col% (1 / exp(-r * .a))) ))
bfi <- (1 / (1+SRBi)) * (1 / sum((sigma) * rowSums(dxM %col% (1 / exp(-r * .a))) +
                            ((1 - sigma)) * rowSums(dxF %col% (1 / exp(-r * .a))) ))
bm <- (SRB / (1+SRB)) * (1 / sum((sigma) * rowSums(dxM %col% (1 / exp(-r * .a))) +
                            ((1 - sigma)) * rowSums(dxF %col% (1 / exp(-r * .a))) ))

bf <- (1 / (1+SRB)) * (1 / sum((sigma) * rowSums(dxM %col% (1 / exp(-r * .a))) +
                            ((1 - sigma)) * rowSums(dxF %col% (1 / exp(-r * .a))) ))
bmi - bm
bfi - bf

(bm + bf) == bTot


bTot <- (1 / sum((sigma) * rowSums(dxM %col% (1 / exp(-r * .a))) +
                            ((1 - sigma)) * rowSums(dxF %col% (1 / exp(-r * .a))) ))
1 / sum(sigma * rowSums(dxM %col% (1 / exp(-r * .a))))
bTot == (bm + bf)

Cy <- bTot * ( rowSums(dxM %col% (1 / exp(-r * .a))) + rowSums(dxF %col% (1 / exp(-r * .a))))
sum(Cy)
sum(bTot * rowSums(dxM %col% (1 / exp(-r * .a))) )+
sum(bTot * rowSums(dxF %col% (1 / exp(-r * .a))) )



sum(bm * rowSums(dxM %col% (1 / exp(-r * .a))))/
sum(bf * rowSums(dxF %col% (1 / exp(-r * .a))))

bT <- bm + bf

sum(bT/2 * rowSums(dxM %col% (1 / exp(-r * .a))))+
        sum(bT/2 * rowSums(dxF %col% (1 / exp(-r * .a))))


