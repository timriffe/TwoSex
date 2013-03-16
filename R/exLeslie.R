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
plot(yearsUS, lambdamUS,  col = "blue", ylim = c(.84,.91), pch = 19)
points(yearsUS, lambdafUS, col = "red", pch = 19)
points(yearsES, lambdamES, col = "blue", pch = 17)
points(yearsES, lambdafES, col = "red", pch = 17)
lines(lowess(lambdafES~yearsES), col = "red", lty = 2)
lines(lowess(lambdamES~yearsES), col = "blue", lty = 2)
lines(lowess(lambdafUS~yearsUS), col = "red")
lines(lowess(lambdamUS~yearsUS), col = "blue")




MakeLExOneSexProjMatrix <- function(Fx, dx, lambda){
    N       <- length(Fx)
    # discount for part of infant mortality not surviving until end of year
    dx[1]   <- dx[1] * (1 - lambda)
    
    # NxN matrix
               # fertility component
    Y       <- outer(dx, Fx, "*") + 
               # add survival element-wise
               rbind(cbind(0,diag(N - 1)),0)
       
    # reduce e0 fertility by 1/2, as only exposed for part of year
    Y[, 1]  <- Y[, 1] / 2
    # do not allow for Inf or NA values: impute 0
    Y       <- Mna0(Minf0(Y))
    # return projection matrix
    Y
}

# 1) sigma, male weight

MakeLExTwoSexProjMatrix <- function(dxm, dxf, FexFF, FexFM, FexMM, FexMF, lambdaM, lambdaF, sigma = .5){
    N <- length(dxm)
    dxm[1] <- dxm[1] * (1 - lambdaM)
    dxf[1] <- dxf[1] * (1 - lambdaF)
    
    # define matrix, then discount first column for year t mortality
    Y <-
    cbind(
    # --------------------------------------------
      # Male side (left half)
      rbind(
        
        # upper left block
        # male-male fert
        sigma * outer(dxm, FexMM, "*") + 
        # add element-wise
        # male survival 
        rbind(cbind(0, diag(N - 1)),0),
        
        # lower left block
        # male- female fert
        sigma * outer(dxf, FexMF, "*") 
      ),
    # --------------------------------------------
      # Female side (right half)
        rbind(
        # upper right block
        # female- male fert
        (1 - sigma) * outer(dxm, FexFM, "*"), 
    
        # lower right block
        # female survival and female-female fert
        (1 - sigma) * outer(dxf, FexFF, "*") + 
        rbind(cbind(0, diag(N - 1)), 0)
      )
    )
    # discount fertility of those dying in year t (column 1) by half.
    Y[,1] <- Y[,1]/2
    
    # do not allow for Inf or NA values: impute 0
    Y     <- Mna0(Minf0(Y))
    # return projection matrix
    Y
}

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

ev <- eigen(Y_1sexES[[1]][[1]])
ev$values
Mod(ev$values)
lmax<-which.max(Re(ev$values))
lmax
Re(ev$values)[lmax]
log(max(Re(ev$values)))
names(ev)
plot(Re(ev$vectors)[,lmax])


rES1sex <- do.call(rbind,lapply(Y_1sexES, function(x){
            c(r.f = log(max(Re(eigen(x[["Yf"]])$values))),
            r.m = log(max(Re(eigen(x[["Ym"]])$values))))
        }))
rUS1sex <- do.call(rbind,lapply(Y_1sexUS, function(x){
                    c(r.f = log(max(Re(eigen(x[["Yf"]])$values))),
                            r.m = log(max(Re(eigen(x[["Ym"]])$values))))
                }))

plot(yearsUS, rUS1sex[,"r.f"], type = 'l', col = "red", ylim = c(-.02,.02))
lines(yearsUS, rUS1sex[,"r.m"], col = "blue")
lines(yearsES, rES1sex[,"r.f"], col = "red",lty=2)
lines(yearsES, rES1sex[,"r.m"], col = "blue",lty=2)

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
        .lambdaf = lambdafUS, .lambdam = lambdafES, .sigma= .5)

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



# exercise: measure time to one sex is double size of other:
names(Y_1sexUS[[1]])
.Y_1sex <- Y_1sexUS
yr <- "1975"
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
mean(doubleES)
mean(doubleUS)




        
        
        
        
        