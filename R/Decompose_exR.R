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

mxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxmUS.Rdata"))) 
mxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxfUS.Rdata"))) 
mxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxmES.Rdata"))) 
mxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxfES.Rdata"))) 


# update function to determine sex gap
# should run on mx and not dx:

#derive dx from mx
mx2dxHMD <- compiler::cmpfun(function(mx, Mna0){
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
            dx
        })

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
#					BexMM     <- rowSums(ExpectedDx(rowSums(.Bxymf[[yr]][["Bxym"]]), dxm1)) 
#					BexMF     <- rowSums(ExpectedDx(rowSums(.Bxymf[[yr]][["Bxyf"]]), dxf1))
#					.Fexm.    <- Mna0(Minf0((BexMM+BexMF) /  rowSums(ExpectedDx(.Ex$Male[.Ex$Year == as.integer(yr)], dxm1))))
#					.sigexm.  <- Mna0(Minf0(BexMM / (BexMM + BexMF)))
#					
#					
#					BexFM     <- rowSums(ExpectedDx(colSums(.Bxymf[[yr]][["Bxym"]]), dxm1)) 
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
#
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
#					BexMF     <- rowSums(ExpectedDx(rowSums(.Bxymf[[yr]][["Bxyf"]]), dxf1))
#					.Fexm.    <- Mna0(Minf0((BexMM+BexMF) /  rowSums(ExpectedDx(.Ex$Male[.Ex$Year == as.integer(yr)], dxm1))))
#					.sigexm.  <- Mna0(Minf0(BexMM / (BexMM + BexMF)))
#					
#					
#					BexFM     <- rowSums(ExpectedDx(colSums(.Bxymf[[yr]][["Bxym"]]), dxm1)) 
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















