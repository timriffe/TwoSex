
# Author: triffe
###############################################################################

# ---------------------------------------------------------
# where Mat is an age x age matrix with age labels in rows and columns
makeBlockAgeGroups <- function(Mat, N = 2){
    # assign grouped rownames and colnames, names must be integers...
    x.col  <- as.integer(colnames(Mat)) - as.integer(colnames(Mat)) %% N
    x.row  <- as.integer(rownames(Mat)) - as.integer(rownames(Mat)) %% N
    colnames(Mat) <- x.col
    rownames(Mat) <- x.row
    # shape to long, then back to wide, summing. Can be done with tapply() but this is cleaner...
    reshape2:::acast(reshape2:::melt(Mat, varnames = c("x.row", "x.col")), x.row ~ x.col, sum, value.var = "value")
}
# ---------------------------------------------------------
# where vec is a vector of data to sum into groups and ages is of the same length
makeVectorAgeGroups <- function(vec, ages = 0:110, N = 2){
    x.new  <- ages - ages %% N
    c(tapply(vec, x.new, sum))
}
# ---------------------------------------------------------
# reshape from cohort to period 
# where lexis is either 1 or 2:
# 1: lower triangle (TL) or Dec 31st populations
# 2: upper triangle (TU) or Jan 1st populations
# watch out, if unspecified and no attribute, then defaults to 2 with no warning.
# rownames must be ages (no '+'!), and colnames cohorts.

AC2AP <- function(ACmatrix, Lexis){
    if (missing(Lexis)){
        Lexis         <- ifelse(is.null(attr(ACmatrix, "Lexis")), 2, attr(ACmatrix, "Lexis"))
    }
    # back to LDB format, but lacking Year
    longform        <- reshape2:::melt(ACmatrix, varnames = c("Age", "Cohort"), value.name = "value")
    # assume Year is t+x+1. where 1 is the year.offset. Could also be 0. Test if not sure
    longform$Year   <- longform$Cohort + longform$Age + ifelse(Lexis == 2, 1, 0)
    longform        <- longform[!is.na(longform$value), ]
    # cast back to Age x Year matrix
    APmatrix        <- reshape2:::acast(longform, Age ~ Year, value.var = "value")
    attr(APmatrix, "Lexis") <- Lexis
    APmatrix
}

# ---------------------------------------------------------
# reshape from period to cohort
# where lexis is either 1 or 2:
# 1: lower triangle (TL) or Dec 31st populations
# 2: upper triangle (TU) or Jan 1st populations
# watch out, if unspecified and no attribute, then defaults to 2 with no warning.
# rownames must be ages (no '+'!), and colnames years.
AP2AC <- function(APmatrix, Lexis){
    if (missing(Lexis)){
        Lexis         <- ifelse(is.null(attr(APmatrix, "Lexis")), 2, attr(APmatrix, "Lexis"))
    }
    # back to LDB format, but lacking Cohort
    longform        <- reshape2:::melt(APmatrix, varnames = c("Age", "Year"), value.name = "value")
    # assume cohort is t-x-1. where -1 is the cohort.offset. Could also be 0. Test if not sure
    longform$Cohort <- longform$Year - longform$Age + ifelse(Lexis == 2, -1, 0)
    # cast back to Age x Cohort matrix
    ACmatrix        <- reshape2:::acast(longform, Age ~ Cohort, value.var = "value")
    attr(ACmatrix, "Lexis") <- Lexis
    ACmatrix
}

# need to do a fair amount of transforming:

logit <- function(x){
    log(x / (1 - x))
}

expit <- function(x){
    exp(x) / (1 + exp(x))
}


# for easier matrix row and column division.
# column division is otherwise a mess with
# transposing. Will be nice to internalize 
# matrix ops, though this might not help all
# that much
'%row%' <- function(M, v){
    M / v
    #diag(v) %*% M    
}
'%col%' <- function(M, v){
    t(t(M) / v)
    # diag(v) %*% M
}

Mna0 <- function(M){
    M[is.na(M)]  <- 0
    M[is.nan(M)] <- 0
    M
}

Minf0 <- function(M){
    M[is.infinite(M)]  <- 0
    M
}
MinfNA <- function(M){
    M[is.infinite(M)]  <- NA
    M
}

# used in rm_rf_divergence.R
dx2pxLes <- function(dx){
    lx <- sum(dx) - cumsum(dx)
    px <- lx[2:length(lx)] / lx[1:(length(lx)-1)]
    px <- px
    Minf0(Mna0(px))
}
Leslie <- function(fx, px){
    fx <- fx * (1 - ((1 - px[1]) / 2)) # simple cheat for infant mort
    cbind(rbind(fx[1:length(px)],diag(px)),0)
}


# this function rounds N-dimensional arrays, while maintaining margin totals
contround2 <- function(origvalue){
    contround <-function(origvalue){
        value   <- origvalue
        newval  <- value
        for(i in 2:length(newval)){
            newval[i]<  -value[i] + (newval[i - 1] - round(newval[i - 1]))
        }
        round(newval)
    }
    # aperm() is used to transpose- it generalizes to arrays.
    # it is however nasty business to transpose back to the 
    # original config, hence this helper function.
    revaperm <- function(.array, perm){
        N <- sum(perm != seq_along(dim(.array))) - 1
        for (i in 1:N){
            .array <- aperm(.array, perm)
        }
        .array
    }
    
    # what are all the unique roundings that countround would give us?
    all.roundings <- lapply(combinat::permn(seq_along(dim(origvalue))), function(dims, .origvalue){
                revaperm(contround(aperm(.origvalue, dims)), dims)
            }, .origvalue = origvalue)
    
    # which is closest to the original cell values?
    Best <- which.min(unlist(lapply(all.roundings, function(x, .origvalue){
                                sum((c(x) - c(.origvalue)) ^ 2)
                            }, .origvalue = origvalue)))
    
    return(all.roundings[[Best]])
}

# used so far in Mitra.R
Rmomentn <- compiler::cmpfun(function(fx,Lx,x,n=0){
            sum((x^n)*fx*Lx)
        })
LotkaRCoale <- compiler::cmpfun(function(fx,Lx,x){
    # from Coale, Ansley J. (1957) A New Method for Calculating Lotka's r- the Intrinsic Rate of Growth in a Stable Population.
    # Population Studies, Vol. 11 no. 1, pp 92-94
    R0 <- Rmomentn(fx,Lx,x,0)
    # first assuming a mean generation time of 29
    ri <- log(R0)/29
    for (i in 1:15){ # 10 is more than enough!
        deltai <- sum(exp(-ri*x)*fx*Lx)-1
        # the mean generation time self-corrects 
        # according to the error produced by the Lotka equation
        ri <- ri+(deltai/(29-(deltai/ri)))
    }
    return(ri)  
})

# find a way to use this is dissertation:
ExpectedDx <- compiler::cmpfun(function(Px, dx){
    dxi      <- dx / sum(dx, na.rm = TRUE)
    N        <- length(dx)
    EDx      <- matrix(0, nrow = N, ncol = N, dimnames = list(Ex = 0:(N-1), Age =  0:(N-1)))
    # Population age loop
    for (i in 1:N){
        # distribute each age of Populatin over death times
        EDx[1:length(dxi), i]    <- Px[i] * dxi
        # remove firs element and rescale
        dxi                      <- dxi[2:length(dxi)] / sum(dxi[2:length(dxi)], na.rm = TRUE)
    }
    EDx[is.na(EDx)] <- 0
    EDx
})
# same as above but for cross-classified matrix
ExpectedDxMxFmatrix <- function(Mat, dxm, dxf){
    # first go over rows, redisitributing mother's age (in columns) by mother's ex.
    # this is a male age x female age matrix still
    FemEx <- t(apply(Mat, 1, function(.bx, .dxf){
                        rowSums(ExpectedDx(Px = .bx, dx = .dxf))
                    },.dxf = dxf))   
    # now we do the same over columns, redisitributing male age (which is in rows)
    MemEx <- t(apply(FemEx, 2, function(.bx, .dxm){
                rowSums(ExpectedDx(Px = .bx, dx = .dxm))
            }, .dxm = dxm))
    dimnames(MemEx) <- dimnames(Mat)
    MemEx
}

PyramidOutline <- function(males, females, scale = sum(c(males, females)), gap = 0, x.shift = 0, y.shift = 0, ...){
    N       <- length(males)
    Total   <- sum(c(males, females), na.rm = TRUE)
    widths  <- rep(1, N)
    age     <- c(0,cumsum(widths)[-N])
    u.age   <- age[N] + widths[N]
   
    males   <- scale * males / Total
    females <- scale * females / Total
    
    
    polygon(x = c(0, rep(females, each = 2) + 0, 0) + gap / 2 + x.shift, 
            y =  c(rep(c(age, u.age), each = 2)) + y.shift, ...)
    polygon(x = c(-0, rep(-males, each = 2) - 0, -0) - gap / 2 + x.shift, 
            y =  c(rep(c(age, u.age), each = 2)) + y.shift, ...)
}

# from my package "decompHoriuchi"
DecompContinuousOrig <- function(func,rates1,rates2,N,...){
    # number of interval jumps   
    y1          <- func(rates1,...)
    y2          <- func(rates2,...)
    d           <- rates2-rates1
    n           <- length(rates1)
    delta       <- d/N
    x <- rates1+d*matrix(rep(.5:(N-.5)/N,length(rates1)),byrow=TRUE,ncol=N)
    cc <- matrix(0,nrow=n,ncol=N)
    for (j in 1:N){
        for (i in 1:n){
            z       <- rep(0,n)
            z[i]    <- delta[i]/2
            cc[i,j] <- func((x[,j]+z),...)-func((x[,j]-z),...)
        }
    }   
    cat("Error of decomposition (in %)\ne =",100*(sum(cc)/(y2-y1)-1),"\n")
    return(rowSums(cc))
}
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
            ax[i.openage]       <- Minf0(Mna0(1 / mx[i.openage]))        
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

mx2dxHMD <- compiler::cmpfun(function(mx){
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

# sum a matrix in blocks in just one margin, not both.
ReduceDimension <- compiler::cmpfun(function(Mat, struct = 0:110, N = 5, margin){
    fac <- struct - struct %% N
    Mat2 <- apply(Mat, margin, function(x, .fac){
                tapply(x, .fac, sum)
            },.fac = fac)
    if (margin == 1){
        return(t(Mat2))
    } else {
        return(Mat2)
    }
})

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
