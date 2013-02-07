
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
'%row%' <- function(M,v){
    diag(v) %*% M    
}
'%col%' <- function(M,v){
    M %*% diag(v)  
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



