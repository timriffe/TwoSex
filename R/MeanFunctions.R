
# Author: triffe
###############################################################################

# let's define a simple generalized mean function

# in no particular order (I did Google searches for many of these and found nothing, 
# so here's to mixing wikipedia and R!

# a caveat, most of these are designed to only take x as a vector with 2 values, 
# *even though* there may be generalizations out there for more values in the vector

# also, only a couple of these accept weights, 
# but they can all be properly weighted with very little tinkering

logorithmic.mean <- function(x){
    stopifnot(all(length(x) == 2, x >= 0))
    if (any(x == 0)){return(0)}
    if (diff(x) == 0) {return(x[1])}
    return(diff(x) / (diff(log(x))))
}

harmonic.mean <- function(x){
    stopifnot(length(x) == 2)
    2 * (prod(x) / sum(x))
}

# arithmetic-harmonic = contraharmonic-arithmetic (flipped about arithmetic mean):
# i.e. harmonic hugs the smaller number and contraharmonic hugs the bigger number
contraharmonic.mean <- function(x){
    mean(x ^ 2) / mean(x)
}

geometric.mean <- function(x){
    stopifnot(length(x) == 2)
    sqrt(prod(x))
}

identric.mean <- function(x){
    stopifnot(length(x) == 2)
    if (diff(x) == 0){
        return(x[1])
    }
    ((1 / (exp(1))) * ((x[1] ^ x[1]) / (x[2] ^ x[2])) ^ (1 / diff(rev(x))))
}

arithmeticgeometric.mean <- function(x, tol = 1e-10){
    stopifnot(length(x) == 2)
    res <- 1
    while(res > tol){
        x <- c(mean(x), sqrt(prod(x)))
        res <- abs(diff(x))
    }
    return(x[1])
}

geometricharmonic.mean <- function(x, tol = 1e-10){
    stopifnot(length(x) == 2)
    res <- 1
    while(res > tol){
        x <- c(sqrt(prod(x)), 2 / (sum(1 / x)))
        res <- abs(diff(x))
    }
    return(x[1])
}

# geometricharmonic.mean(c(10,100))
heronian.mean <- function(x){
    stopifnot(length(x) == 2)
    1 / 3 * (sum(x, sqrt(prod(x))))
}

# root mean square, also called quadratic mean:
rms.mean <- function(x){
    sqrt(mean(x ^ 2))
}

# lehmer mean can get to many of the above depending on p,
# not generalized - only vectors of length 2 are acceptable
lehmer.mean <- function(x, p, w = c(1, 1)){
    stopifnot(length(x) == 2 & length(x) == length(w))
    sum(w * x ^ p, na.rm = TRUE) / sum(w * x ^ (p - 1), na.rm = TRUE)
}

# improved stolarsky mean
stolarsky.mean <- function(x, y, p){
    # where xy must be a vector of length 2...
    xy <- c(x, y)
    # infinite p are either min or max of xy
    ifelse(is.infinite(p), ifelse(p < 0, min(xy), max(xy)),
      # equality
        # p = 0 reduces to logorithmic mean      
        ifelse(p == 0, ifelse(any(xy == 0), 0, (y - x) / (log(y) - log(x))),
          # p = 1 reduces to identric mean
          ifelse(p == 1,  ((1 / (exp(1))) * ((y ^ y) / (x ^ x)) ^ (1 / (y - x))),
            ((x ^ p - y ^ p) / (p*(x - y))) ^ (1 / (p - 1))
          )
        )
    )
}
stolarsky.mean.v <- Vectorize(compiler:::cmpfun(function(x, y, p){
    # infinite p are either min or max of xy
    xy  <- c(x, y)
    ifelse(is.infinite(p), ifelse(p < 0, min(xy), max(xy)),
            # equality
            # p = 0 reduces to logorithmic mean      
            ifelse(p == 0, ifelse(any(xy == 0), 0, (y - x) / (log(y) - log(x))),
                    # p = 1 reduces to identric mean
                    ifelse(p == 1,  ((1 / (exp(1))) * ((y ^ y) / (x ^ x)) ^ (1 / (y - x))),
                            ((x ^ p - y ^ p) / (p*(x - y))) ^ (1 / (p - 1))
                    )
            )
    )
}))
# same as above, though a bit more flexible: this takes either x and y as matrices, with a single
# p; or x and y as single numbers over a vector of p
StolarskyM <- function(x, y, p){
    if (length(x) != length(y)){
        stop("x and y must have the same length / dimensions")
    }
    if (!is.null(dim(x))){
        if (length(p) > 1){
            stop("if x and y are matrices then p must be a single number and not a vector")
        }
        M       <- stolarsky.mean.v(x,y,p)
        dim(M)  <- dim(x)
        if (!is.null(dimnames(x))){
            dimnames(M) <- dimnames(x)
        }
        return(M)
    }
    stolarsky.mean(x, y, p)
    
}
#StolarskyM(10,100,-5)

#A <- matrix(rpois(16,lambda=100),4, dimnames=list(1:4,1:4))
#B <- matrix(rpois(16,lambda=200),4, dimnames=list(1:4,1:4))
#ABm <- StolarskyM(A,B,.5)

# weighted (arithmetic) mean, used i.a. for MAC (mean age at childbearing)
wmean <- function(x, w){
    sum(w * x, na.rm = TRUE) / sum(w, na.rm = TRUE)
}
# counterpart for weighted variance
wvar  <- function(x, w){
    sum(w * ((x -  wmean(x, w)) ^2)) / (sum(w) - 1)    
}

# weighted mean of two matrices:
wmeanM2 <- compiler::cmpfun(function(M1, M2, w.frac){
            M1 * w.frac + M2 * w.frac
        })
#p.i <-  seq(-10,10,by = .1)
#xy <- c(10,100)
#
#stolarsky.mean(c(10,100), 0)
#plot(p.i, stolarsky.mean(c(10,100), p.i), type ='l')

# not certain about my implementation of the gini mean. 
# It captures many of the above according to your choice of r and s.
gini.mean <- function(x, r, s){
    stopifnot(all(length(x) == 2,x >= 0))
    if (r != s){
        return(((sum(x ^ r) / sum(x ^ s))) ^ (1 / (r - s)))
    }
    if (r == s ){
        return(exp(sum((x ^ r) * log(x)) / sum(x ^ r)))
    }
}

