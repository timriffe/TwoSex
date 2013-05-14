# apendage script, not used 
# Author: triffe
###############################################################################


lambda <- .2 # search rate
Kappa <- .5  # energy assessing mates
K      <- 3 # max total energy 
dK     <- 1 # increment in K per iteration

Nm <- 100
Nf <- 100

Hmean <- function(a,b){
    (a * b * 2) / (a + b)
}
Hmean(100,120)
# maybe decreasing returns to increasing lambda?

Qm <- runif(100)
Qf <- runif(100)

IDf <- 1:100
IDm <- 1:100


(1 - abs(Qm[1] - Qf[1])) * Kappa < runif(1)

POP             <- list()
POP$AvailMales  <- data.frame(ID = 1:Nm, Q = runif(Nm), Avail = 1, K = 3)
POP$Females     <- data.frame(ID = 1:Nf, Q = runif(Nf), Avail = 1, K = 3)

maleslooking    <- sample(POP$Females$ID, Nm * lambda, replace = TRUE)
femaleslooking  <- sample(POP$Males$ID, Nf * lambda, replace = TRUE)


lambda <- .05
L <- ceiling(lambda * Nf) # where Nf is the full pool of female candidates
# i.e. Kappa, investment per mate, drops as frequency increases
Kappa <- (1 - L / (1 + L)) ^ 2
Qd <- (1 - abs(Qm[1] - Qf[sample(POP$Females$ID, L, replace = TRUE)]))

sum(replicate(1e5,any(runif(L) < Kappa * ))) / 1e5


# -------------------------------------------------------- #
# 1) standing ovation model for entering the marriage market
# 2) males and females together or in superimposed matrices?
# 2.1) superimposed- easier for competition?
# 2.2) same- less problems with world size?
# 3) in market, works as schelling process
# 3.1) if it's a match, stays put.
# 3.2) if not a match, move.
# 3.3) match determined based on probability and comparison of characteristic(s)
# 4) characteristics can bend to neighbors, partners?
# 5) dissolution also governed by probability from these characteristics?




















