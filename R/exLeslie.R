source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")
yearsUS <- 1969:2009
yearsES <- 1975:2009
# BxES is 0:110, years 1975:2009
BxUS  <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
BxES  <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65
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



# 1975 US
ExmUS1975 <- rowSums(ExpectedDx( with(ExUS, Male[Year == 1975]), dxmUS[, "1975"]))
BxmUS1975 <- rowSums(ExpectedDx( rowSums(BxUS[["1975"]]), dxmUS[, "1975"]))
ExfUS1975 <- rowSums(ExpectedDx( with(ExUS, Female[Year == 1975]), dxfUS[, "1975"]))
BxfUS1975 <- rowSums(ExpectedDx( colSums(BxUS[["1975"]]), dxfUS[, "1975"]))

FxmUS1975 <- BxmUS1975 / ExmUS1975
FxfUS1975 <- BxfUS1975 / ExfUS1975


getFexM <- function(yr){
    yr <- as.character(yr)
    rowSums(ExpectedDx( rowSums(BxUS[[yr]]), dxmUS[, yr])) / 
            rowSums(ExpectedDx( with(ExUS, Male[Year == as.integer(yr)]), dxmUS[, yr]))
}

plot(getFexM(2000))

sum(FxmUS1975)



plot(0:110, FxmUS1975, type = 'l', col = "blue", ylim = c(0,.07))
lines(0:110, FxfUS1975, col = "red")
Fx <- FxmUS1975
dx <- makeVectorAgeGroups(dxmUS[, "1975"], N = 5)
Fx <- makeVectorAgeGroups(FxmUS1975, N = 5)
Px <- makeVectorAgeGroups(ExmUS1975, N = 5)

MakeLeslie <- function(Fx, dx){
    N           <- length(Fx)
    Fert        <- outer(dx, Fx, "*") 
    Fert[, 1]   <- Fert[, 1] / 2 # discount fertility from 0 remaining years by 1/2
    Fert[1, ]   <- Fert[1, ] / 2 # discount fert to 0 by 1/2 (won't survive until interval) 
    Fert + rbind(cbind(0,diag(N - 1)),0)
}
xtable::xtable(L)
# assume 1/2 of fert from age 0 exposed
# and 1/2 of fert to age 0 dies



Lfert <- cbind(0,outer(dx, Fx[2:N], "*"))  
L <- MakeLeslie(Fx,dx)
dim(Lfert)
Births <- c(Lfert %*% Px)
Births

plot(seq(0,110,by=5), Px, type = "l", col = "blue", ylim = c(0,2e7))
for (i in 1:15){
    Px <-c(L %*% Px)
    lines(seq(0,110,by=5), Px , col = gray(.5))    
}
colSums(L) - 1
Fx

M <- cbind(0,outer(dx, Fx[2:N], "*")) + 
        rbind(cbind(0,diag(N -1)),0)
rowSums(M)
colSums(M)

M <- rbind(cbind(0,diag(N -1)),0)
M %*% c(M %*% c(M %*% Px))
Px <- rep(1,5)
Px <- 1:5
Fx <- c(.1,.2,.3,.2,0)

dx <- c(.2,.25,.35,.15,.05)
.2 *.05


MakeLeslie(fx)



#install.packages("popbio")
library(popbio)
EigAn <- eigen.analysis(Lex, zero=TRUE)
names(EigAn)
plot(0:110,EigAn$stable.stage)
stage.vector.plot(pop.projection(Lex, FxmUS1975, 100)$stage.vectors)
SRB <- 1.05
Lex <- MakeLeslie(getFexM(2000) * (SRB / (1+SRB)), dxmUS[, "2000"])
log(popbio::eigen.analysis(Lex, zero=TRUE)[[1]])
eigen.analysis( MakeLeslie(FxfUS1975, dxfUS[, "1975"]), zero=TRUE)[1]

ExmUS1975i <- ExmUS1975
plot(0:110, ExmUS1975i / sum(ExmUS1975i), type = 'l', ylim = c(0,.05))
cols <- paste0(gray(.5),"50")
for (i in 400){
    ExmUS1975i <- c(Lex %*% ExmUS1975i)
    lines(0:110,  ExmUS1975i / sum(ExmUS1975i), col = cols[i])
}
lines()
ExmUS19752 <- c(Lex %*% ExmUS1975i)
log(sum(ExmUS19752)/sum(ExmUS1975i))
EigAn[[1]]


sum(ExmUS1975i) - sum(ExmUS1975)

(B <- sum(ExmUS1975 * FxmUS1975))
(Deaths <- ExmUS1975[1])




(B - Deaths) - (sum(ExmUS1975i) - sum(ExmUS1975))
(ExmUS1975 * FxmUS1975)[1] # infant mortality
ExmUS1975i[1] - ExmUS1975[2]

plot(0:110, ExmUS1975i / sum(ExmUS1975i), type = 'l')
lines(0:110, ExmUS1975i / sum(ExmUS1975i), col = "red")


# 1) sigma, male weight
sig <- .5


SRBm <- 1.05
SRBf <- 1.05
MakeLeslieExTwoSex <- function(Fexm, Fexf, dxm, dxf, sigma = .5, SRBm = 1.05, SRBf = 1.05){
    N <- length(Fexm)
    
    cbind(
    # Male side (left half)
    rbind(
    # upper left block
    # male survival and male-male fert
    cbind(0, outer(dxm, (SRBm / (1 + SRBm)) * sigma * Fexm[2:N], "*")) + 
            rbind(cbind(0,diag(N - 1)),0),
    # lower left block
    # male- female fert
    cbind(0, outer(dxf, (1 / (1 + SRBm)) * sigma * Fexm[2:N], "*")) 
    ),
    
    # Female side (right half)
    rbind(
    # upper right block
    # female- male fert
    cbind(0, outer(dxm, (SRBf / (1 + SRBf)) * (1 - sigma) * Fexf[2:N], "*")), 
    
    # lower right block
    # female survival and female-female fert
    cbind(0, outer(dxf, (1 / (1 + SRBf)) * (1 - sigma) * Fexf[2:N], "*")) + 
            rbind(cbind(0, diag(N - 1)), 0)
    ))
}

Lex2 <- MakeLeslieExTwoSex(FxmUS1975, FxfUS1975, dxmUS[, "1975"], dxfUS[, "1975"], sigma = .5)


#ExLotkaMinM <- function(r, Fx, dx, SRB, a = .5:110.5){
#    (1 - sum(exp(-r * a) * (SRB / (1 + SRB)) * rowSums(outer(dx, Fx, "*"))))^2
#}
#ExLotkaMinF <- function(r, Fx, dx, SRB, a = .5:110.5){
#    (1 - sum(exp(-r * a) * (1 / (1 + SRB)) * rowSums(outer(dx, Fx, "*"))))^2
#}

outer(FxmUS1975, dxmUS[, "1975"], "*")
ExLotkaMinM(.0005,FxmUS1975,  dxmUS[, "1975"], 1.05)

optimize(ExLotkaMinM, c(-.05,.05),  Fx = getFexM(1970), dx = dxmUS[, "1970"], SRB = 1.05)
optimize(ExLotkaMinF, c(-.05,.05),  Fx = FxfUS1975, dx = dxfUS[, "1975"], SRB = 1.05)


TwoSexRit <- function(Fexm, Fexf, dxm, dxf, sigma = .5, 
        SRBm = 1.05, SRBf = 1.05, 
        .a = .5:110.5, maxit = 1e5, tol = 1e-15, T.guess = 60, r.start = .01){  
            for (signs in c(-1,1)){
                r2 <- r.start * signs
                for (i in 1:maxit){ # 15 is more than enough!
                    r1 <- r2
                    deltai <-
                            (1 - sigma) * sum(exp(-r1 * .a) * (1 / (1 + SRBf))   * rowSums(outer(dxf, Fexf, "*"), na.rm = TRUE), na.rm = TRUE) + 
                            (1 - sigma) * sum(exp(-r1 * .a) * (SRBf / (1 + SRBf)) * rowSums(outer(dxm, Fexf, "*"), na.rm = TRUE), na.rm = TRUE) + # Female - Female
                            # Female - Male
                            sigma       * sum(exp(-r1 * .a) * (SRBm / (1 + SRBm)) * rowSums(outer(dxm, Fexm, "*"), na.rm = TRUE), na.rm = TRUE) + 
                            sigma       * sum(exp(-r1 * .a) * (1 / (1 + SRBm))    * rowSums(outer(dxf, Fexm, "*"), na.rm = TRUE), na.rm = TRUE) - 1
                    
                    # the mean generation time self-corrects 
                    # according to the error produced by the Lotka equation
                    r2 <- r1 + (deltai / (T.guess - (deltai / r1)))
                    # in case approaching from wrong side of zero    
                    if (abs(r1 - r2) < tol){
                        break
                    }
                }
       
                if (round(r2,  log10(1/tol)) == 0){
                    next
                } else {
                    break
                }
                
            }
            
            cat(i, "iterations\n")
            return(r2)  
        }

(r <- TwoSexRit(Fexm = FxmUS1975, Fexf = FxfUS1975, 
                dxm = dxmUS[, "1975"], dxf = dxfUS[, "1975"], sigma = .5,
                SRBm = 1.08, SRBf = 1.07,
        maxit = 1e4, tol = 1e-8, T.guess = 60, r.start = .01))


ExLotkaTwoMin <- function(r, Fexm, Fexf, dxm, dxf, sigma = .5, SRBm = 1.05, SRBf = 1.05, .a = .5:110.5){
    (1 - (
       sum(exp(-r * .a) * 
            (
              (1 - sigma) *  # female weight
                ((1 / (1 + SRBf)) * rowSums(outer(dxf, Fexf, "*"), na.rm = TRUE) +       # female -> female
                (SRBf / (1 + SRBf)) * rowSums(outer(dxm, Fexf, "*"), na.rm = TRUE)) +   # female -> male
              sigma * # male weight
                ((SRBm / (1 + SRBm)) * rowSums(outer(dxm, Fexm, "*"), na.rm = TRUE) +    # male -> male
                (1 / (1 + SRBm)) * rowSums(outer(dxf, Fexm, "*"), na.rm = TRUE))         # male -> female
             )
           )
         )) ^ 2
}
