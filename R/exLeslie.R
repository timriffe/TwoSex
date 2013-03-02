source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxES is 0:110, years 1975:2009
BxUS  <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
BxES  <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

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

yearsUS <- 1969:2009
yearsES <- 1975:2009


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
dx <- dxmUS[, "1975"]
MakeLeslie <- function(Fx, dx){
    rowSums(outer(dx,Fx[2:N], "*"))
    rowSums(ExpectedDx(Fx[2:N],dx))
   
    N <- length(Fx)
    rbind(cbind(0,diag(Fx[2:N] + 1)),0)
}
MakeLeslie <- function(Fx, dx){
    N <- length(Fx)
    cbind(0,outer(dx, Fx[2:N], "*")) + 
          rbind(cbind(0,diag(N -1)),0)
}
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

log(popbio::eigen.analysis(Lex2, zero=TRUE)$lambda1)
a <- popbio::eigen.analysis(Lex2, zero=TRUE)
names(a)
image(log(a$elasticities))
males <- a$stable.stage[1:111]
females <- a$stable.stage[112:222]
barplot(-males, space= 0, horiz = TRUE, xlim=c(-.02,.02))
barplot(females, space= 0, horiz = TRUE, add =TRUE)
image(t(log(Lex2))[222:1,])

Pi <- c(ExmUS1975, ExfUS1975)

for (i in 4000){
    Pi <- c(Lex2 %*% Pi)
}
Fst <- Pi[112:222]/ sum(Pi[112:222])
dft <- dxfUS[, "1975"] * exp(-log(a$lambda1) * (.5:110.5))
dft <- dft / sum(dft)

dft - Fst

plot(Fst)

plot(0:110, Pi[1:111] / sum(Pi[1:111]), type = 'l', col="royalblue")
lines(0:110, Pi[112:222]/ sum(Pi[112:222]), col = "orange")
lines(0:110,ExmUS1975 / sum(ExmUS1975), col = "blue")
lines(0:110, ExfUS1975 / sum(ExfUS1975), col = "red")

sum(Pi[1:111]) / sum(Pi[112:222])


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


#ExLotkaTwoMin <- function(r, Fexm, Fexf, dxm, dxf, sigma = .5, SRBm = 1.05, SRBf = 1.05, .a = .5:110.5){
#    (1 - (
#          (1 - sigma) * sum(exp(-r * .a) * (1 / (1 + SRBf))     * rowSums(outer(dxf, Fexf, "*"), na.rm = TRUE), na.rm = TRUE) + 
#          (1 - sigma) * sum(exp(-r * .a) * (SRBf / (1 + SRBf))  * rowSums(outer(dxm, Fexf, "*"), na.rm = TRUE), na.rm = TRUE) + # Female - Female
#                    # Female - Male
#          sigma       * sum(exp(-r * .a) * (SRBm / (1 + SRBm))  * rowSums(outer(dxm, Fexm, "*"), na.rm = TRUE), na.rm = TRUE) + 
#          sigma       * sum(exp(-r * .a) * (1 / (1 + SRBm))     * rowSums(outer(dxf, Fexm, "*"), na.rm = TRUE), na.rm = TRUE)
#          )# Male - Male
#    ) ^ 2
#}
ExLotkaTwoMin(TwoSexRit(Fexm = FxmUS1975, Fexf = FxfUS1975, dxm = dxmUS[, "1975"], dxf = dxfUS[, "1975"], sigma = .5,
                maxit = 1e4, tol = 1e-8, T.guess = 60, r.start = .01), 
        FxmUS1975, FxfUS1975, dxmUS[, "1975"], dxfUS[, "1975"], sigma = .5)
SRBf <- SRBm <- 1.05
r <- .0015
optimize(ExLotkaTwoMin, c(-.1,.1),Fexm = FxmUS1975, Fexf = FxfUS1975, dxm = dxmUS[, "1975"], dxf = dxfUS[, "1975"])

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

getOption("digits")
zapsmall(1e-8,6)


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
ExLotkaTwoMin(r, FxmUS1975, FxfUS1975, dxmUS[, "1975"], dxfUS[, "1975"], sigma = .8)


maxit <- 15

OneSexStructureGivenR <- function(r, dx, .a = .5:110.5){
    gvec <- exp(r * .a)
    dxg <- matrix(0,ncol=111,nrow=111)
    dxi <- rev(dx)
    for (i in 1:111){
        dxg[i:111, i] <- dxi * gvec[i]
        dxi <- dxi[1:(length(dxi)-1)]
    } 
    rev(rowSums(dxg))
}



r <- .008
dx <- dxfUS[, "1975"]
image(t(ExpectedDx(rep(1,111), dxfUS[, "1975"])))

LotkaCoaleR

barplot(-OneSexStructureGivenR(.008, dxmUS[, "1975"]), space = 0, horiz = TRUE, xlim = c(-1,1))
barplot(OneSexStructureGivenR(.008, dxfUS[, "1975"]), space = 0, horiz = TRUE, add = TRUE)



plot(rev(rowSums(rmat * dxg)))
maxit <- 15
signs <- -1
signs <- 1
T.guess <- 45
OneSexRit <- function(Fex, dx, .a = .5:110.5, 
        maxit = 1e5, tol = 1e-15, T.guess = 60, r.start = .01){  
    
    exTFR <- sum(Fex)
   
    # be careful to discount Fex by SRB appropriately for males / females
    # prior to specification
    # Based on Coale (1957)
    for (signs in c(-1,1)){
        r2 <- r.start * signs
        for (i in 1:maxit){ # 15 is more than enough!
            r1 <- r2 
            gvec <- cumsum(rev(exp(-r1 * .a)))
            dxg <- matrix(0,ncol=111,nrow=111)
            dxi <- rev(dx)
            for (i in 1:111){
                dxg[i:111, i] <- dxi * gvec[i]
                dxi <- dxi[1:(length(dxi)-1)]
            }     
           
            deltai <-  1 -  sum( rowSums(dxg,na.rm = TRUE) * Fex) 
            # the mean generation time self-corrects 
            # according to the error produced by the Lotka equation
            r2 <- r1 + (deltai / (T.guess - (deltai / r1)))
            # in case approaching from wrong side of zero    
           cat(r2, "\n")
        }
        if (round(r2,  log10(1 / tol)) == 0){
            next
        } else {
            break
        }
    }
    return(r2)  
}

.a <- .5:110.5

dx   <- dxfUS[,"1975"]
dxM  <- matrix(0,ncol=111,nrow=111)
dxi  <- rev(dx)
for (i in 1:111){
    dxM[i:111, i] <- dxi 
    dxi <- dxi[1:(length(dxi)-1)]
} 
r <- -.004
dxG <- dxM %col% (1 /  exp(-r * .a))

sum(rowSums(dxG) * Fex)

exOneSexMin <- function(r, dx, Fex, .a = .5:110.5){
    # get the overlapped / staggered dx structure
    dxM  <- matrix(0,ncol=111,nrow=111)
    dxi  <- rev(dx)
    for (i in 1:111){
        dxM[i:111, i] <- dxi 
        dxi <- dxi[1:(length(dxi)-1)]
    }     
    (1 - sum(rowSums(dxM %col% (1 /  exp(-r * .a))) * Fex)) ^ 2
}
optimize(exOneSexMin, interval = c(-.2,.2), dx = dxfUS[, "1975"], Fex = FxfUS1975 * (1/2.05))

Fex <- FxfUS1975 * (1/2.05)
OneSexRit(Fex = FxfUS1975 * (1/2.05), dx = dxfUS[, "1975"], T.guess = 68, r.start = -.01, maxit = 1e3)
OneSexRit(Fex = FxmUS1975 * (1.05/2.05), dx = dxmUS[, "1975"], T.guess = 68)

wmean(.5:110.5,sum(FxmUS1975) * dxmUS[, "1975"])

