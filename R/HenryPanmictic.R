# this code was abandoned because I wanted data in single ages, yet
# everyone says the panmictic circles should be done with max 6 circles
# it was unclear how to proceed

setwd("/home/triffe/git/DISS/")
source("R/UtilityFunctions.R")
source("R/MeanFunctions.R")

# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
# BxES is 0:110, years 1975:2009
BxUS <- local(get(load("Data/USbirths/USBxy0_110.Rdata")))
BxES <- local(get(load("Data/ESbirths/ESBxy.Rdata")))

BxymfUS <- local(get(load("Data/USbirths/USBxymf0_110.Rdata")))
BxymfES <- local(get(load("Data/ESbirths/ESBxymf.Rdata")))
# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS <- local(get(load("Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("Data/Exposures/ESexp.Rdata")))
LxmUS <- local(get(load("Data/HMD_Lx/LxmUS.Rdata"))) / 1e5
LxfUS <- local(get(load("Data/HMD_Lx/LxfUS.Rdata"))) / 1e5
LxmES <- local(get(load("Data/HMD_Lx/LxmES.Rdata"))) / 1e5
LxfES <- local(get(load("Data/HMD_Lx/LxfES.Rdata"))) / 1e5
yearsUS <- 1969:2009
yearsES <- 1975:2009

# -------------------------------------------------------


# matrix of marriage counts with female ages in columns and male ages in rows
MAR <- matrix(c(4145,24435,8140,1865,1655,54515,45010,15030,80,6735,20870,19530,5,920,5435,42470),ncol=4)
rownames(MAR) <- colnames(MAR)  <- c("15-19","20-24","25-29","30-60")
# unmarried males and females at start of the year
unMARf  <- c(254876,147705,61804,415497)
unMARm  <- c(265755,199437,114251,429655)
PanmicticRates <- function(MAR,unMARf,unMARm){
    # allocate P, array of stacked component counts
    P <- array(0,dim=c(dim(MAR),min(dim(MAR))))
    Pi <- Ri <- MAR
    
    # scale residuals by proportions of submatrices going iteratively down the diagonal
    # i <- 1
    for (i in 1:(min(dim(MAR))-1)){
        P[,,i]  <-  Pi
        P[(i+1):nrow(P),(i+1):ncol(P),i] <- Pi[(i+1):nrow(P),i]%o%Pi[i,(i+1):ncol(P)]/Pi[i,i]
        Pi <- Ri <- Ri - P[,,i]
    }
    
    # component counts, stacked
    P[,,min(dim(MAR))] <- Pi
    
    # component rates
    Mcomp <- apply(P,3,rowSums)/unMARm
    Fcomp <- apply(P,3,colSums)/unMARf
    rownames(Mcomp) <- rownames(Fcomp) <- rownames(MAR)
    colnames(Mcomp) <- colnames(Fcomp) <- paste(1:ncol(Mcomp),"comp",sep="")
    
    return(list(Mcomp=Mcomp,Fcomp=Fcomp))
}
# remember females columns, males rows
rates <- PanmicticRates(MAR,unMARf,unMARm)
Mr <- rates[["Mcomp"]]
Fr <- rates[["Fcomp"]]

# hmmm
rowSums(Fr * unMARf ) - colSums(MAR)
rowSums(unMARm * Mr) - rowSums(MAR)
176.6 + 75.1+27.7+28.1
170.9+68.3+23.1+27.8

Pma1 <- with(ExUS, Male[Year == 1969])
Pfa1 <- with(ExUS, Female[Year == 1969])

N <- 20
Mat <- BxymfUS[[1]][[1]]
PanmicticDecomp <- function(Mat, Pma1, Pfa1, Pma2, Pma2, Age.step = 3, Circle.width = Age.step * 10, .a = 0:110){
    
    dims <- dim(Mat)
    if(dims[1] != dims[2]){
        stop("need a square matrix, sorry, haven't generalized yet")
    }
    
    chunks <- .a - .a %% Age.step
    u.chunks <- unique(chunks)
    P <- array(0,dim=c(dim(Mat),min(dim(Mat))))
    Ri <- Pi <- Mat
    # i<- 5
    men <- 10
    jump <- 3
    widt <- 30
    #
    for (i in 1:10){
        P[,,i]  <-  Pi
        P[(jump*i+1):(jump*i+widt),(jump*i+1):(jump*i+widt),i] <- 
                Minf0(Mna0(
                                rowSums(Pi[(jump*i+1):(jump*i+widt),(jump*i+1):(jump*i+jump)]) %o% 
                                colSums(Pi[(jump*i+1):(jump*i+jump),(jump*i+1):(jump*i+widt)]) / 
                                   sum(Pi[(jump*i+1):(jump*i+jump),(jump*i+1):(jump*i+jump)])
                          ))
        Pi <- Ri <- Ri - P[,,i]
    }
    for (i in 1:10){
        P[,,i]  <-  Pi
        P[(jump*i+1):(jump*i+widt),(jump*i+1):(jump*i+widt),i] <- 
                Minf0(Mna0(
                                Pi[(jump*i+1):(jump*i+widt),(jump*i+1)] %o% 
                                        Pi[(jump*i+1),(jump*i+1):(jump*i+widt)] / 
                                        sum(Pi[(jump*i+1),(jump*i+1)])
                        ))
        Pi <- Ri <- Ri - P[,,i]
    }
    Mcomp <- Minf0(Mna0(apply(P,3,rowSums)/Pma1))
    Fcomp <- Minf0(Mna0(apply(P,3,colSums)/Pfa1))
}
dim(Pi)
sum(Fcomp)
image(Mcomp)
hist(apply(P,3,rowSums))
image(log(Fcomp))
ExpectedValue <- function(row.mar, col.mar){
    Minf0(Mna0(outer(row.mar, col.mar, "*") / ((sum(row.mar) + sum(col.mar)) / 2)))
}
HenryE <- function(row.mar, col.mar){
    Minf0(Mna0(outer(row.mar, col.mar, "*") / row.mar[1]))
}
HenryE2 <- function(miniMat){
    Minf0(Mna0(outer(miniMat[,1],miniMat[1,], "*") / miniMat[1,1]))
}
HenryE(Mat[20:40,20],Mat[20,20:40])

ExpectedValue(Mat[20:40,20],Mat[20,20:40])

blocks <-


