setwd("/home/triffe/git/DISS/")
source("R/UtilityFunctions.R")
source("R/MeanFunctions.R")

# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
# BxES is 0:110, years 1975:2009
BxUS <- local(get(load("Data/USbirths/USBxy0_110.Rdata")))
BxES <- local(get(load("Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

#BxymfES <- local(get(load("Data/ESbirths/ESBxymf10_65.Rdata")))
#BxymfUS <- local(get(load("Data/USbirths/USBxymf10_65.Rdata")))

# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS <- local(get(load("Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
dxmUS <- local(get(load("Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("Data/HMD_dx/dxfES.Rdata"))) 

dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)

yearsUS <- 1969:2009
yearsES <- 1975:2009

# birthsc crosstabsby ex
ExBxyAallUS <- lapply(as.character(yearsUS), function(yr, .BxUS, .dxmUS, .dxfUS){
            ExpectedDxMxFmatrix( .BxUS[[yr]], .dxmUS[,yr], .dxfUS[,yr])
        }, .BxUS = BxUS, .dxmUS = dxmUS, .dxfUS = dxfUS)
names(ExBxyAallUS) <- yearsUS
ExBxyAallES <- lapply(as.character(yearsES), function(yr, .BxES, .dxmES, .dxfES){
            ExpectedDxMxFmatrix( .BxES[[yr]], .dxmES[,yr], .dxfES[,yr])
        }, .BxES = BxES, .dxmES = dxmES, .dxfES = dxfES)
names(ExBxyAallES) <- yearsES
# exposures by ex
ExExmUS      <-lapply(as.character(yearsUS), function(yr, .ExUS, .dxmUS){
            ExpectedDx( with(.ExUS, Male[Year == as.integer(yr)]), .dxmUS[, yr])
        },.ExUS = ExUS, .dxmUS = dxmUS )
ExExfUS      <-lapply(as.character(yearsUS), function(yr, .ExUS, .dxfUS){
            ExpectedDx( with(.ExUS, Female[Year == as.integer(yr)]), .dxfUS[, yr])
        },.ExUS = ExUS, .dxfUS = dxfUS )
names(ExExfUS) <- names(ExExmUS) <-  yearsUS
ExExmES      <-lapply(as.character(yearsES), function(yr, .ExES, .dxmES){
            ExpectedDx( with(.ExES, Male[Year == as.integer(yr)]), .dxmES[, yr])
        },.ExES = ExES, .dxmES = dxmES )
ExExfES      <-lapply(as.character(yearsES), function(yr, .ExES, .dxfES){
            ExpectedDx( with(.ExES, Female[Year == as.integer(yr)]), .dxfES[, yr])
        },.ExES = ExES, .dxfES = dxfES )
names(ExExfES) <- names(ExExmES) <-  yearsES
# rates, single-sex
FxExmUS <- lapply(as.character(yearsUS), function(yr, .ExBxyAallUS, .ExExmUS){
            .ExBxyAallUS[[yr]] %row% rowSums(.ExExmUS[[yr]], na.rm = TRUE)
        },.ExBxyAallUS = ExBxyAallUS, .ExExmUS = ExExmUS)
FxExfUS <- lapply(as.character(yearsUS), function(yr, .ExBxyAallUS, .ExExfUS){
            .ExBxyAallUS[[yr]] %col% rowSums(.ExExfUS[[yr]], na.rm = TRUE)
        },.ExBxyAallUS = ExBxyAallUS, .ExExfUS = ExExfUS)
names(FxExmUS) <- names(FxExfUS) <- yearsUS
FxExmES <- lapply(as.character(yearsES), function(yr, .ExBxyAallES, .ExExmES){
            .ExBxyAallES[[yr]] %row% rowSums(.ExExmES[[yr]], na.rm = TRUE)
        },.ExBxyAallES = ExBxyAallES, .ExExmES = ExExmES)
FxExfES <- lapply(as.character(yearsES), function(yr, .ExBxyAallES, .ExExfES){
            .ExBxyAallES[[yr]] %col% rowSums(.ExExfES[[yr]], na.rm = TRUE)
        },.ExBxyAallES = ExBxyAallES, .ExExfES = ExExfES)
names(FxExmES) <- names(FxExfES) <- yearsES

# what are the male and female predictions of births at year t+1 given year t rates?
BxPredMUSt1 <- lapply(1:(length(yearsUS)-1), function(i, .FxExmUS, .ExExmUS){
            Minf0(Mna0(.FxExmUS[[i]] %row% (1 / rowSums(.ExExmUS[[i+1]]))))
        },.FxExmUS = FxExmUS, .ExExmUS = ExExmUS)
BxPredFUSt1 <- lapply(1:(length(yearsUS)-1), function(i, .FxExfUS, .ExExfUS){
            Minf0(Mna0(.FxExfUS[[i]] %col% (1 / rowSums(.ExExfUS[[i+1]]))))
        },.FxExfUS = FxExfUS, .ExExfUS = ExExfUS)
names(BxPredMUSt1) <- names(BxPredFUSt1) <- yearsUS[2:length(yearsUS)]
BxPredMESt1 <- lapply(1:(length(yearsES)-1), function(i, .FxExmES, .ExExmES){
            Minf0(Mna0(.FxExmES[[i]] %row% (1 / rowSums(.ExExmES[[i+1]]))))
        },.FxExmES = FxExmES, .ExExmES = ExExmES)
BxPredFESt1 <- lapply(1:(length(yearsES)-1), function(i, .FxExfES, .ExExfES){
            Minf0(Mna0(.FxExfES[[i]] %col% (1 / rowSums(.ExExfES[[i+1]]))))
        },.FxExfES = FxExfES, .ExExfES = ExExfES)
names(BxPredMESt1) <- names(BxPredFESt1) <- yearsES[2:length(yearsES)]

# differences for males and females?
# (these are the year t+1 predictions)
diffsUSexPred <- do.call(rbind,lapply(names(BxPredMUSt1), function(yr, .BxPredMUSt1, .BxPredFUSt1){
            Diff  <- sum(.BxPredMUSt1[[yr]]) - sum(.BxPredFUSt1[[yr]])
            RDiff <- (2 * Diff) / (sum(.BxPredMUSt1[[yr]]) + sum(.BxPredFUSt1[[yr]]))
            RMS   <- sqrt(mean((.BxPredMUSt1[[yr]] - .BxPredFUSt1[[yr]]) ^2, na.rm = TRUE))
            c(Diff = Diff, RDiff = RDiff, RMS = RMS)
        }, .BxPredMUSt1 = BxPredMUSt1, .BxPredFUSt1 = BxPredFUSt1))
diffsESexPred <- do.call(rbind,lapply(names(BxPredMESt1), function(yr, .BxPredMESt1, .BxPredFESt1){
            Diff  <- sum(.BxPredMESt1[[yr]]) - sum(.BxPredFESt1[[yr]])
            RDiff <- (2 * Diff) / (sum(.BxPredMESt1[[yr]]) + sum(.BxPredFESt1[[yr]]))
            RMS   <- sqrt(mean((.BxPredMESt1[[yr]] - .BxPredFESt1[[yr]]) ^2, na.rm = TRUE))
            c(Diff = Diff, RDiff = RDiff, RMS = RMS)
        }, .BxPredMESt1 = BxPredMESt1, .BxPredFESt1 = BxPredFESt1))
diffsUSaPred <- do.call(rbind,lapply(1:(length(yearsUS)-1), function(i, .BxUS, .ExUS){
            BpredM <- Minf0(Mna0(Minf0(Mna0(.BxUS[[i]] %row% 
                            with(.ExUS, Male[Year == yearsUS[i]]))) %row%
                                  (1 / with(.ExUS, Male[Year == yearsUS[i + 1]]))))
            BpredF <- Minf0(Mna0(Minf0(Mna0(.BxUS[[i]] %col% 
                            with(.ExUS, Female[Year == yearsUS[i]]))) %col%
                                  (1 / with(.ExUS, Female[Year == yearsUS[i + 1]]))) )
            Diff  <- sum(BpredM) - sum(BpredF)
            RDiff <- (2 * Diff) / (sum(BpredM) + sum(BpredF))
            RMS   <- sqrt(mean((BpredM - BpredF) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                },  .BxUS = BxUS, .ExUS = ExUS))
diffsESaPred <- do.call(rbind,lapply(1:(length(yearsES)-1), function(i, .BxES, .ExES){
                    BpredM <- Minf0(Mna0(Minf0(Mna0(.BxES[[i]] %row% 
                            with(.ExES, Male[Year == yearsES[i]]))) %row%
                                   (1 / with(.ExES, Male[Year == yearsES[i + 1]]))))
                    BpredF <- Minf0(Mna0(Minf0(Mna0(.BxES[[i]] %col% 
                            with(.ExES, Female[Year == yearsES[i]]))) %col%
                                   (1 / with(.ExES, Female[Year == yearsES[i + 1]]))))
                    Diff  <- sum(BpredM) - sum(BpredF)
                    RDiff <- (2 * Diff) / (sum(BpredM) + sum(BpredF))
                    RMS   <- sqrt(mean((BpredM - BpredF) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                },  .BxES = BxES, .ExES = ExES))

#plot(1970:2009, diffsUSaPred[,2], type = 'l', col = "blue", ylim = c(-.011,.011))
#lines(1970:2009, diffsUSexPred[,2], col = "red")
#lines(1976:2009, diffsESaPred[,2], col = "blue", lty = 5)
#lines(1976:2009, diffsESexPred[,2], col = "red", lty = 5)
#abline(h=0)
#abline(h = mean(diffsUSaPred[,2]),col = "blue")
#abline(h = mean(diffsUSexPred[,2]),col = "red")
#abline(h = mean(diffsESaPred[,2]),col = "blue",lty=2)
#abline(h = mean(diffsESexPred[,2]),col = "red",lty=2)

# for the US, predicting year t+1 births from year t single sex rates 
# reduces the RMS by half, and has a mean relative error closer to 0 than
# age projections. This means less divergence for ex proj than in age proj
# 50% less for the US!
# approx equal for ES.



# percentage differences are miniscule
plot(1970:2009, diffsUS[,3], type = 'l', ylim = c(0,11))
lines(1976:2009, diffsES[,3])
plot(diffsES[,2])

for (yr in 1:(length(yearsUS)-1)){
    BxPredM <- FxExmUS[[i]] %row% (1 / rowSums(ExExmUS[[i+1]]))
    BxPredF <- FxExfUS[[i]] %col% (1 / rowSums(ExExfUS[[i+1]]))
}

fields::image.plot(BxPredM-BxPredF)


head(ExBxmUS[[i+1]])
head(ExBxfUS[[i+1]])

# how do rates from year t predict births in year t+1?

#plot(unlist(lapply(FxExmUS,sum, na.rm = TRUE)))
#plot(unlist(lapply(FxExfUS,sum, na.rm = TRUE)))
#plot(unlist(lapply(FxExmES,sum, na.rm = TRUE)))
#plot(unlist(lapply(FxExfES,sum, na.rm = TRUE)))
#
RMS <-vector(length = length(yearsUS)-1)
for (i in 1:(length(yearsUS) - 1)){
    RMS[i] <- sqrt(mean(((FxExmUS[[i]] * ExBxmUS[[i+1]]) - ExBxyAallUS[[i+1]]) ^ 2))
}



plot(RMS)
Ex  <- with(ExUS, Male[Year == 2000])
A   <- rowSums(ExpectedDx(Ex, dx = dxmUS[,"2000"]))
Ex  <- with(ExUS, Male[Year == 2000])
dxi <- dxmUS[,"2000"]
N   <- length(dxi)
EDx2      <- EDx      <- matrix(0, nrow = N, ncol = N, dimnames = list(Ex = 0:(N-1), Age =  0:(N-1)))
# Population age loop
for (i in 1:N){
    # distribute each age of Populatin over death times
    EDx[1:length(dxi), i]    <- dxi
    # remove firs element and rescale
    dxi                      <- dxi[2:length(dxi)] / sum(dxi[2:length(dxi)], na.rm = TRUE)
}
dxi <- dxmUS[,"2000"]
for (i in 1:N){
    # distribute each age of Populatin over death times
    EDx2[1:length(dxi), i]    <- rev(dxi)
    # remove firs element and rescale
    dxi                       <- dxi[2:length(dxi)] / sum(dxi[2:length(dxi)], na.rm = TRUE)
}
#EDx2 <- EDx2 / rowSums(EDx2, na.rm = TRUE)

colSums(EDx2 * A)

A         <- ExpectedDx(Ex, dx = dxmUS[,"2000"])
RedistMat <- A / rowSums(A)

age <- seq(.5, 110.5)
# average age for each life expectancy:
ax  <- apply(RedistMat,1,function(x, .age){
            wmean(.age, x)}, .age = age)


# 5-year predictions compare:
BxPredMUSt5 <- lapply(1:(length(yearsUS)-5), function(i, .FxExmUS, .ExExmUS){
            Minf0(Mna0(.FxExmUS[[i]] %row% (1 / rowSums(.ExExmUS[[i+5]]))))
        },.FxExmUS = FxExmUS, .ExExmUS = ExExmUS)
BxPredFUSt5 <- lapply(1:(length(yearsUS)-5), function(i, .FxExfUS, .ExExfUS){
            Minf0(Mna0(.FxExfUS[[i]] %col% (1 / rowSums(.ExExfUS[[i+5]]))))
        },.FxExfUS = FxExfUS, .ExExfUS = ExExfUS)

names(BxPredMUSt5) <- names(BxPredFUSt5) <- yearsUS[6:length(yearsUS)]
BxPredMESt5 <- lapply(1:(length(yearsES)-5), function(i, .FxExmES, .ExExmES){
            Minf0(Mna0(.FxExmES[[i]] %row% (1 / rowSums(.ExExmES[[i+5]]))))
        },.FxExmES = FxExmES, .ExExmES = ExExmES)
BxPredFESt5 <- lapply(1:(length(yearsES)-5), function(i, .FxExfES, .ExExfES){
            Minf0(Mna0(.FxExfES[[i]] %col% (1 / rowSums(.ExExfES[[i+5]]))))
        },.FxExfES = FxExfES, .ExExfES = ExExfES)
names(BxPredMESt5) <- names(BxPredFESt5) <- yearsES[6:length(yearsES)]

# differences for males and females?
# (these are the year t+1 predictions)
diffsUSexPred5 <- do.call(rbind,lapply(names(BxPredMUSt5), function(yr, .BxPredMUSt5, .BxPredFUSt5){
                    Diff  <- sum(.BxPredMUSt5[[yr]]) - sum(.BxPredFUSt5[[yr]])
                    RDiff <- (2 * Diff) / (sum(.BxPredMUSt5[[yr]]) + sum(.BxPredFUSt5[[yr]]))
                    RMS   <- sqrt(mean((.BxPredMUSt5[[yr]] - .BxPredFUSt5[[yr]]) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                }, .BxPredMUSt5 = BxPredMUSt5, .BxPredFUSt5 = BxPredFUSt5))
diffsESexPred5 <- do.call(rbind,lapply(names(BxPredMESt5), function(yr, .BxPredMESt5, .BxPredFESt5){
                    Diff  <- sum(.BxPredMESt5[[yr]]) - sum(.BxPredFESt5[[yr]])
                    RDiff <- (2 * Diff) / (sum(.BxPredMESt5[[yr]]) + sum(.BxPredFESt5[[yr]]))
                    RMS   <- sqrt(mean((.BxPredMESt5[[yr]] - .BxPredFESt5[[yr]]) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                }, .BxPredMESt5 = BxPredMESt5, .BxPredFESt5 = BxPredFESt5))
diffsUSaPred5 <- do.call(rbind,lapply(1:(length(yearsUS)-5), function(i, .BxUS, .ExUS){
                    BpredM <- Minf0(Mna0(Minf0(Mna0(.BxUS[[i]] %row% 
                                                                    with(.ExUS, Male[Year == yearsUS[i]]))) %row%
                                            (1 / with(.ExUS, Male[Year == yearsUS[i + 5]]))))
                    BpredF <- Minf0(Mna0(Minf0(Mna0(.BxUS[[i]] %col% 
                                                                    with(.ExUS, Female[Year == yearsUS[i]]))) %col%
                                            (1 / with(.ExUS, Female[Year == yearsUS[i + 5]]))) )
                    Diff  <- sum(BpredM) - sum(BpredF)
                    RDiff <- (2 * Diff) / (sum(BpredM) + sum(BpredF))
                    RMS   <- sqrt(mean((BpredM - BpredF) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                },  .BxUS = BxUS, .ExUS = ExUS))
diffsESaPred5 <- do.call(rbind,lapply(1:(length(yearsES)-5), function(i, .BxES, .ExES){
                    BpredM <- Minf0(Mna0(Minf0(Mna0(.BxES[[i]] %row% 
                                                                    with(.ExES, Male[Year == yearsES[i]]))) %row%
                                            (1 / with(.ExES, Male[Year == yearsES[i + 5]]))))
                    BpredF <- Minf0(Mna0(Minf0(Mna0(.BxES[[i]] %col% 
                                                                    with(.ExES, Female[Year == yearsES[i]]))) %col%
                                            (1 / with(.ExES, Female[Year == yearsES[i +5]]))))
                    Diff  <- sum(BpredM) - sum(BpredF)
                    RDiff <- (2 * Diff) / (sum(BpredM) + sum(BpredF))
                    RMS   <- sqrt(mean((BpredM - BpredF) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                },  .BxES = BxES, .ExES = ExES))

#plot(1974:2009, diffsUSaPred5[,2], type = 'l', col = "blue", ylim = c(-.05,.05))
#lines(1974:2009, diffsUSexPred5[,2], col = "red")
#lines(1980:2009, diffsESaPred5[,2], col = "blue", lty = 5)
#lines(1980:2009, diffsESexPred5[,2], col = "red", lty = 5)
#abline(h=0)
#abline(h = mean(diffsUSaPred5[,2]),col = "blue")
#abline(h = mean(diffsUSexPred5[,2]),col = "red")
#abline(h = mean(diffsESaPred5[,2]),col = "blue",lty=2)
#abline(h = mean(diffsESexPred5[,2]),col = "red",lty=2)



BxPredMUSt10 <- lapply(1:(length(yearsUS)-10), function(i, .FxExmUS, .ExExmUS){
            Minf0(Mna0(.FxExmUS[[i]] %row% (1 / rowSums(.ExExmUS[[i+10]]))))
        },.FxExmUS = FxExmUS, .ExExmUS = ExExmUS)
BxPredFUSt10 <- lapply(1:(length(yearsUS)-10), function(i, .FxExfUS, .ExExfUS){
            Minf0(Mna0(.FxExfUS[[i]] %col% (1 / rowSums(.ExExfUS[[i+10]]))))
        },.FxExfUS = FxExfUS, .ExExfUS = ExExfUS)

names(BxPredMUSt10) <- names(BxPredFUSt10) <- yearsUS[11:length(yearsUS)]
BxPredMESt10 <- lapply(1:(length(yearsES)-10), function(i, .FxExmES, .ExExmES){
            Minf0(Mna0(.FxExmES[[i]] %row% (1 / rowSums(.ExExmES[[i+10]]))))
        },.FxExmES = FxExmES, .ExExmES = ExExmES)
BxPredFESt10 <- lapply(1:(length(yearsES)-10), function(i, .FxExfES, .ExExfES){
            Minf0(Mna0(.FxExfES[[i]] %col% (1 / rowSums(.ExExfES[[i+10]]))))
        },.FxExfES = FxExfES, .ExExfES = ExExfES)
names(BxPredMESt10) <- names(BxPredFESt10) <- yearsES[11:length(yearsES)]

# differences for males and females?
# (these are the year t+1 predictions)
diffsUSexPred10 <- do.call(rbind,lapply(names(BxPredMUSt10), function(yr, .BxPredMUSt10, .BxPredFUSt10){
                    Diff  <- sum(.BxPredMUSt10[[yr]]) - sum(.BxPredFUSt10[[yr]])
                    RDiff <- (2 * Diff) / (sum(.BxPredMUSt10[[yr]]) + sum(.BxPredFUSt10[[yr]]))
                    RMS   <- sqrt(mean((.BxPredMUSt10[[yr]] - .BxPredFUSt10[[yr]]) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                }, .BxPredMUSt10 = BxPredMUSt10, .BxPredFUSt10 = BxPredFUSt10))
diffsESexPred10 <- do.call(rbind,lapply(names(BxPredMESt10), function(yr, .BxPredMESt10, .BxPredFESt10){
                    Diff  <- sum(.BxPredMESt10[[yr]]) - sum(.BxPredFESt10[[yr]])
                    RDiff <- (2 * Diff) / (sum(.BxPredMESt10[[yr]]) + sum(.BxPredFESt10[[yr]]))
                    RMS   <- sqrt(mean((.BxPredMESt10[[yr]] - .BxPredFESt10[[yr]]) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                }, .BxPredMESt10 = BxPredMESt10, .BxPredFESt10 = BxPredFESt10))
diffsUSaPred10 <- do.call(rbind,lapply(1:(length(yearsUS)-10), function(i, .BxUS, .ExUS){
                    BpredM <- Minf0(Mna0(Minf0(Mna0(.BxUS[[i]] %row% 
                                                                    with(.ExUS, Male[Year == yearsUS[i]]))) %row%
                                            (1 / with(.ExUS, Male[Year == yearsUS[i + 10]]))))
                    BpredF <- Minf0(Mna0(Minf0(Mna0(.BxUS[[i]] %col% 
                                                                    with(.ExUS, Female[Year == yearsUS[i]]))) %col%
                                            (1 / with(.ExUS, Female[Year == yearsUS[i + 10]]))) )
                    Diff  <- sum(BpredM) - sum(BpredF)
                    RDiff <- (2 * Diff) / (sum(BpredM) + sum(BpredF))
                    RMS   <- sqrt(mean((BpredM - BpredF) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                },  .BxUS = BxUS, .ExUS = ExUS))
diffsESaPred10 <- do.call(rbind,lapply(1:(length(yearsES)-10), function(i, .BxES, .ExES){
                    BpredM <- Minf0(Mna0(Minf0(Mna0(.BxES[[i]] %row% 
                                                                    with(.ExES, Male[Year == yearsES[i]]))) %row%
                                            (1 / with(.ExES, Male[Year == yearsES[i + 10]]))))
                    BpredF <- Minf0(Mna0(Minf0(Mna0(.BxES[[i]] %col% 
                                                                    with(.ExES, Female[Year == yearsES[i]]))) %col%
                                            (1 / with(.ExES, Female[Year == yearsES[i +10]]))))
                    Diff  <- sum(BpredM) - sum(BpredF)
                    RDiff <- (2 * Diff) / (sum(BpredM) + sum(BpredF))
                    RMS   <- sqrt(mean((BpredM - BpredF) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                },  .BxES = BxES, .ExES = ExES))

#plot(1979:2009, diffsUSaPred10[,2], type = 'l', col = "blue", ylim = c(-.1,.1))
#lines(1979:2009, diffsUSexPred10[,2], col = "red")
#lines(1985:2009, diffsESaPred10[,2], col = "blue", lty = 10)
#lines(1985:2009, diffsESexPred10[,2], col = "red", lty = 10)
#abline(h=0)
#abline(h = mean(diffsUSaPred10[,2]),col = "blue")
#abline(h = mean(diffsUSexPred10[,2]),col = "red")
#abline(h = mean(diffsESaPred10[,2]),col = "blue",lty=2)
#abline(h = mean(diffsESexPred10[,2]),col = "red",lty=2)

# 15-year jumps
BxPredMUSt15 <- lapply(1:(length(yearsUS)-15), function(i, .FxExmUS, .ExExmUS){
            Minf0(Mna0(.FxExmUS[[i]] %row% (1 / rowSums(.ExExmUS[[i+15]]))))
        },.FxExmUS = FxExmUS, .ExExmUS = ExExmUS)
BxPredFUSt15 <- lapply(1:(length(yearsUS)-15), function(i, .FxExfUS, .ExExfUS){
            Minf0(Mna0(.FxExfUS[[i]] %col% (1 / rowSums(.ExExfUS[[i+15]]))))
        },.FxExfUS = FxExfUS, .ExExfUS = ExExfUS)

names(BxPredMUSt15) <- names(BxPredFUSt15) <- yearsUS[16:length(yearsUS)]
BxPredMESt15 <- lapply(1:(length(yearsES)-15), function(i, .FxExmES, .ExExmES){
            Minf0(Mna0(.FxExmES[[i]] %row% (1 / rowSums(.ExExmES[[i+15]]))))
        },.FxExmES = FxExmES, .ExExmES = ExExmES)
BxPredFESt15 <- lapply(1:(length(yearsES)-15), function(i, .FxExfES, .ExExfES){
            Minf0(Mna0(.FxExfES[[i]] %col% (1 / rowSums(.ExExfES[[i+15]]))))
        },.FxExfES = FxExfES, .ExExfES = ExExfES)
names(BxPredMESt15) <- names(BxPredFESt15) <- yearsES[16:length(yearsES)]

# differences for males and females?
# (these are the year t+1 predictions)
diffsUSexPred15 <- do.call(rbind,lapply(names(BxPredMUSt15), function(yr, .BxPredMUSt15, .BxPredFUSt15){
                    Diff  <- sum(.BxPredMUSt15[[yr]]) - sum(.BxPredFUSt15[[yr]])
                    RDiff <- (2 * Diff) / (sum(.BxPredMUSt15[[yr]]) + sum(.BxPredFUSt15[[yr]]))
                    RMS   <- sqrt(mean((.BxPredMUSt15[[yr]] - .BxPredFUSt15[[yr]]) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                }, .BxPredMUSt15 = BxPredMUSt15, .BxPredFUSt15 = BxPredFUSt15))
diffsESexPred15 <- do.call(rbind,lapply(names(BxPredMESt15), function(yr, .BxPredMESt15, .BxPredFESt15){
                    Diff  <- sum(.BxPredMESt15[[yr]]) - sum(.BxPredFESt15[[yr]])
                    RDiff <- (2 * Diff) / (sum(.BxPredMESt15[[yr]]) + sum(.BxPredFESt15[[yr]]))
                    RMS   <- sqrt(mean((.BxPredMESt15[[yr]] - .BxPredFESt15[[yr]]) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                }, .BxPredMESt15 = BxPredMESt15, .BxPredFESt15 = BxPredFESt15))
diffsUSaPred15 <- do.call(rbind,lapply(1:(length(yearsUS)-15), function(i, .BxUS, .ExUS){
                    BpredM <- Minf0(Mna0(Minf0(Mna0(.BxUS[[i]] %row% 
                                                                    with(.ExUS, Male[Year == yearsUS[i]]))) %row%
                                            (1 / with(.ExUS, Male[Year == yearsUS[i + 15]]))))
                    BpredF <- Minf0(Mna0(Minf0(Mna0(.BxUS[[i]] %col% 
                                                                    with(.ExUS, Female[Year == yearsUS[i]]))) %col%
                                            (1 / with(.ExUS, Female[Year == yearsUS[i + 15]]))) )
                    Diff  <- sum(BpredM) - sum(BpredF)
                    RDiff <- (2 * Diff) / (sum(BpredM) + sum(BpredF))
                    RMS   <- sqrt(mean((BpredM - BpredF) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                },  .BxUS = BxUS, .ExUS = ExUS))
diffsESaPred15 <- do.call(rbind,lapply(1:(length(yearsES)-15), function(i, .BxES, .ExES){
                    BpredM <- Minf0(Mna0(Minf0(Mna0(.BxES[[i]] %row% 
                                                                    with(.ExES, Male[Year == yearsES[i]]))) %row%
                                            (1 / with(.ExES, Male[Year == yearsES[i + 15]]))))
                    BpredF <- Minf0(Mna0(Minf0(Mna0(.BxES[[i]] %col% 
                                                                    with(.ExES, Female[Year == yearsES[i]]))) %col%
                                            (1 / with(.ExES, Female[Year == yearsES[i +15]]))))
                    Diff  <- sum(BpredM) - sum(BpredF)
                    RDiff <- (2 * Diff) / (sum(BpredM) + sum(BpredF))
                    RMS   <- sqrt(mean((BpredM - BpredF) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                },  .BxES = BxES, .ExES = ExES))

#plot(1984:2009, diffsUSaPred15[,2], type = 'l', col = "blue", ylim = c(-.15,.15))
#lines(1984:2009, diffsUSexPred15[,2], col = "red")
#lines(1990:2009, diffsESaPred15[,2], col = "blue", lty = 15)
#lines(1990:2009, diffsESexPred15[,2], col = "red", lty = 15)
#abline(h=0)
#abline(h = mean(diffsUSaPred15[,2]),col = "blue")
#abline(h = mean(diffsUSexPred15[,2]),col = "red")
#abline(h = mean(diffsESaPred15[,2]),col = "blue",lty=2)
#abline(h = mean(diffsESexPred15[,2]),col = "red",lty=2)


# -----------------------------------------------------------------------------
# can the male and female exSFR curves be related to one another in such a way
# so that one could be derived from the other?


# can also put full tables in appendix.

# Difference is male - female

MeanTableUS <- matrix(nrow = 4, ncol = 2, 
        dimnames = list(c("1-year", "5-year", "10-year", "15-year"), c("$e_x$", "Age")))
MeanTableES <- matrix(nrow = 4, ncol = 2, 
        dimnames = list(c("1-year", "5-year", "10-year", "15-year"), c("$e_x$", "Age")))
MeanTableUS[,"$e_x$"] <- c(mean(diffsUSexPred[,"RDiff"]),
        mean(diffsUSexPred5[,"RDiff"]),
        mean(diffsUSexPred10[,"RDiff"]),
        mean(diffsUSexPred15[,"RDiff"]))
MeanTableUS[,"Age"] <- c(mean(diffsUSaPred[,"RDiff"]),
        mean(diffsUSaPred5[,"RDiff"]),
        mean(diffsUSaPred10[,"RDiff"]),
        mean(diffsUSaPred15[,"RDiff"]))
MeanTableES[,"$e_x$"] <- c(mean(diffsESexPred[,"RDiff"]),
        mean(diffsESexPred5[,"RDiff"]),
        mean(diffsESexPred10[,"RDiff"]),
        mean(diffsESexPred15[,"RDiff"]))
MeanTableES[,"Age"] <- c(mean(diffsESaPred[,"RDiff"]),
        mean(diffsESaPred5[,"RDiff"]),
        mean(diffsESaPred10[,"RDiff"]),
        mean(diffsESaPred15[,"RDiff"]))

MeanAbsTableUS <- matrix(nrow = 4, ncol = 2, 
        dimnames = list(c("1-year", "5-year", "10-year", "15-year"), c("$e_x$", "Age")))
MeanAbsTableES <- matrix(nrow = 4, ncol = 2, 
        dimnames = list(c("1-year", "5-year", "10-year", "15-year"), c("$e_x$", "Age")))
MeanAbsTableUS[,"$e_x$"] <- c(mean(abs(diffsUSexPred[,"RDiff"])),
        mean(abs(diffsUSexPred5[,"RDiff"])),
        mean(abs(diffsUSexPred10[,"RDiff"])),
        mean(abs(diffsUSexPred15[,"RDiff"])))
MeanAbsTableUS[,"Age"] <- c(mean(abs(diffsUSaPred[,"RDiff"])),
        mean(abs(diffsUSaPred5[,"RDiff"])),
        mean(abs(diffsUSaPred10[,"RDiff"])),
        mean(abs(diffsUSaPred15[,"RDiff"])))
MeanAbsTableES[,"$e_x$"] <- c(mean(abs(diffsESexPred[,"RDiff"])),
        mean(abs(diffsESexPred5[,"RDiff"])),
        mean(abs(diffsESexPred10[,"RDiff"])),
        mean(abs(diffsESexPred15[,"RDiff"])))
MeanAbsTableES[,"Age"] <- c(mean(abs(diffsESaPred[,"RDiff"])),
        mean(abs(diffsESaPred5[,"RDiff"])),
        mean(abs(diffsESaPred10[,"RDiff"])),
        mean(abs(diffsESaPred15[,"RDiff"])))


library(xtable)
print(xtable(MeanTableES, digits = c(0,4,4), align = c("l","|","c","c")),
        sanitize.colnames.function = identity, 
        file = "latex/xtables/exDivergenceMeanTableES.tex",floating=FALSE)
print(xtable(MeanTableUS, digits = c(0,4,4), align = c("l","|","c","c")),
        sanitize.colnames.function = identity, 
        file = "latex/xtables/exDivergenceMeanTableUS.tex",floating=FALSE)
print(xtable(MeanAbsTableES, digits = c(0,4,4)),
        sanitize.colnames.function = identity, 
        file = "latex/xtables/exDivergenceMeanAbsTableES.tex",floating=FALSE,
        include.rownames = FALSE)
print(xtable(MeanAbsTableUS, digits = c(0,4,4)),
        sanitize.colnames.function = identity, 
        file = "latex/xtables/exDivergenceMeanAbsTableUS.tex",floating=FALSE,
        include.rownames = FALSE)

