source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
# BxES is 0:110, years 1975:2009
BxUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
BxES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

#BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf10_65.Rdata")))
#BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf10_65.Rdata")))

# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))


yearsUS <- 1969:2009
yearsES <- 1975:2009



# differences for males and females?
diffsUSaPred <- do.call(rbind,lapply(1:(length(yearsUS)-1), function(i, .BxUS, .ExUS, .yearsUS){
                    yr <- as.character(.yearsUS[i])
                    BpredM <- Minf0(Mna0(Minf0(Mna0(.BxUS[[yr]] %row% 
                                                                    with(.ExUS, Male[Year == .yearsUS[i]]))) %row%
                                            (1 / with(.ExUS, Male[Year == .yearsUS[i + 1]]))))
                    BpredF <- Minf0(Mna0(Minf0(Mna0(.BxUS[[yr]] %col% 
                                                                    with(.ExUS, Female[Year == .yearsUS[i]]))) %col%
                                            (1 / with(.ExUS, Female[Year == .yearsUS[i + 1]]))) )
                    Diff  <- sum(BpredM) - sum(BpredF)
                    RDiff <- (2 * Diff) / (sum(BpredM) + sum(BpredF))
                    RMS   <- sqrt(mean((BpredM - BpredF) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                },  .BxUS = BxUS, .ExUS = ExUS, .yearsUS = yearsUS))
diffsESaPred <- do.call(rbind,lapply(1:(length(yearsES)-1), function(i, .BxES, .ExES, .yearsES){
                    yr <- as.character(.yearsES[i])
                    BpredM <- Minf0(Mna0(Minf0(Mna0(.BxES[[yr]] %row% 
                                                                    with(.ExES, Male[Year == .yearsES[i]]))) %row%
                                            (1 / with(.ExES, Male[Year == .yearsES[i + 1]]))))
                    BpredF <- Minf0(Mna0(Minf0(Mna0(.BxES[[yr]] %col% 
                                                                    with(.ExES, Female[Year == .yearsES[i]]))) %col%
                                            (1 / with(.ExES, Female[Year == .yearsES[i + 1]]))))
                    Diff  <- sum(BpredM) - sum(BpredF)
                    RDiff <- (2 * Diff) / (sum(BpredM) + sum(BpredF))
                    RMS   <- sqrt(mean((BpredM - BpredF) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                },  .BxES = BxES, .ExES = ExES, .yearsES = yearsES))

# (these are the year t+5 predictions)

diffsUSaPred5 <- do.call(rbind,lapply(1:(length(yearsUS)-5), function(i, .BxUS, .ExUS, .yearsUS){
                    yr <- as.character(.yearsUS[i])
                    BpredM <- Minf0(Mna0(Minf0(Mna0(.BxUS[[yr]] %row% 
                                                                    with(.ExUS, Male[Year == .yearsUS[i]]))) %row%
                                            (1 / with(.ExUS, Male[Year == .yearsUS[i + 5]]))))
                    BpredF <- Minf0(Mna0(Minf0(Mna0(.BxUS[[yr]] %col% 
                                                                    with(.ExUS, Female[Year == .yearsUS[i]]))) %col%
                                            (1 / with(.ExUS, Female[Year == .yearsUS[i + 5]]))) )
                    Diff  <- sum(BpredM) - sum(BpredF)
                    RDiff <- (2 * Diff) / (sum(BpredM) + sum(BpredF))
                    RMS   <- sqrt(mean((BpredM - BpredF) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                },  .BxUS = BxUS, .ExUS = ExUS, .yearsUS = yearsUS))
diffsESaPred5 <- do.call(rbind,lapply(1:(length(yearsES)-5), function(i, .BxES, .ExES, .yearsES){
                    yr <- as.character(.yearsES[i])
                    BpredM <- Minf0(Mna0(Minf0(Mna0(.BxES[[yr]] %row% 
                                                                    with(.ExES, Male[Year == .yearsES[i]]))) %row%
                                            (1 / with(.ExES, Male[Year == .yearsES[i + 5]]))))
                    BpredF <- Minf0(Mna0(Minf0(Mna0(.BxES[[yr]] %col% 
                                                                    with(.ExES, Female[Year == .yearsES[i]]))) %col%
                                            (1 / with(.ExES, Female[Year == .yearsES[i + 5]]))))
                    Diff  <- sum(BpredM) - sum(BpredF)
                    RDiff <- (2 * Diff) / (sum(BpredM) + sum(BpredF))
                    RMS   <- sqrt(mean((BpredM - BpredF) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                },  .BxES = BxES, .ExES = ExES, .yearsES = yearsES))

diffsUSaPred10 <- do.call(rbind,lapply(1:(length(yearsUS)-10), function(i, .BxUS, .ExUS, .yearsUS){
                    yr <- as.character(.yearsUS[i])
                    BpredM <- Minf0(Mna0(Minf0(Mna0(.BxUS[[yr]] %row% 
                                                                    with(.ExUS, Male[Year == .yearsUS[i]]))) %row%
                                            (1 / with(.ExUS, Male[Year == .yearsUS[i + 10]]))))
                    BpredF <- Minf0(Mna0(Minf0(Mna0(.BxUS[[yr]] %col% 
                                                                    with(.ExUS, Female[Year == .yearsUS[i]]))) %col%
                                            (1 / with(.ExUS, Female[Year == .yearsUS[i + 10]]))) )
                    Diff  <- sum(BpredM) - sum(BpredF)
                    RDiff <- (2 * Diff) / (sum(BpredM) + sum(BpredF))
                    RMS   <- sqrt(mean((BpredM - BpredF) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                },  .BxUS = BxUS, .ExUS = ExUS, .yearsUS = yearsUS))
diffsESaPred10 <- do.call(rbind,lapply(1:(length(yearsES)-10), function(i, .BxES, .ExES,.yearsES){
                    yr <- as.character(.yearsES[i])
                    BpredM <- Minf0(Mna0(Minf0(Mna0(.BxES[[yr]] %row% 
                                                                    with(.ExES, Male[Year == .yearsES[i]]))) %row%
                                            (1 / with(.ExES, Male[Year == .yearsES[i + 10]]))))
                    BpredF <- Minf0(Mna0(Minf0(Mna0(.BxES[[yr]] %col% 
                                                                    with(.ExES, Female[Year == .yearsES[i]]))) %col%
                                            (1 / with(.ExES, Female[Year == .yearsES[i +10]]))))
                    Diff  <- sum(BpredM) - sum(BpredF)
                    RDiff <- (2 * Diff) / (sum(BpredM) + sum(BpredF))
                    RMS   <- sqrt(mean((BpredM - BpredF) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                },  .BxES = BxES, .ExES = ExES, .yearsES = yearsES))

diffsUSaPred15 <- do.call(rbind,lapply(1:(length(yearsUS)-15), function(i, .BxUS, .ExUS, .yearsUS){
                    yr <- as.character(.yearsUS[i])
                    BpredM <- Minf0(Mna0(Minf0(Mna0(.BxUS[[yr]] %row% 
                                                                    with(.ExUS, Male[Year == .yearsUS[i]]))) %row%
                                            (1 / with(.ExUS, Male[Year == .yearsUS[i + 15]]))))
                    BpredF <- Minf0(Mna0(Minf0(Mna0(.BxUS[[yr]] %col% 
                                                                    with(.ExUS, Female[Year == .yearsUS[i]]))) %col%
                                            (1 / with(.ExUS, Female[Year == .yearsUS[i + 15]]))) )
                    Diff  <- sum(BpredM) - sum(BpredF)
                    RDiff <- (2 * Diff) / (sum(BpredM) + sum(BpredF))
                    RMS   <- sqrt(mean((BpredM - BpredF) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                },  .BxUS = BxUS, .ExUS = ExUS, .yearsUS = yearsUS))
diffsESaPred15 <- do.call(rbind,lapply(1:(length(yearsES)-15), function(i, .BxES, .ExES, .yearsES){
                    yr <- as.character(.yearsES[i])
                    BpredM <- Minf0(Mna0(Minf0(Mna0(.BxES[[yr]] %row% 
                                                                    with(.ExES, Male[Year == .yearsES[i]]))) %row%
                                            (1 / with(.ExES, Male[Year == .yearsES[i + 15]]))))
                    BpredF <- Minf0(Mna0(Minf0(Mna0(.BxES[[yr]] %col% 
                                                                    with(.ExES, Female[Year == .yearsES[i]]))) %col%
                                            (1 / with(.ExES, Female[Year == .yearsES[i + 15]]))))
                    Diff  <- sum(BpredM) - sum(BpredF)
                    RDiff <- (2 * Diff) / (sum(BpredM) + sum(BpredF))
                    RMS   <- sqrt(mean((BpredM - BpredF) ^2, na.rm = TRUE))
                    c(Diff = Diff, RDiff = RDiff, RMS = RMS)
                },  .BxES = BxES, .ExES = ExES, .yearsES = yearsES))

# -----------------------------------------------------------------------------
blues <- RColorBrewer::brewer.pal(6,"Blues")[2:5]
reds <- RColorBrewer::brewer.pal(6,"Reds")[2:5]

pdf("/home/triffe/git/DISS/latex/Figures/BirthCountDivergenceAge.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(1969:2008, diffsUSaPred[, "RDiff"], type = 'n', ylim = c(-.03,.13), xlim = c(1968, 2009),
        axes = FALSE, xlab = "", ylab = "",
        panel.first = list(rect(1968, -.03, 2010, .13, col = gray(.95), border=NA),
                abline(h = seq(-.03,.12,by=.01), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.03,.13,by=.01), seq(-.03,.13,by=.01), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10), -.03, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.04, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1962, .142, "Rel. Diff", cex = 1, xpd = TRUE, pos = 4)))
lines(1969:2008,diffsUSaPred[, "RDiff"], col = blues[1], lwd = 4)
lines(1975:2008,diffsESaPred[, "RDiff"], col = reds[1], lwd = 4)
lines(1969:(2009-5), diffsUSaPred5[, "RDiff"], col = blues[2], lwd = 3)
lines(1975:(2009-5),diffsESaPred5[, "RDiff"], col = reds[2], lwd = 3)
lines(1969:(2009-10), diffsUSaPred10[, "RDiff"], col = blues[3], lwd = 2)
lines(1975:(2009-10),diffsESaPred10[, "RDiff"], col = reds[3], lwd = 2)
lines(1969:(2009-15), diffsUSaPred15[, "RDiff"], col = blues[4], lwd = 1)
lines(1975:(2009-15),diffsESaPred15[, "RDiff"], col = reds[4], lwd = 1)
text(1969, c(diffsUSaPred[1, "RDiff"],diffsUSaPred5[1, "RDiff"],diffsUSaPred10[1, "RDiff"],diffsUSaPred15[1, "RDiff"]) +
                c(.005, 0 , - .0025, 0),
        c("US t+1", "US t+5","US t+10","US t+15"), pos = 4)
text(2009 - c(1,5,10,15) + c(-3,0,0,0), 
        c(diffsESaPred[34, "RDiff"],diffsESaPred5[30, "RDiff"],diffsESaPred10[25, "RDiff"],diffsESaPred15[20, "RDiff"]) +
                c(.01,0,0,0),
        c("ES t+1", "ES t+5","ES t+10","ES t+15"), pos = 4, xpd =TRUE)
dev.off()
