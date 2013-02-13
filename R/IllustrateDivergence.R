
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
# BxES is 0:110, years 1975:2009
BxUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy10_65.Rdata")))
BxES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65
# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))

# cut ES to same dimentions ages 10:65, square (summing boundary ages)
BxES <- lapply(BxES, function(x){
            x[11, ] <- colSums(x[1:11, ])
            x[, 11] <- rowSums(x[, 1:11])
            x[66, ] <- colSums(x[66:110, ])
            x[, 66] <- rowSums(x[, 66:110])
            x[11:66,11:66]
        })

# get years for each pop
yearsES <- 1975:2009
yearsUS <- 1969:2009
ages    <- 10:65
# ASFR sample year 1975:

BxyUS <- BxUS[["1975"]]
BxyES <- BxES[["1975"]]

ExmUS <- with(ExUS, Male[Year == 1975 & Age >= 10 & Age <= 65])
ExfUS <- with(ExUS, Female[Year == 1975 & Age >= 10 & Age <= 65])
ExmES <- with(ExES, Male[Year == 1975 & Age >= 10 & Age <= 65])
ExfES <- with(ExES, Female[Year == 1975 & Age >= 10 & Age <= 65])

ASFRmUS <- rowSums(BxyUS) / ExmUS
ASFRfUS <- colSums(BxyUS) / ExfUS
ASFRmES <- rowSums(BxyES) / ExmES
ASFRfES <- colSums(BxyES) / ExfES

# max(c(ASFRmUS,ASFRfUS,ASFRmES,ASFRfES))

pdf("/home/triffe/git/DISS/latex/Figures/ASFR1975.pdf", height = 5, width = 5)
par(mar = c(3,3,2,2),xaxp = "i", yaxp = "i")
plot(ages, ASFRmUS, type = 'l', ylim = c(0, .21), xlim = c(10,66), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(10,0,66,.21,col = gray(.95), border=NA),
                abline(h = seq(0,.2,by = .02), col = "white"),
                abline(v = seq(10, 65, by = 5), col = "white"),
                text(10, seq(0, .2, by = .02),seq(0,.2, by = .02),pos = 2, cex = .8, xpd = TRUE),
                text(seq(10, 60, by = 10),0, seq(10, 60, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(35, -.01, "Age", cex = 1, pos = 1, xpd = TRUE),
                text(8,.22, "ASFR", cex = 1, xpd = TRUE)))
lines(ages, ASFRfUS, lwd = 2.5, col = gray(.5))
lines(ages, ASFRmES, lwd = 2, col = gray(.2), lty = 5)
lines(ages, ASFRfES, lwd = 2.5, col = gray(.5), lty = 5)

legend(42,.21, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
dev.off()

# --------------------------------------------------------------
# simple series of TFR
TFRES <- as.matrix(do.call(rbind,lapply(as.character(yearsES), function(yr, .BxES, .ExES){
            Exm <- with(.ExES, Male[Year == yr & Age >= 10 & Age <= 65])
            Exf <- with(.ExES, Female[Year == yr & Age >= 10 & Age <= 65])
            Bxy <- .BxES[[yr]]
            c(TFRm = sum(Bxy %row% Exm, na.rm = TRUE), TFRf = sum(Bxy %col% Exf, na.rm = TRUE))
        }, .ExES = ExES, .BxES = BxES)))
TFRUS <- as.matrix(do.call(rbind,lapply(as.character(yearsUS), function(yr, .BxUS, .ExUS){
            Exm <- with(.ExUS, Male[Year == yr & Age >= 10 & Age <= 65])
            Exf <- with(.ExUS, Female[Year == yr & Age >= 10 & Age <= 65])
            Bxy <- .BxUS[[yr]]
            c(TFRm = sum(Bxy %row% Exm, na.rm = TRUE), TFRf = sum(Bxy %col% Exf, na.rm = TRUE))
        }, .ExUS = ExUS, .BxUS = BxUS)))

pdf("/home/triffe/git/DISS/latex/Figures/TFR.pdf", height = 5, width = 5)
par(mar = c(3, 3, 2, 2),xaxp = "i", yaxp = "i")
plot(yearsUS, TFRUS[, 1], type = 'l', ylim = c(1, 2.9), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968,1,2010,2.9,col = gray(.95), border=NA),
                abline(h = seq(0,2.9,by = .2), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(1, 2.9, by = .2),seq(1, 2.9, by = .2), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),1, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.1, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1967,3, "TFR", cex = 1, xpd = TRUE)))
lines(yearsUS, TFRUS[, 2], lwd = 2.5, col = gray(.5))
lines(yearsES, TFRES[, 1], lwd = 2, col = gray(.2), lty = 5)
lines(yearsES, TFRES[, 2], lwd = 2.5, col = gray(.5), lty = 5)

legend(1990,2.85, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
dev.off()

# --------------------------------------------------------------
# simple series of 


