
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

yearsES <- 1975:2009
yearsUS <- 1969:2009
ages    <- 10:65


# cut ES to same dimentions ages 10:65, square (summing boundary ages)
BxES <- lapply(BxES, function(x){
            x[11, ] <- colSums(x[1:11, ])
            x[, 11] <- rowSums(x[, 1:11])
            x[66, ] <- colSums(x[66:110, ])
            x[, 66] <- rowSums(x[, 66:110])
            x[11:66,11:66]
        })

# get years for each pop

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

# plot it, save out
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
# simple series of NRR male and female

BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf10_65.Rdata")))
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf10_65.Rdata")))

# get Lx from HMD, divide by l0 (100000)

# these are age by year matrices:
LxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmUS.Rdata")))[11:66, ] / 1e5
LxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfUS.Rdata")))[11:66, ] / 1e5
LxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmES.Rdata")))[11:66, ] / 1e5
LxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfES.Rdata")))[11:66, ] / 1e5

# calculate male-male and female-female fertility:
FxmfUS <- lapply(as.character(yearsUS), function(yr, .BxymfUS, .ExUS){ 
          
                            Exm <- with(.ExUS, Male[Year == yr & Age >= 10 & Age <= 65])
                            Exf <- with(.ExUS, Female[Year == yr & Age >= 10 & Age <= 65])
                            
                            Bxym <- .BxymfUS[[yr]][["Bxym"]]
                            Bxyf <- .BxymfUS[[yr]][["Bxyf"]]
                            
                            list(ASFRm = rowSums(Bxym %row% Exm),
                                 ASFRf = colSums(Bxyf %col% Exf))
                        }, .ExUS = ExUS, .BxymfUS = BxymfUS)
               
names(FxmfUS) <- yearsUS
FxmfES <- lapply(as.character(yearsES), function(yr, .BxymfES, .ExES){ 
                            Exm <- with(.ExES, Male[Year == yr & Age >= 10 & Age <= 65])
                            Exf <- with(.ExES, Female[Year == yr & Age >= 10 & Age <= 65])
                            
                            Bxym <- .BxymfES[[yr]][["Bxym"]]
                            Bxyf <- .BxymfES[[yr]][["Bxyf"]]
                            
                            list(ASFRm = rowSums(Bxym %row% Exm),
                                    ASFRf = colSums(Bxyf %col% Exf))
                        }, .ExES = ExES, .BxymfES = BxymfES)
names(FxmfES) <- yearsES     
R0mfUS <- as.matrix(do.call(rbind, lapply(as.character(yearsUS), function(yr, .FxmfUS, .LxmUS, .LxfUS){
            ASFRm <- .FxmfUS[[yr]][["ASFRm"]]
            ASFRf <- .FxmfUS[[yr]][["ASFRf"]]
            Lxm   <- .LxmUS[, yr]
            Lxf   <- .LxfUS[, yr]
            c(R0m = sum(ASFRm * Lxm), R0f = sum(ASFRf * Lxf))
        },.FxmfUS = FxmfUS, .LxmUS = LxmUS, .LxfUS = LxfUS)))
rownames(R0mfUS) <- yearsUS

R0mfES <- as.matrix(do.call(rbind, lapply(as.character(yearsES), function(yr, .FxmfES, .LxmES, .LxfES){
            ASFRm <- .FxmfES[[yr]][["ASFRm"]]
            ASFRf <- .FxmfES[[yr]][["ASFRf"]]
            Lxm   <- .LxmES[, yr]
            Lxf   <- .LxfES[, yr]
            c(R0m = sum(ASFRm * Lxm), R0f = sum(ASFRf * Lxf))
        },.FxmfES = FxmfES, .LxmES = LxmES, .LxfES = LxfES)))
rownames(R0mfES) <- yearsES

R0mfUS[,1] < R0mfUS[,2]
R0mfES[,1] < R0mfES[,2]

# plot it, save out:
pdf("/home/triffe/git/DISS/latex/Figures/R0mf.pdf", height = 5, width = 5)
par(mar = c(3, 3, 2, 2),xaxp = "i", yaxp = "i")
plot(yearsUS, R0mfUS[, 1], type = 'l', ylim = c(.5, 1.45), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968,.5,2010,1.45,col = gray(.95), border=NA),
                abline(h = seq(.5,1.4,by = .1), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(.5, 1.4, by = .1),seq(.5, 1.4, by = .1), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),.5, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, .45, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1967,1.5, "TFR", cex = 1, xpd = TRUE)))
lines(yearsUS, R0mfUS[, 2], lwd = 2.5, col = gray(.5))
lines(yearsES, R0mfES[, 1], lwd = 2, col = gray(.2), lty = 5)
lines(yearsES, R0mfES[, 2], lwd = 2.5, col = gray(.5), lty = 5)

legend(1990,1.35, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
dev.off()

# -----------------------------------------------------------
# SRB?

# SRB varies over age: 1975:
wlsSRBmf <- function(Bxyl, sex = "m", ages = 10:65){
    if (sex == "m"){
        SRB <- log(MinfNA(colSums(Bxyl[["Bxym"]]) / colSums(Bxyl[["Bxyf"]])))
        return(lm(SRB ~ ages, weights = sqrt(colSums(Bxyl[["Bxym"]]) + colSums(Bxyl[["Bxyf"]]))))
    } else {
        SRB <- log(MinfNA(rowSums(Bxyl[["Bxym"]]) / rowSums(Bxyl[["Bxyf"]])))
        return(lm(SRB ~ ages, weights = sqrt(rowSums(Bxyl[["Bxym"]]) + rowSums(Bxyl[["Bxyf"]]))))
    }
}


calcSRBmf <- function(Bxyl, sex = "m"){
    if (sex == "m"){
        SRB <- log(MinfNA(colSums(Bxyl[["Bxym"]]) / colSums(Bxyl[["Bxyf"]])))
        SRB[colSums(Bxyl[["Bxym"]]) + colSums(Bxyl[["Bxyf"]]) < 100] <- NA
    } else {
        SRB <- log(MinfNA(rowSums(Bxyl[["Bxym"]]) / rowSums(Bxyl[["Bxyf"]])))
        SRB[rowSums(Bxyl[["Bxym"]]) + rowSums(Bxyl[["Bxyf"]]) < 100] <- NA
    }
    SRB
}

SRBfUS <- calcSRBmf(BxymfUS[["1975"]], "f")
SRBmUS <- calcSRBmf(BxymfUS[["1975"]], "m")
SRBfES <- calcSRBmf(BxymfES[["1975"]], "f")
SRBmES <- calcSRBmf(BxymfES[["1975"]], "m")


# log SRB x age, m, f:
pdf("/home/triffe/git/DISS/latex/Figures/SRBagemf.pdf", height = 5, width = 5)
par(mar = c(5, 2, 2, 2),xaxp = "i", yaxp = "i")
ylim <- c(-.3,.3)
plot(ages, SRBmUS, type = 'l', ylim = ylim, xlim = c(10,65), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(10,ylim[1],65,ylim[2],col = gray(.95), border=NA),
                abline(h = seq(ylim[1], ylim[2], by = .1), col = "white"),
                abline(v = seq(10, 65, by = 5), col = "white"),
                text(10, zapsmall(seq(ylim[1], ylim[2], by = .1)), zapsmall(seq(ylim[1], ylim[2], by = .1)), pos = 2, cex = .8, xpd = TRUE),
                text(seq(10, 60, by = 10),ylim[1], seq(10, 60, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(35, -.65, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(9,.65, "TFR", cex = 1, xpd = TRUE)))
abline(wlsSRBmf(BxymfUS[["1975"]], "m"), col = "red")

lines(ages, SRBfUS, lwd = 2.5, col = gray(.5))
abline(wlsSRBmf(BxymfUS[["1975"]], "f"), col = "blue")

lines(ages, SRBmES, lwd = 2, col = gray(.2), lty = 5)
abline(wlsSRBmf(BxymfES[["1975"]], "m"), col = "orange")

lines(ages, SRBfES, lwd = 2.5, col = gray(.5), lty = 5)
abline(wlsSRBmf(BxymfES[["1975"]], "f"), col = "pink")
abline(h = log(1.05), col = "red")
legend(4,-.7, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE, horiz = TRUE,
        cex = .8)

dev.off()
