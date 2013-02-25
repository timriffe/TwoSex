# ----------------------------------------------
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
xlabs   <- c("1.0%","0.8%","0.6%","0.4%","0.2%","0.0%","0.2%","0.4%","0.6%","0.8%","1.0%")


Mex75US <- rowSums(ExpectedDx(Px = with(ExUS, Male[Year == 1975]), dx = dxmUS[,"1975"]))
Fex75US <- rowSums(ExpectedDx(Px = with(ExUS, Female[Year == 1975]), dx = dxfUS[,"1975"]))
Mex09US <- rowSums(ExpectedDx(Px = with(ExUS, Male[Year == 2009]), dx = dxmUS[,"2009"]))
Fex09US <- rowSums(ExpectedDx(Px = with(ExUS, Female[Year == 2009]), dx = dxfUS[,"2009"]))
pdf("/home/triffe/git/DISS/latex/Figures/exPyramidUS.pdf", height = 5, width = 5)
par(mai = c(.6,.6,.3,.3), xaxs = "i", yaxs = "i")
plot(NULL, type = "n",axes = FALSE, xlab = "",ylab = "", xlim = c(-1, 1), ylim = c(0,111),
        panel.first = list(
                rect(-1, 0, 1, 111, col = gray(.95), border = NA),
                abline(v = seq(-1, 1, by = .2), col = "white"),
                abline(h = seq(0, 110, by = 10), col = "white"),
                text(seq(-1, 1, by = .2),0, xlabs, xpd = TRUE, pos = 1, cex = .7),
                text(-1, seq(0, 110, by = 10), seq(0, 110, by = 10), pos = 2, xpd = TRUE, cex = .7),
                text(-1.17, 116, expression(e[x]), xpd = TRUE, cex = 1),
                text(0, -12, "Percentage", xpd = TRUE, cex = 1),
                text(c(-.5, .5), 115, c("Males", "Females"), cex = .9, xpd = TRUE)
               ))
barplot(-100 * (Mex75US / sum(Mex75US + Fex75US)), border = NA, col = "#44444450", 
        add = TRUE, horiz = TRUE, space = 0, axes = FALSE)
barplot(100 * (Fex75US / sum(Mex75US + Fex75US)), border = NA, col = "#44444450", 
        add = TRUE, horiz = TRUE, space = 0, axes = FALSE)
PyramidOutline(Mex09US, Fex09US, prop = TRUE, border = gray(.2), xpd = TRUE, lwd = 1)    
text(c(-.5,-.55), c(50, 91), c("1975", "2009"), col = c("white", "black"), cex = 1.2, pos = 4)
segments(-.4, 88, -.3, 83)
dev.off()
# for Spain:
Mex75ES <- rowSums(ExpectedDx(Px = with(ExES, Male[Year == 1975]), dx = dxmES[,"1975"]))
Fex75ES <- rowSums(ExpectedDx(Px = with(ExES, Female[Year == 1975]), dx = dxfES[,"1975"]))
Mex09ES <- rowSums(ExpectedDx(Px = with(ExES, Male[Year == 2009]), dx = dxmES[,"2009"]))
Fex09ES <- rowSums(ExpectedDx(Px = with(ExES, Female[Year == 2009]), dx = dxfES[,"2009"]))

pdf("/home/triffe/git/DISS/latex/Figures/exPyramidES.pdf", height = 5, width = 5)
par(mai = c(.5,.3,.3,.3))
plot(NULL, type = "n",axes = FALSE, xlab = "",ylab = "", xlim = c(-1, 1), ylim = c(0,111),
        panel.first = list(
                rect(-1, 0, 1, 111, col = gray(.95), border = NA),
                abline(v = seq(-1, 1, by = .2), col = "white"),
                abline(h = seq(0, 110, by = 10), col = "white"),
                text(seq(-1, 1, by = .2),0, xlabs, xpd = TRUE, pos = 1, cex = .7),
                text(-1, seq(0, 110, by = 10), seq(0, 110, by = 10), pos = 2, xpd = TRUE, cex = .7),
                text(-1.17, 116, expression(e[x]), xpd = TRUE, cex = 1),
                text(0, -12, "Percentage", xpd = TRUE, cex = 1),
                text(c(-.5, .5), 115, c("Males", "Females"), cex = .9, xpd = TRUE)
        ))
barplot(-100*(Mex75ES/sum(Mex75ES+Fex75ES)), border = NA, col = "#44444450", 
        add = TRUE, horiz = TRUE, space = 0, axes = FALSE)
barplot(100*(Fex75ES/sum(Mex75ES+Fex75ES)), border = NA, col = "#44444450", 
        add = TRUE, horiz = TRUE, space = 0, axes = FALSE)
PyramidOutline(Mex09ES, Fex09ES, prop = TRUE, border = gray(.2), xpd = TRUE, lwd = 1)   
text(c(-.5,-.55), c(50, 90), c("1975", "2009"), col = c("white", "black"), cex = 1.2, pos = 4)
segments(-.27,90.5,-.17,87)
dev.off()
    
################## 
# Now ex-specific rates for 2009, male, female US, Spain
x <- 0:110


Bex09USm <- rowSums(ExpectedDx(Px = rowSums(BxUS[["2009"]]), dx = dxmUS[,"2009"]))
Bex09USf <- rowSums(ExpectedDx(Px = colSums(BxUS[["2009"]]), dx = dxfUS[,"2009"]))
Bex09ESm <- rowSums(ExpectedDx(Px = rowSums(BxES[["2009"]]), dx = dxmES[,"2009"]))
Bex09ESf <- rowSums(ExpectedDx(Px = colSums(BxES[["2009"]]), dx = dxfES[,"2009"]))


Fxex09USm <- Bex09USm / Mex09US
Fxex09USf <- Bex09USf / Fex09US
Fxex09ESm <- Bex09ESm / Mex09ES
Fxex09ESf <- Bex09ESf / Fex09ES

pdf("/home/triffe/git/DISS/latex/Figures/eSFR2009.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(x, Fxex09USm, type = 'l', ylim = c(0, .07), xlim = c(0,111), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(0,0,111,0.07,col = gray(.95), border=NA),
                abline(h = seq(0,.07,by = .01), col = "white"),
                abline(v = seq(0, 110, by = 10), col = "white"),
                text(seq(0, 110, by = 10), 0,seq(0, 110, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(0,seq(0, .07, by = .01), seq(0, .07, by = .01), pos = 2, cex = .8, xpd = TRUE),
                text(55, -.004, expression(e[x]), cex = 1, pos = 1, xpd = TRUE),
                text(-15,.076, "Fertility Rate", cex = 1, xpd = TRUE, pos = 4)))
lines(x, Fxex09USf, lwd = 2.5, col = gray(.5))
lines(x, Fxex09ESm, lwd = 2, col = gray(.2), lty = 5)
lines(x, Fxex09ESf, lwd = 2.5, col = gray(.5), lty = 5)

legend(70,.067, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
dev.off()

#################################################
# once again with females 13-50 and males 15-65:
ages <- 0:110
ExUSflim <- with(ExUS,Female[Year == 2009])
ExUSmlim <- with(ExUS, Male[Year == 2009])
ExESflim <- with(ExES, Female[Year == 2009])
ExESmlim <- with(ExES, Male[Year == 2009])

ExUSflim[ages < 13 | ages > 50] <- 0
ExUSmlim[ages < 15 | ages > 65] <- 0
ExESflim[ages < 13 | ages > 50] <- 0
ExESmlim[ages < 15 | ages > 65] <- 0

Mex09USlim <- rowSums(ExpectedDx(Px = ExUSmlim, dx = dxmUS[,"2009"]))
Fex09USlim <- rowSums(ExpectedDx(Px = ExUSflim, dx = dxfUS[,"2009"]))
Mex09ESlim <- rowSums(ExpectedDx(Px = ExESmlim, dx = dxmES[,"2009"]))
Fex09ESlim <- rowSums(ExpectedDx(Px = ExESflim, dx = dxfES[,"2009"]))


Fxex09USm <- Bex09USm / Mex09USlim
Fxex09USf <- Bex09USf / Fex09USlim
Fxex09ESm <- Bex09ESm / Mex09ESlim
Fxex09ESf <- Bex09ESf / Fex09ESlim

pdf("/home/triffe/git/DISS/latex/Figures/eSFR2009limits.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(x, Fxex09USm, type = 'l', ylim = c(0, .09), xlim = c(0,111), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(0,0,111,0.09,col = gray(.95), border=NA),
                abline(h = seq(0,.09,by = .01), col = "white"),
                abline(v = seq(0, 110, by = 10), col = "white"),
                text(seq(0, 110, by = 10), 0,seq(0, 110, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(0,seq(0, .09, by = .01), seq(0, .09, by = .01), pos = 2, cex = .8, xpd = TRUE),
                text(55, -.004, expression(e[x]), cex = 1, pos = 1, xpd = TRUE),
                text(-15,.096, "Fertility Rate", cex = 1, xpd = TRUE, pos = 4)))
lines(x, Fxex09USf, lwd = 2.5, col = gray(.5))
lines(x, Fxex09ESm, lwd = 2, col = gray(.2), lty = 5)
lines(x, Fxex09ESf, lwd = 2.5, col = gray(.5), lty = 5)

legend(70,.087, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
dev.off()

#################################################
# Surface test:

exSFRUSm <- do.call(cbind, lapply(as.character(yearsUS), function(yr, .BxUS, .ExUS, .dxmUS){
            Bex <- rowSums(ExpectedDx(Px = rowSums(Mna0(.BxUS[[yr]])), dx = .dxmUS[,yr]))
            Eex <- rowSums(ExpectedDx(Px = with(.ExUS, Male[Year == as.integer(yr)]), dx = .dxmUS[,yr]))
            Bex / Eex
        }, .BxUS = BxUS, .ExUS = ExUS, .dxmUS = dxmUS))
exSFRUSf <- do.call(cbind, lapply(as.character(yearsUS), function(yr, .BxUS, .ExUS, .dxfUS){
            Bex <- rowSums(ExpectedDx(Px = colSums(Mna0(.BxUS[[yr]])), dx = .dxfUS[,yr]))
            Eex <- rowSums(ExpectedDx(Px = with(.ExUS,Female[Year == as.integer(yr)]), dx = .dxfUS[,yr]))
            Bex / Eex
        }, .BxUS = BxUS, .ExUS = ExUS, .dxfUS = dxfUS))
colnames(exSFRUSm) <- colnames(exSFRUSf) <- yearsUS

exSFRESm <- do.call(cbind, lapply(as.character(yearsES), function(yr, .BxES, .ExES, .dxmES){
                    Bex <- rowSums(ExpectedDx(Px = rowSums(Mna0(.BxES[[yr]])), dx = .dxmES[,yr]))
                    Eex <- rowSums(ExpectedDx(Px = with(.ExES, Male[Year == as.integer(yr)]), dx = .dxmES[,yr]))
                    Bex / Eex
                }, .BxES = BxES, .ExES = ExES, .dxmES = dxmES))
exSFRESf <- do.call(cbind, lapply(as.character(yearsES), function(yr, .BxES, .ExES, .dxfES){
                    Bex <- rowSums(ExpectedDx(Px = colSums(Mna0(.BxES[[yr]])), dx = .dxfES[,yr]))
                    Eex <- rowSums(ExpectedDx(Px = with(.ExES,Female[Year == as.integer(yr)]), dx = .dxfES[,yr]))
                    Bex / Eex
                }, .BxES = BxES, .ExES = ExES, .dxfES = dxfES))
colnames(exSFRESm) <- colnames(exSFRESf) <- yearsES


colramp <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"), space = "Lab")

ages <- 10:65

#dev.new(height = 4, width = 6.5)

brks <- seq(0,.1, length.out = 51)
levs <- seq(0,.12,by = .01)
shift <- 50

# for USA:
pdf("/home/triffe/git/DISS/latex/Figures/eSFRsurfacesUS.pdf", height = 5, width = 5)
par(mai = c(.3,.3,.3,1))
# males
image(x = yearsUS + .5, y = 0:110 + .5, t(exSFRUSm), 
        xlim = c(1969, 2010 + shift), ylim = c(0, 104), zlim = c(0, .12), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1,
        xlab = "", ylab = "")
abline(v=seq(1970,2010,by=10), col = "#FFFFFF50")
abline(h=seq(0,100,by=10), col = "#FFFFFF50")
# contours
contour(x = yearsUS + .5, y = 0:110 + .5, t(exSFRUSm),
        levels = levs, labels = levs, add = TRUE)

# females
image(x = yearsUS + .5 + shift, y = 0:110 + .5, t(exSFRUSf), 
        xlim = c(1969, 2010) + shift, ylim = c(0, 104), zlim = c(0, .12), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1,
        xlab = "", ylab = "", add = TRUE)
abline(v=seq(1970,2010,by=10) + shift, col = "#FFFFFF50")
abline(h=seq(0,100,by=10), col = "#FFFFFF50")
# contours
levs <- seq(0,.1,by = .01)
contour(x = yearsUS + .5 + shift, y = 0:110 + .5, t(exSFRUSf),
        levels = levs, labels = levs, add = TRUE)

# tick labels: 
text(1969, seq(0, 100, by = 10), seq(0, 100, by = 10), pos = 2, cex = .7, xpd = TRUE)
text(1969 + shift, seq(0, 100, by = 10), seq(0, 100, by = 10), pos = 2, cex = .7, xpd = TRUE)
text(seq(1970, 2010, by = 10), 0, seq(1970, 2010, by = 10), pos = 1, cex = .7, xpd = TRUE)
text(seq(1970, 2010, by = 10) + shift, 0, seq(1970, 2010, by = 10), pos = 1, cex = .7, xpd = TRUE)
# axis labels
text(c(1965, 1965 + shift), 108, expression(e[x]), xpd = TRUE)
text(c(1990, 1990 + shift), -7, "Year", xpd = TRUE)
text(c(1990, 1990 + shift), 108, c("Male","Female"), cex = .8, xpd = TRUE)
# legend
.y  <- seq(0,104, length.out = 21)
.y2 <- seq(0,104, length.out = 11)
rect(2064, .y[1:(length(.y)-1)], 2070, .y[2:length(.y)], col = colramp(20), border = NA, xpd = TRUE)
segments(2070, .y2, 2071, .y2, xpd = TRUE)
text(2071, .y2, seq(0, .1, length.out = 11), pos = 4, xpd = TRUE, cex = .7)
dev.off()

# for Spain:
pdf("/home/triffe/git/DISS/latex/Figures/eSFRsurfacesES.pdf", height = 5, width = 5)
par(mai = c(.3,.3,.3,1))
# males
image(x = yearsES + .5, y = 0:110 + .5, t(exSFRESm), 
        xlim = c(1969, 2010 + shift), ylim = c(0, 104), zlim = c(0, .12), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1,
        xlab = "", ylab = "")
abline(v=seq(1970,2010,by=10), col = "#FFFFFF50")
abline(h=seq(0,100,by=10), col = "#FFFFFF50")
# contours
contour(x = yearsES + .5, y = 0:110 + .5, t(exSFRESm),
        levels = levs, labels = levs, add = TRUE)

# females
shift <- 50
image(x = yearsES + .5 + shift, y = 0:110 + .5, t(exSFRESf), 
        xlim = c(1969, 2010) + shift, ylim = c(0, 104), zlim = c(0, .12), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1,
        xlab = "", ylab = "", add = TRUE)
abline(v=seq(1980,2010,by=10) + shift, col = "#FFFFFF50")
abline(h=seq(0,100,by=10), col = "#FFFFFF50")
# contours
levs <- seq(0,.1,by = .01)
contour(x = yearsES + .5 + shift, y = 0:110 + .5, t(exSFRESf),
        levels = levs, labels = levs, add = TRUE)

# tick labels: 
text(1975, seq(0, 100, by = 10), seq(0, 100, by = 10), pos = 2, cex = .7, xpd = TRUE)
text(1975 + shift, seq(0, 100, by = 10), seq(0, 100, by = 10), pos = 2, cex = .7, xpd = TRUE)
text(seq(1980, 2010, by = 10), 0, seq(1980, 2010, by = 10), pos = 1, cex = .7, xpd = TRUE)
text(seq(1980, 2010, by = 10) + shift, 0, seq(1980, 2010, by = 10), pos = 1, cex = .7, xpd = TRUE)
# axis labels
text(c(1971, 1971 + shift), 108, expression(e[x]), xpd = TRUE)
text(c(1990, 1990 + shift), -7, "Year", xpd = TRUE)
text(c(1990, 1990 + shift), 108, c("Male","Female"), cex = .8, xpd = TRUE)
# legend
.y  <- seq(0,104, length.out = 21)
.y2 <- seq(0,104, length.out = 11)
rect(2064, .y[1:(length(.y)-1)], 2070, .y[2:length(.y)], col = colramp(20), border = NA, xpd = TRUE)
segments(2070, .y2, 2071, .y2, xpd = TRUE)
text(2071, .y2, seq(0, .1, length.out = 11), pos = 4, xpd = TRUE, cex = .7)
dev.off()

# ----------------------------------------------
# doesn't look good for Spain, but what about removing non-reproductive ages?
# see appendix

#----------------------------------------------------------------
# exTFR plot
exTFRmUS <- colSums(exSFRUSm, na.rm = TRUE)
exTFRfUS <- colSums(exSFRUSf, na.rm = TRUE)
exTFRmES <- colSums(exSFRESm, na.rm = TRUE)
exTFRfES <- colSums(exSFRESf, na.rm = TRUE)

pdf("/home/triffe/git/DISS/latex/Figures/exTFR.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .3, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, exTFRmUS, type = 'l', ylim = c(1.3, 3), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968, 1.3, 2010, 3,col = gray(.95), border=NA),
                abline(h = seq(1.3, 3, by = .2), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(1.3, 3., by = .2),seq(1.3, 3, by = .2), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10), 1.3, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, 1.2, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1965, 3.05, expression(e[x], "             -TFR"), cex = 1, xpd = TRUE)))
lines(yearsUS, exTFRfUS, lwd = 2.5, col = gray(.5))
lines(yearsES, exTFRmES, lwd = 2, col = gray(.2), lty = 5)
lines(yearsES, exTFRfES, lwd = 2.5, col = gray(.5), lty = 5)

legend(1995, 3, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
dev.off()
