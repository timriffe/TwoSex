# appendix, showing different combinations of reproductive spans:
# source("/home/triffe/git/DISS/R/Parents_exAppendix.R")
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# (1969 - 2010)
# BxES is 0:110, years 1975:2009
BxUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
BxES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

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

# ---------------------------------------
# 15-55 repro spans for both males and females
{
ages <- 0:110
fkeep <- as.integer(ages >= 15 & ages <= 55)
mkeep <- as.integer(ages >= 15 & ages <= 55)
exSFRUSmlim <- Minf0(do.call(cbind, lapply(as.character(yearsUS), function(yr, .BxUS, .ExUS, .dxmUS){
                            Bx <-  rowSums(Mna0(.BxUS[[yr]]))
                            Ex <-  Mna0(with(.ExUS, Male[Year == as.integer(yr)])) * mkeep
                            Bex <- rowSums(ExpectedDx(Px = Bx, dx = .dxmUS[,yr]))
                            Eex <- rowSums(ExpectedDx(Px = Ex, dx = .dxmUS[,yr]))
                            Bex / Eex
                        }, .BxUS = BxUS, .ExUS = ExUS, .dxmUS = dxmUS)))
exSFRUSflim <- Minf0(do.call(cbind, lapply(as.character(yearsUS), function(yr, .BxUS, .ExUS, .dxfUS){
                            Bx  <- colSums(Mna0(.BxUS[[yr]]))
                            Ex  <- Mna0(with(.ExUS, Female[Year == as.integer(yr)])) * fkeep
                            Bex <- rowSums(ExpectedDx(Px = Bx, dx = .dxfUS[,yr]))
                            Eex <- rowSums(ExpectedDx(Px = Ex, dx = .dxfUS[,yr]))
                            Bex / Eex
                        }, .BxUS = BxUS, .ExUS = ExUS, .dxfUS = dxfUS)))
colnames(exSFRUSmlim) <- colnames(exSFRUSflim) <- yearsUS

exSFRESmlim <- Minf0(do.call(cbind, lapply(as.character(yearsES), function(yr, .BxES, .ExES, .dxmES){
                            Bx  <- rowSums(Mna0(.BxES[[yr]]))
                            Ex  <- Mna0(with(.ExES, Male[Year == as.integer(yr)])) * mkeep
                            Bex <- rowSums(ExpectedDx(Px = Bx, dx = .dxmES[,yr]))
                            Eex <- rowSums(ExpectedDx(Px = Ex, dx = .dxmES[,yr]))
                            Bex / Eex
                        }, .BxES = BxES, .ExES = ExES, .dxmES = dxmES)))
exSFRESflim <- Minf0(do.call(cbind, lapply(as.character(yearsES), function(yr, .BxES, .ExES, .dxfES){
                            Bx  <- colSums(Mna0(.BxES[[yr]]))
                            Ex  <- Mna0(with(.ExES, Female[Year == as.integer(yr)])) * fkeep
                            Bex <- rowSums(ExpectedDx(Px = Bx, dx = .dxfES[,yr]))
                            Eex <- rowSums(ExpectedDx(Px = Ex, dx = .dxfES[,yr]))
                            Bex / Eex
                        }, .BxES = BxES, .ExES = ExES, .dxfES = dxfES)))
colnames(exSFRESmlim) <- colnames(exSFRESflim) <- yearsES


colramp <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"), space = "Lab")
brks <- seq(0,.15, length.out = 51)
levs <- seq(0,.15,by = .01)
shift <- 50

# for USA:
pdf("/home/triffe/git/DISS/latex/Figures/eSFRsurfacesUSlim1555.pdf", height = 5, width = 5)
par(mai = c(.3,.3,.3,1))
# males
image(x = yearsUS + .5, y = 0:110 + .5, t(exSFRUSmlim), 
        xlim = c(1969, 2010 + shift), ylim = c(0, 104), zlim = c(0, .15), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1,
        xlab = "", ylab = "")
abline(v=seq(1970,2010,by=10), col = "#FFFFFF50")
abline(h=seq(0,100,by=10), col = "#FFFFFF50")
# contours
contour(x = yearsUS + .5, y = 0:110 + .5, t(exSFRUSmlim),
        levels = levs, labels = levs, add = TRUE)

# females
image(x = yearsUS + .5 + shift, y = 0:110 + .5, t(exSFRUSflim), 
        xlim = c(1969, 2010) + shift, ylim = c(0, 104), zlim = c(0, .15), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1,
        xlab = "", ylab = "", add = TRUE)
abline(v=seq(1970,2010,by=10) + shift, col = "#FFFFFF50")
abline(h=seq(0,100,by=10), col = "#FFFFFF50")
# contours
contour(x = yearsUS + .5 + shift, y = 0:110 + .5, t(exSFRUSflim),
        levels = levs, labels = levs, add = TRUE)

# tick labels: 
text(1969, seq(0, 100, by = 10), seq(0, 100, by = 10), pos = 2, cex = .7, xpd = TRUE)
text(1969 + shift, seq(0, 100, by = 10), seq(0, 100, by = 10), pos = 2, cex = .7, xpd = TRUE)
text(seq(1970, 2010, by = 10), 0, seq(1970, 2010, by = 10), pos = 1, cex = .7, xpd = TRUE)
text(seq(1970, 2010, by = 10) + shift, 0, seq(1970, 2010, by = 10), pos = 1, cex = .7, xpd = TRUE)
# axis labels
text(c(1965, 1965 + shift), 108, expression(e[x]), xpd = TRUE)
text(c(1990, 1990 + shift), -7, "Year", xpd = TRUE)
text(c(1990, 1990 + shift), 105, c("Male","Female"), cex = .8, xpd = TRUE)
# legend
.y  <- seq(0,104, length.out = 31)
.y2 <- seq(0,104, length.out = 16)
rect(2064, .y[1:(length(.y)-1)], 2070, .y[2:length(.y)], col = colramp(30), border = NA, xpd = TRUE)
segments(2070, .y2, 2071, .y2, xpd = TRUE)
text(2071, .y2, seq(0, .15, length.out = 16), pos = 4, xpd = TRUE, cex = .7)
dev.off()

# for Spain:
pdf("/home/triffe/git/DISS/latex/Figures/eSFRsurfacesESlim1555.pdf", height = 5, width = 5)
par(mai = c(.3,.3,.3,1))
# males
image(x = yearsES + .5, y = 0:110 + .5, t(exSFRESmlim), 
        xlim = c(1969, 2010 + shift), ylim = c(0, 104), zlim = c(0, .15), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1,
        xlab = "", ylab = "")
abline(v=seq(1970,2010,by=10), col = "#FFFFFF50")
abline(h=seq(0,100,by=10), col = "#FFFFFF50")
# contours
contour(x = yearsES + .5, y = 0:110 + .5, t(exSFRESmlim),
        levels = levs, labels = levs, add = TRUE)

# females
image(x = yearsES + .5 + shift, y = 0:110 + .5, t(exSFRESflim), 
        xlim = c(1969, 2010) + shift, ylim = c(0, 104), zlim = c(0, .12), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1,
        xlab = "", ylab = "", add = TRUE)
abline(v=seq(1980,2010,by=10) + shift, col = "#FFFFFF50")
abline(h=seq(0,100,by=10), col = "#FFFFFF50")
# contours
contour(x = yearsES + .5 + shift, y = 0:110 + .5, t(exSFRESflim),
        levels = levs, labels = levs, add = TRUE)

# tick labels: 
text(1975, seq(0, 100, by = 10), seq(0, 100, by = 10), pos = 2, cex = .7, xpd = TRUE)
text(1975 + shift, seq(0, 100, by = 10), seq(0, 100, by = 10), pos = 2, cex = .7, xpd = TRUE)
text(seq(1980, 2010, by = 10), 0, seq(1980, 2010, by = 10), pos = 1, cex = .7, xpd = TRUE)
text(seq(1980, 2010, by = 10) + shift, 0, seq(1980, 2010, by = 10), pos = 1, cex = .7, xpd = TRUE)
# axis labels
text(c(1971, 1971 + shift), 108, expression(e[x]), xpd = TRUE)
text(c(1990, 1990 + shift), -7, "Year", xpd = TRUE)
text(c(1990, 1990 + shift), 105, c("Male","Female"), cex = .8, xpd = TRUE)
# legend
.y  <- seq(0,104, length.out = 31)
.y2 <- seq(0,104, length.out = 16)
rect(2064, .y[1:(length(.y)-1)], 2070, .y[2:length(.y)], col = colramp(30), border = NA, xpd = TRUE)
segments(2070, .y2, 2071, .y2, xpd = TRUE)
text(2071, .y2, seq(0, .15, length.out = 16), pos = 4, xpd = TRUE, cex = .7)
dev.off()

#----------------------------------------------------------------
# exTFRlim plot


exTFRmUSlim <- colSums(exSFRUSmlim, na.rm = TRUE)
exTFRfUSlim <- colSums(exSFRUSflim, na.rm = TRUE)
exTFRmESlim <- colSums(exSFRESmlim, na.rm = TRUE)
exTFRfESlim <- colSums(exSFRESflim, na.rm = TRUE)

pdf("/home/triffe/git/DISS/latex/Figures/exTFRlim1555.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .3, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, exTFRmUSlim, type = 'l', ylim = c(2, 7), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968, 2, 2010, 7,col = gray(.95), border=NA),
                abline(h = seq(2, 7, by = .5), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(2, 7, by = .5),seq(2,7, by = .5), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10), 2, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, 1.8, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1965, 7.2, expression(e[x], "             -TFR"), cex = 1, xpd = TRUE)))
lines(yearsUS, exTFRfUSlim, lwd = 2.5, col = gray(.5))
lines(yearsES, exTFRmESlim, lwd = 2, col = gray(.2), lty = 5)
lines(yearsES, exTFRfESlim, lwd = 2.5, col = gray(.5), lty = 5)

legend(1995, 7, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
dev.off()


}


# ---------------------------------------
# 13-49 for females ; 15-64 for males
{
#----------------------------------------------------------------
#
ages <- 0:110
fkeep <- as.integer(ages >= 13 & ages <= 49)
mkeep <- as.integer(ages >= 15 & ages <= 64)
exSFRUSmlim <- Minf0(do.call(cbind, lapply(as.character(yearsUS), function(yr, .BxUS, .ExUS, .dxmUS){
                            Bx <-  rowSums(Mna0(.BxUS[[yr]]))
                            Ex <-  Mna0(with(.ExUS, Male[Year == as.integer(yr)])) * mkeep
                            Bex <- rowSums(ExpectedDx(Px = Bx, dx = .dxmUS[,yr]))
                            Eex <- rowSums(ExpectedDx(Px = Ex, dx = .dxmUS[,yr]))
                            Bex / Eex
                        }, .BxUS = BxUS, .ExUS = ExUS, .dxmUS = dxmUS)))
exSFRUSflim <- Minf0(do.call(cbind, lapply(as.character(yearsUS), function(yr, .BxUS, .ExUS, .dxfUS){
                            Bx  <- colSums(Mna0(.BxUS[[yr]]))
                            Ex  <- Mna0(with(.ExUS, Female[Year == as.integer(yr)])) * fkeep
                            Bex <- rowSums(ExpectedDx(Px = Bx, dx = .dxfUS[,yr]))
                            Eex <- rowSums(ExpectedDx(Px = Ex, dx = .dxfUS[,yr]))
                            Bex / Eex
                        }, .BxUS = BxUS, .ExUS = ExUS, .dxfUS = dxfUS)))
colnames(exSFRUSmlim) <- colnames(exSFRUSflim) <- yearsUS

exSFRESmlim <- Minf0(do.call(cbind, lapply(as.character(yearsES), function(yr, .BxES, .ExES, .dxmES){
                            Bx  <- rowSums(Mna0(.BxES[[yr]]))
                            Ex  <- Mna0(with(.ExES, Male[Year == as.integer(yr)])) * mkeep
                            Bex <- rowSums(ExpectedDx(Px = Bx, dx = .dxmES[,yr]))
                            Eex <- rowSums(ExpectedDx(Px = Ex, dx = .dxmES[,yr]))
                            Bex / Eex
                        }, .BxES = BxES, .ExES = ExES, .dxmES = dxmES)))
exSFRESflim <- Minf0(do.call(cbind, lapply(as.character(yearsES), function(yr, .BxES, .ExES, .dxfES){
                            Bx  <- colSums(Mna0(.BxES[[yr]]))
                            Ex  <- Mna0(with(.ExES, Female[Year == as.integer(yr)])) * fkeep
                            Bex <- rowSums(ExpectedDx(Px = Bx, dx = .dxfES[,yr]))
                            Eex <- rowSums(ExpectedDx(Px = Ex, dx = .dxfES[,yr]))
                            Bex / Eex
                        }, .BxES = BxES, .ExES = ExES, .dxfES = dxfES)))
colnames(exSFRESmlim) <- colnames(exSFRESflim) <- yearsES

colramp <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"), space = "Lab")
brks <- seq(0,.15, length.out = 51)
levs <- seq(0,.15,by = .01)
shift <- 50

# for USA:
pdf("/home/triffe/git/DISS/latex/Figures/eSFRsurfacesUSlim1364mixed.pdf", height = 5, width = 5)
par(mai = c(.3,.3,.3,1))
# males
image(x = yearsUS + .5, y = 0:110 + .5, t(exSFRUSmlim), 
        xlim = c(1969, 2010 + shift), ylim = c(0, 104), zlim = c(0, .15), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1,
        xlab = "", ylab = "")
abline(v=seq(1970,2010,by=10), col = "#FFFFFF50")
abline(h=seq(0,100,by=10), col = "#FFFFFF50")
# contours
contour(x = yearsUS + .5, y = 0:110 + .5, t(exSFRUSmlim),
        levels = levs, labels = levs, add = TRUE)

# females
image(x = yearsUS + .5 + shift, y = 0:110 + .5, t(exSFRUSflim), 
        xlim = c(1969, 2010) + shift, ylim = c(0, 104), zlim = c(0, .15), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1,
        xlab = "", ylab = "", add = TRUE)
abline(v=seq(1970,2010,by=10) + shift, col = "#FFFFFF50")
abline(h=seq(0,100,by=10), col = "#FFFFFF50")
# contours
contour(x = yearsUS + .5 + shift, y = 0:110 + .5, t(exSFRUSflim),
        levels = levs, labels = levs, add = TRUE)

# tick labels: 
text(1969, seq(0, 100, by = 10), seq(0, 100, by = 10), pos = 2, cex = .7, xpd = TRUE)
text(1969 + shift, seq(0, 100, by = 10), seq(0, 100, by = 10), pos = 2, cex = .7, xpd = TRUE)
text(seq(1970, 2010, by = 10), 0, seq(1970, 2010, by = 10), pos = 1, cex = .7, xpd = TRUE)
text(seq(1970, 2010, by = 10) + shift, 0, seq(1970, 2010, by = 10), pos = 1, cex = .7, xpd = TRUE)
# axis labels
text(c(1965, 1965 + shift), 108, expression(e[x]), xpd = TRUE)
text(c(1990, 1990 + shift), -7, "Year", xpd = TRUE)
text(c(1990, 1990 + shift), 105, c("Male","Female"), cex = .8, xpd = TRUE)
# legend
.y  <- seq(0,104, length.out = 31)
.y2 <- seq(0,104, length.out = 16)
rect(2064, .y[1:(length(.y)-1)], 2070, .y[2:length(.y)], col = colramp(30), border = NA, xpd = TRUE)
segments(2070, .y2, 2071, .y2, xpd = TRUE)
text(2071, .y2, seq(0, .15, length.out = 16), pos = 4, xpd = TRUE, cex = .7)
dev.off()

# for Spain:
pdf("/home/triffe/git/DISS/latex/Figures/eSFRsurfacesESlim1364mixed.pdf", height = 5, width = 5)
par(mai = c(.3,.3,.3,1))
# males
image(x = yearsES + .5, y = 0:110 + .5, t(exSFRESmlim), 
        xlim = c(1969, 2010 + shift), ylim = c(0, 104), zlim = c(0, .15), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1,
        xlab = "", ylab = "")
abline(v=seq(1970,2010,by=10), col = "#FFFFFF50")
abline(h=seq(0,100,by=10), col = "#FFFFFF50")
# contours
contour(x = yearsES + .5, y = 0:110 + .5, t(exSFRESmlim),
        levels = levs, labels = levs, add = TRUE)

# females
image(x = yearsES + .5 + shift, y = 0:110 + .5, t(exSFRESflim), 
        xlim = c(1969, 2010) + shift, ylim = c(0, 104), zlim = c(0, .12), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1,
        xlab = "", ylab = "", add = TRUE)
abline(v=seq(1980,2010,by=10) + shift, col = "#FFFFFF50")
abline(h=seq(0,100,by=10), col = "#FFFFFF50")
# contours
contour(x = yearsES + .5 + shift, y = 0:110 + .5, t(exSFRESflim),
        levels = levs, labels = levs, add = TRUE)

# tick labels: 
text(1975, seq(0, 100, by = 10), seq(0, 100, by = 10), pos = 2, cex = .7, xpd = TRUE)
text(1975 + shift, seq(0, 100, by = 10), seq(0, 100, by = 10), pos = 2, cex = .7, xpd = TRUE)
text(seq(1980, 2010, by = 10), 0, seq(1980, 2010, by = 10), pos = 1, cex = .7, xpd = TRUE)
text(seq(1980, 2010, by = 10) + shift, 0, seq(1980, 2010, by = 10), pos = 1, cex = .7, xpd = TRUE)
# axis labels
text(c(1971, 1971 + shift), 108, expression(e[x]), xpd = TRUE)
text(c(1990, 1990 + shift), -7, "Year", xpd = TRUE)
text(c(1990, 1990 + shift), 105, c("Male","Female"), cex = .8, xpd = TRUE)
# legend
.y  <- seq(0,104, length.out = 31)
.y2 <- seq(0,104, length.out = 16)
rect(2064, .y[1:(length(.y)-1)], 2070, .y[2:length(.y)], col = colramp(30), border = NA, xpd = TRUE)
segments(2070, .y2, 2071, .y2, xpd = TRUE)
text(2071, .y2, seq(0, .15, length.out = 16), pos = 4, xpd = TRUE, cex = .7)
dev.off()

#----------------------------------------------------------------
# exTFRlim plot


exTFRmUSlim <- colSums(exSFRUSmlim, na.rm = TRUE)
exTFRfUSlim <- colSums(exSFRUSflim, na.rm = TRUE)
exTFRmESlim <- colSums(exSFRESmlim, na.rm = TRUE)
exTFRfESlim <- colSums(exSFRESflim, na.rm = TRUE)

pdf("/home/triffe/git/DISS/latex/Figures/exTFRlim1364mixed.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .3, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, exTFRmUSlim, type = 'l', ylim = c(1.5, 5.5), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968, 1.5, 2010, 5.5,col = gray(.95), border=NA),
                abline(h = seq(2, 7, by = .5), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(1.5, 5.5, by = .5),seq(1.5,5.5, by = .5), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10), 1.5, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, 1.3, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1965, 5.7, expression(e[x], "             -TFR"), cex = 1, xpd = TRUE)))
lines(yearsUS, exTFRfUSlim, lwd = 2.5, col = gray(.5))
lines(yearsES, exTFRmESlim, lwd = 2, col = gray(.2), lty = 5)
lines(yearsES, exTFRfESlim, lwd = 2.5, col = gray(.5), lty = 5)

legend(1995, 5.5, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
dev.off()
}

#----------------------------------------------------------------

#----------------------------------------------------------------
# 1st-99th quantile of ASFR ages (quantiles of counts = quicker), separate by year = NO!
{
#----------------------------------------------------------------
ages <- 0:110
findKeeperAges <- function(Bx, Ex, q = c(.01,.99), sex = "m"){
    if (sex == "m"){
        marfun <- rowSums
        coln   <- "Male"
    } else {
        marfun <- colSums
        coln   <- "Female"
    }
    years      <- names(Bx)
    fxt        <- rowSums(do.call(cbind,lapply(as.character(years), function(yr, .Bx, .Ex){
                       Minf0(Mna0(marfun(Mna0(.Bx[[yr]])) / Mna0(.Ex[[coln]][.Ex$Year == as.integer(yr)])))
                            },.Bx = Bx, .Ex = Ex)))
    fxt[85:111] <- 0
    fxt         <- cumsum(fxt) / sum(fxt)
    keep.vec    <- as.integer(fxt > min(q) & fxt < max(q))
    keep.vec
}
mkeepUS <- findKeeperAges(Bx = BxUS, Ex = ExUS, sex = "m")
fkeepUS <- findKeeperAges(Bx = BxUS, Ex = ExUS, sex = "f")
mkeepES <- findKeeperAges(Bx = BxES, Ex = ExES, sex = "m")
fkeepES <- findKeeperAges(Bx = BxES, Ex = ExES, sex = "f")

exSFRUSmlim <- Minf0(do.call(cbind, lapply(as.character(yearsUS), function(yr, .BxUS, .ExUS, .dxmUS, .mkeepUS){
                            Bx <-  rowSums(Mna0(.BxUS[[yr]]))
                            Ex <-  Mna0(with(.ExUS, Male[Year == as.integer(yr)])) * .mkeepUS
                            Bex <- rowSums(ExpectedDx(Px = Bx, dx = .dxmUS[,yr]))
                            Eex <- rowSums(ExpectedDx(Px = Ex, dx = .dxmUS[,yr]))
                            Bex / Eex
                        }, .BxUS = BxUS, .ExUS = ExUS, .dxmUS = dxmUS, .mkeepUS = mkeepUS)))
exSFRUSflim <- Minf0(do.call(cbind, lapply(as.character(yearsUS), function(yr, .BxUS, .ExUS, .dxfUS, .fkeepUS){
                            Bx  <- colSums(Mna0(.BxUS[[yr]]))
                            Ex  <- Mna0(with(.ExUS, Female[Year == as.integer(yr)])) * .fkeepUS
                            Bex <- rowSums(ExpectedDx(Px = Bx, dx = .dxfUS[,yr]))
                            Eex <- rowSums(ExpectedDx(Px = Ex, dx = .dxfUS[,yr]))
                            Bex / Eex
                        }, .BxUS = BxUS, .ExUS = ExUS, .dxfUS = dxfUS, .fkeepUS = fkeepUS)))
colnames(exSFRUSmlim) <- colnames(exSFRUSflim) <- yearsUS

exSFRESmlim <- Minf0(do.call(cbind, lapply(as.character(yearsES), function(yr, .BxES, .ExES, .dxmES, .mkeepES){
                            Bx  <- rowSums(Mna0(.BxES[[yr]]))
                            Ex  <- Mna0(with(.ExES, Male[Year == as.integer(yr)])) * .mkeepES
                            Bex <- rowSums(ExpectedDx(Px = Bx, dx = .dxmES[,yr]))
                            Eex <- rowSums(ExpectedDx(Px = Ex, dx = .dxmES[,yr]))
                            Bex / Eex
                        }, .BxES = BxES, .ExES = ExES, .dxmES = dxmES, .mkeepES = mkeepES)))
exSFRESflim <- Minf0(do.call(cbind, lapply(as.character(yearsES), function(yr, .BxES, .ExES, .dxfES, .fkeepES){
                            Bx  <- colSums(Mna0(.BxES[[yr]]))
                            Ex  <- Mna0(with(.ExES, Female[Year == as.integer(yr)])) * .fkeepES
                            Bex <- rowSums(ExpectedDx(Px = Bx, dx = .dxfES[,yr]))
                            Eex <- rowSums(ExpectedDx(Px = Ex, dx = .dxfES[,yr]))
                            Bex / Eex
                        }, .BxES = BxES, .ExES = ExES, .dxfES = dxfES, .fkeepES = fkeepES)))
colnames(exSFRESmlim) <- colnames(exSFRESflim) <- yearsES


colramp <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"), space = "Lab")
brks <- seq(0,.15, length.out = 51)
levs <- seq(0,.15,by = .01)
shift <- 50

# for USA:
pdf("/home/triffe/git/DISS/latex/Figures/eSFRsurfacesUSlimBquant.pdf", height = 5, width = 5)
par(mai = c(.3,.3,.3,1))
# males
image(x = yearsUS + .5, y = 0:110 + .5, t(exSFRUSmlim), 
        xlim = c(1969, 2010 + shift), ylim = c(0, 104), zlim = c(0, .15), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1,
        xlab = "", ylab = "")
abline(v=seq(1970,2010,by=10), col = "#FFFFFF50")
abline(h=seq(0,100,by=10), col = "#FFFFFF50")
# contours
contour(x = yearsUS + .5, y = 0:110 + .5, t(exSFRUSmlim),
        levels = levs, labels = levs, add = TRUE)

# females
image(x = yearsUS + .5 + shift, y = 0:110 + .5, t(exSFRUSflim), 
        xlim = c(1969, 2010) + shift, ylim = c(0, 104), zlim = c(0, .15), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1,
        xlab = "", ylab = "", add = TRUE)
abline(v=seq(1970,2010,by=10) + shift, col = "#FFFFFF50")
abline(h=seq(0,100,by=10), col = "#FFFFFF50")
# contours
contour(x = yearsUS + .5 + shift, y = 0:110 + .5, t(exSFRUSflim),
        levels = levs, labels = levs, add = TRUE)

# tick labels: 
text(1969, seq(0, 100, by = 10), seq(0, 100, by = 10), pos = 2, cex = .7, xpd = TRUE)
text(1969 + shift, seq(0, 100, by = 10), seq(0, 100, by = 10), pos = 2, cex = .7, xpd = TRUE)
text(seq(1970, 2010, by = 10), 0, seq(1970, 2010, by = 10), pos = 1, cex = .7, xpd = TRUE)
text(seq(1970, 2010, by = 10) + shift, 0, seq(1970, 2010, by = 10), pos = 1, cex = .7, xpd = TRUE)
# axis labels
text(c(1965, 1965 + shift), 108, expression(e[x]), xpd = TRUE)
text(c(1990, 1990 + shift), -7, "Year", xpd = TRUE)
text(c(1990, 1990 + shift), 105, c("Male","Female"), cex = .8, xpd = TRUE)
# legend
.y  <- seq(0,104, length.out = 31)
.y2 <- seq(0,104, length.out = 16)
rect(2064, .y[1:(length(.y)-1)], 2070, .y[2:length(.y)], col = colramp(30), border = NA, xpd = TRUE)
segments(2070, .y2, 2071, .y2, xpd = TRUE)
text(2071, .y2, seq(0, .15, length.out = 16), pos = 4, xpd = TRUE, cex = .7)
dev.off()

# for Spain:
pdf("/home/triffe/git/DISS/latex/Figures/eSFRsurfacesESlimBquant.pdf", height = 5, width = 5)
par(mai = c(.3,.3,.3,1))
# males
image(x = yearsES + .5, y = 0:110 + .5, t(exSFRESmlim), 
        xlim = c(1969, 2010 + shift), ylim = c(0, 104), zlim = c(0, .15), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1,
        xlab = "", ylab = "")
abline(v=seq(1970,2010,by=10), col = "#FFFFFF50")
abline(h=seq(0,100,by=10), col = "#FFFFFF50")
# contours
contour(x = yearsES + .5, y = 0:110 + .5, t(exSFRESmlim),
        levels = levs, labels = levs, add = TRUE)

# females
image(x = yearsES + .5 + shift, y = 0:110 + .5, t(exSFRESflim), 
        xlim = c(1969, 2010) + shift, ylim = c(0, 104), zlim = c(0, .12), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1,
        xlab = "", ylab = "", add = TRUE)
abline(v=seq(1980,2010,by=10) + shift, col = "#FFFFFF50")
abline(h=seq(0,100,by=10), col = "#FFFFFF50")
# contours
contour(x = yearsES + .5 + shift, y = 0:110 + .5, t(exSFRESflim),
        levels = levs, labels = levs, add = TRUE)

# tick labels: 
text(1975, seq(0, 100, by = 10), seq(0, 100, by = 10), pos = 2, cex = .7, xpd = TRUE)
text(1975 + shift, seq(0, 100, by = 10), seq(0, 100, by = 10), pos = 2, cex = .7, xpd = TRUE)
text(seq(1980, 2010, by = 10), 0, seq(1980, 2010, by = 10), pos = 1, cex = .7, xpd = TRUE)
text(seq(1980, 2010, by = 10) + shift, 0, seq(1980, 2010, by = 10), pos = 1, cex = .7, xpd = TRUE)
# axis labels
text(c(1971, 1971 + shift), 108, expression(e[x]), xpd = TRUE)
text(c(1990, 1990 + shift), -7, "Year", xpd = TRUE)
text(c(1990, 1990 + shift), 105, c("Male","Female"), cex = .8, xpd = TRUE)
# legend
.y  <- seq(0,104, length.out = 31)
.y2 <- seq(0,104, length.out = 16)
rect(2064, .y[1:(length(.y)-1)], 2070, .y[2:length(.y)], col = colramp(30), border = NA, xpd = TRUE)
segments(2070, .y2, 2071, .y2, xpd = TRUE)
text(2071, .y2, seq(0, .15, length.out = 16), pos = 4, xpd = TRUE, cex = .7)
dev.off()

#----------------------------------------------------------------
# exTFRlim plot


exTFRmUSlim <- colSums(exSFRUSmlim, na.rm = TRUE)
exTFRfUSlim <- colSums(exSFRUSflim, na.rm = TRUE)
exTFRmESlim <- colSums(exSFRESmlim, na.rm = TRUE)
exTFRfESlim <- colSums(exSFRESflim, na.rm = TRUE)

pdf("/home/triffe/git/DISS/latex/Figures/exTFRlimBquant.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .3, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, exTFRmUSlim, type = 'l', ylim = c(2.5, 9.5), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968, 2, 2010, 9.5,col = gray(.95), border=NA),
                abline(h = seq(2.5, 9.5, by = .5), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(3, 9, by = 1),seq(3, 9, by = 1), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10), 2.5, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, 2.2, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1965, 9.85, expression(e[x], "             -TFR"), cex = 1, xpd = TRUE)))
lines(yearsUS, exTFRfUSlim, lwd = 2.5, col = gray(.5))
lines(yearsES, exTFRmESlim, lwd = 2, col = gray(.2), lty = 5)
lines(yearsES, exTFRfESlim, lwd = 2.5, col = gray(.5), lty = 5)

legend(1994, 9.5, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
dev.off()

}
#----------------------------------------------------------------
# finally the 1st -99th percentile on a year-by year basis, exTFR only
{
exTFRUSmlim <- colSums(Minf0(do.call(cbind, lapply(as.character(yearsUS), function(yr, .BxUS, .ExUS, .dxmUS){
                                    Bx    <- rowSums(Mna0(.BxUS[[yr]]))
                                    Ex    <- Mna0(with(.ExUS, Male[Year == as.integer(yr)])) 
                                    
                                    Fx    <- Bx / Ex
                                    Fx[85:111]  <- 0
                                    Fx          <- cumsum(Fx) / sum(Fx)
                                    keep.vec    <- as.integer(Fx > .01 & Fx < .99)
                                    
                                    Bex   <- rowSums(ExpectedDx(Px = Bx, dx = .dxmUS[,yr]))
                                    Eex   <- rowSums(ExpectedDx(Px = Ex * keep.vec, dx = .dxmUS[,yr]))
                                    Bex / Eex
                                }, .BxUS = BxUS, .ExUS = ExUS, .dxmUS = dxmUS))), na.rm = TRUE)
exTFRUSflim <- colSums(Minf0(do.call(cbind, lapply(as.character(yearsUS), function(yr, .BxUS, .ExUS, .dxfUS){
                                    Bx    <- colSums(Mna0(.BxUS[[yr]]))
                                    Ex    <- Mna0(with(.ExUS, Female[Year == as.integer(yr)])) 
                                    
                                    Fx    <- Bx / Ex
                                    Fx[85:111]  <- 0
                                    Fx          <- cumsum(Fx) / sum(Fx)
                                    keep.vec    <- as.integer(Fx > .01 & Fx < .99)
                                    
                                    Bex   <- rowSums(ExpectedDx(Px = Bx, dx = .dxfUS[,yr]))
                                    Eex   <- rowSums(ExpectedDx(Px = Ex * keep.vec, dx = .dxfUS[,yr]))
                                    Bex / Eex
                                }, .BxUS = BxUS, .ExUS = ExUS, .dxfUS = dxfUS))), na.rm = TRUE)
exTFRESmlim <- colSums(Minf0(do.call(cbind, lapply(as.character(yearsES), function(yr, .BxES, .ExES, .dxmES){
                                    Bx    <- rowSums(Mna0(.BxES[[yr]]))
                                    Ex    <- Mna0(with(.ExES, Male[Year == as.integer(yr)])) * mkeepES
                                    
                                    Fx    <- Bx / Ex
                                    Fx[85:111]  <- 0
                                    Fx          <- cumsum(Fx) / sum(Fx)
                                    keep.vec    <- as.integer(Fx > .01 & Fx < .99)
                                    
                                    Bex   <- rowSums(ExpectedDx(Px = Bx, dx = .dxmES[,yr]))
                                    Eex   <- rowSums(ExpectedDx(Px = Ex * keep.vec, dx = .dxmES[,yr]))
                                    Bex / Eex
                                }, .BxES = BxES, .ExES = ExES, .dxmES = dxmES))), na.rm = TRUE)
exTFRESflim <- colSums(Minf0(do.call(cbind, lapply(as.character(yearsES), function(yr, .BxES, .ExES, .dxfES){
                                    Bx    <- colSums(Mna0(.BxES[[yr]]))
                                    Ex    <- Mna0(with(.ExES, Female[Year == as.integer(yr)])) * fkeepES
                                    
                                    Fx    <- Bx / Ex
                                    Fx[85:111]  <- 0
                                    Fx          <- cumsum(Fx) / sum(Fx)
                                    keep.vec    <- as.integer(Fx > .01 & Fx < .99)
                                    
                                    Bex   <- rowSums(ExpectedDx(Px = Bx, dx = .dxfES[,yr]))
                                    Eex   <- rowSums(ExpectedDx(Px = Ex * keep.vec, dx = .dxfES[,yr]))
                                    Bex / Eex
                                }, .BxES = BxES, .ExES = ExES, .dxfES = dxfES))), na.rm = TRUE)

pdf("/home/triffe/git/DISS/latex/Figures/exTFRlimBquantyr.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .3, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, exTFRmUSlim, type = 'l', ylim = c(2.5, 9.5), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968, 2, 2010, 9.5,col = gray(.95), border=NA),
                abline(h = seq(2.5, 9.5, by = .5), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(3, 9, by = 1),seq(3, 9, by = 1), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10), 2.5, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, 2.2, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1965, 9.85, expression(e[x], "             -TFR"), cex = 1, xpd = TRUE)))
lines(yearsUS, exTFRfUSlim, lwd = 2.5, col = gray(.5))
lines(yearsES, exTFRmESlim, lwd = 2, col = gray(.2), lty = 5)
lines(yearsES, exTFRfESlim, lwd = 2.5, col = gray(.5), lty = 5)

legend(1994, 9.5, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
dev.off()


}