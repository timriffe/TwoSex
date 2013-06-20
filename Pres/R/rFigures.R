# figure of r for both sexes and standard Lotka r
setwd("/home/triffe/git/DISS/")

rLotkaES    <- local(get(load("Data/results/agerSRB/rLotkaES.Rdata")))
rLotkaUS    <- local(get(load("Data/results/agerSRB/rLotkaUS.Rdata")))
rfES        <- local(get(load("/home/triffe/git/DISS/Data/results/exSingleSex/rfES.Rdata")))
rmES        <- local(get(load("/home/triffe/git/DISS/Data/results/exSingleSex/rmES.Rdata")))
rfUS        <- local(get(load("/home/triffe/git/DISS/Data/results/exSingleSex/rfUS.Rdata")))
rmUS        <- local(get(load("/home/triffe/git/DISS/Data/results/exSingleSex/rmUS.Rdata")))
yearsUS     <- 1969:2009
yearsES     <- 1975:2009
Cols        <- RColorBrewer::brewer.pal(9,"Set1")

# toward intro, figure showing differences in sing-sex lotka
pdf("Pres/FiguresStatic/rSingleSexLotka.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .3, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rmUS[, 1], type = 'n', ylim = c(-.02,.011),xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,-.02,2010,.011,col = gray(.95), border=NA),
                abline(h = seq(-.02,.011,by = .0025), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.02,.011,by = .005),seq(-.02,.011,by = .005), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.02, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.022, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,.0125, "r", cex = 1, xpd = TRUE)))

lines(yearsES, rLotkaES[,1], col = Cols[2], lwd = 2.5, lty = 5)
lines(yearsES, rLotkaES[,2], col = Cols[8], lwd = 2.5, lty = 5)
lines(yearsUS, rLotkaUS[,1], col = Cols[2], lwd = 2)
lines(yearsUS, rLotkaUS[,2], col = Cols[8], lwd = 2)
text(c(1976, 1975, 1992, 1985),c(0.0021, -0.008387479, -0.008657443, -0.015028590),
        c("US males","US females","ES males","ES females"))
dev.off()


# for age US
pdf("Pres/FiguresStatic/rSingleSex1.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rmUS[, 1], type = 'n', ylim = c(-.02,.011),xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,-.02,2010,.011,col = gray(.95), border=NA),
                abline(h = seq(-.02,.011,by = .0025), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.02,.011,by = .005),seq(-.02,.011,by = .005), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.02, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.022, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,.0125, "r", cex = 1, xpd = TRUE)))

lines(yearsES, rmES[, 1],col = grey(.9), lwd = 2.5)
lines(yearsES, rfES[, 1],col = grey(.9), lwd = 3)
lines(yearsES, rLotkaES[,1], col = grey(.9), lwd = 1)
lines(yearsES, rLotkaES[,2], col = grey(.9), lwd = 1)
lines(yearsUS, rmUS[, 1],col = grey(.9), lwd = 2.5)
lines(yearsUS, rfUS[, 1],col = grey(.9), lwd = 3)
lines(yearsUS, rLotkaUS[,1], col = Cols[2], lwd = 1)
lines(yearsUS, rLotkaUS[,2], col = Cols[8], lwd = 1)
text(c(1974, 1975),c(0.005, -0.008387479),
        c("US males\nchrono","US females\nchrono"))
dev.off()

# for age and ey US
pdf("Pres/FiguresStatic/rSingleSex2.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rmUS[, 1], type = 'n', ylim = c(-.02,.011),xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,-.02,2010,.011,col = gray(.95), border=NA),
                abline(h = seq(-.02,.011,by = .0025), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.02,.011,by = .005),seq(-.02,.011,by = .005), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.02, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.022, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,.0125, "r", cex = 1, xpd = TRUE)))

lines(yearsES, rmES[, 1],col = grey(.9), lwd = 2.5)
lines(yearsES, rfES[, 1],col = grey(.9), lwd = 3)
lines(yearsES, rLotkaES[,1], col = grey(.9), lwd = 1)
lines(yearsES, rLotkaES[,2], col = grey(.9), lwd = 1)
lines(yearsUS, rmUS[, 1],col = Cols[2], lwd = 2.5)
lines(yearsUS, rfUS[, 1],col = Cols[8], lwd = 3)
lines(yearsUS, rLotkaUS[,1], col = Cols[2], lwd = 1)
lines(yearsUS, rLotkaUS[,2], col = Cols[8], lwd = 1)
text(c(1974, 1975, 1986, 1986),c(0.005, -0.008387479, .004, -.0002),
        c("US males\nchrono","US females\nchrono","US males\nthano","US females\nthano"))
dev.off()

pdf("Pres/FiguresStatic/rSingleSex3.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rmUS[, 1], type = 'n', ylim = c(-.02,.011),xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,-.02,2010,.011,col = gray(.95), border=NA),
                abline(h = seq(-.02,.011,by = .0025), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.02,.011,by = .005),seq(-.02,.011,by = .005), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.02, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.022, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,.0125, "r", cex = 1, xpd = TRUE)))
lines(yearsES, rmES[, 1],col = grey(.9), lwd = 2.5)
lines(yearsES, rfES[, 1],col = grey(.9), lwd = 3)
lines(yearsUS, rmUS[, 1],col = grey(.9), lwd = 2.5)
lines(yearsUS, rfUS[, 1],col = grey(.9), lwd = 3)
lines(yearsUS, rLotkaUS[,1], col = grey(.9), lwd = 1)
lines(yearsUS, rLotkaUS[,2], col = grey(.9), lwd = 1)
lines(yearsES, rLotkaES[,1], col = Cols[2], lwd = 1)
lines(yearsES, rLotkaES[,2], col =  Cols[8], lwd = 1)
text(c(1995.5, 1985),c( -0.016, -0.015028590),
        c("ES males\nchrono","ES females\nchrono"))
dev.off()

pdf("Pres/FiguresStatic/rSingleSex4.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rmUS[, 1], type = 'n', ylim = c(-.02,.011),xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,-.02,2010,.011,col = gray(.95), border=NA),
                abline(h = seq(-.02,.011,by = .0025), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.02,.011,by = .005),seq(-.02,.011,by = .005), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.02, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.022, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,.0125, "r", cex = 1, xpd = TRUE)))

lines(yearsUS, rmUS[, 1],col = grey(.9), lwd = 2.5)
lines(yearsUS, rfUS[, 1],col = grey(.9), lwd = 3)
lines(yearsUS, rLotkaUS[,1], col = grey(.9), lwd = 1)
lines(yearsUS, rLotkaUS[,2], col = grey(.9), lwd = 1)
lines(yearsES, rLotkaES[,1], col = Cols[2], lwd = 1)
lines(yearsES, rLotkaES[,2], col =  Cols[8], lwd = 1)
lines(yearsES, rmES[, 1],col = Cols[2], lwd = 2.5)
lines(yearsES, rfES[, 1],col = Cols[8], lwd = 3)
text(c(1995.5, 1985,2000,2004.5),c( -0.016, -0.015028590,-0.007145646, -0.010169241),
        c("ES males\nchrono","ES females\nchrono","ES males\nthano","ES females\nthano"))
dev.off()

# ----------------------------------------------------------------------------- ## ----------------------------------------------------------------------------- #

