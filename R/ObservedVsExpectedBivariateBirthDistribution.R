
# only one plot is needed to demonstrate this, since the basic observation will hold
# for any human population.

source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

BxUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy10_65.Rdata")))
# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))

Bxy         <- BxUS[["1970"]]
rowTotals   <- rowSums(Bxy)
colTotals   <- colSums(Bxy)
nOfCases    <- sum(rowTotals)
expected    <- outer(rowSums(Bxy), colSums(Bxy), "*") / nOfCases

colramp <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"), space = "Lab")

Zmax <- max(pretty(max(c(Bxy, expected))))
ages <- 10:65

#outer(10^c(0:-10))

logBxy <- MinfNA(log(Bxy))
logExy <- MinfNA(log(expected))

ticks <- c(.00001, .00002, .00005, .0001, .0002, .0005, 0.001, .002, .005, .01,.02, .05, .1)
brks  <- seq(-11.5, -2.5, by = .05)

image.plot(x = ages + .5, y = ages + .5, t(log(Bxy)), 
        xlim = c(10, 66), ylim = c(10, 66), zlim = c(0, 11), 
        col = colramp(50), axis.args = list(at = log(ticks), labels = ticks))




















