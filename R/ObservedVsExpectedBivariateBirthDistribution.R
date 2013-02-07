
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

logBxy      <- MinfNA(log(Bxy))
logExy      <- MinfNA(log(expected))

TickMaj     <- rev(c(10^(1:-6)))
TickMajLab  <- c(".000001",".00001",".0001",".001",".01",".1","1","10")
TickMin     <- c(t(outer(TickMaj, 2:9)))



g.xy <- seq(10, 65, by = 5)
gb.xy <- seq(10, 65, by = 10)

graphics.off()
dev.new(height = 4, width = 6.5)
par(mfrow=c(1,2), mar = c(1,1,1,1))
image(x = ages + .5, y = ages + .5, t(logBxy), 
        xlim = c(10, 66), ylim = c(10, 66), zlim = c(0, 11), 
        col = colramp(50), axes = FALSE, asp = 1, 
        panel.first = list(rect(10, 10, 66, 66, col = "#EEEEEE", xpd = TRUE, border = NA), 
                           segments(g.xy, 10, g.xy, 66, col = "white", lwd = .5),
                           segments(10, g.xy, 66, g.xy, col = "white", lwd = .5),
                           text(10, g.xy, g.xy, pos = 2, cex = .5, xpd = TRUE),
                           text(g.xy, 10, g.xy, pos = 1, cex = .5, xpd = TRUE)))
leg.y <- seq(0,11,length = 20) * (56 / 11) + 10
rect(10,leg.y[1:(length(leg.y) - 1)] )





















