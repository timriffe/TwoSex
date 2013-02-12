
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

sum(expected[lower.tri(expected, TRUE)]) /
sum(expected[upper.tri(expected, TRUE)])

sum(Bxy[lower.tri(Bxy, TRUE)]) /
sum(Bxy[upper.tri(Bxy, TRUE)])

upper.tri(matrix(0,ncol=3,nrow = 3))

colramp <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"), space = "Lab")

Zmax <- max(pretty(max(c(Bxy, expected))))
ages <- 10:65

#outer(10^c(0:-10))

logBxy      <- MinfNA(log(Bxy))
logExy      <- MinfNA(log(expected))
logBxy[nrow(logBxy), ] <- NA # remove open age males
logExy[nrow(logExy), ] <- NA
TickMaj     <- rev(c(10^(0:4)))
TickMajLab  <- c("10000","1000","100","10","1")
TickMin     <- t(outer(TickMaj[2:length(TickMaj)], 2:9))



g.xy <- seq(10, 65, by = 5)
gb.xy <- seq(10, 65, by = 10)

graphics.off()
dev.new(height = 4, width = 6.5)

brks <- seq(min(logBxy, na.rm = TRUE),max(logBxy, na.rm = TRUE), length.out = 51)

pdf("/home/triffe/git/DISS/DiagnosticPlots/ObsvsExpectedBxy/ObservedvsExpectedBxy.pdf", 
        height = 4, width = 6.5)
par(mfrow=c(1,2), mar = c(2,1,2,3))
image(x = ages + .5, y = ages + .5, logBxy, 
        xlim = c(10, 66), ylim = c(10, 66), zlim = c(0, 11), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1, 
        panel.first = list(rect(10, 10, 66, 66, col = "#EEEEEE", xpd = TRUE, border = NA), 
                           segments(g.xy, 10, g.xy, 66, col = "white", lwd = .5),
                           segments(10, g.xy, 66, g.xy, col = "white", lwd = .5),
                           text(10, g.xy, g.xy, pos = 2, cex = .5, xpd = TRUE),
                           text(g.xy, 10, g.xy, pos = 1, cex = .5, xpd = TRUE)),
        xlab = "", ylab = "")
# contours
levs <- c(1,10,100,1000,10000)         
contour(x = ages[1:55] + .5, y = ages[1:55] + .5, logBxy[1:55,1:55], 
        levels = log(levs), labels = levs, add = TRUE)

# legend
leg.y <- seq(0,11,length = 50) * (56 / 11) + 10
rect(68,leg.y[1:(length(leg.y) - 1)],72, leg.y[2:length(leg.y)],
        col = colramp(50), border = NA, xpd = TRUE)
minmax <- range(brks)

TickMajsc <- ((log(TickMaj)) / diff(range(minmax))) * diff(range(leg.y)) + 10
TickMinsc <- ((log(TickMin)) / diff(range(minmax))) * diff(range(leg.y)) + 10

segments(72,TickMajsc,72.5,TickMajsc,col = "#444444", xpd = TRUE)
text(72,TickMajsc, TickMajLab, pos = 4, xpd = TRUE, cex = .5)
segments(72,TickMinsc,72.3,TickMinsc,col = "#444444", xpd = TRUE)
# line of homogamy:
segments(10,10,66,66,col = "#50505050")

# axis labels
text(40,4,"Age of father", xpd = TRUE, cex = .7)
text(4,68,"Age of Mother", xpd = TRUE, pos = 4, cex = .7)
# expected bivariate distribution:
par(mar=c(2,2,2,2))
image(x = ages + .5, y = ages + .5, logExy, 
        xlim = c(10, 66), ylim = c(10, 66), zlim = c(0, 11), 
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1, 
        panel.first = list(rect(10, 10, 66, 66, col = "#EEEEEE", xpd = TRUE, border = NA), 
                segments(g.xy, 10, g.xy, 66, col = "white", lwd = .5),
                segments(10, g.xy, 66, g.xy, col = "white", lwd = .5),
                text(10, g.xy, g.xy, pos = 2, cex = .5, xpd = TRUE),
                text(g.xy, 10, g.xy, pos = 1, cex = .5, xpd = TRUE)),
        xlab = "", ylab = "")
# axis labels
text(40,4,"Age of father", xpd = TRUE, cex = .7)
text(4,68,"Age of Mother", xpd = TRUE, pos = 4, cex = .7)
# contours
levs <- c(1,10,100,1000,10000)         
contour(x = ages[1:55] + .5, y = ages[1:55] + .5, logExy[1:55,1:55], 
        levels = log(levs), labels = levs, add = TRUE)
# line of homogamy:
segments(10,10,66,66,col = "#50505050")
dev.off()


# total variation distance:
(ToalVar <- sum(abs(Bxy / sum(Bxy) - expected / sum(expected))) / 2 )

# overlap 
overlap     <- pmin(Bxy / sum(Bxy), expected / sum(expected))
(coefOverlap <- sum(overlap))

# two indices sum to 1
(coefOverlap + ToalVar)


# TFR that would result from no honogamy bias:
Exm <- with(ExUS, Male[Age >= 10 & Age <= 65 & Year == 1970])
Exf <- with(ExUS, Female[Age >= 10 & Age <= 65 & Year == 1970])

(TFRfexpected <- sum(expected %col% (1/Exf)))
(TFRfobserved <- sum(Bxy %col% (1/Exf))) # identical


TotalVar <- unlist(lapply(BxUS, function(x){
            x[is.na(x)] <- 0
            rowTotals   <- rowSums(x)
            colTotals   <- colSums(x)
            nOfCases    <- sum(rowTotals)
            # expected distr
            E           <- outer(rowSums(x), colSums(x), "*") / nOfCases
            # coef of diff
            sum(abs(x / sum(x) - E / sum(E))) / 2
 
        })  )
# i.e. convergence as fertility decreased.
plot(as.integer(names(BxUS)), TotalVar, type = 'l')
head(ExUS)
ExfList <- tapply(ExUS$Female, ExUS$Year, function(x){
            x[11:66]
        })
TFRf <- unlist(mapply(function(x,y){
                    x[is.na(x)] <- 0
                    
                    sum(colSums(x) / y)
                    
                },BxUS,ExfList, )  )









