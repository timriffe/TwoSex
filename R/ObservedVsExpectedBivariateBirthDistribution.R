
# only one plot is needed to demonstrate this, since the basic observation will hold
# for any human population.

source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

BxES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata")))
BxUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy10_65.Rdata")))
# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))

Bxy         <- BxUS[["1970"]]
rowTotals   <- rowSums(Bxy)
colTotals   <- colSums(Bxy)
nOfCases    <- sum(rowTotals)
expected    <- outer(rowSums(Bxy), colSums(Bxy), "*") / nOfCases

ks.test(Bxy, expected)$statistic

# structural hypogamy due only to single-sex spans and distributions.

# -----------------------------------------------
# including diagonal
# lower = hyergamy, upper = hypogamy
sum(expected[lower.tri(expected, TRUE)]) /
sum(expected[upper.tri(expected, TRUE)])

# lower = hyergamy, upper = hypogamy
sum(Bxy[lower.tri(Bxy, TRUE)]) /
sum(Bxy[upper.tri(Bxy, TRUE)])
# -----------------------------------------------
# percent on the diagonal
sum(diag(Bxy)) / sum(Bxy)
sum(diag(expected)) / sum(expected)
# -----------------------------------------------
# excluding diagonal
sum(expected[lower.tri(expected)]) /
sum(expected[upper.tri(expected)])

# lower = hyergamy, upper = hypogamy
sum(Bxy[lower.tri(Bxy)]) /
sum(Bxy[upper.tri(Bxy)])
#------------------------------
# excess hypergamy:
(sum(Bxy[lower.tri(Bxy)]) / sum(Bxy[upper.tri(Bxy)])) / 
(sum(expected[lower.tri(expected)]) / sum(expected[upper.tri(expected)]))

# --------------------------------
# surfaces for observed and expected Bxy
{
colramp <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"), space = "Lab")

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
#dev.new(height = 4, width = 6.5)

brks <- seq(min(logBxy, na.rm = TRUE),max(logBxy, na.rm = TRUE), length.out = 51)

pdf("/home/triffe/git/DISS/latex/Figures/ObservedvsExpectedBxy.pdf", 
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
}

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

# ------------------------------------------------------------------
# calculate total difference (theta) between observed and expected:
# add confidence bands:
BxyES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata")))
BxyUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
yearsES <- 1975:2009
yearsUS <- 1969:2009
names(BxyES) <- yearsES
TotalVarUS <- unlist(lapply(as.character(yearsUS), function(yr, .BxyUS){
                    
            x           <-Mna0(.BxyUS[[yr]])
            rowTotals   <- rowSums(x)
            colTotals   <- colSums(x)
            nOfCases    <- sum(rowTotals)
            # expected distr
            E           <- outer(rowSums(x), colSums(x), "*") / nOfCases
            # coef of diff
            sum(abs(x / sum(x) - E / sum(E))) / 2
 
        }, .BxyUS=BxyUS)  )

# i.e. convergence as fertility decreased.

TotalVarES <- unlist(lapply(as.character(yearsES), function(yr, .BxyES){
                    
                    x           <-Mna0(.BxyES[[yr]])
                    rowTotals   <- rowSums(x)
                    colTotals   <- colSums(x)
                    nOfCases    <- sum(rowTotals)
                    # expected distr
                    E           <- outer(rowSums(x), colSums(x), "*") / nOfCases
                    # coef of diff
                    sum(abs(x / sum(x) - E / sum(E))) / 2
                    
                }, .BxyES=BxyES)  )
# these take about 3 minutes each to run...
# curse my not having an analytic solution!
TotalVarUS95 <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .BxyUS){
                    yr <- "1975"
                 x           <- Mna0(.BxyUS[[yr]])
               
                 quantile(replicate(1000, {
                             X           <- matrix(rpois(n = length(x), lambda = x), ncol = ncol(x))
                             rowTotals   <- rowSums(X)
                             colTotals   <- colSums(X)
                             # expected distr
                             E           <- outer(rowSums(X), colSums(X), "*") / sum(X)
                             1 - sum(pmin(X / sum(X), E / sum(E)))
                         }) , probs = c(.025, .975)) 
                }, .BxyUS=BxyUS)  )
TotalVarES95 <- do.call(rbind,lapply(as.character(yearsES), function(yr, .BxyES){
                  x           <- Mna0(.BxyES[[yr]])
                  quantile(replicate(1000, {
                             X <- matrix(rpois(n = length(x), lambda = x), ncol=ncol(x))
                             rowTotals   <- rowSums(X)
                             colTotals   <- colSums(X)
                             # expected distr
                             E           <- outer(rowSums(X), colSums(X), "*") / sum(X)
                             1 - sum(pmin(X / sum(X), E / sum(E)))
                                    }) , probs = c(.025, .975)) 
                }, .BxyES = BxyES)  )
# --------------------------------------------
# plot theta:
yearsES <- 1975:2009
yearsUS <- 1969:2009

pdf("/home/triffe/git/DISS/latex/Figures/TotalVariationObsvsExpectedUSES.pdf", height = 4.5, width = 4.5)
par(mai = c(.4,.4,.4,.2), xaxs = "i", yaxs = "i")
plot(yearsUS, TotalVarUS, type = 'n', ylim = c(.33,.48), xlim = c(1967, 2012), 
        axes = FALSE, xlab = "", ylab = "",
        panel.first = list(rect(1960, .3, 2012, .5, col = gray(.95), border = NA),
                           abline(v = seq(1970, 2010, by = 5),col = "white"),
                           abline(h = seq(.3, .45, by = .025),col = "white"),
                           text(1990,.318, "Year", xpd = TRUE),
                           text(1965,.49, expression(theta), xpd = TRUE),
                           text(seq(1970, 2010, by = 5), .33, seq(1970, 2010, by = 5), pos = 1 ,cex = .7, xpd = TRUE),
                           text(1967, seq(.35, .45, by = .05), seq(.35, .45, by = .05), pos = 2, cex = .7, xpd = TRUE)
                           ))
polygon(c(yearsUS,rev(yearsUS)), c(TotalVarUS95[,1],rev(TotalVarUS95[,2])), border = NA, col = gray(.6))
polygon(c(yearsES,rev(yearsES)), c(TotalVarES95[,1],rev(TotalVarES95[,2])), border = NA, col = gray(.6))         
lines(yearsUS, TotalVarUS, col = gray(.2), lwd = 2, lty = 1)                
lines(yearsES, TotalVarES, col = gray(.4), lwd = 3, lty = 5)
legend("bottomleft", col = gray(c(.2,.4)), lwd = c(2,3), lty = c(1,5),
    legend = c(expression(paste(theta," USA")), expression(paste(theta," ES"))), bty = "n")
dev.off()

# --------------------------------------------
# calculate strength of hypergamy, H, (ratio of total to structural hypergamy)
      
HypergamyES <- as.matrix(do.call(rbind,lapply(BxES, function(x){
                    x[is.na(x)] <- 0
                    rowTotals   <- rowSums(x)
                    colTotals   <- colSums(x)
                    nOfCases    <- sum(rowTotals)
                    # expected distr
                    E           <- outer(rowSums(x), colSums(x), "*") / nOfCases
                    # structural hypergamy:
                    sH <- sum(E[lower.tri(E)]) /
                            sum(E[upper.tri(E)])
                    # observed hypergamy:
                    tH <- sum(x[lower.tri(x)]) /
                            sum(x[upper.tri(x)])
                    # excess hypergamy: 
                    list(structural = sH, total = tH, excess = tH / sH, homogamy = sum(diag(x)) / sum(x))
                })  ))
HypergamyUS <- as.matrix(do.call(rbind,lapply(BxUS, function(x){
                    x[is.na(x)] <- 0
                    rowTotals   <- rowSums(x)
                    colTotals   <- colSums(x)
                    nOfCases    <- sum(rowTotals)
                    # expected distr
                    E           <- outer(rowSums(x), colSums(x), "*") / nOfCases
                    # structural hypergamy:
                    sH <- sum(E[lower.tri(E)]) /
                            sum(E[upper.tri(E)])
                    # observed hypergamy:
                    tH <- sum(x[lower.tri(x)]) /
                            sum(x[upper.tri(x)])
                    # excess hypergamy: 
                    list(structural = sH, total = tH, excess = tH / sH, homogamy = sum(diag(x)) / sum(x))
                })  ))
# --------------------------------------------
# plot H:
pdf("/home/triffe/git/DISS/latex/Figures/StrengthHypergamy.pdf", height = 4.5, width = 4.5)
par(mai = c(.4,.6,.4,.2), xaxs = "i", yaxs = "i")
USyrs <- as.integer(rownames(HypergamyUS))
ESyrs <- as.integer(rownames(HypergamyES))
plot(USyrs, HypergamyUS[,"structural"], type = 'l', ylim = c(0,8), xlim = c(1967,2012), 
        col = gray(.2), lty = 1, lwd = 2, axes = FALSE, xlab = "", ylab = "",
        panel.first = list(rect(1960,0,2012,9, col = gray(.95), border = NA),
                abline(v = seq(1970,2010,by = 5),col = "white"),
                abline(h = seq(0,8,by = 1),col = "white"),
                text(seq(1970,2010,by = 5),0,seq(1970,2010,by = 5),pos = 1,cex = .7, xpd = TRUE),
                text(1967,seq(0,8,by = 1),seq(0,8,by = 1),pos = 2,cex = .7, xpd = TRUE)
                
        ))
lines(USyrs, HypergamyUS[,"total"], col = gray(.1), lwd = 1, lty = 1)
lines(USyrs, HypergamyUS[,"excess"], col = gray(.2), lwd = 3, lty = 1)

lines(ESyrs, HypergamyES[,"structural"], col = gray(.2), lwd = 2, lty = 5)
lines(ESyrs, HypergamyES[,"total"], col = gray(.1), lwd = 1, lty = 5)
lines(ESyrs, HypergamyES[,"excess"], col = gray(.2), lwd = 3, lty = 5)

text(1966,1,"structural hypergamy", cex = .8, pos = 4)
text(1966,4.5,"excess hypergamy", cex = .8, pos = 4)
text(1972,7.7,"total observed hypergamy", cex = .8, pos = 4)

text(1962,4.5,expression(frac(B["x>y"], B["x<y"])),xpd = TRUE)

text(c(1970.984, 1968.685, 1969.160),c(6.955055, 4.072075, 1.549466), c("US","US","US"), cex = .7)
text(c(1979.548, 1981.213, 1986.447),c(6.868566, 3.754947, 1.95), c("ES","ES","ES"), cex = .7)
dev.off()
# --------------------------------------------










