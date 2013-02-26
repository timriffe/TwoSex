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

ages <- 0:110
# example for 1975 US:

Mex75US <- rowSums(ExpectedDx(Px = with(ExUS, Male[Year == 1970]), dx = dxmUS[,"1970"]))
Fex75US <- rowSums(ExpectedDx(Px = with(ExUS, Female[Year == 1970]), dx = dxfUS[,"1970"]))

# tricky. Matrix must have male ages in rows, female ages in columns

ExBxy <- ExpectedDxMxFmatrix( BxUS[["1970"]], dxmUS[,"1970"], dxfUS[,"1970"])

expected    <- outer(rowSums(ExBxy), colSums(ExBxy), "*") / sum(ExBxy)

Fxym <- ExBxy %row% Mex75US
Fxyf <- ExBxy %col% Fex75US

# let's do an image, same 1975
colramp <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"), space = "Lab")


# plotting pars
brks <- seq(min(ExBxy, na.rm = TRUE), max(ExBxy, na.rm = TRUE), length.out = 51)
#brks <- seq(min(lExBxy, na.rm = TRUE), max(lExBxy, na.rm = TRUE), length.out = 51)
g.xy <- seq(0, 100, by = 5)
gb.xy <- seq(0, 100, by = 10)
levs <- c(100,seq(500,3000,by=500))     # for contour plot

ExBxy[ExBxy == 0]       <- NA
expected[expected == 0] <- NA
pdf("/home/triffe/git/DISS/latex/Figures/ObservedvsExpectedBexey.pdf", 
        height = 4, width = 6.5)

par(mfrow=c(1,2), mar = c(3,1,2,3))
image(x = ages + .5, y = ages + .5, t(ExBxy), 
        xlim = c(0, 101), ylim = c(0, 101), zlim = c(0,3100),
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1, 
        panel.first = list(rect(0, 0, 101, 101, col = "#EEEEEE", xpd = TRUE, border = NA), 
                abline(h = g.xy, col = "white", lwd = .5),
                abline(v = g.xy, col = "white", lwd = .5),
                text(2, gb.xy, gb.xy, pos = 2, cex = .5, xpd = TRUE),
                text(gb.xy, 0, gb.xy, pos = 1, cex = .5, xpd = TRUE),
                segments(0, gb.xy, -1, gb.xy, xpd = TRUE),
                segments(gb.xy, 0, gb.xy, -1, xpd = TRUE)),
        
        xlab = "", ylab = "")
# contours
contour(x = ages + .5, y = ages + .5, t(ExBxy), 
        levels = levs, labels = levs, add = TRUE)
# legend
leg.y <- seq(0,11,length = 50) * (101 / 11) 
rect(105,leg.y[1:(length(leg.y) - 1)],113, leg.y[2:length(leg.y)],
        col = colramp(50), border = NA, xpd = TRUE)
minmax <- range(brks)
TickMaj <- seq(0,3000,by = 200)
TickMajsc <- (TickMaj / diff(range(minmax))) * diff(range(leg.y)) 
segments(113,TickMajsc,113.5,TickMajsc,col = "#444444", xpd = TRUE)
text(113,TickMajsc, TickMaj, pos = 4, xpd = TRUE, cex = .5)

# line of homogamy:
segments(0,0,111,111,col = "#50505050")

# axis labels
fath <- "Father"
moth <- "Mother"
text(50,-10, bquote(.(fath) ~ e[x]), xpd = TRUE, cex = .7, pos =1)
text(-10,110,bquote(.(moth) ~ e[x]), xpd = TRUE, pos = 4, cex = .7)

# expected bivariate distribution:

par(mar=c(3,2,2,2))
image(x = ages + .5, y = ages + .5, t(expected), 
        xlim = c(0, 101), ylim = c(0, 101), zlim = c(0,3100),
        col = colramp(50), breaks = brks, axes = FALSE, asp = 1, 
        panel.first = list(rect(0, 0, 101, 101, col = "#EEEEEE", xpd = TRUE, border = NA), 
                abline(h = g.xy, col = "white", lwd = .5),
                abline(v = g.xy, col = "white", lwd = .5),
                text(2, gb.xy, gb.xy, pos = 2, cex = .5, xpd = TRUE),
                text(gb.xy, 0, gb.xy, pos = 1, cex = .5, xpd = TRUE),
                segments(0, gb.xy, -1, gb.xy, xpd = TRUE),
                segments(gb.xy, 0, gb.xy, -1, xpd = TRUE)),
        xlab = "", ylab = "")
# axis labels
text(50,-10, bquote(.(fath) ~ e[x]), xpd = TRUE, cex = .7, pos = 1)
text(-15,110,bquote(.(moth) ~ e[x]), xpd = TRUE, pos = 4, cex = .7)

# contours
 
contour(x = ages + .5, y = ages + .5, t(expected), 
        levels = levs, labels = levs, add = TRUE)
# line of homogamy:
segments(0,0,101,101,col = "#50505050")
dev.off()

# observed vs expected:


#image(x = ages + .5, y = ages + .5, t(expected), ylim = c(0,111),xlim = c(0,111))
#image(x = ages + .5, y = ages + .5, t( expected / sum(expected)) - t(ExBxy / sum(ExBxy)) , ylim = c(0,111),xlim = c(0,111))
# total variation distance:

ExBxy <- ExpectedDxMxFmatrix( BxUS[["1970"]], dxmUS[,"1970"], dxfUS[,"1970"])

expected    <- outer(rowSums(ExBxy), colSums(ExBxy), "*") / sum(ExBxy)

# takes a minute or two to run...
ExBxyAallUS <- lapply(as.character(yearsUS), function(yr, .BxUS, .dxmUS, .dxfUS){
            ExpectedDxMxFmatrix( .BxUS[[yr]], .dxmUS[,yr], .dxfUS[,yr])
        }, .BxUS = BxUS, .dxmUS = dxmUS, .dxfUS = dxfUS)
ExBxyAallES <- lapply(as.character(yearsES), function(yr, .BxES, .dxmES, .dxfES){
            ExpectedDxMxFmatrix( .BxES[[yr]], .dxmES[,yr], .dxfES[,yr])
        }, .BxES = BxES, .dxmES = dxmES, .dxfES = dxfES)
# animate run-through
#lapply(ExBxyAall, function(ExBxy){
#            brks <- seq(min(ExBxy, na.rm = TRUE), max(ExBxy, na.rm = TRUE), length.out = 51)
##brks <- seq(min(lExBxy, na.rm = TRUE), max(lExBxy, na.rm = TRUE), length.out = 51)
#            g.xy <- seq(0, 100, by = 5)
#            gb.xy <- seq(0, 100, by = 10)
#            levs <- c(100,seq(500,3500,by=500))     # for contour plot
#            
#            ExBxy[ExBxy == 0]       <- NA
#        image(x = ages + .5, y = ages + .5, t(ExBxy), 
#        xlim = c(0, 101), ylim = c(0, 101), zlim = c(0,3500),
#        col = colramp(50), breaks = brks, axes = FALSE, asp = 1, 
#        panel.first = list(rect(0, 0, 101, 101, col = "#EEEEEE", xpd = TRUE, border = NA), 
#                abline(h = g.xy, col = "white", lwd = .5),
#                abline(v = g.xy, col = "white", lwd = .5),
#                text(2, gb.xy, gb.xy, pos = 2, cex = .5, xpd = TRUE),
#                text(gb.xy, 0, gb.xy, pos = 1, cex = .5, xpd = TRUE),
#                segments(0, gb.xy, -1, gb.xy, xpd = TRUE),
#                segments(gb.xy, 0, gb.xy, -1, xpd = TRUE)),
#        
#        xlab = "", ylab = "")
## contours
#       contour(x = ages + .5, y = ages + .5, t(ExBxy), 
#        levels = levs, labels = levs, add = TRUE)
#Sys.sleep(1)
#        })

TotalVarUS <- unlist(lapply(ExBxyAallUS, function(.ExBxy){
            expected    <- outer(rowSums(.ExBxy), colSums(.ExBxy), "*") / sum(.ExBxy)
            sum(abs(.ExBxy / sum(.ExBxy) - expected / sum(expected))) / 2
        }))
TotalVarES <- unlist(lapply(ExBxyAallES, function(.ExBxy){
                    expected    <- outer(rowSums(.ExBxy), colSums(.ExBxy), "*") / sum(.ExBxy)
                    sum(abs(.ExBxy / sum(.ExBxy) - expected / sum(expected))) / 2
                }))

# plot it
pdf("/home/triffe/git/DISS/latex/Figures/TotalVariationObsvsExpectedexUSES.pdf", height = 4.5, width = 4.5)
par(mai = c(.5,.4,.4,.2), xaxs = "i", yaxs = "i")
plot(yearsUS, TotalVarUS, type = 'l', ylim = c(.04,.07), xlim = c(1968,2010), 
        col = gray(.2), lwd = 2, axes = FALSE, xlab = "", ylab = "",
        panel.first = list(rect(1968,.04,2010,.07, col = gray(.95), border = NA),
                abline(v = seq(1970,2010,by = 5),col = "white"),
                abline(h = seq(.04,.07,by = .005),col = "white"),
                text(1990, .037, "Year", xpd = TRUE),
                text(1966, .072, expression(theta), xpd = TRUE),
                text(seq(1970,2010,by = 5),.04,seq(1970,2010,by = 5),pos = 1,cex = .7, xpd = TRUE),
                text(1968,seq(.04,.07,by = .005),seq(.04,.07,by = .005), pos = 2,cex = .7, xpd = TRUE)
        ))
lines(yearsES, TotalVarES, col = gray(.4), lwd = 3, lty = 5)
legend("topright", col = gray(c(.2,.4)), lwd = c(2,3), lty = c(1,5),
        legend = c(expression(paste(theta," USA")), expression(paste(theta," ES"))), bty = "n")
dev.off()










(ToalVar <- sum(abs(ExBxy / sum(ExBxy) - expected / sum(expected))) / 2 )

# overlap 
(coefOverlap     <- sum(pmin(ExBxy / sum(ExBxy), expected / sum(expected))))

# much more perfect distribution.


















