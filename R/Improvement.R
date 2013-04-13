
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")
yearsUS <- 1969:2009
yearsES <- 1975:2009

dxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfES.Rdata"))) 
# make sum to 1
dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)

mxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxmUS.Rdata"))) 
mxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxfUS.Rdata"))) 
mxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxmES.Rdata"))) 
mxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxfES.Rdata"))) 

PxUS  <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxUS.Rdata")))
PxES  <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxES.Rdata")))

iota <- .997
iota ^ c(1:110)
yr <- "1975"
wmean(.5:110.5,mx2dxHMD(mxmUS[,yr]*iota ^ c(1:111)))
wmean(.5:110.5,mx2dxHMD(mxmUS[,yr]))
mx <- mxmUS[,"1975"]
Px <- with(PxUS, Male[Year == 1975])
ExpectedDx2 <- function(Px, mx, iota = 1){
    impr     <- iota ^ c(1:111)
    
    N        <- length(mx)
    EDx      <- matrix(0, nrow = N, ncol = N, dimnames = list(Ex = 0:(N-1), Age =  0:(N-1)))
    # Population age loop
    
    for (i in 1:110){
        
        dxn  <- mx2dxHMD(mx * c(rep(1,(i)),impr[1:(N-i)]))[i:N]
        dxn  <- dxn / sum(dxn)
        # distribute each age of Populatin over death times
        EDx[1:length(dxn), i]    <- Px[i] * dxn
    }
    EDx[1,N] <- Px[N]
    EDx[is.na(EDx)] <- 0
    EDx
}
Pxm <- with(PxUS, Male[Year == 2009])
Pxf <- with(PxUS, Female[Year == 2009])
xlabs   <- c("1.0%","0.8%","0.6%","0.4%","0.2%","0.0%","0.2%","0.4%","0.6%","0.8%","1.0%")

Males2009 <- rowSums(ExpectedDx2(Pxm, mxmUS[,"2009"]))
Females2009 <- rowSums(ExpectedDx2(Pxf, mxfUS[,"2009"]))
pdf("/home/triffe/git/DISS/latex/Figures/exPyramidUSimpr.pdf", height = 5, width = 5)
par(mai = c(.6,.6,.3,.3), xaxs = "i", yaxs = "i")
plot(NULL, type = "n",axes = FALSE, xlab = "",ylab = "", xlim = c(-1, 1), ylim = c(0,111),
        panel.first = list(
                rect(-1, 0, 1, 111, col = gray(.95), border = NA),
                abline(v = seq(-1, 1, by = .2), col = "white"),
                abline(h = seq(0, 110, by = 10), col = "white"),
                text(seq(-1, 1, by = .2),0, xlabs, xpd = TRUE, pos = 1, cex = .7),
                text(-1, seq(0, 110, by = 10), seq(0, 110, by = 10), pos = 2, xpd = TRUE, cex = .7),
                text(-1.17, 116, expression(e[y]), xpd = TRUE, cex = 1),
                text(0, -12, "Percentage", xpd = TRUE, cex = 1)
        ))
barplot(-100 * (Males2009 / sum(Males2009 + Females2009)), border = NA, col = "#44444450", 
        add = TRUE, horiz = TRUE, space = 0, axes = FALSE,axisnames=FALSE)
barplot(100 * (Females2009 / sum(Males2009 + Females2009)), border = NA, col = "#44444450", 
        add = TRUE, horiz = TRUE, space = 0, axes = FALSE,axisnames=FALSE) 
PyramidOutline(rowSums(ExpectedDx2(Pxm, mxmUS[,"2009"],.995)), 
        rowSums(ExpectedDx2(Pxf, mxfUS[,"2009"],.993)), scale =100, border = gray(.2), xpd = TRUE, lwd = .5) 
text(c(0.27, 0.44), c(50, 100), c("2009 fixed", expression(iota==0.995)), 
        col = c("white", "black"), cex = 1.2)
segments(0.2439776, 98.40316, 0.1937150, 95.93115)
dev.off()

# now spain 2009
Pxm <- with(PxES, Male[Year == 2009])
Pxf <- with(PxES, Female[Year == 2009])
Males2009 <- rowSums(ExpectedDx2(Pxm, mxmES[,"2009"]))
Females2009 <- rowSums(ExpectedDx2(Pxf, mxfES[,"2009"]))
pdf("/home/triffe/git/DISS/latex/Figures/exPyramidESimpr.pdf", height = 5, width = 5)
par(mai = c(.6,.6,.3,.3), xaxs = "i", yaxs = "i")
plot(NULL, type = "n",axes = FALSE, xlab = "",ylab = "", xlim = c(-1, 1), ylim = c(0,111),
        panel.first = list(
                rect(-1, 0, 1, 111, col = gray(.95), border = NA),
                abline(v = seq(-1, 1, by = .2), col = "white"),
                abline(h = seq(0, 110, by = 10), col = "white"),
                text(seq(-1, 1, by = .2),0, xlabs, xpd = TRUE, pos = 1, cex = .7),
                text(-1, seq(0, 110, by = 10), seq(0, 110, by = 10), pos = 2, xpd = TRUE, cex = .7),
                text(-1.17, 116, expression(e[y]), xpd = TRUE, cex = 1),
                text(0, -12, "Percentage", xpd = TRUE, cex = 1)
        ))
barplot(-100 * (Males2009 / sum(Males2009 + Females2009)), border = NA, col = "#44444450", 
        add = TRUE, horiz = TRUE, space = 0, axes = FALSE,axisnames=FALSE)
barplot(100 * (Females2009 / sum(Males2009 + Females2009)), border = NA, col = "#44444450", 
        add = TRUE, horiz = TRUE, space = 0, axes = FALSE,axisnames=FALSE) 
PyramidOutline(rowSums(ExpectedDx2(Pxm, mxmES[,"2009"],.995)), 
        rowSums(ExpectedDx2(Pxf, mxfES[,"2009"],.993)), scale =100, border = gray(.2), xpd = TRUE, lwd = .5) 
text(c(0.27, 0.45), c(50, 100), c("2009 fixed", expression(iota==0.995)), 
        col = c("white", "black"), cex = 1.2)
segments(0.2774860, 97.64254, 0.2037675, 93.83945)
dev.off()