

# from Mitra (1976)
# Mitra used 1966 US data to compare with Das Gupta. Can't reproduce his results
# but I can apply the method.
# ----------------------------------------------
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")
yearsES <- 1975:2009
yearsUS <- 1969:2009

BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata")))
names(BxymfES) <- yearsES
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))

# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
LxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmUS.Rdata"))) / 1e5
LxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfUS.Rdata"))) / 1e5
LxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmES.Rdata"))) / 1e5
LxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfES.Rdata"))) / 1e5


yearsES <- 1975:2009
yearsUS <- 1969:2009
ages    <- 0:110


# example year to get started:

# calculate ES and US single sex r estimates:

rmfUS <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .BxymfUS, .ExUS, .LxmUS, .LxfUS, age = 0:110){
                    
                    Bxym <- Mna0(.BxymfUS[[yr]][["Bxym"]])
                    Bxyf <- Mna0(.BxymfUS[[yr]][["Bxyf"]])
                    Exm  <- with(.ExUS, Male[Year == as.integer(yr)])
                    Eyf  <- with(.ExUS, Female[Year == as.integer(yr)])
                    Lxm  <- .LxmUS[, yr]
                    Lyf  <- .LxfUS[, yr]
                    mm   <- Mna0(rowSums(Bxym) / Exm)
                    mf   <- Mna0(colSums(Bxyf) / Eyf)
                    c(r.m = LotkaRCoale(mm, Lxm, age + .5),
                      r.f = LotkaRCoale(mf, Lyf, age + .5))
                }, .BxymfUS = BxymfUS, .ExUS = ExUS, .LxmUS = LxmUS, .LxfUS = LxfUS))
rownames(rmfUS) <- yearsUS

rmfES <- do.call(rbind,lapply(as.character(yearsES), function(yr, .BxymfES, .ExES, .LxmES, .LxfES, age = 0:110){        
                    Bxym <- Mna0(.BxymfES[[yr]][["Bxym"]])
                    Bxyf <- Mna0(.BxymfES[[yr]][["Bxyf"]])
                    Exm  <- with(.ExES, Male[Year == as.integer(yr)])
                    Eyf  <- with(.ExES, Female[Year == as.integer(yr)])
                    Lxm  <- .LxmES[, yr]
                    Lyf  <- .LxfES[, yr]
                    mm   <- Mna0(rowSums(Bxym) / Exm)
                    mf   <- Mna0(colSums(Bxyf) / Eyf)
                    c(r.m = LotkaRCoale(mm, Lxm, age + .5),
                      r.f = LotkaRCoale(mf, Lyf, age))
                }, .BxymfES = BxymfES, .ExES = ExES, .LxmES = LxmES, .LxfES = LxfES))
rownames(rmfES) <- yearsES

# plot it:
pdf("/home/triffe/git/DISS/latex/Figures/rmf.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .3, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rmfUS[, 1], type = 'l', ylim = c(-.02, .011), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968,-.02,2010,.011,col = gray(.95), border=NA),
                abline(h = seq(-.02,.011,by = .005), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.02,.011,by = .005),seq(-.02,.011,by = .005), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.02, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.0225, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,.012, "r", cex = 1, xpd = TRUE)))
lines(yearsUS, rmfUS[, 2], lwd = 2.5, col = gray(.5))
lines(yearsES, rmfES[, 1], lwd = 2, col = gray(.2), lty = 5)
lines(yearsES, rmfES[, 2], lwd = 2.5, col = gray(.5), lty = 5)

legend(1993,.009, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
dev.off()

fem.ind <- rmfUS[,1] < rmfUS[,2]
# 1994 - 1996 ; 2001
all(rmfES[,1] > rmfES[,2])

yearsUS[apply(rmfUS, 1, function(x){sign(x[1]) != sign(x[2])})]
yearsES[apply(rmfES, 1, function(x){sign(x[1]) != sign(x[2])})]

# ------------------------------------------------------------
# years to double?


dxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata"))) / 1e5
dxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata"))) / 1e5
dxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmES.Rdata"))) / 1e5
dxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfES.Rdata"))) / 1e5

pxmUS <- apply(dxmUS, 2, dx2pxLes)
pxfUS <- apply(dxfUS, 2, dx2pxLes)
pxmES <- apply(dxmES, 2, dx2pxLes)
pxfES <- apply(dxfES, 2, dx2pxLes)

# spit back years to observed doubling
yrs2double <- compiler::cmpfun(function(Exm, Eyf, Lf, Lm, maxit = 5000){
    Exmi  <- Exm
    Eyfi  <- Eyf
    SR    <- sum(Exmi) / sum(Eyfi) 
    yrt   <- 0
    while (SR < 2 & SR > .5 & yrt < maxit){
        Exmi <- c(Lm %*% Exmi)
        Eyfi <- c(Lf %*% Eyfi)
        SR   <- sum(Exmi) / sum(Eyfi) 
        yrt  <- yrt + 1
    }
    yrt
})

doubleES <- unlist(lapply(as.character(yearsES), function(yr, .BxymfES, .ExES, .pxmES, .pxfES, age = 0:110){        
                    Bxym <- Mna0(.BxymfES[[yr]][["Bxym"]])
                    Bxyf <- Mna0(.BxymfES[[yr]][["Bxyf"]])
                    Exm  <- with(.ExES, Male[Year == as.integer(yr)])
                    Eyf  <- with(.ExES, Female[Year == as.integer(yr)])

                    mm   <- Mna0(rowSums(Bxym) / Exm)
                    mf   <- Mna0(colSums(Bxyf) / Eyf)
                    
                    yrs2double(Exm, Eyf, Lf = Leslie(mf, pxfUS[,yr]), Lm = Leslie(mm, pxmUS[,yr]), maxit = 500000)
                }, .BxymfES = BxymfES, .ExES = ExES, .pxmES = pxmES, .pxfES = pxfES))
doubleUS <- unlist(lapply(as.character(yearsUS), function(yr, .BxymfUS, .ExUS, .pxmUS, .pxfUS, age = 0:110){        
                    Bxym <- Mna0(.BxymfUS[[yr]][["Bxym"]])
                    Bxyf <- Mna0(.BxymfUS[[yr]][["Bxyf"]])
                    Exm  <- with(.ExUS, Male[Year == as.integer(yr)])
                    Eyf  <- with(.ExUS, Female[Year == as.integer(yr)])
                    
                    mm   <- Mna0(rowSums(Bxym) / Exm)
                    mf   <- Mna0(colSums(Bxyf) / Eyf)
                    
                    yrs2double(Exm, Eyf, Lf = Leslie(mf, pxfUS[,yr]), Lm = Leslie(mm, pxmUS[,yr]), maxit = 500000)
                }, .BxymfUS = BxymfUS, .ExUS = ExUS, .pxmUS = pxmUS, .pxfUS = pxfUS))

# plot it (log scale)
TicksMaj <- 10 ^ (2:5)
TicksMagLab <- c("100","1000","10000","10000")
TickMin     <- log(c(t(outer(TicksMaj, 2:9))))

pdf("/home/triffe/git/DISS/latex/Figures/rSRdoubling.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, log(doubleUS), type = 'l', ylim = log(c(100,150000)), xlim = c(1968, 2010),
        axes = FALSE, xlab = "", ylab = "", col = gray(.2), lwd = 2, 
        panel.first = list(rect(1968, log(100), 2010, log(150000), col = gray(.95), border=NA),
                abline(h = log(TicksMaj), col = "white"),
                segments(1968, TickMin, 1969, TickMin,col = "white"),
                segments(2009, TickMin, 2010, TickMin,  col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, log(TicksMaj), TicksMagLab, pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10), log(100), seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, 4.1, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1962, 12.4, "years to\nSR > 2 or < .5", cex = 1, xpd = TRUE, pos = 4)))
lines(yearsES, log(doubleES), lty = 1, col = gray(.5), lwd = 2.5)
points(yearsUS[fem.ind], log(doubleUS)[fem.ind], pch = 9, col = "black", xpd = TRUE)
legend(1970,11.7, lty = c(1,1,NA), col = gray(c(.2,.5,0)), lwd = c(2,2.5, NA), pch = c(NA,NA,9),bty = "n",
        legend = c("US", "ES", expression(r^f > r^m)), xpd = TRUE)
dev.off()
