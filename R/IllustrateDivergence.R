
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
# BxES is 0:110, years 1975:2009
BxUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy10_65.Rdata")))
BxES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65
# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))

yearsES <- 1975:2009
yearsUS <- 1969:2009
ages    <- 10:65


# cut ES to same dimentions ages 10:65, square (summing boundary ages)
BxES <- lapply(BxES, function(x){
            x[11, ] <- colSums(x[1:11, ])
            x[, 11] <- rowSums(x[, 1:11])
            x[66, ] <- colSums(x[66:110, ])
            x[, 66] <- rowSums(x[, 66:110])
            x[11:66,11:66]
        })

# get years for each pop

# ASFR sample year 1975:

BxyUS <- BxUS[["1975"]]
BxyES <- BxES[["1975"]]

ExmUS <- with(ExUS, Male[Year == 1975 & Age >= 10 & Age <= 65])
ExfUS <- with(ExUS, Female[Year == 1975 & Age >= 10 & Age <= 65])
ExmES <- with(ExES, Male[Year == 1975 & Age >= 10 & Age <= 65])
ExfES <- with(ExES, Female[Year == 1975 & Age >= 10 & Age <= 65])

ASFRmUS <- rowSums(BxyUS) / ExmUS
ASFRfUS <- colSums(BxyUS) / ExfUS
ASFRmES <- rowSums(BxyES) / ExmES
ASFRfES <- colSums(BxyES) / ExfES

# max(c(ASFRmUS,ASFRfUS,ASFRmES,ASFRfES))

pdf("/home/triffe/git/DISS/latex/Figures/ASFR1975.pdf", height = 5, width = 5)
par(mai = c(.5,.5,.3,.3),xaxs = "i", yaxs = "i")
plot(ages, ASFRmUS, type = 'l', ylim = c(0, .21), xlim = c(10,66), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(10,0,66,.21,col = gray(.95), border=NA),
                abline(h = seq(0,.2,by = .02), col = "white"),
                abline(v = seq(10, 65, by = 5), col = "white"),
                text(10, seq(0, .2, by = .02),seq(0,.2, by = .02),pos = 2, cex = .8, xpd = TRUE),
                text(seq(10, 60, by = 10),0, seq(10, 60, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(35, -.01, "Age", cex = 1, pos = 1, xpd = TRUE),
                text(8,.22, "ASFR", cex = 1, xpd = TRUE)))
lines(ages, ASFRfUS, lwd = 2.5, col = gray(.5))
lines(ages, ASFRmES, lwd = 2, col = gray(.2), lty = 5)
lines(ages, ASFRfES, lwd = 2.5, col = gray(.5), lty = 5)

legend(42,.21, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
dev.off()

# --------------------------------------------------------------
# simple series of TFR
TFRES <- as.matrix(do.call(rbind,lapply(as.character(yearsES), function(yr, .BxES, .ExES){
            Exm <- with(.ExES, Male[Year == yr & Age >= 10 & Age <= 65])
            Exf <- with(.ExES, Female[Year == yr & Age >= 10 & Age <= 65])
            Bxy <- .BxES[[yr]]
            c(TFRm = sum(Bxy %row% Exm, na.rm = TRUE), TFRf = sum(Bxy %col% Exf, na.rm = TRUE))
        }, .ExES = ExES, .BxES = BxES)))
TFRUS <- as.matrix(do.call(rbind,lapply(as.character(yearsUS), function(yr, .BxUS, .ExUS){
            Exm <- with(.ExUS, Male[Year == yr & Age >= 10 & Age <= 65])
            Exf <- with(.ExUS, Female[Year == yr & Age >= 10 & Age <= 65])
            Bxy <- .BxUS[[yr]]
            c(TFRm = sum(Bxy %row% Exm, na.rm = TRUE), TFRf = sum(Bxy %col% Exf, na.rm = TRUE))
        }, .ExUS = ExUS, .BxUS = BxUS)))

# plot it, save out
pdf("/home/triffe/git/DISS/latex/Figures/TFR.pdf", height = 5, width = 5)
par(mai = c(.5, .3, .3, .3),xaxs = "i", yaxs = "i")
plot(yearsUS, TFRUS[, 1], type = 'l', ylim = c(1, 2.9), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968,1,2010,2.9,col = gray(.95), border=NA),
                abline(h = seq(0,2.9,by = .2), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(1, 2.9, by = .2),seq(1, 2.9, by = .2), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),1, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, .88, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1967,3, "TFR", cex = 1, xpd = TRUE)))
lines(yearsUS, TFRUS[, 2], lwd = 2.5, col = gray(.5))
lines(yearsES, TFRES[, 1], lwd = 2, col = gray(.2), lty = 5)
lines(yearsES, TFRES[, 2], lwd = 2.5, col = gray(.5), lty = 5)

legend(1993,2.85, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
dev.off()

# --------------------------------------------------------------
# simple series of NRR male and female

BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf10_65.Rdata")))
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf10_65.Rdata")))
names(BxymfES) <- yearsES
# get Lx from HMD, divide by l0 (100000)

# these are age by year matrices:
LxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmUS.Rdata")))[11:66, ] / 1e5
LxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfUS.Rdata")))[11:66, ] / 1e5
LxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmES.Rdata")))[11:66, ] / 1e5
LxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfES.Rdata")))[11:66, ] / 1e5

# calculate male-male and female-female fertility:
FxmfUS <- lapply(as.character(yearsUS), function(yr, .BxymfUS, .ExUS){ 
          
                            Exm <- with(.ExUS, Male[Year == yr & Age >= 10 & Age <= 65])
                            Exf <- with(.ExUS, Female[Year == yr & Age >= 10 & Age <= 65])
                            
                            Bxym <- .BxymfUS[[yr]][["Bxym"]]
                            Bxyf <- .BxymfUS[[yr]][["Bxyf"]]
                            
                            list(ASFRm = rowSums(Bxym %row% Exm),
                                 ASFRf = colSums(Bxyf %col% Exf))
                        }, .ExUS = ExUS, .BxymfUS = BxymfUS)
               
names(FxmfUS) <- yearsUS
FxmfES <- lapply(as.character(yearsES), function(yr, .BxymfES, .ExES){ 
                            Exm <- with(.ExES, Male[Year == yr & Age >= 10 & Age <= 65])
                            Exf <- with(.ExES, Female[Year == yr & Age >= 10 & Age <= 65])
                            
                            Bxym <- Mna0(.BxymfES[[yr]][["Bxym"]])
                            Bxyf <- Mna0(.BxymfES[[yr]][["Bxyf"]])
                            
                            list(ASFRm = rowSums(Bxym %row% Exm),
                                    ASFRf = colSums(Bxyf %col% Exf))
                        }, .ExES = ExES, .BxymfES = BxymfES)
names(FxmfES) <- yearsES     
R0mfUS <- as.matrix(do.call(rbind, lapply(as.character(yearsUS), function(yr, .FxmfUS, .LxmUS, .LxfUS){
            ASFRm <- .FxmfUS[[yr]][["ASFRm"]]
            ASFRf <- .FxmfUS[[yr]][["ASFRf"]]
            Lxm   <- .LxmUS[, yr]
            Lxf   <- .LxfUS[, yr]
            c(R0m = sum(ASFRm * Lxm), R0f = sum(ASFRf * Lxf))
        },.FxmfUS = FxmfUS, .LxmUS = LxmUS, .LxfUS = LxfUS)))
rownames(R0mfUS) <- yearsUS

R0mfES <- as.matrix(do.call(rbind, lapply(as.character(yearsES), function(yr, .FxmfES, .LxmES, .LxfES){
            ASFRm <- .FxmfES[[yr]][["ASFRm"]]
            ASFRf <- .FxmfES[[yr]][["ASFRf"]]
            Lxm   <- .LxmES[, yr]
            Lxf   <- .LxfES[, yr]
            c(R0m = sum(ASFRm * Lxm), R0f = sum(ASFRf * Lxf))
        },.FxmfES = FxmfES, .LxmES = LxmES, .LxfES = LxfES)))
rownames(R0mfES) <- yearsES

# plot it, save out:
pdf("/home/triffe/git/DISS/latex/Figures/R0mf.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .3, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, R0mfUS[, 1], type = 'l', ylim = c(.5, 1.45), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968,.5,2010,1.45,col = gray(.95), border=NA),
                abline(h = seq(.5,1.4,by = .1), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(.5, 1.4, by = .1),seq(.5, 1.4, by = .1), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),.5, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, .45, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,1.5, expression(R[0]), cex = 1, xpd = TRUE)))
lines(yearsUS, R0mfUS[, 2], lwd = 2.5, col = gray(.5))
lines(yearsES, R0mfES[, 1], lwd = 2, col = gray(.2), lty = 5)
lines(yearsES, R0mfES[, 2], lwd = 2.5, col = gray(.5), lty = 5)

legend(1993,1.45, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
dev.off()

# -----------------------------------------------------------
# SRB

calcSRBmf <- function(Bxy, sex = "m"){
    if (sex == "m"){
        SRB <- MinfNA(rowSums(Bxy[["Bxym"]]) / rowSums(Bxy[["Bxyf"]]))
    } else {
        SRB <- MinfNA(colSums(Bxy[["Bxym"]]) / colSums(Bxy[["Bxyf"]]))
    }
    SRB
}



BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))
names(BxymfES) <- yearsES
SRBfUS <- calcSRBmf(BxymfUS[["1975"]], sex ="f")
SRBmUS <- calcSRBmf(BxymfUS[["1975"]], sex ="m")
SRBfES <- calcSRBmf(BxymfES[["1975"]], sex ="f")
SRBmES <- calcSRBmf(BxymfES[["1975"]], sex ="m")

pdf("/home/triffe/git/DISS/latex/Figures/SRB1975.pdf", height = 5, width = 5)
age <- .5:110.5
ind <- age >14 & age < 50
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
ylim <- c(.8,1.25)
plot(NULL, type = 'n', ylim = ylim, xlim = range(age[ind]), axes = FALSE,
       xlab = "", ylab = "",
        panel.first = list(rect(10,ylim[1],65,ylim[2],col = gray(.95), border=NA),
                abline(h = seq(ylim[1], ylim[2], by = .05), col = "white"),
                abline(v = seq(15, 50, by = 5), col = "white"),
                text(min(age[ind]), seq(ylim[1], ylim[2], by = .05), seq(ylim[1], ylim[2], by = .05), pos = 2, cex = .8, xpd = TRUE),
                text(seq(15, 50, by = 5),ylim[1], seq(15, 50, by = 5), pos = 1, cex = .8, xpd = TRUE),
                text(32, .77, "Age", cex = 1, pos = 1, xpd = TRUE),
                text(13,1.29, "SRB", cex = 1, xpd = TRUE)))
lines(age[ind], SRBmUS[ind], lwd = 2, col = gray(.2))
lines(age[ind], SRBfUS[ind], lwd = 2.5, col = gray(.5))
lines(age[ind], SRBmES[ind], lwd = 2, col = gray(.2), lty = 5)
lines(age[ind], SRBfES[ind], lwd = 2.5, col = gray(.5), lty = 5)

legend(15,.9, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
dev.off()

# total SRB over time

SRBUS <- unlist(lapply(as.character(yearsUS), function(yr, .Bxymf){
             sum(.Bxymf[[yr]][["Bxym"]]) / sum(.Bxymf[[yr]][["Bxyf"]]) 
        }, .Bxymf = BxymfUS))
SRBES <- unlist(lapply(as.character(yearsES), function(yr, .Bxymf){
                    sum(.Bxymf[[yr]][["Bxym"]]) / sum(.Bxymf[[yr]][["Bxyf"]]) 
                }, .Bxymf = BxymfES))

pdf("/home/triffe/git/DISS/latex/Figures/SRByear.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
ylim <- c(1.04,1.1)
plot(NULL, type = 'n', ylim = ylim, xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,ylim[1],2010,ylim[2],col = gray(.95), border=NA),
                abline(h = seq(ylim[1], ylim[2], by = .01), col = "white"),
                abline(v = seq(1970, 2005, by = 5), col = "white"),
                text(1968, seq(ylim[1], ylim[2], by = .01), seq(ylim[1], ylim[2], by = .01), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2005, by = 5),ylim[1], seq(1970, 2005, by = 5), pos = 1, cex = .8, xpd = TRUE),
                text(1990, 1.036, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,1.105, "SRB", cex = 1, xpd = TRUE)))
lines(yearsUS, SRBUS, lwd = 2, col = gray(.2))
lines(yearsES, SRBES, lwd = 2.5, col = gray(.5), lty=5)
legend(1995,1.1, lty = c(1,5), col = gray(c(.2,.5)), lwd = c(2,2.5),bty = "n",
        legend = c("US", "Spain"), xpd = TRUE)
dev.off()


# --------------------------------------------------------------------
# reproductive spans, self-sufficient code
yearsES <- 1975:2009
yearsUS <- 1969:2009
BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))
names(BxymfES) <- yearsES
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))



Bounds99US <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxy, .Ex, .a = .5:110.5){
           
            Fxm <- Mna0(Minf0(rowSums(.Bxy[[yr]][[1]] + .Bxy[[yr]][[1]] ) / 
                                    with(.Ex, Male[Year == as.integer(yr)])))
            Fxf <- Mna0(Minf0(colSums(.Bxy[[yr]][[1]] + .Bxy[[yr]][[1]] ) / 
                                    with(.Ex, Female[Year == as.integer(yr)])))
            
            keep <- .a <= 85        
            
            .anew     <- seq(min(.a[keep]),max(.a[keep]),by = .01)
            Fxminterp <- approx(x=.a[keep], y=Fxm[keep], xout =  .anew)
            Fxfinterp <- approx(x=.a[keep], Fxf[keep], xout =  .anew)

            csfFxm    <- cumsum(Fxminterp$y / sum(Fxminterp$y))
            csfFxf    <- cumsum(Fxfinterp$y / sum(Fxfinterp$y))
            
            c(m.lower = rev(.anew[csfFxm <= .005])[1],
            m.upper = rev(.anew[csfFxm <= .995])[1],
            m.med = rev(.anew[csfFxm <= .5])[1],
            f.lower = rev(.anew[csfFxf <= .005])[1],
            f.upper = rev(.anew[csfFxf <= .995])[1],
            f.med = rev(.anew[csfFxf <= .5])[1])
            
        },  .Bxy = BxymfUS, .Ex = ExUS))
Bounds99ES <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxy, .Ex, .a = .5:110.5){
            
            Fxm <- Mna0(Minf0(rowSums(.Bxy[[yr]][[1]] + .Bxy[[yr]][[1]] ) / 
                                            with(.Ex, Male[Year == as.integer(yr)])))
            Fxf <- Mna0(Minf0(colSums(.Bxy[[yr]][[1]] + .Bxy[[yr]][[1]] ) / 
                                            with(.Ex, Female[Year == as.integer(yr)])))
            keep <- .a <= 85        
                    
            .anew     <- seq(min(.a[keep]),max(.a[keep]),by = .01)
            Fxminterp <- approx(x=.a[keep], y=Fxm[keep], xout =  .anew)
            Fxfinterp <- approx(x=.a[keep], Fxf[keep], xout =  .anew)
                    
            csfFxm    <- cumsum(Fxminterp$y / sum(Fxminterp$y))
            csfFxf    <- cumsum(Fxfinterp$y / sum(Fxfinterp$y))
                    
            c(m.lower = rev(.anew[csfFxm <= .005])[1],
                    m.upper = rev(.anew[csfFxm <= .995])[1],
                    m.med = rev(.anew[csfFxm <= .5])[1],
                    f.lower = rev(.anew[csfFxf <= .005])[1],
                    f.upper = rev(.anew[csfFxf <= .995])[1],
                    f.med = rev(.anew[csfFxf <= .5])[1])
                    
                }, .Bxy = BxymfES, .Ex = ExES))


pdf("/home/triffe/git/DISS/latex/Figures/ASFRbounds.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .5), xaxs = "i", yaxs = "i")
ylim <- c(10,67)
plot(NULL, type = 'n', ylim = ylim, xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,ylim[1],2010,ylim[2],col = gray(.95), border=NA),
                abline(h = seq(ylim[1], ylim[2], by = 5), col = "white"),
                abline(v = seq(1970, 2005, by = 5), col = "white"),
                text(1968, seq(ylim[1], ylim[2], by = 5), seq(ylim[1], ylim[2], by = 5), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2005, by = 5),ylim[1], seq(1970, 2005, by = 5), pos = 1, cex = .8, xpd = TRUE),
                text(1988, 7, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,71, "Age", cex = 1, xpd = TRUE)))
lines(yearsUS, Bounds99US[, "f.lower"], lwd = 1.5, col = gray(.1))
lines(yearsUS, Bounds99US[, "f.upper"], lwd = 1.5, col = gray(.1))
lines(yearsUS, Bounds99US[, "f.med"], lwd = 1.5, col = gray(.1))
lines(yearsUS, Bounds99US[, "m.lower"], lwd = 2.9, col = gray(.5))
lines(yearsUS, Bounds99US[, "m.upper"], lwd = 2.9, col = gray(.5))
lines(yearsUS, Bounds99US[, "m.med"], lwd = 2.9, col = gray(.5))

lines(yearsES, Bounds99ES[, "f.lower"], lwd = 1.5, col = gray(.1), lty=4)
lines(yearsES, Bounds99ES[, "f.upper"], lwd = 1.5, col = gray(.1), lty=4)
lines(yearsES, Bounds99ES[, "f.med"], lwd = 1.5, col = gray(.1), lty=4)
lines(yearsES, Bounds99ES[, "m.lower"], lwd = 2.9, col = gray(.5), lty=4)
lines(yearsES, Bounds99ES[, "m.upper"], lwd = 2.9, col = gray(.5), lty=4)
lines(yearsES, Bounds99ES[, "m.med"], lwd = 2.9, col = gray(.5), lty=4)


text(2009.2,c(16.5,30,48),c("0.5%","50%","99.5%"), pos=4,xpd=TRUE)
segments(c(2010,2010),c(48,48),c(2006,2006),c(45,52))

legend(1968,69, lty = c(1,1,4,4), col = gray(c(.1,.5,.1,.5)), lwd = c(1.5,2.9,1.5,2.9),bty = "n",
        legend = c("US females", "US males", "ES females","ES males"), xpd = TRUE)
dev.off()

# ----------------------------------------------
# self-sufficient dissimilarity code:
yearsES <- 1975:2009
yearsUS <- 1969:2009
BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))
names(BxymfES) <- yearsES
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))



TotASFRUS <- unlist(lapply(as.character(yearsUS), function(yr, .Bxy, .Ex, .a = .5:110.5){
                    Bxy <- .Bxy[[yr]][[1]] + .Bxy[[yr]][[1]]
                    Fxm <- Mna0(Minf0(rowSums(Bxy) / 
                                            with(.Ex, Male[Year == as.integer(yr)])))
                    Fxf <- Mna0(Minf0(colSums(Bxy) / 
                                            with(.Ex, Female[Year == as.integer(yr)])))
                    
 
                    keep <- .a <= 85        
                    
                    .anew     <- seq(min(.a[keep]),max(.a[keep]),by = .01)
                    Fxminterp <- approx(x=.a[keep], y=Fxm[keep], xout =  .anew)
                    Fxfinterp <- approx(x=.a[keep], Fxf[keep], xout =  .anew)
                    
                    fFxm    <- Fxminterp$y / sum(Fxminterp$y)
                    fFxf    <- Fxfinterp$y / sum(Fxfinterp$y)
                    
                    theta <- 1 - sum(pmin(fFxm,fFxf))
                    
                },  .Bxy = BxymfUS, .Ex = ExUS))



TotASFRES <- unlist(lapply(as.character(yearsES), function(yr, .Bxy, .Ex, .a = .5:110.5){
                    Fxm <- Mna0(Minf0(rowSums(.Bxy[[yr]][[1]] + .Bxy[[yr]][[1]] ) / 
                                            with(.Ex, Male[Year == as.integer(yr)])))
                    Fxf <- Mna0(Minf0(colSums(.Bxy[[yr]][[1]] + .Bxy[[yr]][[1]] ) / 
                                            with(.Ex, Female[Year == as.integer(yr)])))
                    
                    keep <- .a <= 85        
                    
                    .anew     <- seq(min(.a[keep]),max(.a[keep]),by = .01)
                    Fxminterp <- approx(x=.a[keep], y=Fxm[keep], xout =  .anew)
                    Fxfinterp <- approx(x=.a[keep], Fxf[keep], xout =  .anew)
                    
                    fFxm    <- Fxminterp$y / sum(Fxminterp$y)
                    fFxf    <- Fxfinterp$y / sum(Fxfinterp$y)
                    
                    1 - sum(pmin(fFxm,fFxf))
                    
                },  .Bxy = BxymfES, .Ex = ExES))

# same thing with monte-carlo confidence bands...
# takes about 5 minutes to run. need to be described
TotASFRUS95 <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxy, .Ex, .a = .5:110.5){
                    Bxy <- .Bxy[[yr]][[1]] + .Bxy[[yr]][[1]]
                    Exf <- with(.Ex, Female[Year == as.integer(yr)])
                    Exm <-  with(.Ex, Male[Year == as.integer(yr)])
                    keep <- .a <= 85 
                   
                    quantile(replicate(1000,{
                               Bxyi <- matrix(rpois(n=length(Bxy), lambda = Bxy), ncol = ncol (Bxy))
                               Fxmi <- Mna0(Minf0(rowSums(Bxyi) / Exm))
                               Fxfi <- Mna0(Minf0(colSums(Bxyi) / Exf))
                               .anew     <- seq(min(.a[keep]),max(.a[keep]),by = .01)
                               Fxminterp <- approx(x=.a[keep], y=Fxmi[keep], xout =  .anew)
                               Fxfinterp <- approx(x=.a[keep], Fxfi[keep], xout =  .anew)
                               1 - sum(pmin( Fxminterp$y / sum(Fxminterp$y),Fxfinterp$y / sum(Fxfinterp$y)))
                           }), probs = c(.025, .975))          
                },  .Bxy = BxymfUS, .Ex = ExUS))

TotASFRES95 <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxy, .Ex, .a = .5:110.5){
                    Bxy <- .Bxy[[yr]][[1]] + .Bxy[[yr]][[1]]
                    Exf <- with(.Ex, Female[Year == as.integer(yr)])
                    Exm <-  with(.Ex, Male[Year == as.integer(yr)])
                    keep <- .a <= 85 
                    
                    quantile(replicate(1000,{
                                        Bxyi <- matrix(rpois(n=length(Bxy), lambda = Bxy), ncol = ncol (Bxy))
                                        Fxmi <- Mna0(Minf0(rowSums(Bxyi) / Exm))
                                        Fxfi <- Mna0(Minf0(colSums(Bxyi) / Exf))
                                        .anew     <- seq(min(.a[keep]),max(.a[keep]),by = .01)
                                        Fxminterp <- approx(x=.a[keep], y=Fxmi[keep], xout =  .anew)
                                        Fxfinterp <- approx(x=.a[keep], Fxfi[keep], xout =  .anew)
                                        1 - sum(pmin( Fxminterp$y / sum(Fxminterp$y),Fxfinterp$y / sum(Fxfinterp$y)))
                                    }), probs = c(.025, .975))          
                },  .Bxy = BxymfES, .Ex = ExES))

pdf("/home/triffe/git/DISS/latex/Figures/ASFRdissimilarity.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
ylim <- c(.12,.22)
plot(NULL, type = 'n', ylim = ylim, xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,ylim[1],2010,ylim[2],col = gray(.95), border=NA),
                abline(h = seq(ylim[1], ylim[2], by = .02), col = "white"),
                abline(v = seq(1970, 2005, by = 5), col = "white"),
                text(1968, seq(ylim[1], ylim[2], by = .02), seq(ylim[1], ylim[2], by = .02), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2005, by = 5),ylim[1], seq(1970, 2005, by = 5), pos = 1, cex = .8, xpd = TRUE),
                text(1988, .113, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,.229, expression(theta), cex = 1, xpd = TRUE)))
polygon(c(yearsUS,rev(yearsUS)), c(TotASFRUS95[,1],rev(TotASFRUS95[,2])), 
        border = gray(.2), col ="#BBBBBB40", lty = 1, lwd = .5)
polygon(c(yearsES,rev(yearsES)), c(TotASFRES95[,1],rev(TotASFRES95[,2])), 
        border = gray(.2), col = "#BBBBBB40", lty = 1, lwd = .5)
lines(yearsUS, TotASFRUS, lwd = 1, col = gray(.2))
lines(yearsES, TotASFRES, lwd = 1, col = gray(.2), lty=4)
text(c(1971, 1977),
        c(0.1647749, 0.1980414),
        c(expression(paste(theta," USA")), expression(paste(theta," ES"))))
#legend(1968,.14, lty = c(1,5), col = gray(c(.2,.4)), lwd = c(2,3),bty = "n",
#        legend = c("US", "Spain"), xpd = TRUE)
dev.off()