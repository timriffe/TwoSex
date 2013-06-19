setwd("/home/triffe/git/DISS/")

# get some functions
source("R/UtilityFunctions.R")

dxmUS <- local(get(load("Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("Data/HMD_dx/dxfES.Rdata"))) 

dxmUS <- t(t(dxmUS) / colSums(dxmUS))
dxfUS <- t(t(dxfUS) / colSums(dxfUS))
dxmES <- t(t(dxmES) / colSums(dxmES))
dxfES <- t(t(dxfES) / colSums(dxfES))
BxUS  <- local(get(load("Data/USbirths/USBxy0_110.Rdata")))
BxES  <- local(get(load("Data/ESbirths/ESBxy.Rdata")))
ExUS <- local(get(load("Data/Exposures/USexp.Rdata"))) 
ExES <- local(get(load("Data/Exposures/ESexp.Rdata"))) 

# 2 years, 1975, 2009
FyUS   <- do.call(cbind,lapply(c(1975,2009), function(yr, .Bx,.Ex,.dxm, .dxf){
                    yrc <- as.character(yr)
                    Exm <- with(.Ex, Male[Year == yr])
                    Exf <- with(.Ex, Female[Year == yr])
                    Bxy <- .Bx[[yrc]]
                    Bym <- rowSums(ExpectedDx(rowSums(Bxy),.dxm[,yrc]))
                    Byf <- rowSums(ExpectedDx(colSums(Bxy),.dxf[,yrc]))
                    Eym <- rowSums(ExpectedDx(Exm,.dxm[,yrc]))
                    Eyf <- rowSums(ExpectedDx(Exf,.dxf[,yrc]))
                    Fym <- Minf0(Mna0(Bym / Eym))
                    Fyf <- Minf0(Mna0(Byf / Eyf))
                    cbind(Fym, Fyf)
                }, .Bx = BxUS,.Ex = ExUS,.dxm = dxmUS, .dxf = dxfUS))
FyES   <- do.call(cbind,lapply(c(1975,2009), function(yr, .Bx,.Ex,.dxm, .dxf){
                    yrc <- as.character(yr)
                    Exm <- with(.Ex, Male[Year == yr])
                    Exf <- with(.Ex, Female[Year == yr])
                    Bxy <- .Bx[[yrc]]
                    Bym <- rowSums(ExpectedDx(rowSums(Bxy),.dxm[,yrc]))
                    Byf <- rowSums(ExpectedDx(colSums(Bxy),.dxf[,yrc]))
                    Eym <- rowSums(ExpectedDx(Exm,.dxm[,yrc]))
                    Eyf <- rowSums(ExpectedDx(Exf,.dxf[,yrc]))
                    Fym <- Minf0(Mna0(Bym / Eym))
                    Fyf <- Minf0(Mna0(Byf / Eyf))
                    cbind(Fym, Fyf)
                }, .Bx = BxES,.Ex = ExES,.dxm = dxmES, .dxf = dxfES))
FxUS   <- do.call(cbind,lapply(c(1975,2009), function(yr, .Bx,.Ex){
                    yrc <- as.character(yr)
                    Exm <- with(.Ex, Male[Year == yr])
                    Exf <- with(.Ex, Female[Year == yr])
                    Bxy <- .Bx[[yrc]]
                    Fxm <- Minf0(Mna0(rowSums(Bxy) / Exm))
                    Fxf <- Minf0(Mna0(colSums(Bxy) / Exf))
                    cbind(Fxm, Fxf)
                }, .Bx = BxUS,.Ex = ExUS))
FxES   <- do.call(cbind,lapply(c(1975,2009), function(yr, .Bx,.Ex,.dxm, .dxf){
                    yrc <- as.character(yr)
                    Exm <- with(.Ex, Male[Year == yr])
                    Exf <- with(.Ex, Female[Year == yr])
                    Bxy <- .Bx[[yrc]]
                    Fxm <- Minf0(Mna0(rowSums(Bxy) / Exm))
                    Fxf <- Minf0(Mna0(colSums(Bxy) / Exf))
                    cbind(Fxm, Fxf)
                }, .Bx = BxES,.Ex = ExES))



eys <- 0:110 # x values
LC <- RColorBrewer::brewer.pal(3,"Dark2")

# remaining years plots
{
#png("/home/triffe/git/DISS/Pres/FiguresStatic/Fym.png",height=600,width=600) 
    # US males
    pdf("/home/triffe/git/DISS/Pres/FiguresStatic/Fy1.pdf")
    par(mai=c(.5,.5,.5,.5),xaxs = "i", yaxs = "i",xaxs="i")
    plot(NULL, type = "n", xlim = c(0,111),ylim = c(0,.22), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(0,0,111,.22,col = gray(.95), border = NA),
                    abline(v = seq(0,110,by=10), col = "white"),   
                    abline(h = seq(0, .22, by = .02), col = "white"),
                    text(0, seq(0, .22, by = .02), seq(0, .22, by = .02), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(0,110,by=10), 0, seq(0,110,by=10), pos = 1, cex = .8, xpd = TRUE),
                    text(55,-.015,expression(e[y]),xpd=TRUE,cex=1.5),
                    text(-5,.23,expression(F[y]),xpd=TRUE,cex=1.5)
            ))
    lines(eys, FyES[,2], lwd = 3, col = gray(.9))
    lines(eys, FyES[,4], lwd = 3, col = gray(.9))
    points(eys, FyES[,4],pch=19,cex = .6, col = gray(.9))
    lines(eys, FyUS[,2],col = gray(.9), lwd = 3)
    lines(eys, FyUS[,4],col = gray(.9), lwd = 3)
    points(eys, FyUS[,4],col = gray(.9),pch=19,cex = .6)
    lines(eys, FyES[,1], lwd = 3, col = gray(.9))
    lines(eys, FyES[,3], lwd = 3, col = gray(.9))
    points(eys, FyES[,3],pch=19,cex = .6, col = gray(.9))
    lines(eys, FyUS[,1],col = LC[3], lwd = 3)
    lines(eys, FyUS[,3],col = LC[3], lwd = 3)
    points(eys, FyUS[,3],col = LC[3],pch=19,cex = .6)
    text(c(22, 72),c(0.04, 0.05),c("1975","2009"),cex=1.5)
    dev.off()
    
# US females
    pdf("/home/triffe/git/DISS/Pres/FiguresStatic/Fy2.pdf")
    par(mai=c(.5,.5,.5,.5),xaxs = "i", yaxs = "i",xaxs="i")
    plot(NULL, type = "n", xlim = c(0,111),ylim = c(0,.22), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(0,0,111,.22,col = gray(.95), border = NA),
                    abline(v = seq(0,110,by=10), col = "white"),   
                    abline(h = seq(0, .22, by = .02), col = "white"),
                    text(0, seq(0, .22, by = .02), seq(0, .22, by = .02), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(0,110,by=10), 0, seq(0,110,by=10), pos = 1, cex = .8, xpd = TRUE),
                    text(55,-.015,expression(e[y]),xpd=TRUE,cex=1.5),
                    text(-5,.23,expression(F[y]),xpd=TRUE,cex=1.5)
            ))
    lines(eys, FyES[,2], lwd = 3, col = gray(.9))
    lines(eys, FyES[,4], lwd = 3, col = gray(.9))
    points(eys, FyES[,4],pch=19,cex = .6, col = gray(.9))
    lines(eys, FyES[,2], lwd = 3, col = gray(.9))
    lines(eys, FyES[,4], lwd = 3, col = gray(.9))
    points(eys, FyES[,4],pch=19,cex = .6, col = gray(.9))
    lines(eys, FyUS[,1],col = gray(.9), lwd = 3)
    lines(eys, FyUS[,3],col = gray(.9), lwd = 3)
    points(eys, FyUS[,3],col = gray(.9),pch=19,cex = .6)
    lines(eys, FyUS[,2],col = LC[3], lwd = 3)
    lines(eys, FyUS[,4],col = LC[3], lwd = 3)
    points(eys, FyUS[,4],col = LC[3],pch=19,cex = .6)
    text(c(37, 77),c(0.04, 0.05),c("1975","2009"),cex=1.5)
    dev.off()
    
# ES males
    pdf("/home/triffe/git/DISS/Pres/FiguresStatic/Fy3.pdf")
    par(mai=c(.5,.5,.5,.5),xaxs = "i", yaxs = "i",xaxs="i")
    plot(NULL, type = "n", xlim = c(0,111),ylim = c(0,.22), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(0,0,111,.22,col = gray(.95), border = NA),
                    abline(v = seq(0,110,by=10), col = "white"),   
                    abline(h = seq(0, .22, by = .02), col = "white"),
                    text(0, seq(0, .22, by = .02), seq(0, .22, by = .02), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(0,110,by=10), 0, seq(0,110,by=10), pos = 1, cex = .8, xpd = TRUE),
                    text(55,-.015,expression(e[y]),xpd=TRUE,cex=1.5),
                    text(-5,.23,expression(F[y]),xpd=TRUE,cex=1.5)
            ))
    lines(eys, FyES[,2],col = gray(.9), lwd = 3)
    lines(eys, FyES[,4],col = gray(.9), lwd = 3)
    points(eys, FyES[,4],col = gray(.9),pch=19,cex = .6)
    lines(eys, FyUS[,2],col = gray(.9), lwd = 3)
    lines(eys, FyUS[,4],col = gray(.9), lwd = 3)
    points(eys, FyUS[,4],col = gray(.9),pch=19,cex = .6)
    lines(eys, FyUS[,1],col = gray(.9), lwd = 3)
    lines(eys, FyUS[,3],col = gray(.9), lwd = 3)
    points(eys, FyUS[,3],col = gray(.9),pch=19,cex = .6)
    lines(eys, FyES[,1], lwd = 3, col = LC[2])
    lines(eys, FyES[,3], lwd = 3, col = LC[2])
    points(eys, FyES[,3],pch=19,cex = .6, col = LC[2])
    text(c(27,48),c(.065,.03),c("1975","2009"),cex=1.5)
    dev.off()
    
# ES females
    pdf("/home/triffe/git/DISS/Pres/FiguresStatic/Fy4.pdf")
    par(mai=c(.5,.5,.5,.5),xaxs = "i", yaxs = "i",xaxs="i")
    plot(NULL, type = "n", xlim = c(0,111),ylim = c(0,.22), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(0,0,111,.22,col = gray(.95), border = NA),
                    abline(v = seq(0,110,by=10), col = "white"),   
                    abline(h = seq(0, .22, by = .02), col = "white"),
                    text(0, seq(0, .22, by = .02), seq(0, .22, by = .02), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(0,110,by=10), 0, seq(0,110,by=10), pos = 1, cex = .8, xpd = TRUE),
                    text(55,-.015,expression(e[y]),xpd=TRUE,cex=1.5),
                    text(-5,.23,expression(F[y]),xpd=TRUE,cex=1.5)
            ))
    lines(eys, FyUS[,2],col = gray(.9), lwd = 3)
    lines(eys, FyUS[,4],col = gray(.9), lwd = 3)
    points(eys, FyUS[,4],col = gray(.9),pch=19,cex = .6)
    lines(eys, FyES[,1], lwd = 3, col = gray(.9))
    lines(eys, FyES[,3], lwd = 3, col = gray(.9))
    points(eys, FyES[,3],pch=19,cex = .6, col = gray(.9))
    lines(eys, FyUS[,1],col = gray(.9), lwd = 3)
    lines(eys, FyUS[,3],col = gray(.9), lwd = 3)
    points(eys, FyUS[,3],col = gray(.9),pch=19,cex = .6)
    lines(eys, FyES[,2],col = LC[2], lwd = 3)
    lines(eys, FyES[,4],col = LC[2], lwd = 3)
    points(eys, FyES[,4],col = LC[2],pch=19,cex = .6)
    text(c(30, 55),c(0.045, 0.03),c("1975","2009"),cex=1.5)
    dev.off()
    
    
# ------------------------------------ #
    
# ------------------------------------ #
    
# same curves by age
    pdf("/home/triffe/git/DISS/Pres/FiguresStatic/Fx1.pdf")
    par(mai=c(.5,.5,.5,.5),xaxs = "i", yaxs = "i",xaxs="i")
    plot(NULL, type = "n", xlim = c(0,111),ylim = c(0,.22), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(0,0,111,.22,col = gray(.95), border = NA),
                    abline(v = seq(0,110,by=10), col = "white"),   
                    abline(h = seq(0, .22, by = .02), col = "white"),
                    text(0, seq(0, .22, by = .02), seq(0, .22, by = .02), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(0,110,by=10), 0, seq(0,110,by=10), pos = 1, cex = .8, xpd = TRUE),
                    text(55,-.015,"Age",xpd=TRUE,cex=1.5),
                    text(-4,.23,expression(F[x]),xpd=TRUE,cex=1.5)
            ))
    lines(eys, FxUS[,2], lwd = 3, col = gray(.9))
    lines(eys, FxUS[,4], lwd = 3, col = gray(.9))
    points(eys, FxUS[,4],pch=19,cex = .6, col = gray(.9))
    lines(eys, FxES[,2], lwd = 3, col = gray(.9))
    lines(eys, FxES[,4], lwd = 3, col = gray(.9))
    points(eys, FxES[,4],pch=19,cex = .6, col = gray(.9))
    lines(eys, FxES[,1], lwd = 3, col = gray(.9))
    lines(eys, FxES[,3], lwd = 3, col = gray(.9))
    points(eys, FxES[,3],pch=19,cex = .6, col = gray(.9))
    lines(eys, FxUS[,1],col = LC[3], lwd = 3)
    lines(eys, FxUS[,3],col = LC[3], lwd = 3)
    points(eys, FxUS[,3],col = LC[3],pch=19,cex = .6)
    text(c(13, 42),c(0.1, 0.1),c("1975","2009"),cex=1.5)
    dev.off()
    
    pdf("/home/triffe/git/DISS/Pres/FiguresStatic/Fx2.pdf")
    par(mai=c(.5,.5,.5,.5),xaxs = "i", yaxs = "i",xaxs="i")
    plot(NULL, type = "n", xlim = c(0,111),ylim = c(0,.22), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(0,0,111,.22,col = gray(.95), border = NA),
                    abline(v = seq(0,110,by=10), col = "white"),   
                    abline(h = seq(0, .22, by = .02), col = "white"),
                    text(0, seq(0, .22, by = .02), seq(0, .22, by = .02), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(0,110,by=10), 0, seq(0,110,by=10), pos = 1, cex = .8, xpd = TRUE),
                    text(55,-.015,"Age",xpd=TRUE,cex=1.5),
                    text(-4,.23,expression(F[x]),xpd=TRUE,cex=1.5)
            ))
    lines(eys, FxES[,2], lwd = 3, col = gray(.9))
    lines(eys, FxES[,4], lwd = 3, col = gray(.9))
    points(eys, FxES[,4],pch=19,cex = .6, col = gray(.9))
    lines(eys, FxES[,1], lwd = 3, col = gray(.9))
    lines(eys, FxES[,3], lwd = 3, col = gray(.9))
    points(eys, FxES[,3],pch=19,cex = .6, col = gray(.9))
    lines(eys, FxUS[,1],col = gray(.9), lwd = 3)
    lines(eys, FxUS[,3],col = gray(.9), lwd = 3)
    points(eys, FxUS[,3],col = gray(.9),pch=19,cex = .6)
    lines(eys, FxUS[,2], lwd = 3, col = LC[3])
    lines(eys, FxUS[,4], lwd = 3, col = LC[3])
    points(eys, FxUS[,4],pch=19,cex = .6, col = LC[3])
    text(c(11, 40),c(0.1, 0.1),c("1975","2009"),cex=1.5)
    dev.off()
    
    
    pdf("/home/triffe/git/DISS/Pres/FiguresStatic/Fx3.pdf")
    par(mai=c(.5,.5,.5,.5),xaxs = "i", yaxs = "i",xaxs="i")
    plot(NULL, type = "n", xlim = c(0,111),ylim = c(0,.22), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(0,0,111,.22,col = gray(.95), border = NA),
                    abline(v = seq(0,110,by=10), col = "white"),   
                    abline(h = seq(0, .22, by = .02), col = "white"),
                    text(0, seq(0, .22, by = .02), seq(0, .22, by = .02), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(0,110,by=10), 0, seq(0,110,by=10), pos = 1, cex = .8, xpd = TRUE),
                    text(55,-.015,"Age",xpd=TRUE,cex=1.5),
                    text(-4,.23,expression(F[x]),xpd=TRUE,cex=1.5)
            ))
    lines(eys, FxES[,2], lwd = 3, col = gray(.9))
    lines(eys, FxES[,4], lwd = 3, col = gray(.9))
    points(eys, FxES[,4],pch=19,cex = .6, col = gray(.9))
    lines(eys, FxUS[,1],col = gray(.9), lwd = 3)
    lines(eys, FxUS[,3],col = gray(.9), lwd = 3)
    points(eys, FxUS[,3],col = gray(.9),pch=19,cex = .6)
    lines(eys, FxUS[,2], lwd = 3, col =  gray(.9))
    lines(eys, FxUS[,4], lwd = 3, col =  gray(.9))
    points(eys, FxUS[,4],pch=19,cex = .6, col =  gray(.9))
    lines(eys, FxES[,1], lwd = 3, col =  LC[2])
    lines(eys, FxES[,3], lwd = 3, col = LC[2])
    points(eys, FxES[,3],pch=19,cex = .6, col = LC[2])
    text(c(15, 32),c(0.1, 0.02),c("1975","2009"),cex=1.5)
    dev.off()
    
    pdf("/home/triffe/git/DISS/Pres/FiguresStatic/Fx4.pdf")
    par(mai=c(.5,.5,.5,.5),xaxs = "i", yaxs = "i",xaxs="i")
    plot(NULL, type = "n", xlim = c(0,111),ylim = c(0,.22), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(0,0,111,.22,col = gray(.95), border = NA),
                    abline(v = seq(0,110,by=10), col = "white"),   
                    abline(h = seq(0, .22, by = .02), col = "white"),
                    text(0, seq(0, .22, by = .02), seq(0, .22, by = .02), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(0,110,by=10), 0, seq(0,110,by=10), pos = 1, cex = .8, xpd = TRUE),
                    text(55,-.015,"Age",xpd=TRUE,cex=1.5),
                    text(-4,.23,expression(F[x]),xpd=TRUE,cex=1.5)
            ))
    
    lines(eys, FxUS[,1],col = gray(.9), lwd = 3)
    lines(eys, FxUS[,3],col = gray(.9), lwd = 3)
    points(eys, FxUS[,3],col = gray(.9),pch=19,cex = .6)
    lines(eys, FxUS[,2], lwd = 3, col =  gray(.9))
    lines(eys, FxUS[,4], lwd = 3, col =  gray(.9))
    points(eys, FxUS[,4],pch=19,cex = .6, col = gray(.9))
    lines(eys, FxES[,1], lwd = 3, col = gray(.9))
    lines(eys, FxES[,3], lwd = 3, col = gray(.9))
    points(eys, FxES[,3],pch=19,cex = .6, col = gray(.9))
    lines(eys, FxES[,2], lwd = 3, col = LC[2])
    lines(eys, FxES[,4], lwd = 3, col = LC[2])
    points(eys, FxES[,4],pch=19,cex = .6, col = LC[2])
    text(c(13, 26),c(0.1, 0.02),c("1975","2009"),cex=1.5)
    dev.off()
}

