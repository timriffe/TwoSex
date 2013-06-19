
setwd("/home/triffe/git/DISS/")

# get some functions
source("R/UtilityFunctions.R")
source("latex/Pres/R/UtilitiesPresentation.R")

# load data
dxmUS <- local(get(load("Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("Data/HMD_dx/dxfES.Rdata"))) 

dxmUS <- t(t(dxmUS) / colSums(dxmUS))
dxfUS <- t(t(dxfUS) / colSums(dxfUS))
dxmES <- t(t(dxmES) / colSums(dxmES))
dxfES <- t(t(dxfES) / colSums(dxfES))

PxUS <- local(get(load("Data/HMD_Px/PxUS.Rdata"))) 
PxES <- local(get(load("Data/HMD_Px/PxES.Rdata"))) 

# -----------------------------------------
# shape data for pyramid animations
# errr 1975 US
Males   <- ExpectedDx(with(PxUS, Male[Year == 1975]),dxmUS[,"1975"])
Females <- ExpectedDx(with(PxUS, Female[Year == 1975]),dxfUS[,"1975"])

# proportions:
Tot     <- sum(Males) + sum(Females)
Males   <- -Males / Tot
Females <- Females / Tot

# no need to draw rectangles for 0s
Males[Males == 0] <- NA
Females[Females == 0] <- NA
MindNA <- !is.na(Males)
FindNA <- !is.na(Females)

# first column = age 0
#take cumsum down the columns
MalesCa     <- apply(Males,2,cumsum)
FemalesCa   <- apply(Females,2,cumsum)
MalesCy     <- t(apply(Males,1,cumsum))
FemalesCy   <- t(apply(Females,1,cumsum))

# centroid x vals will be this minus 1/2 or orig values
MalesCnta   <- MalesCa - Males / 2
FemalesCnta <- FemalesCa - Females / 2
MalesCnty   <- MalesCy - Males / 2
FemalesCnty <- FemalesCy - Females / 2

heights <- 1
# the prop values are now the widths

x1m <- MalesCnta[MindNA]
y1m <- col(Males)[MindNA] - .5 # midpoints

x2m <- MalesCnty[MindNA]
y2m <- row(Males)[MindNA] - .5

x1f <- FemalesCnta[FindNA]
y1f <- col(Females)[FindNA] - .5 # midpoints

x2f <- FemalesCnty[FindNA]
y2f <- row(Females)[FindNA] - .5
#RColorBrewer::display.brewer.all()

colR <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,"YlGnBu"),space="Lab")

colsm <- c(rep(colR(11),each=10),rev(colR(11))[1])[row(Males)] 
dim(colsm) <- dim(Males)
colsm <- colsm[MindNA]
#colsf <- gray(seq(.85,.05,length.out=111))[row(Females)] 
colsf <- c(rep(colR(11),each=10),rev(colR(11))[1])[row(Females)] 
dim(colsf) <- dim(Females)
colsf <- colsf[FindNA]

# colR2 is for coloring age within an ey population
# actually it's easier than it seems:
colR2 <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd"),space="Lab")
# just need new colors to go in the other direction
colsm2 <-c(rep(colR2(11),each=10),rev(colR2(11))[1])[col(Males)] 
dim(colsm2) <- dim(Males)
colsm2 <- colsm2[MindNA]
#colsf <- gray(seq(.85,.05,length.out=111))[row(Females)] 
colsf2 <- c(rep(colR2(11),each=10),rev(colR2(11))[1])[col(Females)] 
dim(colsf2) <- dim(Females)
colsf2 <- colsf2[FindNA]
# widths
wm <- Males[MindNA]
wf <- Females[FindNA]

HM <- t(mapply(makearcv, x1m, y1m, x2m, y2m, MoreArgs = list(bc = .00015, nvert = 200, female = FALSE)))
HF <- t(mapply(makearcv, x1f, y1f, x2f, y2f, MoreArgs = list(bc = .00015, nvert = 200, female = TRUE)))
HM100 <- t(mapply(makearcv, x1m, y1m, x2m, y2m, MoreArgs = list(bc = .00015, nvert = 100, female = FALSE)))
HF100 <- t(mapply(makearcv, x1f, y1f, x2f, y2f, MoreArgs = list(bc = .00015, nvert = 100, female = TRUE)))

# x vals are in first 100 cols, y vals in cols 101-200
# xy vals for each rect arc
xm <- HM[,1:(ncol(HM)/2)]
ym <- HM[,(ncol(HM)/2+1):ncol(HM)]
xf <- HF[,1:(ncol(HF)/2)]
yf <- HF[,(ncol(HF)/2+1):ncol(HF)]

xm100 <- HM100[,1:(ncol(HM100)/2)]
ym100 <- HM100[,(ncol(HM100)/2+1):ncol(HM100)]
xf100 <- HF100[,1:(ncol(HF100)/2)]
yf100 <- HF100[,(ncol(HF100)/2+1):ncol(HF100)]
# save for presentation animation

# margins different than in gif animation!
xlabs <- c("1.0%","0.8%","0.6%","0.4%","0.2%","0%","0.2%","0.4%","0.6%","0.8%","1.0%")
#i<-0

# ----------------------------------------------------------------------------- #
# first plot, a gray age-sex-structured population
{
    
#png("/home/triffe/git/DISS/Pres/FiguresStatic/AgeSexGray.png",height=600,width=600) 
    pdf("/home/triffe/git/DISS/Pres/FiguresStatic/AgeSexGray.pdf")
    par(mai=c(.5,.5,.5,.3),xaxs = "i", yaxs = "i",xaxs="i")
    plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(-.011,0,.011,111,col = gray(.95), border = NA),
                    abline(h = seq(0,110,by=10), col = "white"),   
                    abline(v = seq(-.01, .01, by = .002), col = "white"),
                    text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                    text(0,-6,"Percent",xpd=TRUE,cex=1.3),
                    text(-.0115,114,"Age",xpd=TRUE,cex=1.3)))
    PyramidOutline2(with(PxUS, Male[Year == 1975]),with(PxUS, Female[Year == 1975]),scale=1,border=gray(.1),col = gray(.5))
#makeRect(xm100[,100],ym100[,100],wm,1,col = gray(.5), border = NA, xpd = TRUE)
#makeRect(xf100[,100],yf100[,100],wf,1,col = gray(.5), border = NA, xpd = TRUE)
    segments(0,0,0,111,col="white")
    dev.off()
}
# ----------------------------------------------------------------------------- #
# second plot, age-sex with ey heterogeneity
{
    
#png("/home/triffe/git/DISS/Pres/FiguresStatic/AgeSexEyHetero.png",height=600,width=600) 
    pdf("/home/triffe/git/DISS/Pres/FiguresStatic/AgeSexEyHetero.pdf")
    par(mai=c(.5,.5,.5,.3),xaxs = "i", yaxs = "i",xaxs="i")
    plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(-.011,0,.011,111,col = gray(.95), border = NA),
                    abline(h = seq(0,110,by=10), col = "white"),   
                    abline(v = seq(-.01, .01, by = .002), col = "white"),
                    text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                    text(0,-6,"Percent",xpd=TRUE,cex=1.3),
                    text(-.0115,114,"Age",xpd=TRUE,cex=1.3)))
    makeRect(xm100[,100],ym100[,100],wm,1,col = colsm, border = NA, xpd = TRUE)
    makeRect(xf100[,100],yf100[,100],wf,1,col = colsf, border = NA, xpd = TRUE)
    segments(0,0,0,111,col="white")
    text(-.011,80,"few remaining years",pos = 4,xpd=TRUE,cex=1.3)
    text(-.0114,40,"many remaining \n               years",pos = 4,xpd=TRUE,cex=1.3)
    segments(-0.0045,78, -0.00065579,69.61301)
    segments(-0.006239436,35.200460, -0.007,5.041144)
    PyramidOutline2(with(PxUS, Male[Year == 1975]),with(PxUS, Female[Year == 1975]),scale=1,border=gray(.1))
    segments(0,0,0,111,col="white")
    dev.off()
}
# ----------------------------------------------------------------------------- #
# Age -> ey animation figures (png)
{
    par(mai=c(0,0,0,0),xaxs = "i", yaxs = "i",xaxs="i")
    png("/home/triffe/git/DISS/Pres/Age2eyAnimation/frame000.png",height=600,width=600) 
    plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(-.011,0,.011,111,col = gray(.95), border = NA),
                    abline(h = seq(0,110,by=10), col = "white"),   
                    abline(v = seq(-.01, .01, by = .002), col = "white"),
                    text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                    text(0,-6,"Percent",xpd=TRUE,cex=1.3),
                    text(-.0115,114,"Age",xpd=TRUE,cex=1.3)))
    makeRect(xm100[,100],ym100[,100],wm,1,col = colsm, border = NA, xpd = TRUE)
    makeRect(xf100[,100],yf100[,100],wf,1,col = colsf, border = NA, xpd = TRUE)
    segments(0,0,0,111,col="white")
    text(-.011,80,"few remaining years",pos = 4,xpd=TRUE,cex=1.3)
    text(-.011,40,"many remaining \n                  years",pos = 4,xpd=TRUE,cex=1.3)
    segments(-0.0045,78, -0.00065579,69.61301)
    segments(-0.006239436,35.200460, -0.007,5.041144)
    PyramidOutline2(with(PxUS, Male[Year == 1975]),with(PxUS, Female[Year == 1975]),scale=1,border=gray(.1))
    segments(0,0,0,111,col="white")
    dev.off()
# first frame same as last (frames 001 - 100)
    for (i in (ncol(yf100)-1):0){
        out.path <- file.path("/home/triffe/git/DISS/Pres/Age2eyAnimation",
                paste0("frame",sprintf("%03d",ncol(yf100)-i),".png"))
        par(mai=c(0,0,0,0),xaxs = "i", yaxs = "i",xaxs="i")
        png(out.path,height=600,width=600) 
        plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
                axes = FALSE, xlab = "",ylab="", 
                panel.first = list(
                        rect(-.011,0,.011,111,col = gray(.95), border = NA),
                        abline(h = seq(0,110,by=10), col = "white"),   
                        abline(v = seq(-.01, .01, by = .002), col = "white"),
                        text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = 1.1, xpd = TRUE),
                        text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = 1.1, xpd = TRUE),
                        #text(0,118,"We can re-stack this pyramid",xpd=TRUE,cex=1.3),
                        #text(0,113,"by remaining years of life",xpd=TRUE,cex=1.3),
                        text(0,-8,"Percent",xpd=TRUE,cex=1.4)
                ))
        makeRect(xm100[,i+1],ym100[,i+1],wm,1,col = colsm, border = NA)
        makeRect(xf100[,i+1],yf100[,i+1],wf,1,col = colsf, border = NA)
        dev.off()
    }
# last frame touched-up (frame 101)
    par(mai=c(0,0,0,0),xaxs = "i", yaxs = "i",xaxs="i")
    png("/home/triffe/git/DISS/Pres/Age2eyAnimation/frame101.png",height=600,width=600) 
    plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(-.011,0,.011,111,col = gray(.95), border = NA),
                    abline(h = seq(0,110,by=10), col = "white"),   
                    abline(v = seq(-.01, .01, by = .002), col = "white"),
                    text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                    text(0,-6,"Percent",xpd=TRUE,cex=1.3),
                    text(-.012,114,expression(e[y]),xpd=TRUE,cex=1.3)))
    makeRect(xm100[,1],ym100[,1],wm,1,col = colsm, border = NA, xpd = TRUE)
    makeRect(xf100[,1],yf100[,1],wf,1,col = colsf, border = NA, xpd = TRUE)
    PyramidOutline2(-rowSums(Males,na.rm=TRUE),rowSums(Females,na.rm=TRUE),scale=1,border=gray(.1))
    segments(0,0,0,111,col="white")
    dev.off()
}
# ----------------------------------------------------------------------------- #
# ey stripes. repeat of last frame from previous
{
    par(mai=c(0,0,0,0),xaxs = "i", yaxs = "i",xaxs="i")
    png("/home/triffe/git/DISS/Pres/FiguresStatic/eyPyramid.png",height=600,width=600) 
    plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(-.011,0,.011,111,col = gray(.95), border = NA),
                    abline(h = seq(0,110,by=10), col = "white"),   
                    abline(v = seq(-.01, .01, by = .002), col = "white"),
                    text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                    text(0,-6,"Percent",xpd=TRUE,cex=1.3),
                    text(-.012,114,expression(e[y]),xpd=TRUE,cex=1.3)))
    makeRect(xm100[,1],ym100[,1],wm,1,col = colsm, border = NA, xpd = TRUE)
    makeRect(xf100[,1],yf100[,1],wf,1,col = colsf, border = NA, xpd = TRUE)
    PyramidOutline2(-rowSums(Males,na.rm=TRUE),rowSums(Females,na.rm=TRUE),scale=1,border=gray(.1))
    segments(0,0,0,111,col="white")
    dev.off()
}
# ----------------------------------------------------------------------------- #
# a gray ey-sex-structured population
{
#png("/home/triffe/git/DISS/Pres/FiguresStatic/eyPyramidGray.png",height=600,width=600) 
    pdf("/home/triffe/git/DISS/Pres/FiguresStatic/eyPyramidGray.pdf") 
    par(mai=c(.5,.5,.5,.3),xaxs = "i", yaxs = "i",xaxs="i")
    plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(-.011,0,.011,111,col = gray(.95), border = NA),
                    abline(h = seq(0,110,by=10), col = "white"),   
                    abline(v = seq(-.01, .01, by = .002), col = "white"),
                    text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                    text(0,-6,"Percent",xpd=TRUE,cex=1.3),
                    text(-.012,114,expression(e[y]),xpd=TRUE,cex=1.3)))
    PyramidOutline2(-rowSums(Males,na.rm=TRUE),rowSums(Females,na.rm=TRUE),scale=1,border=gray(.1),col = gray(.5))
#makeRect(xm100[,1],ym100[,1],wm,1,col = gray(.5), border = NA, xpd = TRUE)
#makeRect(xf100[,1],yf100[,1],wf,1,col = gray(.5), border = NA, xpd = TRUE)
    segments(0,0,0,111,col="white")
    dev.off()
}
# ----------------------------------------------------------------------------- #
# ey-sex-structured population with age-heterogeneity
{
    
    #png("/home/triffe/git/DISS/Pres/FiguresStatic/eyPyramidAgeHet.png",height=600,width=600) 
    pdf("/home/triffe/git/DISS/Pres/FiguresStatic/eyPyramidAgeHet.pdf")
    par(mai=c(.5,.5,.5,.3),xaxs = "i", yaxs = "i",xaxs="i")
    plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(-.011,0,.011,111,col = gray(.95), border = NA),
                    abline(h = seq(0,110,by=10), col = "white"),   
                    abline(v = seq(-.01, .01, by = .002), col = "white"),
                    text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                    text(0,-6,"Percent",xpd=TRUE,cex=1.3),
                    text(-.012,114,expression(e[y]),xpd=TRUE,cex=1.3)))
    makeRect(xm100[, 1], ym100[, 1], wm, 1, col = colsm2, border = NA, xpd = TRUE)
    makeRect(xf100[, 1], yf100[, 1], wf, 1, col = colsf2, border = NA, xpd = TRUE)
    PyramidOutline2(-rowSums(Males,na.rm=TRUE),rowSums(Females,na.rm=TRUE),scale=1,border=gray(.1))
    text(-.008, 80, "Young", pos = 4,cex=1.3)
    text(-.008, 10, "Old", pos = 4,cex=1.3)
    segments(-0.0054899542,79.85945, -0.0005433743,69.61301)
    segments(-0.006351858,9.874367, -0.003953517,5.041144)
    segments(0,0,0,111,col="white")
    dev.off()
}
# ----------------------------------------------------------------------------- #
# ey --> age animation
{
# first frame = previous plot
    par(mai=c(0,0,0,0),xaxs = "i", yaxs = "i",xaxs="i")
    png("/home/triffe/git/DISS/Pres/ey2ageAnimation/frame000.png",height=600,width=600) 
    plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(-.011,0,.011,111,col = gray(.95), border = NA),
                    abline(h = seq(0,110,by=10), col = "white"),   
                    abline(v = seq(-.01, .01, by = .002), col = "white"),
                    text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                    text(0,-6,"Percent",xpd=TRUE,cex=1.3),
                    text(-.012,114,expression(e[y]),xpd=TRUE,cex=1.3)))
    makeRect(xm100[, 1], ym100[, 1], wm, 1, col = colsm2, border = NA, xpd = TRUE)
    makeRect(xf100[, 1], yf100[, 1], wf, 1, col = colsf2, border = NA, xpd = TRUE)
    PyramidOutline2(-rowSums(Males,na.rm=TRUE),rowSums(Females,na.rm=TRUE),scale=1,border=gray(.1))
    text(-.008, 80, "Young", pos = 4,cex=1.3)
    text(-.008, 10, "Old", pos = 4,cex=1.3)
    segments(-0.0054899542,79.85945, -0.0005433743,69.61301)
    segments(-0.006351858,9.874367, -0.003953517,5.041144)
    segments(0,0,0,111,col="white")
    dev.off()
# now do animation ey -> age
    for (i in 1:100){
        out.path <- file.path("/home/triffe/git/DISS/Pres/ey2ageAnimation",
                paste0("frame",sprintf("%03d",i),".png"))
        par(mai=c(0,0,0,0),xaxs = "i", yaxs = "i",xaxs="i")
        png(out.path,height=600,width=600) 
        plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
                axes = FALSE, xlab = "",ylab="", 
                panel.first = list(
                        rect(-.011,0,.011,111,col = gray(.95), border = NA),
                        abline(h = seq(0,110,by=10), col = "white"),   
                        abline(v = seq(-.01, .01, by = .002), col = "white"),
                        text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                        text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                        text(0,-6,"Percent",xpd=TRUE,cex=1.3)
                ))
        makeRect(xm100[,i],ym100[,i],wm,1,col = colsm2, border = NA, xpd = TRUE)
        makeRect(xf100[,i],yf100[,i],wf,1,col = colsf2, border = NA, xpd = TRUE)
        dev.off()
    }
# last frame repeated
    par(mai=c(0,0,0,0),xaxs = "i", yaxs = "i",xaxs="i")
    png( "/home/triffe/git/DISS/Pres/ey2ageAnimation/frame101.png",height=600,width=600) 
    plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(-.011,0,.011,111,col = gray(.95), border = NA),
                    abline(h = seq(0,110,by=10), col = "white"),   
                    abline(v = seq(-.01, .01, by = .002), col = "white"),
                    text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                    text(0,-6,"Percent",xpd=TRUE,cex=1.3)
            ))
    makeRect(xm100[,100],ym100[,100],wm,1,col = colsm2, border = NA, xpd = TRUE)
    makeRect(xf100[,100],yf100[,100],wf,1,col = colsf2, border = NA, xpd = TRUE)
    segments(0,0,0,111,col="white")
    dev.off()
}
# ----------------------------------------------------------------------------- #
# animation Bx -> By
{
    BxUS  <- local(get(load("Data/USbirths/USBxy0_110.Rdata")))
    Bxf   <- colSums(BxUS[["1975"]])
    Byf   <- ExpectedDx(Bxf, dxfUS[,"1975"])
    BCa     <- apply(Byf,2,cumsum)
    BCy     <- t(apply(Byf,1,cumsum))
    BCnta   <- BCa - Byf / 2
    BCnty   <- BCy - Byf / 2
    BindNA  <- Byf > 0
    
    y1b <- BCnta[BindNA]
    x1b <- col(Byf)[BindNA] - .5 # midpoints
    y2b <- BCnty[BindNA]
    x2b <- row(Byf)[BindNA] - .5
# get colors
    colB <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,"RdPu"),space="Lab")
    colsB <- c(rep(colB(11),each=10),rev(colB(11))[1])[row(BCnta)] 
    dim(colsB) <- dim(BCnta)
    colsB <- colsB[BindNA]
    wB <- Byf[BindNA]
    
    HB100 <- t(mapply(makearcv, x1b, y1b, x2b, y2b, MoreArgs = list(bc = .00035, 
                            nvert = 100, female = TRUE)))
    xb100 <- HB100[,1:(ncol(HB100)/2)]
    yb100 <- HB100[,(ncol(HB100)/2+1):ncol(HB100)]
    
    par(mai=c(0,0,0,0),xaxs = "i", yaxs = "i",xaxs="i")
    png("/home/triffe/git/DISS/Pres/Ba2ByAnimation/frame000.png",height=600,width=600) 
    plot(NULL, type = "n", xlim = c(0,111),ylim = c(0,250000), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(0,0,111,250000,col = gray(.95), border = NA),
                    abline(v = seq(0,110,by=10), col = "white"),   
                    abline(h = seq(0, 250000, by = 25000), col = "white"),
                    text(0, seq(0,250000,by=50000),  c("0","50000","100000","150000","200000","250000"), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(0, 110, by = 10), 0, seq(0, 110, by = 10), pos = 1, cex = .8, xpd = TRUE),
                    text(55,-10000,"Age",xpd=TRUE,cex=1.3,pos=1),
                    text(-8,265000,"Births",xpd=TRUE,cex=1.3,pos=4)
            ))
    rect(x1b-.5, y1b - wB/2, x1b+.5, y1b+wB/2, col = colsB, border = NA)
    segments(0,0,0,111,col="white")
    dev.off()
    for (i in 99:0){
        out.path <- file.path("/home/triffe/git/DISS/Pres/Ba2ByAnimation",
                paste0("frame", sprintf("%03d", ncol(yb100) - i), ".png"))
        par(mai=c(0,0,0,0),xaxs = "i", yaxs = "i",xaxs="i")
        png(out.path, height = 600, width = 600) 
        plot(NULL, type = "n", xlim = c(0,111),ylim = c(0,250000), 
                axes = FALSE, xlab = "",ylab="", 
                panel.first = list(
                        rect(0,0,111,250000,col = gray(.95), border = NA),
                        abline(v = seq(0,110,by=10), col = "white"),   
                        abline(h = seq(0, 250000, by = 25000), col = "white"),
                        text(0, seq(0,250000,by=50000),  c("0","50000","100000","150000","200000","250000"), pos = 2, cex = .8, xpd = TRUE),
                        text(seq(0, 110, by = 10), 0, seq(0, 110, by = 10), pos = 1, cex = .8, xpd = TRUE),
                        text(55,-10000,"Age",xpd=TRUE,cex=1.3,pos=1),
                        text(-8,265000,"Births",xpd=TRUE,cex=1.3,pos=4)
                ))
        rect(xb100[, i + 1] - .5, yb100[, i + 1] - wB / 2, xb100[, i + 1] + .5, yb100[, i + 1] + wB / 2, 
                col = colsB, border = NA)
        dev.off()
    }
    par(mai=c(0,0,0,0),xaxs = "i", yaxs = "i",xaxs="i")
    png("/home/triffe/git/DISS/Pres/Ba2ByAnimation/frame101.png",height=600,width=600)  
    plot(NULL, type = "n", xlim = c(0,111),ylim = c(0,250000), 
            axes = FALSE, xlab = "",ylab="", 
            panel.first = list(
                    rect(0,0,111,250000,col = gray(.95), border = NA),
                    abline(v = seq(0,110,by=10), col = "white"),   
                    abline(h = seq(0, 250000, by = 25000), col = "white"),
                    text(0, seq(0,250000,by=50000),  c("0","50000","100000","150000","200000","250000"), pos = 2, cex = .8, xpd = TRUE),
                    text(seq(0, 110, by = 10), 0, seq(0, 110, by = 10), pos = 1, cex = .8, xpd = TRUE),
                    text(55,-10000,expression(e[y]),xpd=TRUE,cex=1.3,pos=1),
                    text(-8,265000,"Births",xpd=TRUE,cex=1.3,pos=4)
            )) 
    rect(x2b-.5, y2b - wB/2, x2b+.5, y2b+wB/2, col = colsB, border = NA)
    segments(0,0,0,111,col="white")
    dev.off()
}
# ----------------------------------------------------------------------------- #









