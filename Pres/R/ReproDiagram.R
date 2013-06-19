# Figure for ey Reproduction, like the diagram in dissertation, grays

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






colR3 <- grDevices::colorRampPalette(gray(c(.3,.8)),space="Lab")

# just need new colors to go in the other direction
colsm3 <-c(rep(colR3(11),each=10),rev(colR3(11))[1])[col(Males)] 
dim(colsm3) <- dim(Males)
colsm3 <- colsm3[MindNA]
#colsf <- gray(seq(.85,.05,length.out=111))[row(Females)] 
colsf3 <- c(rep(colR3(11),each=10),rev(colR3(11))[1])[col(Females)] 
dim(colsf3) <- dim(Females)
colsf3 <- colsf3[FindNA]

colR4 <- grDevices::colorRampPalette(rev(gray(c(.3,.8))),space="Lab")
# just need new colors to go in the other direction
colsm4 <-c(rep(colR4(11),each=10),rev(colR4(11))[1])[row(Males)] 
dim(colsm4) <- dim(Males)
colsm4 <- 
        colsm4 <- colsm4[MindNA]
#colsf <- gray(seq(.85,.05,length.out=111))[row(Females)] 
colsf4 <- c(rep(colR4(11),each=10),rev(colR4(11))[1])[row(Females)] 
dim(colsf4) <- dim(Females)
colsf4 <- colsf4[FindNA]
# flat gray, no indications
pdf("/home/triffe/git/DISS/Pres/FiguresStatic/eyRepro1.pdf") 
par(mai=c(.6,.5,.5,.3),xaxs = "i", yaxs = "i",xaxs="i")
plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
        axes = FALSE, xlab = "",ylab="", 
        panel.first = list(
                rect(-.011,0,.011,111,col = gray(.95), border = NA),
                abline(h = seq(0,110,by=10), col = "white"),   
                abline(v = seq(-.01, .01, by = .002), col = "white"),
                text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                text(0,-9,"Percent",xpd=TRUE,cex=1.5),
                text(-.012,116,expression(e[y]),xpd=TRUE,cex=1.5)))
points(0,1,col="#FFFF0050")
PyramidOutline2(-rowSums(Males,na.rm=TRUE),rowSums(Females,na.rm=TRUE),scale=1,border=gray(.1),col = gray(.6))
segments(0,0,0,111,col="white")
dev.off()

pdf("/home/triffe/git/DISS/Pres/FiguresStatic/exRepro1.pdf") 
par(mai=c(.6,.5,.5,.3),xaxs = "i", yaxs = "i",xaxs="i")
plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
        axes = FALSE, xlab = "",ylab="", 
        panel.first = list(
                rect(-.011,0,.011,111,col = gray(.95), border = NA),
                abline(h = seq(0,110,by=10), col = "white"),   
                abline(v = seq(-.01, .01, by = .002), col = "white"),
                text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                text(0,-9,"Percent",xpd=TRUE,cex=1.5),
                text(-.012,116,"Age",xpd=TRUE,cex=1.5)))
points(0,1,col="#FFFF0050")
PyramidOutline2(-colSums(Males,na.rm=TRUE),colSums(Females,na.rm=TRUE),scale=1,border=gray(.1),col = gray(.6))
segments(0,0,0,111,col="white")
dev.off()

##### flat gray. Down Arrow only
pdf("/home/triffe/git/DISS/Pres/FiguresStatic/eyRepro2.pdf")
par(mai=c(.6,.5,.5,.3),xaxs = "i", yaxs = "i",xaxs="i")
plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
        axes = FALSE, xlab = "",ylab="", 
        panel.first = list(
                rect(-.011,0,.011,111,col = gray(.95), border = NA),
                abline(h = seq(0,110,by=10), col = "white"),   
                abline(v = seq(-.01, .01, by = .002), col = "white"),
                text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                text(0,-9,"Percent",xpd=TRUE,cex=1.5),
                text(-.012,116,expression(e[y]),xpd=TRUE,cex=1.5)))
points(0,1,col="#FFFF0050")
PyramidOutline2(-rowSums(Males,na.rm=TRUE),rowSums(Females,na.rm=TRUE),scale=1,border=gray(.1),col = gray(.6))
segments(0,0,0,111,col="white")
arrows(-.008,80,-.008,40,lwd=2)
text(-.008,80,"Progression",cex=1.5,pos=3)
dev.off()
# Age (up arrow for prgression)
pdf("/home/triffe/git/DISS/Pres/FiguresStatic/exRepro2.pdf")
par(mai=c(.6,.5,.5,.3),xaxs = "i", yaxs = "i",xaxs="i")
plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
        axes = FALSE, xlab = "",ylab="", 
        panel.first = list(
                rect(-.011,0,.011,111,col = gray(.95), border = NA),
                abline(h = seq(0,110,by=10), col = "white"),   
                abline(v = seq(-.01, .01, by = .002), col = "white"),
                text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                text(0,-9,"Percent",xpd=TRUE,cex=1.5),
                text(-.012,116,"Age",xpd=TRUE,cex=1.5)))
points(0,1,col="#FFFF0050")
PyramidOutline2(-colSums(Males,na.rm=TRUE),colSums(Females,na.rm=TRUE),scale=1,border=gray(.1),col = gray(.6))
segments(0,0,0,111,col="white")
arrows(-.008,40,-.008,80,lwd=2)
text(-.008,40,"Progression",cex=1.5,pos=1)
dev.off()
# increment
grDevices::colorRampPalette(c("#00FFAAFF","#00FFAA00"),space="Lab")



pdf("/home/triffe/git/DISS/Pres/FiguresStatic/eyRepro3.pdf")
par(mai=c(.6,.5,.5,.3),xaxs = "i", yaxs = "i",xaxs="i")
plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
        axes = FALSE, xlab = "",ylab="", 
        panel.first = list(
                rect(-.011,0,.011,111,col = gray(.95), border = NA),
                abline(h = seq(0,110,by=10), col = "white"),   
                abline(v = seq(-.01, .01, by = .002), col = "white"),
                text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                text(0,-9,"Percent",xpd=TRUE,cex=1.5),
                text(-.012,116,expression(e[y]),xpd=TRUE,cex=1.5)))
points(0,1,col="#FFFF0050")

makeRect(xm100[, 1], ym100[, 1], wm, 1, col = colsm3, border = NA, xpd = TRUE)
makeRect(xf100[, 1], yf100[, 1], wf, 1, col = colsf3, border = NA, xpd = TRUE)
MalesInc   <- Males[,1:10]
FemalesInc <-Females[,1:10]
Mcol  <- Fcol  <- matrix(rgb(0, 1, .4, seq(1,0,length=10), names = NULL, maxColorValue = 1),
        ncol=ncol(MalesInc),nrow=nrow(MalesInc),byrow=TRUE)
MyMat <- FyMat <- row(MalesInc) - .5
Find       <- is.na(FemalesInc)
Mind       <- is.na(MalesInc)
MalesInc[Mind] <- 0
FemalesInc[Find] <- 0
MalesIncx    <- t(apply(MalesInc,1,cumsum)) - MalesInc / 2
FemalesIncx  <- t(apply(FemalesInc,1,cumsum)) - FemalesInc / 2
MalesInc     <- MalesInc[!Mind] ;  FemalesInc     <- FemalesInc[!Find]
MalesIncx    <- MalesIncx[!Mind] ; FemalesIncx    <- FemalesIncx[!Find]
Mcol         <- Mcol[!Mind] ;      Fcol           <- Fcol[!Find]
MyMat        <- MyMat[!Mind] ;     FyMat          <- FyMat[!Find]
makeRect(x=MalesIncx / 2,y=MyMat, 
        w=MalesInc,h=1,col = Mcol, border = NA)
makeRect(x=FemalesIncx / 2,y=FyMat, 
        w=FemalesInc,h=1,col = Fcol, border = NA)
text(.0003,81,"Increment",srt=270,col = "black",cex=1.5)
arrows(-.008,80,-.008,40,lwd=2)
text(-.008,80,"Progression",cex=1.5,pos=3)
PyramidOutline2(-rowSums(Males,na.rm=TRUE),rowSums(Females,na.rm=TRUE),scale=1,border=gray(.1))
dev.off()

pdf("/home/triffe/git/DISS/Pres/FiguresStatic/exRepro3.pdf")
par(mai=c(.6,.5,.5,.3),xaxs = "i", yaxs = "i",xaxs="i")
plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
        axes = FALSE, xlab = "",ylab="", 
        panel.first = list(
                rect(-.011,0,.011,111,col = gray(.95), border = NA),
                abline(h = seq(0,110,by=10), col = "white"),   
                abline(v = seq(-.01, .01, by = .002), col = "white"),
                text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                text(0,-9,"Percent",xpd=TRUE,cex=1.5),
                text(-.012,116,"Age",xpd=TRUE,cex=1.5)))
points(0,1,col="#FFFF0050")
PyramidOutline2(-colSums(Males,na.rm=TRUE),colSums(Females,na.rm=TRUE),scale=1,border=gray(.1),col = gray(.6))
segments(0,0,0,111,col="white")
arrows(-.008,40,-.008,80,lwd=2)
text(-.008,40,"Progression",cex=1.5,pos=1)
rect(colSums(Males,na.rm=TRUE)[1:10],0:9,colSums(Females,na.rm=TRUE)[1:10],1:10,
        col=rgb(0, 1, .4, seq(1,0,length=10), names = NULL, maxColorValue = 1),border = NA)
text(.0015,2.5,"Increment",cex = 1.5,pos=4)
PyramidOutline2(-colSums(Males,na.rm=TRUE),colSums(Females,na.rm=TRUE),scale=1,border=gray(.1),col = NA)
dev.off()
# decrement

pdf("/home/triffe/git/DISS/Pres/FiguresStatic/eyRepro4.pdf")
par(mai=c(.6,.5,.5,.3),xaxs = "i", yaxs = "i",xaxs="i")
plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
        axes = FALSE, xlab = "",ylab="", 
        panel.first = list(
                rect(-.011,0,.011,111,col = gray(.95), border = NA),
                abline(h = seq(0,110,by=10), col = "white"),   
                abline(v = seq(-.01, .01, by = .002), col = "white"),
                text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                text(0,-9,"Percent",xpd=TRUE,cex=1.5),
                text(-.012,116,expression(e[y]),xpd=TRUE,cex=1.5)))
makeRect(xm100[, 1], ym100[, 1], wm, 1, col = colsm3, border = NA, xpd = TRUE)
makeRect(xf100[, 1], yf100[, 1], wf, 1, col = colsf3, border = NA, xpd = TRUE)
makeRect(x=MalesIncx / 2,y=MyMat, 
        w=MalesInc,h=1,col = Mcol, border = NA)
makeRect(x=FemalesIncx / 2,y=FyMat, 
        w=FemalesInc,h=1,col = Fcol, border = NA)
rect(rowSums(Males,na.rm=TRUE)[1:10],0:9,rowSums(Females,na.rm=TRUE)[1:10],1:10,
        col=rgb(1, 1, 0, seq(1,0,length=10), names = NULL, maxColorValue = 1),border = NA)
text(.005,2.5,"Decrement",cex = 1.5,pos=4)
arrows(.005,2.5,0.004,1,length=.05)
text(.0003,81,"Increment",srt=270,col = "black",cex=1.5)
arrows(-.008,80,-.008,40,lwd=2)
text(-.008,80,"Advance",cex=1.5,pos=3)
PyramidOutline2(-rowSums(Males,na.rm=TRUE),rowSums(Females,na.rm=TRUE),scale=1,border=gray(.1))
dev.off()

pdf("/home/triffe/git/DISS/Pres/FiguresStatic/exRepro4.pdf")
par(mai=c(.6,.5,.5,.3),xaxs = "i", yaxs = "i",xaxs="i")
plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
        axes = FALSE, xlab = "",ylab="", 
        panel.first = list(
                rect(-.011,0,.011,111,col = gray(.95), border = NA),
                abline(h = seq(0,110,by=10), col = "white"),   
                abline(v = seq(-.01, .01, by = .002), col = "white"),
                text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                text(0,-9,"Percent",xpd=TRUE,cex=1.5),
                text(-.012,116,"Age",xpd=TRUE,cex=1.5)))
points(0,1,col="#FFFF0050")

makeRect(xm100[, 100], ym100[, 100], wm, 1, col = colsm4, border = NA, xpd = TRUE)
makeRect(xf100[, 100], yf100[, 100], wf, 1, col = colsf4, border = NA, xpd = TRUE)
MalesIncA   <- Males[1:10,]
FemalesIncA <- Females[1:10,]
McolA  <- FcolA  <- matrix(rgb(1, 1, 0, seq(1,0,length=10), names = NULL, maxColorValue = 1),
        ncol=ncol(MalesIncA),nrow=nrow(MalesIncA))
MyMatA <- FyMatA <- col(MalesIncA) - .5
FindA       <- is.na(FemalesIncA)
MindA       <- is.na(MalesIncA)
MalesIncA[MindA] <- 0
FemalesIncA[FindA] <- 0
MalesIncxA    <- apply(MalesIncA,2,cumsum) - MalesIncA / 2
FemalesIncxA  <- apply(FemalesIncA,2,cumsum) - FemalesIncA / 2
MalesIncA     <- MalesIncA[!MindA] ;  FemalesIncA     <- FemalesIncA[!FindA]
MalesIncxA    <- MalesIncxA[!MindA] ; FemalesIncxA    <- FemalesIncxA[!FindA]
makeRect(x=MalesIncxA / 2,y=MyMatA, 
        w=MalesIncA,h=1,col = McolA, border = NA)
makeRect(x=FemalesIncxA / 2,y=FyMatA, 
        w=FemalesIncA,h=1,col = FcolA, border = NA)

text(.0002,78,"Decrement",srt=270,col = "black",cex=1.5)
arrows(-.008,40,-.008,80,lwd=2)
text(-.008,40,"Progression",cex=1.5,pos=1)
rect(colSums(Males,na.rm=TRUE)[1:10],0:9,colSums(Females,na.rm=TRUE)[1:10],1:10,
        col=rgb(0, 1, .4, seq(1,0,length=10), names = NULL, maxColorValue = 1),border = NA)
text(.0015,2.5,"Increment",cex = 1.5,pos=4)
PyramidOutline2(-colSums(Males,na.rm=TRUE),colSums(Females,na.rm=TRUE),scale=1,border=gray(.1),col = NA)
dev.off()
