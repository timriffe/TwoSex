# make a pyramid rearrange itself!

makeRect <- function(x, y, w, h, ...){
    w2 <- w / 2
    h2 <- h / 2
    rect(x - w2, y - h2, x + w2, y + h2, ...)
}
rotationmat <- compiler::cmpfun(function(theta){
    matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),ncol=2)
})
rotatevec <- compiler::cmpfun(function(xy,Rmat){
    Rmat %*% c(xy)
})
# gold here. work on this one
makearc <- compiler::cmpfun(function(x1,y1,x2,y2,bc=.05,nvert=100){
            
    if (x1==x2 & y1 == y2){
        return(cbind(x=rep(x1,nvert),y = rep(y1,nvert)))
    }
    rise    <- y2 - y1 #y2 - y1
    run     <- x2 - x1 #x2-x1
    a       <- .5 * (rise^2 + run^2)^.5 # 1/2 hypotenuse
    b       <- a * bc # how bowed?
    add     <- ifelse(x1 < x2, 0, pi) # direction
    theta   <- add + atan(rise/run)
    Rmat    <- rotationmat(theta)
    pivec   <- seq(from = 0, to = pi, length.out = nvert)
    x       <- a * cos(pivec) + a
    y       <- b * sin(pivec)
    t(apply(cbind(x, y), 1, rotatevec, Rmat = Rmat) + c(x1, y1))   
})

makearcv <- compiler::cmpfun(function(x1,y1,x2,y2,bc=.05,nvert=100,female = TRUE){
            
            if (x1 == x2 & y1 == y2){
                return(c(rep(x1,nvert),rep(y1,nvert)))
            }
            rise    <- y2 - y1 #y2 - y1
            run     <- x2 - x1 #x2-x1
            a       <- .5 * (rise^2 + run^2)^.5 # 1/2 hypotenuse
            b       <- a * bc # how bowed?
            #add     <- ifelse(female, ifelse(x1<x2,-2*pi,-pi), ifelse(x1<x2,0,pi))
            add     <- ifelse(x1 < x2, 0, pi) # direction
            #add     <- pi * female
            theta   <- add + atan(rise/run) #* ifelse(female,-1,1)
            theta   <- ifelse(female,{-2*pi+theta},theta)
            Rmat    <- rotationmat(theta)
            pivec   <- seq(from = 0, to = pi, length.out = nvert)
            x       <- a * cos(pivec) + a
            y       <- b * sin(pivec)
            c(t(apply(cbind(x, y), 1, rotatevec, Rmat = Rmat) + c(x1, y1)))   
        })
#plot(NULL,xlim=c(-5,5),ylim=c(-5,5),asp=1)
#arc1 <- makearc(0,0,2,2,.3,100)
#lines(arc1[,1],arc1[,2],col="red")
#points(x=c(0,2),y=c(0,2))
#
#points(c(x1,x2),c(y1,y2),pch=19,cex=.5)
#for (i in 1:nrow(arc1)){
#   # plot(NULL,xlim=c(-10,10),ylim=c(-10,10),asp=1)
#    draw.rect(arc1[i,1],arc1[i,2],.5,.5)
#    Sys.sleep(.1)
#}


# step 2 get ll and ur coords for rectangles
setwd("/home/triffe/git/DISS/")
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

source("R/UtilityFunctions.R")

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
# let's look at some fertility rates:
{
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

}
# ----------------------------------------------------------------------------- #
# Figure for ey Reproduction, like the diagram in dissertation, grays
{
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
}
# ----------------------------------------------------------------------------- #
# figure of r for both sexes and standard Lotka r
{
rLotkaES    <- local(get(load("Data/results/agerSRB/rLotkaES.Rdata")))
rLotkaUS    <- local(get(load("Data/results/agerSRB/rLotkaUS.Rdata")))
rfES        <- local(get(load("/home/triffe/git/DISS/Data/results/exSingleSex/rfES.Rdata")))
rmES        <- local(get(load("/home/triffe/git/DISS/Data/results/exSingleSex/rmES.Rdata")))
rfUS        <- local(get(load("/home/triffe/git/DISS/Data/results/exSingleSex/rfUS.Rdata")))
rmUS        <- local(get(load("/home/triffe/git/DISS/Data/results/exSingleSex/rmUS.Rdata")))
yearsUS <- 1969:2009
yearsES <- 1975:2009
Cols <- RColorBrewer::brewer.pal(9,"Set1")
# for age US
pdf("Pres/FiguresStatic/rSingleSex1.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rmUS[, 1], type = 'n', ylim = c(-.02,.011),xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,-.02,2010,.011,col = gray(.95), border=NA),
                abline(h = seq(-.02,.011,by = .0025), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.02,.011,by = .005),seq(-.02,.011,by = .005), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.02, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.022, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,.0125, "r", cex = 1, xpd = TRUE)))

lines(yearsES, rmES[, 1],col = grey(.9), lwd = 2.5)
lines(yearsES, rfES[, 1],col = grey(.9), lwd = 3)
lines(yearsES, rLotkaES[,1], col = grey(.9), lwd = 1)
lines(yearsES, rLotkaES[,2], col = grey(.9), lwd = 1)
lines(yearsUS, rmUS[, 1],col = grey(.9), lwd = 2.5)
lines(yearsUS, rfUS[, 1],col = grey(.9), lwd = 3)
lines(yearsUS, rLotkaUS[,1], col = Cols[2], lwd = 1)
lines(yearsUS, rLotkaUS[,2], col = Cols[8], lwd = 1)
dev.off()

# for age and ey US
pdf("Pres/FiguresStatic/rSingleSex2.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rmUS[, 1], type = 'n', ylim = c(-.02,.011),xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,-.02,2010,.011,col = gray(.95), border=NA),
                abline(h = seq(-.02,.011,by = .0025), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.02,.011,by = .005),seq(-.02,.011,by = .005), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.02, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.022, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,.0125, "r", cex = 1, xpd = TRUE)))

lines(yearsES, rmES[, 1],col = grey(.9), lwd = 2.5)
lines(yearsES, rfES[, 1],col = grey(.9), lwd = 3)
lines(yearsES, rLotkaES[,1], col = grey(.9), lwd = 1)
lines(yearsES, rLotkaES[,2], col = grey(.9), lwd = 1)
lines(yearsUS, rmUS[, 1],col = Cols[2], lwd = 2.5)
lines(yearsUS, rfUS[, 1],col = Cols[8], lwd = 3)
lines(yearsUS, rLotkaUS[,1], col = Cols[2], lwd = 1)
lines(yearsUS, rLotkaUS[,2], col = Cols[8], lwd = 1)
dev.off()

pdf("Pres/FiguresStatic/rSingleSex3.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rmUS[, 1], type = 'n', ylim = c(-.02,.011),xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,-.02,2010,.011,col = gray(.95), border=NA),
                abline(h = seq(-.02,.011,by = .0025), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.02,.011,by = .005),seq(-.02,.011,by = .005), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.02, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.022, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,.0125, "r", cex = 1, xpd = TRUE)))

lines(yearsES, rmES[, 1],col = grey(.9), lwd = 2.5)
lines(yearsES, rfES[, 1],col = grey(.9), lwd = 3)
lines(yearsUS, rmUS[, 1],col = grey(.9), lwd = 2.5)
lines(yearsUS, rfUS[, 1],col = grey(.9), lwd = 3)
lines(yearsUS, rLotkaUS[,1], col = grey(.9), lwd = 1)
lines(yearsUS, rLotkaUS[,2], col = grey(.9), lwd = 1)
lines(yearsES, rLotkaES[,1], col = Cols[2], lwd = 1)
lines(yearsES, rLotkaES[,2], col =  Cols[8], lwd = 1)
dev.off()

pdf("Pres/FiguresStatic/rSingleSex4.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, rmUS[, 1], type = 'n', ylim = c(-.02,.011),xlim = c(1968,2010), axes = FALSE,
        xlab = "", ylab = "",
        panel.first = list(rect(1968,-.02,2010,.011,col = gray(.95), border=NA),
                abline(h = seq(-.02,.011,by = .0025), col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(-.02,.011,by = .005),seq(-.02,.011,by = .005), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),-.02, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.022, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1966,.0125, "r", cex = 1, xpd = TRUE)))

lines(yearsUS, rmUS[, 1],col = grey(.9), lwd = 2.5)
lines(yearsUS, rfUS[, 1],col = grey(.9), lwd = 3)
lines(yearsUS, rLotkaUS[,1], col = grey(.9), lwd = 1)
lines(yearsUS, rLotkaUS[,2], col = grey(.9), lwd = 1)
lines(yearsES, rLotkaES[,1], col = Cols[2], lwd = 1)
lines(yearsES, rLotkaES[,2], col =  Cols[8], lwd = 1)
lines(yearsES, rmES[, 1],col = Cols[2], lwd = 2.5)
lines(yearsES, rfES[, 1],col = Cols[8], lwd = 3)
dev.off()
}
# ----------------------------------------------------------------------------- ## ----------------------------------------------------------------------------- #





# awesome
#i<- 1
#for (i in 1:200){
#    plot(NULL, type = "n", xlim = c(-.01,.01),ylim = c(0,111))
#    makeRect(xm[,i],ym[,i],wm,1,col = colsm2, border = NA)
#    makeRect(xf[,i],yf[,i],wf,1,col = colsf2, border = NA)
#}
# ----------------------------------
library(animation)
xlabs <- c("1.0%","0.8%","0.6%","0.4%","0.2%","0%","0.2%","0.4%","0.6%","0.8%","1.0%")
#  /tmp/RtmpNwGXh2/age2eyPyramid.gif
saveGIF({
                par(mai=c(.5,.5,.6,.5),xaxs = "i", yaxs = "i")
            # frame 1 age pyramid, flat, starting image
                plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
                        axes = FALSE, xlab = "",ylab="", 
                        panel.first = list(
                           rect(-.011,0,.011,111,col = gray(.95), border = NA),
                           abline(h = seq(0,110,by=10), col = "white"),   
                           abline(v = seq(-.01, .01, by = .002), col = "white"),
                           text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                           text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                           text(0,-6,"Percent",xpd=TRUE,cex=1.3),
                           text(-.0115,114,"Age",xpd=TRUE,cex=1.3),
                           text(0,118,"An age-structured population",xpd=TRUE,cex=1.3)))
                makeRect(xm[,200],ym[,200],wm,1,col = gray(.5), border = NA, xpd = TRUE)
                makeRect(xf[,200],yf[,200],wf,1,col = gray(.5), border = NA, xpd = TRUE)
                segments(0,0,0,111,col="white")
            # frame 2  age pyramid, with ey colored
                plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
                    axes = FALSE, xlab = "",ylab="", 
                    panel.first = list(
                            rect(-.011,0,.011,111,col = gray(.95), border = NA),
                            abline(h = seq(0,110,by=10), col = "white"),   
                            abline(v = seq(-.01, .01, by = .002), col = "white"),
                            text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                            text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                            text(0,-6,"Percent",xpd=TRUE,cex=1.3),
                            text(-.0115,114,"Age",xpd=TRUE,cex=1.3),
                            text(0,118,"An age-structured population",xpd=TRUE,cex=1.3),
                            text(0,113,"is heterogeneous with respect to remaining years of life",xpd=TRUE,cex=1.3)))
                makeRect(xm[,200],ym[,200],wm,1,col = colsm, border = NA, xpd = TRUE)
                makeRect(xf[,200],yf[,200],wf,1,col = colsf, border = NA, xpd = TRUE)
                segments(0,0,0,111,col="white")
                text(-.011,80,"few remaining years",pos = 4,xpd=TRUE,cex=1.3)
                text(-.011,40,"many remaining \n                  years",pos = 4,xpd=TRUE,cex=1.3)
                segments(-0.0045,78, -0.00065579,69.61301)
                segments(-0.006239436,35.200460, -0.007,5.041144)
                PyramidOutline2(with(PxUS, Male[Year == 1975]),with(PxUS, Female[Year == 1975]),scale=1,border=gray(.1))
                segments(0,0,0,111,col="white")
            # frames 3-202 age-> ey
            for (i in 200:1){
                plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
                        axes = FALSE, xlab = "",ylab="", 
                        panel.first = list(
                                rect(-.011,0,.011,111,col = gray(.95), border = NA),
                                abline(h = seq(0,110,by=10), col = "white"),   
                                abline(v = seq(-.01, .01, by = .002), col = "white"),
                                text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                                text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                                text(0,118,"We can re-stack this pyramid",xpd=TRUE,cex=1.3),
                                text(0,113,"by remaining years of life",xpd=TRUE,cex=1.3),
                                text(0,-6,"Percent",xpd=TRUE,cex=1.3)
                                ))
                makeRect(xm[,i],ym[,i],wm,1,col = colsm, border = NA, xpd = TRUE)
                makeRect(xf[,i],yf[,i],wf,1,col = colsf, border = NA, xpd = TRUE)
            }
            # frame 203 ey pyramid with ey color
                plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
                        axes = FALSE, xlab = "",ylab="", 
                        panel.first = list(
                                rect(-.011,0,.011,111,col = gray(.95), border = NA),
                                abline(h = seq(0,110,by=10), col = "white"),   
                                abline(v = seq(-.01, .01, by = .002), col = "white"),
                                text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                                text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                                text(0,118,"and now we have a population",xpd=TRUE,cex=1.3),
                                text(0,113,"structured by remaining years of life",xpd=TRUE,cex=1.3),
                                text(0,-6,"Percent",xpd=TRUE,cex=1.3),
                                text(-.012,114,expression(e[y]),xpd=TRUE,cex=1.3)
                                ))
                makeRect(xm[,1],ym[,1],wm,1,col = colsm, border = NA, xpd = TRUE)
                makeRect(xf[,1],yf[,1],wf,1,col = colsf, border = NA, xpd = TRUE)
                PyramidOutline2(-rowSums(Males,na.rm=TRUE),rowSums(Females,na.rm=TRUE),scale=1,border=gray(.1))
                segments(0,0,0,111,col="white")
            # frame 204 ey flat
                plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
                        axes = FALSE, xlab = "",ylab="", 
                        panel.first = list(
                                rect(-.011,0,.011,111,col = gray(.95), border = NA),
                                abline(h = seq(0,110,by=10), col = "white"),   
                                abline(v = seq(-.01, .01, by = .002), col = "white"),
                                text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                                text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                                text(0,118,"and now we have a population",xpd=TRUE,cex=1.3),
                                text(0,113,"structured by remaining years of life",xpd=TRUE,cex=1.3),
                                text(0,-6,"Percent",xpd=TRUE,cex=1.3),
                                text(-.012,114,expression(e[y]),xpd=TRUE,cex=1.3)
                                 ))
                makeRect(xm[,1],ym[,1],wm,1,col = gray(.5), border = NA, xpd = TRUE)
                makeRect(xf[,1],yf[,1],wf,1,col = gray(.5), border = NA, xpd = TRUE)
                segments(0,0,0,111,col="white")
            # frame 205 ey pyramid with age color
                plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
                        axes = FALSE, xlab = "",ylab="", 
                        panel.first = list(
                                rect(-.011,0,.011,111,col = gray(.95), border = NA),
                                abline(h = seq(0,110,by=10), col = "white"),   
                                abline(v = seq(-.01, .01, by = .002), col = "white"),
                                text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                                text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                                text(0,-6,"Percent",xpd=TRUE,cex=1.3),
                                text(-.012,114,expression(e[y]),xpd=TRUE,cex=1.3),
                                text(0,118,"A remaining-years structured population",xpd=TRUE,cex=1.3),
                                text(0,113,"is heterogeneous with respect to age",xpd=TRUE,cex=1.3)))
                makeRect(xm[, 1], ym[, 1], wm, 1, col = colsm2, border = NA, xpd = TRUE)
                makeRect(xf[, 1], yf[, 1], wf, 1, col = colsf2, border = NA, xpd = TRUE)
                PyramidOutline2(-rowSums(Males,na.rm=TRUE),rowSums(Females,na.rm=TRUE),scale=1,border=gray(.1))
                text(-.008, 80, "Young", pos = 4,cex=1.3)
                text(-.008, 10, "Old", pos = 4,cex=1.3)
                segments(-0.0054899542,79.85945, -0.0005433743,69.61301)
                segments(-0.006351858,9.874367, -0.003953517,5.041144)
                segments(0,0,0,111,col="white")
            # frames 206-406 ey to age
            for (i in 1:200){
                plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
                        axes = FALSE, xlab = "",ylab="", 
                        panel.first = list(
                                rect(-.011,0,.011,111,col = gray(.95), border = NA),
                                abline(h = seq(0,110,by=10), col = "white"),   
                                abline(v = seq(-.01, .01, by = .002), col = "white"),
                                text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                                text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                                text(0,118,"We can un-stack this pyramid",xpd=TRUE,cex=1.3),
                                text(0,113,"back to age",xpd=TRUE,cex=1.3),
                                text(0,-6,"Percent",xpd=TRUE,cex=1.3)
                        ))
                    makeRect(xm[,i],ym[,i],wm,1,col = colsm2, border = NA, xpd = TRUE)
                    makeRect(xf[,i],yf[,i],wf,1,col = colsf2, border = NA, xpd = TRUE)
            }
            # frame 407 age with age color
             plot(NULL, type = "n", xlim = c(-.011,.011),ylim = c(0,111), 
                        axes = FALSE, xlab = "",ylab="", 
                        panel.first = list(
                            rect(-.011,0,.011,111,col = gray(.95), border = NA),
                            abline(h = seq(0,110,by=10), col = "white"),   
                            abline(v = seq(-.01, .01, by = .002), col = "white"),
                            text(-.011, seq(0,110,by=10),  seq(0,110,by=10), pos = 2, cex = .8, xpd = TRUE),
                            text(seq(-.01, .01, by = .002), 0, xlabs, pos = 1, cex = .8, xpd = TRUE),
                            text(0,118,"and now we have a population",xpd=TRUE,cex=1.3),
                            text(0,113,"structured by age",xpd=TRUE,cex=1.3),
                            text(0,-6,"Percent",xpd=TRUE,cex=1.3),
                            text(-.0115,114,"Age",xpd=TRUE,cex=1.3)
                            ))
            makeRect(xm[,200],ym[,200],wm,1,col = colsm2, border = NA, xpd = TRUE)
            makeRect(xf[,200],yf[,200],wf,1,col = colsf2, border = NA, xpd = TRUE)
            PyramidOutline2(with(PxUS, Male[Year == 1975]),with(PxUS, Female[Year == 1975]),scale=1,border=gray(.1))
            segments(0,0,0,111,col="white")
        }, movie.name = "age2eyPyramid.gif", interval = c(4,4,rep(.02,200),4,4,4,rep(.02,200),4), ani.width = 600,
        ani.height = 600, clean = FALSE)
     
list.files("/home/triffe/git/DISS/Pres/Age2eyAnimation")
        
        
        
        