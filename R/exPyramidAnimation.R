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
plot(NULL,xlim=c(-5,5),ylim=c(-5,5),asp=1)
arc1 <- makearc(0,0,2,2,.3,100)
lines(arc1[,1],arc1[,2],col="red")
points(x=c(0,2),y=c(0,2))

points(c(x1,x2),c(y1,y2),pch=19,cex=.5)
for (i in 1:nrow(arc1)){
   # plot(NULL,xlim=c(-10,10),ylim=c(-10,10),asp=1)
    draw.rect(arc1[i,1],arc1[i,2],.5,.5)
    Sys.sleep(.1)
}


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
PxES <- local(get(load("Data/HMD_Px/PxUS.Rdata"))) 

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
#colR <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,"Set3"),space="Lab")
# get colors
#colsm <- gray(seq(.85,.05,length.out=111))[row(Males)] 
#colsm <- colR(111)[row(Males)] 
#dim(colsm) <- dim(Males)
#colsm <- colsm[MindNA]
##colsf <- gray(seq(.85,.05,length.out=111))[row(Females)] 
#colsf <- colR(111)[row(Females)] 
#dim(colsf) <- dim(Females)
#colsf <- colsf[FindNA]

colsm <- c(rep(colR(11),each=10),rev(colR(11))[1])[row(Males)] 
dim(colsm) <- dim(Males)
colsm <- colsm[MindNA]
#colsf <- gray(seq(.85,.05,length.out=111))[row(Females)] 
colsf <- c(rep(colR(11),each=10),rev(colR(11))[1])[row(Females)] 
dim(colsf) <- dim(Females)
colsf <- colsf[FindNA]


# widths
wm <- Males[MindNA]
wf <- Females[FindNA]

HM <- t(mapply(makearcv, x1m, y1m, x2m, y2m, MoreArgs = list(bc = .00015, nvert = 200, female = FALSE)))
HF <- t(mapply(makearcv, x1f, y1f, x2f, y2f, MoreArgs = list(bc = .00015, nvert = 200, female = TRUE)))

# x vals are in first 100 cols, y vals in cols 101-200
# xy vals for each rect arc
xm <- HM[,1:(ncol(HM)/2)]
ym <- HM[,(ncol(HM)/2+1):ncol(HM)]
xf <- HF[,1:(ncol(HF)/2)]
yf <- HF[,(ncol(HF)/2+1):ncol(HF)]


for (i in 200:1){
plot(NULL, type = "n", xlim = c(-.01,.01),ylim = c(0,111))
makeRect(xm[,i],ym[,i],wm,1,col = colsm, border = NA)
makeRect(xf[,i],yf[,i],wf,1,col = colsf, border = NA)
}


#library(animation)
#saveGIF({
#            for (i in 200:1){
#                plot(NULL, type = "n", xlim = c(-.01,.01),ylim = c(0,111))
#                makeRect(xm[,i],ym[,i],wm,1,col = colsm, border = NA)
#                makeRect(xf[,i],yf[,i],wf,1,col = colsf, border = NA)
#            }
#        }, movie.name = "age2eyPyramid.gif", interval = c(2,rep(.02,198),2), ani.width = 600,
#        ani.height = 600)
#getwd()

# Think of making one that goes back and forth. 
# 1) age pyramid, no color (pause 3 sec)
# 2) color by remaining years of life (pause 5 sec)
# 3) move to remaining-years pyramid (200 frames)
# 4) remove color (pause 3 sec)
# 5) recolor by age (pause 5 sec)
# 6) move back to age pyramid (200 frames)
# 

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


# awesome
i<- 1
for (i in 1:200){
    plot(NULL, type = "n", xlim = c(-.01,.01),ylim = c(0,111))
    makeRect(xm[,i],ym[,i],wm,1,col = colsm2, border = NA)
    makeRect(xf[,i],yf[,i],wf,1,col = colsf2, border = NA)
}
# ----------------------------------
library(animation)
seq(-.01, .01, by = .002)
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
                PyramidOutline(with(PxUS, Male[Year == 1975]),with(PxUS, Female[Year == 1975]),scale=1,border=gray(.1))
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
                PyramidOutline(-rowSums(Males,na.rm=TRUE),rowSums(Females,na.rm=TRUE),scale=1,border=gray(.1))
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
                PyramidOutline(-rowSums(Males,na.rm=TRUE),rowSums(Females,na.rm=TRUE),scale=1,border=gray(.1))
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
            PyramidOutline(with(PxUS, Male[Year == 1975]),with(PxUS, Female[Year == 1975]),scale=1,border=gray(.1))
            segments(0,0,0,111,col="white")
        }, movie.name = "age2eyPyramid.gif", interval = c(4,4,rep(.02,200),4,4,4,rep(.02,200),4), ani.width = 600,
        ani.height = 600)
     
        