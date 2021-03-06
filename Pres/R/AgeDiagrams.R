# ----------------------------------------------
# linear age diagram
setwd("/home/triffe/git/DISS/")

pdf("Pres/FiguresStatic/AgeChrono.pdf", height = 2, width = 5)
par(mai=c(0,0,0,0),xaxs="i",yaxs="i")
plot(NULL, xlim = c(-.05,1.05),ylim = c(0,1),axes = FALSE, xlab = "", ylab = "")
arrows(.2,.7,.8,.7, lwd = 4)
text(.5,.92,"time",cex=2)
points(.2,.7,pch=19,cex=2)
text(.08,.7,"Birth",cex=2)
points(.9,.7,cex=2)
segments(seq(.2,.7,by=.1),.7,seq(.2,.7,by=.1),.65)
text(seq(.2,.7,by=.1),.6,0:5)
dev.off()

pdf("Pres/FiguresStatic/AgeThano.pdf", height = 2, width = 5)
par(mai=c(0,0,0,0),xaxs="i",yaxs="i")
plot(NULL, xlim = c(-.05,1.05),ylim = c(0,1),axes = FALSE, xlab = "", ylab = "")
arrows(.2,.7,.8,.7, lwd = 4, col = gray(.8))
text(.5,.92,"time",cex=2)
points(.2,.7,pch=19,cex=2, col = gray(.8))
text(.08,.7,"Birth",cex=2, col = gray(.8))
points(.9,.7,cex=2, col = gray(.8))
segments(seq(.2,.7,by=.1),.7,seq(.2,.7,by=.1),.65, col = gray(.8))
text(seq(.2,.7,by=.1),.6,0:5, col = gray(.8))
arrows(.2,.3,.8,.3, lwd = 4)
points(.8,.3,pch=19,cex=2)
text(.93,.3,"Death",cex=2)
points(.1,.3,cex=2)
segments(seq(.3,.8,by=.1),.3,seq(.3,.8,by=.1),.25)
text(seq(.3,.8,by=.1),.2,5:0)
dev.off()

pdf("Pres/FiguresStatic/AgeFiller.pdf", height = 2, width = 5)
par(mai=c(0,0,0,0),xaxs="i",yaxs="i")
plot(NULL, xlim = c(-.05,1.05),ylim = c(0,1),axes = FALSE, xlab = "", ylab = "")
dev.off()