setwd("/home/triffe/git/DISS/")

# simply load decomposition results:
USdecompR <- local(get(load("Data/rDecompResults/USdecompmxR2.Rdata")))
ESdecompR <- local(get(load("Data/rDecompResults/ESdecompmxR2.Rdata")))

# determine axes compatible with output from both countries
Neg <- USdecompR
Neg[Neg > 0] <- 0
Pos <- USdecompR
Pos[Pos < 0] <- 0
2010-1968

cols <- RColorBrewer::brewer.pal(4,"Set1")[c(3,2,1,4)]
pdf("Pres/FiguresStatic/DecomprUS.pdf", height = 5, width = 5)
ymin <- -.003
ymax <- .007
par(mai = c(.5,.5,.5,.2), xaxs = "i", yaxs = "i")
plot(NULL, type = "n", xlab = "", ylab = "", xlim = c(-1, 42), ylim = c(ymin,ymax), 
        axes = FALSE)
rect(-1,ymin,42,ymax,col = gray(.95),border = NA)
abline(h = seq(ymin,ymax,by=.001), col = "white")
abline(v = seq(1,41,by=5), col = "white")
barplot(t(Neg), add = TRUE, space = 0, border = NA, col = cols,
        width = 1,axes = FALSE,axisnames=FALSE)
barplot(t(Pos), add = TRUE, space = 0, border = NA, col = cols,
        width = 1,axes = FALSE,axisnames=FALSE)
text(seq(1,41,by=5),-.003,seq(1970,2010,by=5),pos=1,cex=.8,xpd=TRUE)
text(-1,seq(ymin,ymax,by=.001),seq(ymin,ymax,by=.001),cex=.8,pos=2,xpd=TRUE)
lines(.5:40.5, rowSums(USdecompR), col = "green",lwd=2)
text(c(7.5, 5, 7.5, 7.5, 29.6246),
        c(-0.0004035491, -0.0019188301,  0.0010072297,  0.0034978640,-0.001187315), 
        c("Mortality", "Fertility shape", "SRB", "TFR", "TFR"), cex = 1.3,pos = 4, 
        col = c("white", gray(.2), "white", "white", "white"))
segments(5.559413,-0.0016575748, 2.001951,-0.0006996385)
text(20,-.004,"Year",xpd=TRUE)
text(-2,.0075,"Contribution\nto difference in r", pos = 4, xpd = TRUE)
text(22.5,.001,"Gap",col="green",cex=1.2)
dev.off()


# Spain
Neg <- ESdecompR
Neg[Neg > 0] <- 0
Pos <- ESdecompR
Pos[Pos < 0] <- 0

pdf("Pres/FiguresStatic/DecomprES.pdf", height = 5, width = 5)
ymin <- -.003
ymax <- .007
par(mai = c(.5,.5,.5,.2), xaxs = "i", yaxs = "i")
plot(NULL, type = "n", xlab = "", ylab = "", xlim = c(-7, 36), ylim = c(ymin,ymax), 
        axes = FALSE)
rect(-7,ymin,36,ymax,col = gray(.95),border = NA)
abline(h = seq(ymin,ymax,by=.001), col = "white")
abline(v = seq(-5,35,by=5), col = "white")
barplot(t(Neg), add = TRUE, space = 0, border = NA, col = cols,
        width = 1,axes = FALSE,axisnames=FALSE)
barplot(t(Pos), add = TRUE, space = 0, border = NA, col = cols,
        width = 1,axes = FALSE,axisnames=FALSE)
text(seq(-5,35,by=5),-.003,seq(1970,2010,by=5),pos=1,cex=.8,xpd=TRUE)
text(-7,seq(ymin,ymax,by=.001),seq(ymin,ymax,by=.001),cex=.8,pos=2,xpd=TRUE)
lines(.5:34.5, rowSums(ESdecompR), col = "green",lwd=2)
text(c(12,12,12,12,28),
        c(-0.0003338810,  0.0003976339,  0.0021393362,  0.004, -0.0011873151), 
        c("Mortality","Fertility shape","SRB", "TFR", "TFR"), cex = 1.3,pos = 4, col = "white")
text(15,-.004,"Year",xpd=TRUE)
text(-9,.0075,"Contribution\nto difference in r", pos = 4, xpd = TRUE)
text(30,.002,"Gap",col="green",cex=1.2)
dev.off()


