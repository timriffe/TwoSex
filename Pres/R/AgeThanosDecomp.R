setwd("/home/triffe/git/DISS/")

USdecompR <- local(get(load("Data/rDecompResults/USdecompExR2.Rdata")))
ESdecompR <- local(get(load("Data/rDecompResults/ESdecompExR2.Rdata")))

cols <- RColorBrewer::brewer.pal(4,"Set1")[c(3,2,1,4)]

# determine axes compatible with output from both countries
Neg <- USdecompR
Neg[Neg > 0] <- 0
Pos <- USdecompR
Pos[Pos < 0] <- 0

pdf("Pres/FiguresStatic/DecomprThanosUS.pdf", height = 5, width = 5)
par(mai = c(.5,.5,.5,.2), xaxs = "i", yaxs = "i")
plot(NULL, type = "n", xlab = "", ylab = "", xlim = c(-1, 42), ylim = c(-.004,.007), 
        axes = FALSE)
rect(-1,-.004,42,.007,col = gray(.95),border = NA)
abline(h = seq(-.004,.007,by=.001), col = "white")
abline(v = seq(1,41,by=5), col = "white")
barplot(t(Neg), add = TRUE, space = 0, border = NA, col = cols,width = 1,axes = FALSE)
barplot(t(Pos), add = TRUE, space = 0, border = NA, col = cols,width = 1,axes = FALSE)
lines(.5:40.5, rowSums(USdecompR), col = "green",lwd=2)
text(seq(1,41,by=5),-.004,seq(1970,2010,by=5),pos=1,cex=.8,xpd=TRUE)
text(-1,seq(-.004,.007,by=.001),seq(-.004,.007,by=.001),cex=.8,pos=2,xpd=TRUE)
text(c(11, 11, 11, 12,32),c(0.000522905,  0.002898587, -0.001163063,  0.005350904,-.0027), 
        c("Fertility shape", "SRB", "Mortality", "eTFR","eTFR"), cex = 1.3,pos = 4, col = c("white","white","white",gray(.2),gray(.2)))
segments(12.604584, 0.005159316, 9.953926, 0.004431285)
segments(32.97280, -0.002465856, 31.50796, -0.001833618)
text(20,-.005,"Year",xpd=TRUE)
text(-2,.0075,"Contribution\nto difference in r", pos = 4, xpd = TRUE)
text(28,.0007,"Gap",col = "green",cex=1.3)
dev.off()

# Spain
Neg <- ESdecompR
Neg[Neg > 0] <- 0
Pos <- ESdecompR
Pos[Pos < 0] <- 0

pdf("Pres/FiguresStatic/DecomprThanosES.pdf", height = 5, width = 5)
par(mai = c(.5,.5,.5,.2), xaxs = "i", yaxs = "i")
plot(NULL, type = "n", xlab = "", ylab = "", xlim = c(-7, 36), ylim = c(-.004,.007), 
        axes = FALSE)
rect(-7,-.004,36,.007,col = gray(.95),border = NA)
abline(h = seq(-.004,.007,by=.001), col = "white")
abline(v = seq(-5,35,by=5), col = "white")
barplot(t(Neg), add = TRUE, space = 0, border = NA, col = cols, width = 1,axes = FALSE)
barplot(t(Pos), add = TRUE, space = 0, border = NA, col = cols, width = 1,axes = FALSE)
lines(.5:34.5, rowSums(ESdecompR), col = "green",lwd=2)
text(seq(-5,35,by=5),-.004,seq(1970,2010,by=5),pos=1,cex=.8,xpd=TRUE)
text(-7,seq(-.004,.007,by=.001),seq(-.004,.007,by=.001),cex=.8,pos=2,xpd=TRUE)
text(c(12,12,12,  6, 29),c(0.0009635557,  0.0040289517, -0.0015079198,  0.0058107131, -0.0033), 
        c("Fertility shape","SRB","Mortality","eTFR","eTFR"), cex = 1.3,pos = 4, 
        col = c("white", "white", "white", gray(.2), gray(.2)))
segments(9.185489,0.005580808 ,10.441064,0.005044364)
segments(32.20436,-0.003002300,32.90191,-0.002389221)
text(15,-.005,"Year",xpd=TRUE)
text(-8,.0075,"Contribution\nto difference in r", pos = 4, xpd = TRUE)
text(18.5,.0025,"Gap",col = "green",cex=1.3)
dev.off()







