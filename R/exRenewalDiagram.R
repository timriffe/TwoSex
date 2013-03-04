# well-drawn diagram required for exPopulation Pyramid.

source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")


BxUS  <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))

# exposures, as such, straight from HMD, all ages 0-110, long form
ExUS  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
# get Lx estimates for R0, r
dxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata"))) 

dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)

ExmUS1975 <- rowSums(ExpectedDx( with(ExUS, Male[Year == 1975]), dxmUS[, "1975"]))
BxmUS1975 <- rowSums(ExpectedDx( rowSums(BxUS[["1975"]]), dxmUS[, "1975"]))
ExfUS1975 <- rowSums(ExpectedDx( with(ExUS, Female[Year == 1975]), dxfUS[, "1975"]))
BxfUS1975 <- rowSums(ExpectedDx( colSums(BxUS[["1975"]]), dxfUS[, "1975"]))

FxmUS1975 <- BxmUS1975 / ExmUS1975
FxfUS1975 <- BxfUS1975 / ExfUS1975


ExM     <- ExpectedDx( with(ExUS, Male[Year == 1975]), dxmUS[, "1975"])
ExF     <- ExpectedDx( with(ExUS, Female[Year == 1975]), dxfUS[, "1975"])
EM      <- rowSums(ExM)
EF      <- rowSums(ExF)
Total   <- sum(EM) + sum(EF)
Scale   <- 1000
ExM     <- Scale * ExM / Total
ExF     <- Scale * ExF / Total
EM      <- Scale * EM / Total
EF      <- Scale * EF / Total



pdf("/home/triffe/git/DISS/latex/Figures/exRenovationDiagram.pdf",height = 6.5, width = 8.5)
par(mai = c(.5,.4,.5,.2))
plot(NULL, type = "n", xlim = c(-10, 50), ylim = c(-50,120), axes = FALSE, xlab = "", ylab = "")

ages<- seq(5,110,by=5)
cols <- gray(seq(.2,.8,length = length(ages)))
for (i in  length(ages):1){
PyramidOutline(males = rowSums(ExM[,1:ages[i]]), 
        females = rowSums(ExF[,1:ages[i]]), 
        scale = sum(c(ExM[,1:ages[i]],ExF[,1:ages[i]])), 
        gap = 3, col = cols[i], border = NA, xpd =TRUE, 
        x.shift = -3, y.shift = 25)
}
PyramidOutline(males = rowSums(ExM[,1:5]), 
        females = rowSums(ExF[,1:5]), 
        scale = sum(c(ExM[,1:5],ExF[,1:5])), 
        gap = 3, col = "royalblue", xpd =TRUE, 
        x.shift = -3, y.shift = 25)
PyramidOutline(males = EM, 
        females = EF, 
        scale = sum(c(EM, EF)), 
        gap = 3, col = NA, xpd = TRUE,
        x.shift = -3, y.shift = 25)
segments(-4.5,seq(0,100,by = 10) + 25, -4, seq(0,100,by = 10) + 25)
segments(-1.5,seq(0,100,by = 10) + 25, -2, seq(0,100,by = 10) + 25)
text(-3, seq(0,100,by = 10) + 25,seq(0,100,by = 10),cex = .5)
# ex
text(-3,115 + 25, expression(e[x]), xpd =TRUE)

rect(-EM[1]-4.5,25,-4.5,26,col="yellow")
rect(EF[1]-1.5,25,-1.5,26,col="yellow")

# eSFR
ASP <- 4
polygon(c(c(0,0:110,110) / ASP) + 12, c(0,FxfUS1975 * 500,0) + 70, col = "lightblue")
segments(seq(0,110,by = 10)/ASP + 12, 70,seq(0,110,by = 10)/ASP + 12,68)
text(seq(0,110,by = 10)/ASP + 12, 70,seq(0,110,by = 10), pos = 1, cex = .5, xpd = TRUE)
text(27,80,expression(F["y,t"]))

# times
text(13.5,52,expression(NULL %*% NULL), cex = 2, lwd = 2)

# exposure
polygon(c(c(0,0:110,110) / ASP) + 12, c(0,EF*ASP,0) + 20, col = "lightblue")
segments(seq(0,110,by = 10)/ASP + 12, 20,seq(0,110,by = 10)/ASP + 12,18)
text(seq(0,110,by = 10)/ASP + 12, 20,seq(0,110,by = 10), pos = 1, cex = .5, xpd = TRUE)
text(27,30,expression(E["y,t"]))
# equals
text(13.5,4.2,expression(NULL == NULL), cex = 2, lwd = 2)

# births
polygon(c(c(0,0:110,110) / ASP) + 12, c(0,EF*FxfUS1975,0) * 30 - 10, col = "lightblue")
segments(seq(0,110,by = 10)/ASP + 12, -10,seq(0,110,by = 10)/ASP + 12,-12)
text(seq(0,110,by = 10)/ASP + 12, -10,seq(0,110,by = 10), pos = 1, cex = .5, xpd = TRUE)
text(27,-3,expression(B["y,t"]))

# sum Byt
text(45,1.5,expression(sum(B["y,t"],y==0,omega) == B[t]), cex = 1.5)

# pointer
lines(c(-1, 2,  8.5,  8.5), 
        c(107, 124, 124, -25))
rect(5,-18,41,-68, col = gray(.95), xpd = TRUE)
# Bt times
text(8,-48, expression(B[t]%*%NULL), cex = 1.5)

polygon(c(c(0,0:110,110) / ASP) + 12, c(0,dxfUS[, "1975"] * 30,0) * 30 - 55, col = "lightblue")
segments(seq(0,110,by = 10)/ASP + 12, -55,seq(0,110,by = 10)/ASP + 12,-57)
text(27,-40,expression(d["a,t"]))
text(seq(0,110,by = 10)/ASP + 12, -55,seq(0,110,by = 10), pos = 1, cex = .5, xpd = TRUE)

# deaths
#text(-12,1,"Year t deaths", pos = 4, xpd =TRUE)
segments(-6.4,6.6,-6.4,25)
segments(-6.4,6.6,.4,25)

# increment
text(27,120,"Increment", cex = 2)
# decrement
text(-13,2,"decrement", cex = 2,pos= 4)

text(16,-35,"(redistribute)", xpd = TRUE)

dev.off()











