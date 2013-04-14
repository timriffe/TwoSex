source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")
yearsUS <- 1969:2009
yearsES <- 1975:2009

dxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfES.Rdata"))) 
# make sum to 1
dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)

PxUS  <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxUS.Rdata")))
PxES  <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxES.Rdata")))




Sums <- outer(.5:110.5,.5:110.5,"+")
Props <- replicate(111,(.5:110.5))/Sums



pUS <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Px, .dxm, .dxf, .Props){
            c(pm = wmean(.Props,ExpectedDx(with(.Px, Male[Year == as.integer(yr)]), .dxm[,yr])),
              pf = wmean(.Props,ExpectedDx(with(.Px, Female[Year == as.integer(yr)]), .dxf[,yr])))    
            
        }, .Px = PxUS, .dxm = dxmUS, .dxf = dxfUS, .Props = Props))
pES <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Px, .dxm, .dxf, .Props){
                    c(pm = wmean(.Props,ExpectedDx(with(.Px, Male[Year == as.integer(yr)]), .dxm[,yr])),
                            pf = wmean(.Props,ExpectedDx(with(.Px, Female[Year == as.integer(yr)]), .dxf[,yr])))    
                    
                }, .Px = PxES, .dxm = dxmES, .dxf = dxfES, .Props = Props))

?cor
cor(yearsUS,pUS[,1])
cor(yearsUS,pUS[,2])
summary(lm(pUS[,1]~yearsUS))
plot(yearsUS, pUS[,1], type = 'l', col = "blue", ylim = c(.5,.58))
lines(yearsUS, pUS[,2], col ="red")

lines(yearsES, pES[,1], col = "blue", lty = 2)
lines(yearsES, pES[,2], col ="red", lty = 2)

pdf("/home/triffe/git/DISS/latex/Figures/PLL.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(yearsUS, pUS[,1], type = 'l', ylim = c(.5, .58), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968,.5,2010,.58,col = gray(.95), border=NA),
                abline(h = seq(.5,.58,by = .005), col = "white"),
                abline(v = seq(1970,2010,by = 5), col = "white"),
                text(seq(1970,2010,by=10), .5,seq(1970,2010,by=10), pos = 1, cex = .8, xpd = TRUE),
                text(1968,seq(.5, .58, by = .01), seq(.5, .58, by = .01), pos = 2, cex = .8, xpd = TRUE),
                text(1963,.585, "Proportion", cex = 1, xpd = TRUE, pos = 4),
                text(1990, .492,"Year", xpd =TRUE, cex = 1)))
lines(yearsUS, pUS[,2], lwd = 2.5, col = gray(.5))
lines(yearsES, pES[,1], lwd = 2, col = gray(.2), lty = 5)
lines(yearsES, pES[,2], lwd = 2.5, col = gray(.5), lty = 5)
legend(1995,.58, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males","US females","ES males","ES females"), xpd = TRUE)
dev.off()


r <- -.003
exp(r)
log(.997)

plot(exp(r*(1:110)))

plot(cumprod(rep(.997,100)))


