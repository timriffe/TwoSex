
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")


# exposures, as such, straiht from HMD, all ages 0-110, long form
PxUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxUS.Rdata")))
PxES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxES.Rdata")))
dxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfES.Rdata"))) 
# make sum to 1
dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)

yearsUS <- 1969:2009
yearsES <- 1975:2009

StabUS <- do.call(rbind,lapply(1969:2008, function(yr, .PxUS, .dxmUS, .dxfUS){
            .yr1 <- as.character(yr)
            .yr2 <- as.character(yr+1)
            Pxm <- .PxUS$Male[.PxUS$Year == yr]
            Pxf <- .PxUS$Female[.PxUS$Year == yr]
            Mex <- rowSums(ExpectedDx(Pxm, .dxmUS[,.yr1]))
            Fex <- rowSums(ExpectedDx(Pxf, .dxfUS[,.yr1]))
            
            cx1 <- c(Pxm,Pxf) / sum( c(Pxm,Pxf))
            cEx1 <- c(Mex,Fex) / sum( c(Mex,Fex))
            
            Pxm <- .PxUS$Male[.PxUS$Year == (yr+1)]
            Pxf <- .PxUS$Female[.PxUS$Year == (yr+1)]
            Mex <- rowSums(ExpectedDx(Pxm, .dxmUS[,.yr2]))
            Fex <- rowSums(ExpectedDx(Pxf, .dxfUS[,.yr2]))
            
            cx2 <- c(Pxm,Pxf) / sum( c(Pxm,Pxf))
            cEx2 <- c(Mex,Fex) / sum( c(Mex,Fex))
            
           1-c(Age = sum(pmin(cx1, cx2)), Ex = sum(pmin(cEx1, cEx2)))
            
        }, .PxUS = PxUS, .dxmUS = dxmUS, .dxfUS = dxfUS))
StabES <- do.call(rbind,lapply(1975:2008, function(yr, .PxES, .dxmES, .dxfES){
                    .yr1 <- as.character(yr)
                    .yr2 <- as.character(yr+1)
                    Pxm <- .PxES$Male[.PxES$Year == yr]
                    Pxf <- .PxES$Female[.PxES$Year == yr]
                    Mex <- rowSums(ExpectedDx(Pxm, .dxmES[,.yr1]))
                    Fex <- rowSums(ExpectedDx(Pxf, .dxfES[,.yr1]))
                    
                    cx1 <- c(Pxm,Pxf) / sum( c(Pxm,Pxf))
                    cEx1 <- c(Mex,Fex) / sum( c(Mex,Fex))
                    
                    Pxm <- .PxES$Male[.PxES$Year == (yr+1)]
                    Pxf <- .PxES$Female[.PxES$Year == (yr+1)]
                    Mex <- rowSums(ExpectedDx(Pxm, .dxmES[,.yr2]))
                    Fex <- rowSums(ExpectedDx(Pxf, .dxfES[,.yr2]))
                    
                    cx2 <- c(Pxm,Pxf) / sum( c(Pxm,Pxf))
                    cEx2 <- c(Mex,Fex) / sum( c(Mex,Fex))
                    
                    1-c(Age = sum(pmin(cx1, cx2)), Ex = sum(pmin(cEx1, cEx2)))
                    
                }, .PxES = PxES, .dxmES = dxmES, .dxfES = dxfES))

#
#plot(1969:2008, StabUS[,1], type = 'l', ylim = c(0,.02), col = "red")
#lines(1969:2008, StabUS[,2], col = "blue")
#lines(1975:2008, StabES[,1], col = "red", lty = 2)
#lines(1975:2008, StabES[,2], col = "blue", lty = 2)

pdf("/home/triffe/git/DISS/latex/Figures/PyramidStabilityThetaRatio.pdf", height = 5, width = 5)
par(mai = c(.5, .3, .6, .2),xaxs = "i", yaxs = "i")
plot(1969:2008, StabUS[,2]/StabUS[,1], type = 'l', ylim = c(0,.8), xlim = c(1968,2009), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968,0,2009,1,col = gray(.95), border=NA),
                abline(h = seq(0,.8,by = .1), col = "white"),
                abline(v = seq(1970, 2009, by = 5), col = "white"),
                text(1968, seq(0, .8, by = .1),seq(0, .8, by = .1), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2009, by = 10),0, seq(1970, 2009, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, -.05, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1968.5,.87, expression(frac(theta~e[x],theta~age)), cex = 1, xpd = TRUE)))

lines(1975:2008, StabES[,2]/StabES[,1], lwd = 2.5, col = gray(.5), lty = 5)

legend(1968,.78, lty = c(1,5), col = gray(c(.2,.5)), lwd = c(2,2.5),bty = "n",
        legend = c("US","Spain"), xpd = TRUE)
dev.off()




















