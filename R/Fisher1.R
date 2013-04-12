source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")


yearsUS <- 1969:2009
yearsES <- 1975:2009

# Bxy totals  (not used()
BxUS  <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
BxES  <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

# repeated below for better SRB assumptions
BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))
names(BxymfES) <- yearsES
names(BxymfUS) <- yearsUS
# exposures, as such, straight from HMD, all ages 0-110, long form
ExUS  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
dxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfES.Rdata"))) 
# make sum to 1
dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)


LxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmUS.Rdata"))) 
LxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfUS.Rdata"))) 
LxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmES.Rdata"))) 
LxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfES.Rdata"))) 


#################################################
        

Fishervx <- function(Fx,Lx){
    R0 <- sum(Fx * Lx)    
    rev(cumsum(rev(Fx * Lx))) / R0
}

VxUS <- do.call(rbind,lapply(as.character(yearsUS),function(yr,.Bxy,.Ex,.Lxm,.Lxf){
            Fxm <- Minf0(Mna0(rowSums(.Bxy[[yr]][["Bxym"]]) / with(.Ex, Male[Year == as.integer(yr)])))
            Fxf <- Minf0(Mna0(colSums(.Bxy[[yr]][["Bxyf"]]) / with(.Ex, Female[Year == as.integer(yr)])))
            cbind(Year = as.integer(yr),Age = 0:110,Mvx = Fishervx(Fxm,.Lxm[,yr]), Fvx = Fishervx(Fxf,.Lxf[,yr]))
           
        }, .Bxy = BxymfUS, .Ex = ExUS, .Lxm = LxmUS, .Lxf = LxfUS))
MvxUS <- reshape2::acast(as.data.frame(VxUS), Age~Year, sum, value.var="Mvx")
FvxUS <- reshape2::acast(as.data.frame(VxUS), Age~Year, sum, value.var="Fvx")


VxexES <- do.call(rbind,lapply(as.character(yearsES), function(yr, .Bxymf, .Ex, .dxm, .dxf){
            N <- 111
            # 1) get mx for year
            dxf     <- .dxf[,yr]
            dxm     <- .dxm[,yr]

            BexMM     <- rowSums(ExpectedDx(rowSums(.Bxymf[[yr]][["Bxym"]]), dxm)) 
            Fexm    <- Mna0(Minf0(BexMM / rowSums(ExpectedDx(.Ex$Male[.Ex$Year == as.integer(yr)], dxm))))
            eTFRm     <- sum(Fexm)

            BexFF     <- rowSums(ExpectedDx(colSums(.Bxymf[[yr]][["Bxyf"]]), dxf))
            Fexf    <- Mna0(Minf0(BexFF /  rowSums(ExpectedDx(.Ex$Female[.Ex$Year == as.integer(yr)], dxf))))
            eTFRf     <- sum(Fexf)
   
            cbind(Year = as.integer(yr),y = 0:110,Vxm = cumsum(Fexm / eTFRm), Vxf = cumsum(Fexf / eTFRf))
        }, .Bxymf = BxymfES, .Ex = ExES, .dxm = dxmES, .dxf = dxfES))
VxexUS <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .Bxymf, .Ex, .dxm, .dxf){
                    N <- 111
                    # 1) get mx for year
                    dxf     <- .dxf[,yr]
                    dxm     <- .dxm[,yr]
                    
                    BexMM     <- rowSums(ExpectedDx(rowSums(.Bxymf[[yr]][["Bxym"]]), dxm)) 
                    Fexm      <- Mna0(Minf0(BexMM / rowSums(ExpectedDx(.Ex$Male[.Ex$Year == as.integer(yr)], dxm))))
                    eTFRm     <- sum(Fexm)
                    
                    BexFF     <- rowSums(ExpectedDx(colSums(.Bxymf[[yr]][["Bxyf"]]), dxf))
                    Fexf      <- Mna0(Minf0(BexFF / rowSums(ExpectedDx(.Ex$Female[.Ex$Year == as.integer(yr)], dxf))))
                    eTFRf     <- sum(Fexf)
                    
                    cbind(Year = as.integer(yr),y = 0:110,Vxm = cumsum(Fexm / eTFRm), Vxf = cumsum(Fexf / eTFRf))
                }, .Bxymf = BxymfUS, .Ex = ExUS, .dxm = dxmUS, .dxf = dxfUS))


MvxexUS <- reshape2::acast(as.data.frame(VxexUS), y~Year, sum, value.var="Vxm")
FvxexUS <- reshape2::acast(as.data.frame(VxexUS), y~Year, sum, value.var="Vxf")

image(t(Mvxex/Fvxex))

plot(0:110, MvxexUS[,1], type = 'l',col = "blue")
lines(0:110, MvxUS[,1], col = "red")
lines(0:110, FvxexUS[,1], col = "blue",lty=2)
lines(0:110, FvxUS[,1], col = "red",lty=2)
image(FvxUS[,1:40] - FvxUS[,2:41])
image(FvxUS[1:110,1:40] - FvxUS[2:111,2:41])

.a <- 0:110
pdf("/home/triffe/git/DISS/latex/Figures/vyUS1990.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .5, .3), xaxs = "i", yaxs = "i")
plot(.a, MvxexUS[,"1990"], type = 'l', ylim = c(-.03, 1.03), xlim = c(-1,111), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(-1,-.03,111,1.03,col = gray(.95), border=NA),
                abline(h = seq(0,1,by = .1), col = "white"),
                abline(v = seq(0, 110, by = 10), col = "white"),
                text(seq(0, 110, by = 10), -.03,seq(0, 110, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(-1,seq(0, 1, by = .1), seq(0, 1, by = .1), pos = 2, cex = .8, xpd = TRUE),
                text(-15,1.1, "reproductive value", cex = 1, xpd = TRUE, pos = 4)))
lines(.a, FvxexUS[,"1990"], lwd = 2.5, col = gray(.5))
lines(.a, MvxUS[,"1990"], lwd = 2, col = gray(.2), lty = 5)
lines(.a, FvxUS[,"1990"], lwd = 2.5, col = gray(.5), lty = 5)

legend(70,.3, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c(expression(v[y]~males),expression(v[y]~females),expression(v[x]~males),expression(v[x]~females)), xpd = TRUE)
text(110, -.14, "remaining years", cex = 1, pos = 2, xpd = TRUE)
text(0, -.14, "Age", cex = 1, pos = 4, xpd = TRUE)
arrows(0,-.1,40,-.1,xpd=TRUE,length=.1)
arrows(111,-.1,71,-.1,xpd=TRUE,length=.1)
dev.off()





