
# periodizing Henry 1965:
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")
yearsES <- 1975:2009
yearsUS <- 1969:2009
# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
# BxES is 0:110, years 1975:2009
BxymfES <-  local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <-  local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))
names(BxymfES) <- yearsES
#BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf10_65.Rdata")))
#BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf10_65.Rdata")))

# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS    <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES    <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
# life expectancies
E0US    <- local(get(load("/home/triffe/git/DISS/Data/HMD_e0period/e0perUS.Rdata")))
E0ES    <- local(get(load("/home/triffe/git/DISS/Data/HMD_e0period/e0perES.Rdata")))
# lifetable Lx:
LxmUS   <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmUS.Rdata"))) / 1e5
LxfUS   <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfUS.Rdata"))) / 1e5
LxmES   <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmES.Rdata"))) / 1e5
LxfES   <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfES.Rdata"))) / 1e5




# Henry R0star
head(E0ES)
R0starUS <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .BxymfUS, .ExUS, .E0US, .LxmUS, .LxfUS){
            Bxm          <- rowSums(Mna0(.BxymfUS[[yr]][["Bxym"]]))
            Bxf          <- colSums(Mna0(.BxymfUS[[yr]][["Bxyf"]]))
            Exm          <- with(.ExUS, Male[Year == as.integer(yr)])
            Exf          <- with(.ExUS, Female[Year == as.integer(yr)])
            Fxm          <- Bxm / Exm
            Fxf          <- Bxf / Exf
            # get e0t
            e0m          <- with(.E0US, Male[Year == as.integer(yr)])
            e0f          <- with(.E0US, Female[Year == as.integer(yr)])
            # get prev e0s
            e0vecm       <- rev(with(.E0US, Male[Year <= as.integer(yr)]))
            e0vecf       <- rev(with(.E0US, Female[Year <= as.integer(yr)]))
         
            N <- length(e0vecm)
            if (N < length(Fxm)){
                e0vecm[(N + 1):length(Fx)] <- e0vecm[N]
                e0vecf[(N + 1):length(Fx)] <- e0vecf[N]
            }
            Lxm <- .LxmUS[,yr]
            Lxf <- .LxfUS[,yr]
            c(R0starm = sum(Lxm * Fxm * (e0m / e0vecm), na.rm = TRUE),
                    R0starf = sum(Lxf * Fxf * (e0f / e0vecf), na.rm = TRUE))
        },.BxymfUS = BxymfUS, .ExUS = ExUS, .E0US = E0US, .LxmUS = LxmUS, .LxfUS = LxfUS))

#plot(yearsES, R0starES[,"R0starm"], type = 'l',col = "blue", lty = 2, ylim = c(.5,2))
#lines(yearsES, R0starES[,"R0starf"], col = "red", lty = 2)
#lines(yearsUS, R0starUS[,"R0starm"], col = "blue")
#lines(yearsUS, R0starUS[,"R0starf"], col = "red")


pdf("/home/triffe/git/DISS/latex/Figures/R0perHenry.pdf", height = 5, width = 5)

par(mar = c(3, 3, 2, 2),xaxp = "i", yaxp = "i")
plot(yearsUS, R0starUS[, "R0starm"], type = 'l', ylim = c(.5, 2), xlim = c(1968,2010), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(1968,.5,2010,2,col = gray(.95), border=NA),
                abline(h = seq(.6,2,by = .2),col = "white"),
                abline(v = seq(1970, 2010, by = 5), col = "white"),
                text(1968, seq(.6, 2, by = .2),seq(.6, 2, by = .2), pos = 2, cex = .8, xpd = TRUE),
                text(seq(1970, 2010, by = 10),.5, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(1990, .4, "Year", cex = 1, pos = 1, xpd = TRUE),
                text(1967,2.1, "R0*", cex = 1, xpd = TRUE)))
lines(yearsUS, R0starUS[, "R0starf"], lwd = 2.5, col = gray(.5))
lines(yearsES, R0starES[, "R0starm"], lwd = 2, col = gray(.2), lty = 5)
lines(yearsES, R0starES[, "R0starf"], lwd = 2.5, col = gray(.5), lty = 5)

legend(1990,2, lty = c(1,1,5,5), col = gray(c(.2,.5,.2,.5)), lwd = c(2,2.5,2,2.5),bty = "n",
        legend = c("US males", "US females", "ES males", "ES females"), xpd = TRUE)
dev.off()