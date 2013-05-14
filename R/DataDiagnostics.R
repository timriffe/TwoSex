setwd("/home/triffe/git/DISS/")
# 1) missingness in US data:
yearsUS <- 1969:2009
USNA <- local(get(load("Data/results/USmissings/USNAt.Rdata")))
names(USNA) <- yearsUS
USNA <- USNA[as.character(yearsUS)] # remove 2010
USNAt <- unlist(lapply(USNA, "[[",1))
pdf("latex/Figures/USmissingAge.pdf", height = 5, width = 5)
par(mai = c(.5,.5,.5,.3), xaxs = "i", yaxs = "i")
plot(yearsUS, USNAt, type = 'l', ylim = c(.07,.18),xlim = c(1968,2010), axes = FALSE,
                panel.first = list(rect(1968, 0.07, 2010, .18,col = gray(.95), border=NA),
                        abline(h = seq(.08,0.18,by = .02), col = "white"),
                        abline(v = seq(1970, 2010, by = 5), col = "white"),
                        text(1968, seq(.08,0.18,by = .02),seq(.08,0.18,by = .02),pos = 2, cex = .8, xpd = TRUE),
                        text(seq(1970, 2010, by = 10),0.07, seq(1970, 2010, by = 10), pos = 1, cex = .8, xpd = TRUE),
                        text(1990, .065, "Year", cex = 1, pos = 1, xpd = TRUE),
                        text(1969,.188, "Prop. Missing", cex = 1, xpd = TRUE)))
dev.off()
ages<- 10:69
allNA <-rep(0,length(ages))
names(allNA) <- 10:69
for (i in 1:length(USNA)){
           NAx <- Minf0(Mna0(USNA[[i]][[2]]))
    
           allNA[names(NAx)] <- allNA[names(NAx)] + NAx
        }
allNA <- allNA[as.character(ages)]
avgprop <- allNA / length(USNA)
ages[avgprop > .2]
ages[avgprop > .4]
ages[avgprop > .6]
ages[avgprop > .8]

