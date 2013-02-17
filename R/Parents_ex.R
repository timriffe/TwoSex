

# ----------------------------------------------
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
# BxES is 0:110, years 1975:2009
BxUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
BxES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

#BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf10_65.Rdata")))
#BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf10_65.Rdata")))

# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
dxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfES.Rdata"))) 

dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)

ParMUS <- rowSums(BxUS[["2009"]])
ParFUS <- colSums(BxUS[["2009"]])


library(Pyramid)

for (yr in as.character(yearsUS)){
    ParMUS <- rowSums(BxUS[[yr]])
    ParFUS <- colSums(BxUS[[yr]])
    
Pyramid(males = ExpectedDx(Px = ParMUS, dx = dxmUS[,1]),
        females = ExpectedDx(Px = ParFUS, dx = dxfUS[,1]),
        widths = rep(1, nrow(dxmUS)),
        fill.males = gray(seq(.1, .9, length = nrow(dxmUS))),
        fill.females = gray(seq(.1, .9, length = nrow(dxmUS))),
        border.males = NA, border.females = NA, grid = FALSE,
        main = yr, ylab.left = "(1933 + t) , Age",cex.axis = .7,
        xlim = c(-1.7,1.7), prop = TRUE, box = FALSE, mar = c(4,4,4,4)
)
Sys.sleep(1)
}


par(mfrow =c(1,2))
ParMUS <- rowSums(BxUS[["1970"]])
ParFUS <- colSums(BxUS[["1970"]])

Pyramid(males = ExpectedDx(Px = ParMUS, dx = dxmUS[,1]),
        females = ExpectedDx(Px = ParFUS, dx = dxfUS[,1]),
        widths = rep(1, nrow(dxmUS)),
        fill.males = gray(seq(.1, .9, length = nrow(dxmUS))),
        fill.females = gray(seq(.1, .9, length = nrow(dxmUS))),
        border.males = NA, border.females = NA, grid = FALSE,
        main = "1970", ylab.left = "(1933 + t) , Age",cex.axis = .7,
        xlim = c(-1.7,1.7), prop = TRUE, box = FALSE, mar = c(4,4,4,4)
)
ParMUS <- rowSums(BxUS[["2009"]])
ParFUS <- colSums(BxUS[["2009"]])

Pyramid(males = ExpectedDx(Px = ParMUS, dx = dxmUS[,"2009"]),
        females = ExpectedDx(Px = ParFUS, dx = dxfUS[,"2009"]),
        widths = rep(1, nrow(dxmUS)),
        fill.males = gray(seq(.1, .9, length = nrow(dxmUS))),
        fill.females = gray(seq(.1, .9, length = nrow(dxmUS))),
        border.males = NA, border.females = NA, grid = FALSE,
        main = "2009", ylab.left = "(1933 + t) , Age",cex.axis = .7,
        xlim = c(-1.7,1.7), prop = TRUE, box = FALSE, mar = c(4,4,4,4)
)


for (yr in c("1975","2009")){
    ParM <- rowSums(BxUS[[yr]])
    ParF <- colSums(BxUS[[yr]])
    Pyramid(males = ExpectedDx(Px = ParM, dx = dxmUS[,yr]),
            females = ExpectedDx(Px = ParF, dx = dxfUS[,yr]),
            widths = rep(1, nrow(dxmUS)),
            fill.males = gray(seq(.1, .9, length = nrow(dxmUS))),
            fill.females = gray(seq(.1, .9, length = nrow(dxmUS))),
            border.males = NA, border.females = NA, grid = FALSE,
            main = yr, ylab.left = "(1933 + t) , Age",cex.axis = .7,
            xlim = c(-2,2), prop = TRUE, box = FALSE, mar = c(4,4,4,4)
    )
}

for (yr in c("1975","2009")){
    ParM <- rowSums(BxES[[yr]])
    ParF <- colSums(BxES[[yr]])
    Pyramid(males = ExpectedDx(Px = ParM, dx = dxmES[,yr]),
            females = ExpectedDx(Px = ParF, dx = dxfES[,yr]),
            widths = rep(1, nrow(dxmUS)),
            fill.males = gray(seq(.1, .9, length = nrow(dxmUS))),
            fill.females = gray(seq(.1, .9, length = nrow(dxmUS))),
            border.males = NA, border.females = NA, grid = FALSE,
            main = yr, ylab.left = "(1933 + t) , Age",cex.axis = .7,
            xlim = c(-2,2), prop = TRUE, box = FALSE, mar = c(4,4,4,4)
    )
}


ityrs <- c("1975","2009")
cols <- c("red","blue")
plot(NULL, type = "n", xlim = c(-2,2), ylim = c(0,111))
for (i in 1:2){
    ParM <- rowSums(BxUS[[ityrs[i]]])
    ParF <- colSums(BxUS[[ityrs[i]]])
    
    PyramidOutline(
            rowSums(ExpectedDx(Px = ParM, dx = dxmUS[,ityrs[i]])), 
            rowSums(ExpectedDx(Px = ParF, dx = dxfUS[,ityrs[i]])), prop = ,
            TRUE,border = cols[i], xpd = TRUE, lwd = 1)    
}

plot(NULL, type = "n", xlim = c(-2,2), ylim = c(0,111))
for (i in 1:2){
    ParM <- rowSums(BxES[[ityrs[i]]])
    ParF <- colSums(BxES[[ityrs[i]]])
    
    PyramidOutline(
            rowSums(ExpectedDx(Px = ParM, dx = dxmES[,ityrs[i]])), 
            rowSums(ExpectedDx(Px = ParF, dx = dxfES[,ityrs[i]])), prop = ,
            TRUE,border = cols[i], xpd = TRUE, lwd = 1)    
}

# Spain
Bxm <- rowSums(BxES[["1975"]][11:66,])
Bxf <- colSums(BxES[["1975"]][,11:66])
plot(10:65, Bxm / sum(Bxm), type = 'l', col = "blue", xlab = "Age",ylab = "Birth density",main = "Spain")
lines(10:65, Bxf / sum(Bxf), col = "red")

Bxm <- rowSums(BxES[["2009"]][11:66,])
Bxf <- colSums(BxES[["2009"]][,11:66])
lines(10:65, Bxm / sum(Bxm), col = "blue", lty = 2)
lines(10:65, Bxf / sum(Bxf), col = "red", lty = 2)
legend("topright", lty = c(1,1,2,2), col = c("blue","red","blue","red"), 
        legend = c("fathers 1975", "mothers 1975","fathers 2009", "mothers 2009"))


# US
Bxm <- rowSums(BxUS[["1975"]][11:66,])
Bxf <- colSums(BxUS[["1975"]][,11:66])
plot(10:65, Bxm / sum(Bxm), type = 'l', col = "blue",xlab = "Age",ylab = "Birth density",main = "USA")
lines(10:65, Bxf / sum(Bxf), col = "red")

Bxm <- rowSums(BxUS[["2009"]][11:66,])
Bxf <- colSums(BxUS[["2009"]][,11:66])
lines(10:65, Bxm / sum(Bxm), col = "blue", lty = 2)
lines(10:65, Bxf / sum(Bxf), col = "red", lty = 2)

legend("topright", lty = c(1,1,2,2), col = c("blue","red","blue","red"), 
        legend = c("fathers 1975", "mothers 1975","fathers 2009", "mothers 2009"))

library(DecompHoriuchi)
