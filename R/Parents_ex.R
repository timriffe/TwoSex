

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

yearsUS <- 1969:2009
yearsES <- 1975:2009
xlabs   <- c("1.0%","0.8%","0.6%","0.4%","0.2%","0.0%","0.2%","0.4%","0.6%","0.8%","1.0%")


Mex75 <- rowSums(ExpectedDx(Px = with(ExUS, Male[Year == 1975]), dx = dxmUS[,"1975"]))
Fex75 <- rowSums(ExpectedDx(Px = with(ExUS, Female[Year == 1975]), dx = dxfUS[,"1975"]))
Mex09 <- rowSums(ExpectedDx(Px = with(ExUS, Male[Year == 2009]), dx = dxmUS[,"2009"]))
Fex09 <- rowSums(ExpectedDx(Px = with(ExUS, Female[Year == 2009]), dx = dxfUS[,"2009"]))
pdf("/home/triffe/git/DISS/latex/Figures/exPyramidUS.pdf", height = 5, width = 5)
par(mai = c(.5,.3,.3,.3))
plot(NULL, type = "n",axes = FALSE, xlab = "",ylab = "", xlim = c(-1, 1), ylim = c(0,111),
        panel.first = list(
                rect(-1, 0, 1, 111, col = gray(.95), border = NA),
                abline(v = seq(-1, 1, by = .2), col = "white"),
                abline(h = seq(0, 110, by = 10), col = "white"),
                text(seq(-1, 1, by = .2),0, xlabs, xpd = TRUE, pos = 1, cex = .7),
                text(-1, seq(0, 110, by = 10), seq(0, 110, by = 10), pos = 2, xpd = TRUE, cex = .7),
                text(-1.15, 118, expression(e[x]), xpd = TRUE, cex = 1),
                text(0, -15, "Percentage", xpd = TRUE, cex = 1),
                text(c(-.5, .5), 115, c("Males", "Females"), cex = .9, xpd = TRUE)
               ))
barplot(-100*(Mex75/sum(Mex75+Fex75)), border = NA, col = "#44444450", 
        add = TRUE, horiz = TRUE, space = 0, axes = FALSE)
barplot(100*(Fex75/sum(Mex75+Fex75)), border = NA, col = "#44444450", 
        add = TRUE, horiz = TRUE, space = 0, axes = FALSE)
PyramidOutline(Mex09, Fex09, prop = TRUE, border = gray(.2), xpd = TRUE, lwd = 1)    
text(c(-.5,-.55), c(50, 91), c("1975", "2009"), col = c("white", "black"), cex = 1.2, pos = 4)
segments(-.4,88,-.3,83)
dev.off()
# for Spain:
Mex75ES <- rowSums(ExpectedDx(Px = with(ExES, Male[Year == 1975]), dx = dxmES[,"1975"]))
Fex75ES <- rowSums(ExpectedDx(Px = with(ExES, Female[Year == 1975]), dx = dxfES[,"1975"]))
Mex09ES <- rowSums(ExpectedDx(Px = with(ExES, Male[Year == 2009]), dx = dxmES[,"2009"]))
Fex09ES <- rowSums(ExpectedDx(Px = with(ExES, Female[Year == 2009]), dx = dxfES[,"2009"]))

pdf("/home/triffe/git/DISS/latex/Figures/exPyramidES.pdf", height = 5, width = 5)
par(mai = c(.5,.3,.3,.3))
plot(NULL, type = "n",axes = FALSE, xlab = "",ylab = "", xlim = c(-1, 1), ylim = c(0,111),
        panel.first = list(
                rect(-1, 0, 1, 111, col = gray(.95), border = NA),
                abline(v = seq(-1, 1, by = .2), col = "white"),
                abline(h = seq(0, 110, by = 10), col = "white"),
                text(seq(-1, 1, by = .2),0, xlabs, xpd = TRUE, pos = 1, cex = .7),
                text(-1, seq(0, 110, by = 10), seq(0, 110, by = 10), pos = 2, xpd = TRUE, cex = .7),
                text(-1.15, 118, expression(e[x]), xpd = TRUE, cex = 1),
                text(0, -15, "Percentage", xpd = TRUE, cex = 1),
                text(c(-.5, .5), 115, c("Males", "Females"), cex = .9, xpd = TRUE)
        ))
barplot(-100*(Mex75ES/sum(Mex75ES+Fex75ES)), border = NA, col = "#44444450", 
        add = TRUE, horiz = TRUE, space = 0, axes = FALSE)
barplot(100*(Fex75ES/sum(Mex75ES+Fex75ES)), border = NA, col = "#44444450", 
        add = TRUE, horiz = TRUE, space = 0, axes = FALSE)
PyramidOutline(Mex09ES, Fex09ES, prop = TRUE, border = gray(.2), xpd = TRUE, lwd = 1)   
text(c(-.5,-.55), c(50, 90), c("1975", "2009"), col = c("white", "black"), cex = 1.2, pos = 4)
segments(-.27,90.5,-.17,87)
dev.off()
    

    
    
    
    
    


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




Meanex <- function(Bx, .dx, x = 0:110 + .5){
    ex <- rowSums(ExpectedDx(Px = Bx, dx = .dx))
    wmean(x, ex)
}

exES <- do.call(rbind,lapply(as.character(yearsES), function(yr, .BxES, .dxf, .dxm){
            c(Father_ex = Meanex(rowSums(Mna0(.BxES[[yr]])), .dxm[, yr]),
              Mother_ex = Meanex(colSums(Mna0(.BxES[[yr]])), .dxf[, yr]))
        }, .BxES = BxES, .dxf = dxfES, .dxm = dxmES))

plot(yearsES, exES[, "Father_ex"] - exES[, "Father_ex"][1], type = 'l', col = "blue", ylim = c(0,5))
lines(yearsES, exES[, "Mother_ex"] - exES[, "Mother_ex"][1], col = "red")

plot(yearsES, exES[, "Mother_ex"] - exES[, "Father_ex"])


exUS <- do.call(rbind,lapply(as.character(yearsUS), function(yr, .BxUS, .dxf, .dxm){
                    c(Father_ex = Meanex(rowSums(Mna0(.BxUS[[yr]])), .dxm[, yr]),
                            Mother_ex = Meanex(colSums(Mna0(.BxUS[[yr]])), .dxf[, yr]))
                }, .BxUS = BxUS, .dxf = dxfUS, .dxm = dxmUS))

plot(yearsUS, exUS[, "Father_ex"], type = 'l', col = "blue", ylim = c(40,55))
lines(yearsUS, exUS[, "Mother_ex"], col = "red")

# compare female advantage tiems series, ES, US
plot(yearsUS, exUS[, "Mother_ex"] - exUS[, "Father_ex"], col = "blue",type= 'l')
lines(yearsES, exES[, "Mother_ex"] - exES[, "Father_ex"], col = "red")

# what about rates by ex instead of Exp?
colnames(dxmUS)
FxExUS <- lapply(as.character(yearsUS), function(yr, .ExpUS, .BxUS, .dxfUS, .dxmUS){
            Bexm <- rowSums(ExpectedDx(Px = rowSums(Mna0(.BxUS[[yr]])), dx =  .dxmUS[,yr]))
            Bexf <- rowSums(ExpectedDx(Px = colSums(Mna0(.BxUS[[yr]])), dx =  .dxfUS[,yr]))
            Exp.exm <- rowSums(ExpectedDx(Px = with(.ExpUS, Male[Year == as.integer(yr)]), dx =  .dxmUS[,yr]))
            Exp.exf <- rowSums(ExpectedDx(Px = with(.ExpUS, Female[Year == as.integer(yr)]), dx =  .dxfUS[,yr]))
            list(F.exf = Bexf / Exp.exf,
            M.exf = Bexm / Exp.exm)
        }, .ExpUS = ExUS, .BxUS = BxUS, .dxfUS = dxfUS, .dxmUS = dxmUS)

FxExUSf <- do.call(cbind, lapply(FxExUS, "[[", 1))
FxExUSm <- do.call(cbind, lapply(FxExUS, "[[", 2))
colnames(FxExUSf) <-colnames(FxExUSm) <- yearsUS
FxExES <- lapply(as.character(yearsES), function(yr, .ExpES, .BxES, .dxfES, .dxmES){
            yr <- "1975"
            Bexm <- rowSums(ExpectedDx(Px = rowSums(Mna0(.BxES[[yr]])), dx =  .dxmES[,yr]))
            Bexf <- rowSums(ExpectedDx(Px = colSums(Mna0(.BxES[[yr]])), dx =  .dxfES[,yr]))
            Exp.exm <- rowSums(ExpectedDx(Px = with(.ExpES, Male[Year == as.integer(yr)]), dx =  .dxmES[,yr]))
            Exp.exf <- rowSums(ExpectedDx(Px = with(.ExpES, Female[Year == as.integer(yr)]), dx =  .dxfES[,yr]))
            list(F.exf = Bexf / Exp.exf,
                    M.exf = Bexm / Exp.exm)
        }, .ExpES = ExES, .BxES = BxES, .dxfES = dxfES, .dxmES = dxmES)
FxExESf <- do.call(cbind, lapply(FxExES, "[[", 1))
FxExESm <- do.call(cbind, lapply(FxExES, "[[", 2))
colnames(FxExESf) <-colnames(FxExESm) <- yearsES

#image(t(FxExUSf))
#image(t(FxExUSm))
#image(t(FxExESf))
#image(t(FxExESm))
plot(0:110, FxExUSf[,1], type = 'l')

plot(yearsUS, colSums(FxExUSf, na.rm = TRUE), type = 'l', col = "red", ylim = c(1,3))
lines(yearsUS, colSums(FxExUSm, na.rm = TRUE), col = "blue")
lines(yearsES, colSums(FxExESf, na.rm = TRUE), col = "red",lty=2)
lines(yearsES, colSums(FxExESm, na.rm = TRUE), col = "blue",lty=2)





