
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")
yearsUS <- 1969:2009
yearsES <- 1975:2009

#PxUS  <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxUS.Rdata")))
#PxES  <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxES.Rdata")))

ExUS  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))

LxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmUS.Rdata"))) / 1e5
LxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfUS.Rdata"))) / 1e5
LxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmES.Rdata"))) / 1e5
LxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfES.Rdata"))) / 1e5


BxymfUS  <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))
BxymfES  <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata"))) 

FxFMUS <- do.call(cbind, lapply(as.character(yearsUS), function(yr, .Bxymf, .Ex){
                    Exf <- with(.Ex, Female[Year == as.integer(yr)])
                    BFM  <- colSums(.Bxymf[[yr]][["Bxym"]])
                    Mna0(Minf0(BFM / Exf))
                }, .Bxymf = BxymfUS, .Ex = ExUS))
FxMFUS <- do.call(cbind, lapply(as.character(yearsUS), function(yr, .Bxymf, .Ex){
                    Exm <- with(.Ex, Male[Year == as.integer(yr)])
                    BMF  <- rowSums(.Bxymf[[yr]][["Bxyf"]])
                    Mna0(Minf0(BMF / Exm))
                }, .Bxymf = BxymfUS, .Ex = ExUS))
FxFMES <- do.call(cbind, lapply(as.character(yearsES), function(yr, .Bxymf, .Ex){
                    Exf <- with(.Ex, Female[Year == as.integer(yr)])
                    BFM  <- colSums(.Bxymf[[yr]][["Bxym"]])
                    Mna0(Minf0(BFM / Exf))
                }, .Bxymf = BxymfES, .Ex = ExES))
FxMFES <- do.call(cbind, lapply(as.character(yearsES), function(yr, .Bxymf, .Ex){
                    Exm <- with(.Ex, Male[Year == as.integer(yr)])
                    BMF  <- rowSums(.Bxymf[[yr]][["Bxyf"]])
                    Mna0(Minf0(BMF / Exm))
                }, .Bxymf = BxymfES, .Ex = ExES))
# equation is interesting
colnames(FxFMUS) <- colnames(FxMFUS) <- yearsUS
colnames(FxFMES) <- colnames(FxMFES) <- yearsES

PollardMin <- compiler::cmpfun(function(r, thetam, thetaf, .a = .5:110.5){
    (1-sum(exp(-outer(.a,.a,"+")*r)*outer(thetam, thetaf)))^2
})

# Niiiiice
#optimize(PollardMin, interval = c(-.02,.02), tol = 1e-15, thetam = thetam,thetaf=thetaf)
# yr <- "1975"
USrPollard <- unlist(lapply(as.character(yearsUS), function(yr, .FxMF, .FxFM, .Lxm, .Lxf){
                    thetam <- .FxMF[, yr] * .Lxm[,yr]
                    thetaf <- .FxFM[, yr] * .Lxf[,yr]
                    optimize(PollardMin, interval = c(-.03,.03), tol = 1e-15, thetam = thetam, thetaf=thetaf)$minimum
                }, .FxMF = FxMFUS, .FxFM = FxFMUS, .Lxm = LxmUS, .Lxf = LxfUS))


ESrPollard <- unlist(lapply(as.character(yearsES), function(yr, .FxMF, .FxFM, .Lxm, .Lxf){
                    thetam <- .FxMF[, yr] * .Lxm[,yr]
                    thetaf <- .FxFM[, yr] * .Lxf[,yr]
                    optimize(PollardMin, interval = c(-.03,.03), tol = 1e-15, thetam = thetam, thetaf=thetaf)$minimum
                }, .FxMF = FxMFES, .FxFM = FxFMES, .Lxm = LxmES, .Lxf = LxfES))


#plot(yearsUS, USr,type = 'l', ylim = c(-.02,.015))
#lines(yearsES, ESr, col = "red")
#abline(h=0)
# confirmed that second unity equation yields same results:
#
#PollardMin2 <- function(r, thetam, thetaf, .a = .5:110.5){
#    (1-(sum(exp(-r*.a)*thetam) * sum(exp(-r*.a)*thetaf)))^2
#}
#
#USr2 <- unlist(lapply(as.character(yearsUS), function(yr, .FxMF, .FxFM, .Lxm, .Lxf){
#                    thetam <- .FxMF[, yr] * .Lxm[,yr]
#                    thetaf <- .FxFM[, yr] * .Lxf[,yr]
#                    optimize(PollardMin2, interval = c(-.03,.03), tol = 1e-15, thetam = thetam, thetaf=thetaf)$minimum
#                }, .FxMF = FxMFUS, .FxFM = FxFMUS, .Lxm = LxmUS, .Lxf = LxfUS))
#plot(yearsUS, USr,type = 'l', ylim = c(-.02,.015))
#lines(yearsUS, USr2, col = "red", lty = 2, lwd =2)
#
#yr <- "1975"
#USr <- unlist(lapply(as.character(yearsUS), function(yr, .FxMF, .FxFM, .Lxm, .Lxf){
#                    thetam <- .FxMF[, yr] * .Lxm[,yr]
#                    thetaf <- .FxFM[, yr] * .Lxf[,yr]
#                    optimize(PollardMin, interval = c(-.03,.03), tol = 1e-15, thetam = thetam, thetaf=thetaf)$minimum
#                }, .FxMF = FxMFUS, .FxFM = FxFMUS, .Lxm = LxmUS, .Lxf = LxfUS))
#names(USr) <- yearsUS
#sum(thetam*exp(-USr[yr]*.a))/
#sum(thetaf*exp(-USr[yr]*.a))
#
#sum(thetam*thetaf*.Lxf[,yr]*exp(-USr[yr]*.a))
#sum(thetam*thetaf*.Lxm[,yr]*exp(-USr[yr]*.a))
#
#(1-sum(exp(-outer(.a,.a,"+")*r)*outer(thetam, thetaf)))^2
#b <- 1 / sum(exp(-outer(.a,.a,"+")*USr[yr])*outer(.Lxm[,yr], .Lxf[,yr]))
#b * .Lxm[,yr]*exp(-.a*USr[yr])



