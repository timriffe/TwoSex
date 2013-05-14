setwd("/home/triffe/git/DISS/")
source("R/UtilityFunctions.R")
# need:
# Lotka           X
# Pollard         X
# Mitra           X 
# Goodman         X (1,.5,0), with SRB
# Das Gupta       X
# HM, GM, LM, Min X with SRB
# IPF             X with SRB

# get observed SRB

rMitraOLSUS <- local(get(load("Data/results/agerSRB/rMitraOLSUS.Rdata")))
rGoodmanUS  <- local(get(load("Data/results/agerSRB/rGoodmanUS.Rdata")))
rGuptaUS    <- local(get(load("Data/results/agerSRB/rGuptaUS.Rdata")))
rHMUS       <- local(get(load("Data/results/agerSRB/rHMUS.Rdata")))
rGMUS       <- local(get(load("Data/results/agerSRB/rGMUS.Rdata")))
rLMUS       <- local(get(load("Data/results/agerSRB/rLMUS.Rdata")))
rMinUS      <- local(get(load("Data/results/agerSRB/rMinUS.Rdata")))
rIPFhmUS      <- local(get(load("Data/results/agerSRB/rIPFhmUS.Rdata")))
rLotkaUS    <- local(get(load("Data/results/agerSRB/rLotkaUS.Rdata")))
rPollardUS  <- local(get(load("Data/results/agerSRB/rPollardUS.Rdata")))


rMitraOLSES <- local(get(load("Data/results/agerSRB/rMitraOLSES.Rdata")))
rGoodmanES  <- local(get(load("Data/results/agerSRB/rGoodmanES.Rdata")))
rGuptaES    <- local(get(load("Data/results/agerSRB/rGuptaES.Rdata")))
rHMES       <- local(get(load("Data/results/agerSRB/rHMES.Rdata")))
rGMES       <- local(get(load("Data/results/agerSRB/rGMES.Rdata")))
rLMES       <- local(get(load("Data/results/agerSRB/rLMES.Rdata")))
rMinES      <- local(get(load("Data/results/agerSRB/rMinES.Rdata")))
rIPFhmES      <- local(get(load("Data/results/agerSRB/rIPFhmES.Rdata")))
rLotkaES    <- local(get(load("Data/results/agerSRB/rLotkaES.Rdata")))
rPollardES  <- local(get(load("Data/results/agerSRB/rPollardES.Rdata")))

yearsES <- 1975:2009
yearsUS <- 1969:2009
BxymfUS <- local(get(load("Data/USbirths/USBxymf0_110.Rdata")))
BxymfES <- local(get(load("Data/ESbirths/ESBxymf.Rdata")))

# get observed SRB
SRBUS <- unlist(lapply(BxymfUS, function(X){
            sum(X[["Bxym"]], na.rm = TRUE) / sum(X[["Bxyf"]], na.rm = TRUE)
        }))
SRBES <- unlist(lapply(BxymfES, function(X){
                    sum(X[["Bxym"]], na.rm = TRUE) / sum(X[["Bxyf"]], na.rm = TRUE)
                }))
SRBUS <- SRBUS[as.character(yearsUS)]
SRBES <- SRBES[as.character(yearsES)]


rUS <- cbind(rLotkaUS,rPollardUS,rMitraOLSUS[,"r"],rGoodmanUS[,c("r.1.r","r.0.r","r0.5.r")],
        rGuptaUS, rHMUS[,3], rGMUS[,3], rLMUS[,3], rMinUS[,3], rIPFhmUS[,1])
rES <- cbind(rLotkaES,rPollardES,rMitraOLSES[,"r"],rGoodmanES[,c("r.1.r","r.0.r","r0.5.r")],
        rGuptaES, rHMES[,3], rGMES[,3], rLMES[,3], rMinES[,3], rIPFhmES[,1])



dimnames(rUS) <- list(yearsUS, c("$r^m$", "$r^f$", "$r^{Pollard}$","$r^{Mitra}$",
                "$r^{(\\sigma = 1)}$", "$r^{(\\sigma = 0)}$", "$r^{(\\sigma = 0.5)}$",
                "$r^{Gupta}$", "$r^{HM}$", "$r^{GM}$", "$r^{LM}$", "$r^{min}$", "$r^{IPF-HM}$"))
dimnames(rES) <- list(yearsES, c("$r^m$", "$r^f$", "$r^{Pollard}$","$r^{Mitra}$",
                "$r^{(\\sigma = 1)}$", "$r^{(\\sigma = 0)}$", "$r^{(\\sigma = 0.5)}$",
                "$r^{Gupta}$", "$r^{HM}$", "$r^{GM}$", "$r^{LM}$", "$r^{min}$", "$r^{IPF-HM}$"))
dim(rUS)
library(xtable)
print(xtable(rES, digits = c(0,rep(4,13)), align = rep("c",14)),
        sanitize.colnames.function = identity, 
        file = "latex/xtables/ESageALL.tex",floating=FALSE)
print(xtable(rUS, digits = c(0,rep(4,13)), align = rep("c",14)),
        sanitize.colnames.function = identity, 
        file = "latex/xtables/USageALL.tex",floating=FALSE)
# -----------------------------------------------
# sex ratio at birth, initial and stable
SUS <- cbind(SRBUS,rGoodmanUS[,c("SRB1.SRB","SRB0.SRB","SRB0.5.SRB")],
        rHMUS[,4], rGMUS[,4], rLMUS[,4], rMinUS[,4], rIPFhmUS[,2])
SES <- cbind(SRBES,rGoodmanES[,c("SRB1.SRB","SRB0.SRB","SRB0.5.SRB")],
        rHMES[,4], rGMES[,4], rLMES[,4], rMinES[,4], rIPFhmES[,2])

dimnames(SES) <- list(yearsES, c("$S(t)$", 
                "$S^{(\\sigma = 1)}$", "$S^{(\\sigma = 0)}$", "$S^{(\\sigma = 0.5)}$",
                "$S^{HM}$", "$S^{GM}$", "$S^{LM}$", "$S^{min}$", "$S^{IPF-HM}$"))
dimnames(SUS) <- list(yearsUS, c("$S(t)$", 
                "$S^{(\\sigma = 1)}$", "$S^{(\\sigma = 0)}$", "$S^{(\\sigma = 0.5)}$",
                "$S^{HM}$", "$S^{GM}$", "$S^{LM}$", "$S^{min}$", "$S^{IPF-HM}$"))

print(xtable(SES, digits = c(0,rep(5,9)), align = rep("c",10)),
        sanitize.colnames.function = identity, 
        file = "latex/xtables/ESageSRBALL.tex",floating=FALSE)
print(xtable(SUS, digits = c(0,rep(5,9)), align = rep("c",10)),
        sanitize.colnames.function = identity, 
        file = "latex/xtables/USageSRBALL.tex",floating=FALSE)

