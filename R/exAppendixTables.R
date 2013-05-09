
rMeanhmUS   <- local(get(load("/home/triffe/git/DISS/Data/results/exMeanr/rUShm.Rdata")))[,1]
rMeangmUS   <- local(get(load("/home/triffe/git/DISS/Data/results/exMeanr/rUSgm.Rdata")))[,1]
rMeanlmUS   <- local(get(load("/home/triffe/git/DISS/Data/results/exMeanr/rUSlm.Rdata")))[,1]
rIPFhmUS    <- local(get(load("/home/triffe/git/DISS/Data/results/exIPFr/rUSipfhm.Rdata")))[,1]
rSigmaUS    <- local(get(load("/home/triffe/git/DISS/Data/results/exGoodmanr/rSigmaUS.Rdata")))
rSigma0US   <- rSigmaUS[,"r.0.r"]
rSigma.5US  <- rSigmaUS[,"r.5.r"]
rSigma1US   <- rSigmaUS[,"r.1.r"]
rCRhmUS     <- local(get(load("/home/triffe/git/DISS/Data/results/exCPr/rUScp.Rdata")))[,1]
rMUS        <- local(get(load("/home/triffe/git/DISS/Data/results/exSingleSex/rmUS.Rdata")))[,1]
rFUS        <- local(get(load("/home/triffe/git/DISS/Data/results/exSingleSex/rfUS.Rdata")))[,1]

US <- cbind(rMUS,rFUS,rSigma1US,rSigma0US,rSigma.5US,rMeanhmUS,rMeangmUS,rMeanlmUS,rCRhmUS,rIPFhmUS)

rMeanhmES   <- local(get(load("/home/triffe/git/DISS/Data/results/exMeanr/rEShm.Rdata")))[,1]
rMeangmES   <- local(get(load("/home/triffe/git/DISS/Data/results/exMeanr/rESgm.Rdata")))[,1]
rMeanlmES   <- local(get(load("/home/triffe/git/DISS/Data/results/exMeanr/rESlm.Rdata")))[,1]
rIPFhmES    <- local(get(load("/home/triffe/git/DISS/Data/results/exIPFr/rESipfhm.Rdata")))[,1]
rSigmaES    <- local(get(load("/home/triffe/git/DISS/Data/results/exGoodmanr/rSigmaES.Rdata")))
rSigma0ES   <- rSigmaES[,"r.0.r"]
rSigma.5ES  <- rSigmaES[,"r.5.r"]
rSigma1ES   <- rSigmaES[,"r.1.r"]
rCRhmES     <- local(get(load("/home/triffe/git/DISS/Data/results/exCPr/rEScp.Rdata")))[,1]
rMES        <- local(get(load("/home/triffe/git/DISS/Data/results/exSingleSex/rmES.Rdata")))[,1]
rFES        <- local(get(load("/home/triffe/git/DISS/Data/results/exSingleSex/rfES.Rdata")))[,1]

ES <- cbind(rMES,rFES,rSigma1ES,rSigma0ES,rSigma.5ES,rMeanhmES,rMeangmES,rMeanlmES,rCRhmES,rIPFhmES)
dim(US)
dimnames(ES) <- list(yearsES, c("$r^m$", "$r^f$", 
                "$r^{(\\sigma = 1)}$", "$r^{(\\sigma = 0)}$", "$r^{(\\sigma = 0.5)}$",
                "$r^{HM}$", "$r^{GM}$", "$r^{LM}$", "$r^{RAdj-HM}$", "$r^{IPF-HM}$"))
dimnames(US) <- list(yearsUS, c("$r^m$", "$r^f$", 
                "$r^{(\\sigma = 1)}$", "$r^{(\\sigma = 0)}$", "$r^{(\\sigma = 0.5)}$",
                "$r^{HM}$", "$r^{GM}$", "$r^{LM}$", "$r^{RAdj-HM}$", "$r^{IPF-HM}$"))
library(xtable)
print(xtable(ES, digits = c(0,rep(4,10)), align = rep("c",11)),
        sanitize.colnames.function = identity, 
        file = "/home/triffe/git/DISS/latex/xtables/ESexALL.tex",floating=FALSE)
print(xtable(US, digits = c(0,rep(4,10)), align = rep("c",11)),
        sanitize.colnames.function = identity, 
        file = "/home/triffe/git/DISS/latex/xtables/USexALL.tex",floating=FALSE)


# How's about SRB too?

SMeanhmUS   <- local(get(load("/home/triffe/git/DISS/Data/results/exMeanr/rUShm.Rdata")))[,2]
SMeangmUS   <- local(get(load("/home/triffe/git/DISS/Data/results/exMeanr/rUSgm.Rdata")))[,2]
SMeanlmUS   <- local(get(load("/home/triffe/git/DISS/Data/results/exMeanr/rUSlm.Rdata")))[,2]
SIPFhmUS    <- local(get(load("/home/triffe/git/DISS/Data/results/exIPFr/rUSipfhm.Rdata")))[,2]
rSigmaUS    <- local(get(load("/home/triffe/git/DISS/Data/results/exGoodmanr/rSigmaUS.Rdata")))
SSigma0US   <- rSigmaUS[,"S.0.SRB"]
SSigma.5US  <- rSigmaUS[,"S.5.SRB"]
SSigma1US   <- rSigmaUS[,"S.1.SRB"]
SCRhmUS     <- local(get(load("/home/triffe/git/DISS/Data/results/exCPr/rUScp.Rdata")))[,2]
# Spain
SMeanhmES   <- local(get(load("/home/triffe/git/DISS/Data/results/exMeanr/rEShm.Rdata")))[,2]
SMeangmES   <- local(get(load("/home/triffe/git/DISS/Data/results/exMeanr/rESgm.Rdata")))[,2]
SMeanlmES   <- local(get(load("/home/triffe/git/DISS/Data/results/exMeanr/rESlm.Rdata")))[,2]
SIPFhmES    <- local(get(load("/home/triffe/git/DISS/Data/results/exIPFr/rESipfhm.Rdata")))[,2]
rSigmaES    <- local(get(load("/home/triffe/git/DISS/Data/results/exGoodmanr/rSigmaES.Rdata")))
SSigma0ES   <- rSigmaES[,"S.0.SRB"]
SSigma.5ES  <- rSigmaES[,"S.5.SRB"]
SSigma1ES   <- rSigmaES[,"S.1.SRB"]
SCRhmES     <- local(get(load("/home/triffe/git/DISS/Data/results/exCPr/rEScp.Rdata")))[,2]

# get observed SRB:
BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))
names(BxymfES) <- yearsES
names(BxymfUS) <- yearsUS

SRBUS <- unlist(lapply(BxymfUS, function(X){
                    sum(X[["Bxym"]], na.rm = TRUE) / sum(X[["Bxyf"]], na.rm = TRUE)
                }))[as.character(yearsUS)]
SRBES <- unlist(lapply(BxymfES, function(X){
                    sum(X[["Bxym"]], na.rm = TRUE) / sum(X[["Bxyf"]], na.rm = TRUE)
                }))[as.character(yearsES)]

USS <- cbind(SRBUS,SSigma1US,SSigma0US,SSigma.5US,SMeanhmUS,SMeangmUS,SMeanlmUS,SCRhmUS,SIPFhmUS)
ESS <- cbind(SRBES,SSigma1ES,SSigma0ES,SSigma.5ES,SMeanhmES,SMeangmES,SMeanlmES,SCRhmES,SIPFhmES)
dimnames(ESS) <- list(yearsES, c("$S(t)$", 
                "$S^{(\\sigma = 1)}$", "$S^{(\\sigma = 0)}$", "$S^{(\\sigma = 0.5)}$",
                "$S^{HM}$", "$S^{GM}$", "$S^{LM}$", "$S^{RAdj-HM}$", "$S^{IPF-HM}$"))
dimnames(USS) <- list(yearsUS, c("$S(t)$", 
                "$S^{(\\sigma = 1)}$", "$S^{(\\sigma = 0)}$", "$S^{(\\sigma = 0.5)}$",
                "$S^{HM}$", "$S^{GM}$", "$S^{LM}$", "$S^{RAdj-HM}$", "$S^{IPF-HM}$"))
print(xtable(ESS, digits = c(0,rep(5,9)), align = rep("c",10)),
        sanitize.colnames.function = identity, 
        file = "/home/triffe/git/DISS/latex/xtables/ESexSRBALL.tex",floating=FALSE)
print(xtable(USS, digits = c(0,rep(5,9)), align = rep("c",10)),
        sanitize.colnames.function = identity, 
        file = "/home/triffe/git/DISS/latex/xtables/USexSRBALL.tex",floating=FALSE)





