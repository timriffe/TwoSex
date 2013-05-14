
# based on Garenne (2013) symmetric IPF 
setwd("/home/triffe/git/DISS/")
source("R/UtilityFunctions.R")
source("R/MeanFunctions.R")

# BxES is 0:110, years 1975:2009
BxUS  <- local(get(load("Data/USbirths/USBxy0_110.Rdata")))
BxES  <- local(get(load("Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

yearsUS <- 1969:2009
yearsES <- 1975:2009

BxymfES <- local(get(load("Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("Data/USbirths/USBxymf0_110.Rdata")))
names(BxymfES) <- yearsES
names(BxymfUS) <- yearsUS

# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS  <- local(get(load("Data/Exposures/USexp.Rdata")))
ExES  <- local(get(load("Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
dxmUS <- local(get(load("Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("Data/HMD_dx/dxfES.Rdata"))) 

dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)

# complement to make difference coef
UScompare <- 1-do.call(rbind,lapply(1969:2008, function(yr, .Ex, .dxm, .dxf, .Bxy){
        yrc1 <- as.character(yr)
        yrc2 <- as.character(yr+1)
        Bxymex1 <- ExpectedDxMxFmatrix(.Bxy[[yrc1]][["Bxym"]], .dxm[, yrc1], .dxf[, yrc1])
        Bxyfex1 <- ExpectedDxMxFmatrix(.Bxy[[yrc1]][["Bxyf"]], .dxm[, yrc1], .dxf[, yrc1])
        BxyTex1 <- Bxymex1 + Bxyfex1
        Bxymex2 <- ExpectedDxMxFmatrix(.Bxy[[yrc2]][["Bxym"]], .dxm[, yrc2], .dxf[, yrc2])
        Bxyfex2 <- ExpectedDxMxFmatrix(.Bxy[[yrc2]][["Bxyf"]], .dxm[, yrc2], .dxf[, yrc2])
        BxyTex2 <- Bxymex2 + Bxyfex2
        Exm1 <- rowSums(ExpectedDx(.Ex$Male[.Ex$Year == yr], .dxm[, yrc1]))
        Exf1 <- rowSums(ExpectedDx(.Ex$Female[.Ex$Year == yr], .dxf[, yrc1]))
        Exm2 <- rowSums(ExpectedDx(.Ex$Male[.Ex$Year == yr+1], .dxm[, yrc2]))
        Exf2 <- rowSums(ExpectedDx(.Ex$Female[.Ex$Year == yr+1], .dxf[, yrc2]))
        
        Bxymexi <- Bxymex1
        Bxyfexi <- Bxyfex1
        BxyTexi <- BxyTex1
        
        # get Ratios of observed to expected:
        RatioM <-  Minf0(Mna0(Bxymex1 / (outer(rowSums(Bxymex1), colSums(Bxymex1)) / sum(Bxymex1))))
        RatioF <-  Minf0(Mna0(Bxyfex1 / (outer(rowSums(Bxyfex1), colSums(Bxyfex1)) / sum(Bxyfex1))))
        RatioT <-  Minf0(Mna0(BxyTex1 / (outer(rowSums(BxyTex1), colSums(BxyTex1)) / sum(BxyTex1))))
        
        # marginal predictions t+1
        BxmMi <- (rowSums(Bxymex1) / Exm1) * Exm2
        BxmFi <- (colSums(Bxymex1) / Exf1) * Exf2
        BxfMi <- (rowSums(Bxyfex1) / Exm1) * Exm2
        BxfFi <- (colSums(Bxyfex1) / Exf1) * Exf2
        BxTMi <- (rowSums(BxyTexi) / Exm1) * Exm2
        BxTFi <- (colSums(BxyTexi) / Exf1) * Exf2
       
        # and corresponding expected bivariate distributions from marginals:
        ExpectedMp <- outer(BxmMi,BxmFi,"*") / ((sum(BxmMi)+sum(BxmFi)) / 2)
        ExpectedFp <- outer(BxfMi,BxfFi,"*") / ((sum(BxfMi)+sum(BxfFi)) / 2)
        ExpectedTp <- outer(BxTMi,BxTFi,"*") / ((sum(BxTMi)+sum(BxTFi)) / 2)
        
        # adjusted with Ratio (scaled back to orig)
        PredM <- (RatioM * ExpectedMp) * (sum(ExpectedMp) / sum(RatioM * ExpectedMp))
        PredF <- (RatioF * ExpectedFp) * (sum(ExpectedFp) / sum(RatioF * ExpectedFp))
        PredT <- (RatioT * ExpectedTp) * (sum(ExpectedTp) / sum(RatioT * ExpectedTp))
        
        # IPF symm 3 ways:
        for (i in 1:10){
            # male births
            # male offer:
            BxymexiM1 <- Bxymexi * Minf0(Mna0(BxmMi / rowSums(Bxymexi)))
            BxymexiM2 <- t(t(BxymexiM1) * Minf0(Mna0((BxmFi / colSums(BxymexiM1)))))
            # female offer
            BxymexiF1 <- t(t(Bxymexi) * Minf0(Mna0((BxmFi / colSums(Bxymexi)))))
            BxymexiF2 <- BxymexiF1 * Minf0(Mna0(BxmMi / rowSums(BxymexiF1)))
            # avg male births iteration
            Bxymexi <- (BxymexiM2 + BxymexiF2) / 2
            
            # female births
            # male offer
            BxyfexiM1 <- Bxyfexi * Minf0(Mna0(BxfMi / rowSums(Bxyfexi)))
            BxyfexiM2 <- t(t(BxyfexiM1) * Minf0(Mna0((BxfFi / colSums(BxyfexiM1)))))
            # female offer
            BxyfexiF1 <- t(t(Bxyfexi) * Minf0(Mna0((BxfFi / colSums(Bxyfexi)))))
            BxyfexiF2 <- BxyfexiF1 * Minf0(Mna0(BxfMi / rowSums(BxyfexiF1)))
            # avg male births iteration
            Bxyfexi <- (BxyfexiM2 + BxyfexiF2) / 2
            
            # Total births
            # male offer
            BxyTexiM1 <- BxyTexi * Minf0(Mna0(BxTMi / rowSums(BxyTexi)))
            BxyTexiM2 <- t(t(BxyTexiM1) * Minf0(Mna0((BxTFi / colSums(BxyTexiM1)))))
            # female offer
            BxyTexiF1 <- t(t(BxyTexi) * Minf0(Mna0((BxTFi / colSums(BxyTexi)))))
            BxyTexiF2 <- BxyTexiF1 * Minf0(Mna0(BxTMi / rowSums(BxyTexiF1)))
            # avg male births iteration
            BxyTexi <- (BxyTexiM2 + BxyTexiF2) / 2
        }
        # now get pdfs
        BxymexPDF <- Bxymexi / sum(Bxymexi)
        BxyfexPDF <- Bxyfexi / sum(Bxyfexi)
        BxyTexPDF <- BxyTexi / sum(BxyTexi)
        
        # year t+1 to test against:
        mt1PDF <- Bxymex2 / sum(Bxymex2)
        ft1PDF <- Bxyfex2 / sum(Bxyfex2)
        Tt1PDF <- BxyTex2 / sum(BxyTex2)
        
        OverlapIt <- c(ipfMal = sum(pmin(BxymexPDF, mt1PDF)),
                ipfFem = sum(pmin(BxyfexPDF, ft1PDF)),
                ipfTot = sum(pmin(BxyTexPDF, Tt1PDF)))
        t0vst1 <- c(pdfMal = sum(pmin(Bxymex1/sum(Bxymex1),mt1PDF)),
                pdfFem = sum(pmin(Bxyfex1/sum(Bxyfex1),ft1PDF)),
                pdfTot = sum(pmin(BxyTex1/sum(BxyTex1),Tt1PDF)))
        # my system
        # pdfs from these can be compared with year t+1
        OverlapNoIt <- c(
                expratMal = sum(pmin(PredM / sum(PredM), mt1PDF)),
                expratFem = sum(pmin(PredF / sum(PredF), ft1PDF)),
                expratTot = sum(pmin(PredT / sum(PredT), Tt1PDF)))
        c(t0vst1, OverlapIt, OverlapNoIt)
        }, .Ex = ExUS, .dxm = dxmUS, .dxf = dxfUS, .Bxy = BxymfUS))
     
 
years <- 1969:2008
plot(years,UScompare[,"pdfTot"], type = 'l', col = "blue", ylim = c(0,.02))
lines(years,UScompare[,"ipfTot"], col = "blue", lty = 5)
lines(years,UScompare[,"expratTot"], col = "blue", lty = 3)

plot(years,UScompare[,"expratTot"] / UScompare[,"pdfTot"], type = 'l', col = "blue")
lines(years,UScompare[,"expratTot"] / UScompare[,"ipfTot"])

plot(years, UScompare[,"expratTot"] / UScompare[,"ipfTot"], type = 'l')
lines(years,UScompare[,"pdfTot"] / UScompare[,"ipfTot"])

plot(years, UScompare[,"expratTot"] , type= 'l')
lines(years,  UScompare[,"pdfTot"], col = "red")

sum(UScompare[,"expratTot"] > UScompare[,"pdfTot"])
sum(UScompare[,"ipfTot"] > UScompare[,"pdfTot"])


var(UScompare[,"expratTot"]) / var(UScompare[,"pdfTot"])
var(UScompare[,"ipfTot"]) / var(UScompare[,"pdfTot"])
var(UScompare[,"expratTot"]) / var(UScompare[,"ipfTot"])

# conclusions:
# 1) my method is easier, no iterations
# 2) my method produces slightly lower variance than ipf in predicted distribution
# 3) my method did better than simple pdf predictions in more years than ipf did = it predicts better
# 4) ipf allows for competition but my method doesn't (directly)

# 5) my method 'might' allow for competition if translated back to age..


McFarlandMarPredict <- function(M,uf,um,nf,nm,tol=1e-6) {
    # stick marriage mat together with those remaining unmarried
    Mmc <- rbind(cbind(M,um),c(uf,sum(uf,um)))
    # rescale rows then columns until margins add up to new margins
    for(i in 1:25){
        Mmc[-nrow(Mmc),] <- Mmc[-nrow(Mmc),]*(nm/(rowSums(Mmc)[-nrow(Mmc)]))
        Mmc[,-ncol(Mmc)] <- t(t(Mmc)[-ncol(Mmc),]*(nf/(colSums(Mmc)[-ncol(Mmc)])))
        if (sum(abs(rowSums(Mmc[-nrow(Mmc),])-nm)+abs(colSums(Mmc[,-ncol(Mmc)])-nf)) < tol) {break}
    }
    # as with M, male age rows, female age cols
    Mmc <- Mmc[-nrow(Mmc),-ncol(Mmc)]
    # male age in rows, female age in cols:
    maleRates <- Mmc/nm
    # female age in rows, male age in cols: (switched, so that all output same)
    femaleRates <- t(Mmc)/nf
    return(list(Marriages=Mmc,maleRates=maleRates,femaleRates=femaleRates))
}
# (same example data as for Henry (SE 1957 *5), from McFarland, 1972
# matrix of marriage counts with female ages in columns and male ages in rows
MAR <- matrix(c(4145,24435,8140,1865,1655,54515,45010,15030,80,6735,20870,19530,5,920,5435,42470),ncol=4)
rownames(MAR) <- colnames(MAR)  <- c("15-19","20-24","25-29","30-60")
# unmarried males and females at start of the year
unMARf  <- c(254876,147705,61804,415497)
unMARm  <- c(265755,199437,114251,429655)
# need knowledge of count of intially unmarried people whose marriage are to be predicted- 
# for example just jitter data a bit
set.seed(1)
initUNMARf2 <- colSums(MAR)+unMARf+runif(4,min=-40000,max=80000)
initUNMARm2 <- rowSums(MAR)+unMARm+runif(4,min=-20000,max=90000)

McFarlandMarPredict(MAR,unMARf,unMARm,initUNMARf2,initUNMARm2)

# try with births:

BxyPred <- McFarlandMarPredict(BxUS[["1969"]],
        with(ExUS, Female[Year == 1969]),
        with(ExUS, Male[Year == 1969]),
        with(ExUS, Female[Year == 1970]),
        with(ExUS, Male[Year == 1970]))

image(log(BxyPred[[1]]))

sum(BxyPred[[1]])
sum(BxUS[["1970"]])

sum(rowSums(BxUS[["1969"]] / with(ExUS, Male[Year == 1969])) * with(ExUS, Male[Year == 1970]))
sum((colSums(BxUS[["1969"]]) / with(ExUS, Female[Year == 1969])) * with(ExUS, Female[Year == 1970]))


# TODO: modify Garenne iterative procedure to accept Stolarsky mean.
# find which stolarsky parameter give harmonic mean. Make a function to do that.

GarenneGeneral <- compiler::cmpfun(function(Bxy1, Exm1, Exf1, Exm2, Exf2, p = -1){
    Mmarg <- (rowSums(Bxy1) / Exm1) * Exm2
    Fmarg <- (colSums(Bxy1) / Exf1) * Exf2
    Bxyi <- Bxy1
    for (i in 1:10){
        
        BxyiM1 <- Bxyi * Minf0(Mna0(Mmarg / rowSums(Bxyi)))
        BxyiM2 <- t(t(BxyiM1) * Minf0(Mna0((Fmarg / colSums(BxyiM1)))))
    # female offer
        BxyiF1 <- t(t(Bxyi) * Minf0(Mna0((Fmarg / colSums(Bxyi)))))
        BxyiF2 <- BxyiF1 * Minf0(Mna0(Mmarg / rowSums(BxyiF1)))
    # avg male births iteration
   
        Bxyi <- stolarsky.mean.v(BxyiM2, BxyiF2, p)
    } 
    Bxyi
})

BxyPred <- GarenneGeneral(BxUS[["1969"]], with(ExUS, Male[Year == 1969]),with(ExUS, Female[Year == 1969]),
        with(ExUS, Male[Year == 1970]),with(ExUS, Female[Year == 1970]))
