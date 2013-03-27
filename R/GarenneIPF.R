
# based on Garenne (2013) symmetric IPF 

source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxES is 0:110, years 1975:2009
BxUS  <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
BxES  <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65

yearsUS <- 1969:2009
yearsES <- 1975:2009

BxymfES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxymf.Rdata")))
BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))
names(BxymfES) <- yearsES
names(BxymfUS) <- yearsUS

# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
# get Lx estimates for R0, r
dxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata"))) 
dxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata"))) 
dxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmES.Rdata"))) 
dxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfES.Rdata"))) 

dxmUS <- dxmUS %col% colSums(dxmUS)
dxfUS <- dxfUS %col% colSums(dxfUS)
dxmES <- dxmES %col% colSums(dxmES)
dxfES <- dxfES %col% colSums(dxfES)


Bxym <- BxymfUS[["1969"]][["Bxym"]]
Bxyf <- BxymfUS[["1969"]][["Bxyf"]]

Bxymex <- ExpectedDxMxFmatrix(Bxym, dxmUS[,"1969"], dxfUS[,"1969"])
Bxyfex <- ExpectedDxMxFmatrix(Bxyf, dxmUS[,"1969"], dxfUS[,"1969"])

Exm <- rowSums(ExpectedDx(ExUS$Male[ExUS$Year == 1969], dxmUS[,"1969"]))
Exf <- rowSums(ExpectedDx(ExUS$Female[ExUS$Year == 1969], dxfUS[,"1969"]))

Bxyex <- Bxymex + Bxyfex 

Fxym <- rowSums(Bxyex) / Exm
Fxyf <- colSums(Bxyex) / Exf

Exm2 <- rowSums(ExpectedDx(ExUS$Male[ExUS$Year == 1970], dxmUS[,"1970"]))
Exf2 <- rowSums(ExpectedDx(ExUS$Female[ExUS$Year == 1970], dxfUS[,"1970"]))

Bxmi <- Fxym * Exm2
Bxfi <- Fxyf * Exf2
sum(Bxmi)
sum(Bxfi)
Bxyexi <- Bxyex
for (i in 1:10){
    # male offer:
    BxyexiM1 <- Bxyexi * Minf0(Mna0(Bxmi / rowSums(Bxyexi)))
    BxyexiM2 <- t(t(BxyexiM1) * Minf0(Mna0((Bxfi / colSums(BxyexiM1)))))
    BxyexiF1 <- t(t(Bxyexi) * Minf0(Mna0((Bxfi / colSums(Bxyexi)))))
    BxyexiF2 <- BxyexiF1 * Minf0(Mna0(Bxmi / rowSums(BxyexiF1)))
    Bxyexi <- (BxyexiM2 + BxyexiF2) / 2
}

Bxyexipdf <- Bxyexi / sum(Bxyexi)

Bxyex2 <- ExpectedDxMxFmatrix( BxymfUS[["1970"]][["Bxym"]], dxmUS[,"1970"], dxfUS[,"1970"]) + 
        ExpectedDxMxFmatrix( BxymfUS[["1970"]][["Bxyf"]], dxmUS[,"1970"], dxfUS[,"1970"])

Bxyex2pdf <- Bxyex2 / sum(Bxyex2)

sum(pmin(Bxyexipdf, Bxyex2pdf))

Bxyexpdf <- Bxyex / sum(Bxyex)
sum(pmin(Bxyexpdf, Bxyex2pdf))

image(Bxyexipdf - Bxyex2pdf)

Bxm2 <- rowSums(ExpectedDx(rowSums(BxymfUS[["1970"]][["Bxym"]]), dxmUS[,"1970"]) + 
                  ExpectedDx(rowSums(BxymfUS[["1970"]][["Bxyf"]]), dxmUS[,"1970"]))
Bxf2 <- rowSums(ExpectedDx(colSums(BxymfUS[["1970"]][["Bxym"]]), dxfUS[,"1970"]) + 
                  ExpectedDx(colSums(BxymfUS[["1970"]][["Bxyf"]]), dxfUS[,"1970"]))
  
Bxyex1970 <- outer(Bxm2,Bxf2,"*") / sum(Bxm2)
Bxyexex <- outer(Bxmi,Bxfi,"*") / ((sum(Bxmi)+sum(Bxfi)) / 2)
Bxyex
BxyexEX <- outer(rowSums(Bxyex),colSums(Bxyex),"*") / sum(Bxyex)

# relates observed vs expected original (static)
Ratio <- Minf0(Mna0(Bxyex / BxyexEX))
# scale IPF output by ratio to adjust shape, then scale to proper total
Pred <- (Ratio * Bxyexex) * (sum(Bxyexex) / sum(Ratio * Bxyexex))

# compare pdfs
sum(pmin(Pred / sum(Pred),Bxyex2pdf)) # pdf predicted from 1969 vs 1970 pdf
sum(pmin(Bxyexpdf, Bxyex2pdf)) # 1969 vs 1970 fit

yr <- 1969
UScompare <- do.call(rbind,lapply(1969:2008, function(yr, .Ex, .dxm, .dxf, .Bxy){
        yrc1 <- as.character(yr)
        yrc2 <- as.character(yr+1)
        Bxymex1 <- ExpectedDxMxFmatrix(.Bxy[[yrc1]][["Bxym"]], .dxm[, yrc1], .dxf[, yrc1])
        Bxyfex1 <- ExpectedDxMxFmatrix(.Bxy[[yrc1]][["Bxyf"]], .dxm[, yrc1], .dxf[, yrc1])
        Bxymex2 <- ExpectedDxMxFmatrix(.Bxy[[yrc2]][["Bxym"]], .dxm[, yrc2], .dxf[, yrc2])
        Bxyfex2 <- ExpectedDxMxFmatrix(.Bxy[[yrc2]][["Bxyf"]], .dxm[, yrc2], .dxf[, yrc2])
        BxyTex2 <- Bxymex2 + Bxyfex2
        Exm1 <- rowSums(ExpectedDx(.Ex$Male[.Ex$Year == yr], .dxm[, yrc1]))
        Exf1 <- rowSums(ExpectedDx(.Ex$Female[.Ex$Year == yr], .dxf[, yrc1]))
        Exm2 <- rowSums(ExpectedDx(.Ex$Male[.Ex$Year == yr+1], .dxm[, yrc2]))
        Exf2 <- rowSums(ExpectedDx(.Ex$Female[.Ex$Year == yr+1], .dxf[, yrc2]))
        
        Bxymexi <- Bxymex1
        Bxyfexi <- Bxyfex1
        BxyTexi <- Bxymex1 + Bxyfex1
        
        # get Ratios of observed to expected:
        RatioM <-  Minf0(Mna0(Bxymex1 / (outer(rowSums(Bxymex1), colSums(Bxymex1)) / sum(Bxymex1))))
        RatioF <-  Minf0(Mna0(Bxyfex1 / (outer(rowSums(Bxyfex1), colSums(Bxyfex1)) / sum(Bxyfex1))))
        RatioT <-  Minf0(Mna0(BxyTexi / (outer(rowSums(BxyTexi), colSums(BxyTexi)) / sum(BxyTexi))))
        
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
        
        # adjusted with Ratio:
        PredM <- (RatioM * ExpectedMp) * (sum(ExpectedMp) / sum(RatioM * ExpectedMp))
        PredF <- (RatioF * ExpectedFp) * (sum(ExpectedFp) / sum(RatioF * ExpectedFp))
        PredT <- (RatioT * ExpectedTp) * (sum(ExpectedTp) / sum(RatioT * ExpectedTp))
        # pdfs from these can be compared with year t+1
        OverlapNoIt <- c(
        Mal = sum(pmin(PredM / sum(PredM), Bxymex2 / sum(Bxymex2))),
        Fem = sum(pmin(PredF / sum(PredF), Bxyfex2 / sum(Bxyfex2))),
        Tot = sum(pmin(PredT / sum(PredT), BxyTex2 / sum(BxyTex2))))
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
        
        OverlapIt <- c(Mal = sum(pmin(BxymexPDF, Bxymex2 / sum(Bxymex2))),
                Fem = sum(pmin(BxyfexPDF, Bxyfex2 / sum(Bxyfex2))),
                Tot = sum(pmin(BxyTexPDF, BxyTex2 / sum(BxyTex2))))
        c(OverlapIt, OverlapNoIt)
        }, .Ex = ExUS, .dxm = dxmUS, .dxf = dxfUS, .Bxy = BxymfUS))

years <- 1969:2008
plot(years,UScompare[,1] / UScompare[,4], type = 'l', col = "blue")
lines(years,UScompare[,2] / UScompare[,5], type = 'l', col = "red")
abline(h=1)


