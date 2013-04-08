
# Author: triffe
###############################################################################
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
yearsUS <- 1969:2009

BxymfUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxymf0_110.Rdata")))
names(BxymfUS)      <- yearsUS
DT                  <- local(get(load("/home/triffe/git/DISS/Data/RopeData/DTLTUUS.Rdata")))
Px                  <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxUS.Rdata")))
Ex                  <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))

dxm                 <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmUS.Rdata")))
dxf                 <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfUS.Rdata")))
dxm <- dxm %col% colSums(dxm)
dxf <- dxf %col% colSums(dxf)

# mx to acct for mort improvement:
mxm                 <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxmUS.Rdata")))[,"2009"]
mxf                 <- local(get(load("/home/triffe/git/DISS/Data/HMD_mux/muxfUS.Rdata")))[,"2009"]

# 2009 US
.dxm <- dxm; .dxf <- dxf; .Bxymf <- BxymfUS; .Ex <- Ex
yr <- "2009"
yri     <- 2009
.dxm.   <- .dxm[, yr]
.dxf.   <- .dxf[, yr]
ExM     <- rowSums(ExpectedDx( with(.Ex, Male[Year == yri]), .dxm.))
BxMM    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxm.))
BxMF    <- rowSums(ExpectedDx( rowSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxm.))

ExF     <- rowSums(ExpectedDx( with(.Ex, Female[Year == yri]), .dxf.))
BxFF    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxyf"]], na.rm = TRUE), .dxf.))
BxFM    <- rowSums(ExpectedDx( colSums(.Bxymf[[yr]][["Bxym"]], na.rm = TRUE), .dxf.))

# sex-sex-ex- specific rates:
FxMM    <- Mna0(Minf0(BxMM / ExM))
FxFF    <- Mna0(Minf0(BxFF / ExF))
FxMF    <- Mna0(Minf0(BxMF / ExM))
FxFM    <- Mna0(Minf0(BxFM / ExF))

#
Males   <- Px$Male[Px$Year == 2009]
Females <- Px$Female[Px$Year == 2009]

# redistribute by future deaths:
MxMp     <- MxM     <- ExpectedDx(Males, dxm[, "2009"] )
MxFp     <- MxF     <- ExpectedDx(Females, dxf[, "2009"] )
#tail(MxM)
sigma   <- .5
N <- 300
ProjMatM <- ProjMatF <- matrix(nrow = N, ncol = 111)

K <- .03
facF <- .003
FacM <- .01

r <- .05
Mdecs <- facs <- vector(length = N)
Mdec <- 0


for (i in 1:N){
    Mdec <- 1 - ((1 - Mdec) * (1 - FacM / (1 + FacM * 1.5 * i)))
    facs[i] <- facF
    Bxm         <- sum(sigma * rowSums(MxMp[i:(110+i),1:111]) * FxMM * (facF + 1)+ 
                    (1 - sigma) * rowSums(MxFp[i:(110+i),1:111]) * FxFM * (facF + 1))
    Bxf         <- sum(sigma * rowSums(MxMp[i:(110+i),1:111]) * FxMF * (facF + 1)+ 
                    (1 - sigma) * rowSums(MxFp[i:(110+i),1:111]) * FxFF * (facF + 1))
    
    Bym         <- Bxm * mx2dxHMD(mxm * (1 - Mdec))
    Bym         <- c(rep(0,i), Bym)
    #Bym[1:2]    <- mean(Bym[1:2])
    MxMp        <- cbind(Bym, rbind(MxMp,0))
    Byf         <- Bxf *  mx2dxHMD(mxf * (1 - Mdec))
    Byf         <- c(rep(0,i), Byf)
    #Byf[1:2]    <- mean(Byf[1:2])
    MxFp        <- cbind(Byf, rbind(MxFp,0))
    
    # logistic increase to fert
    facs[i] <- facF
    Mdecs[i] <- Mdec
    facF <- facF + r * facF * (1 - facF / K)
    # 
} 

#graphics.off()
#plot(1:N,Mdecs,type = 'l')
CohsNew         <- ncol(MxMp):1 + 1897
colnames(MxMp)  <- CohsNew
colnames(MxFp)  <- CohsNew
MxMpC           <- ReduceDimension(MxMp,CohsNew,20,1)
MxFpC           <- ReduceDimension(MxFp,CohsNew,20,1)
Coh20           <- CohsNew - CohsNew %% 20
Coh20N          <- unique(Coh20)
ProjcolFun      <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "BuPu"),space = "Lab")
Projcol         <- rev(ProjcolFun(45)[1:22])
ProjcolM        <- matrix(Projcol, ncol = length(Coh20N), nrow = nrow(MxMpC), byrow = TRUE)

# now prepare presently living population
CohortsL        <- 2008 - 0:110 
MxMn            <- ReduceDimension(MxM,CohortsL,20,1)
FxMn            <- ReduceDimension(MxF,CohortsL,20,1)
CohL            <- CohortsL - CohortsL %% 20
CohLN           <- unique(CohL)
Livcol          <- RColorBrewer::brewer.pal(9, "Greens")[3:9] # must be same length as CohLN
LivcolM          <- matrix(Livcol, ncol = length(unique(Coh2)), nrow = 111, byrow=TRUE)


# get cumsums to stack:
LivCM           <- t(apply(MxMn, 1, cumsum))
LivCF           <- t(apply(FxMn, 1, cumsum))
ProjCM          <- t(apply(MxMpC, 1, cumsum))
ProjCF          <- t(apply(MxFpC, 1, cumsum))

# get historical ready
DT$Age              <- as.integer(sub(DT$Age, pattern = "\\+", replacement = ""))
NAind               <- which(is.na(DT$Cohort))
DT$Cohort[NAind]    <- DT$Cohort[(NAind - 1)]
MalesD              <- reshape2::acast(DT, Year~Cohort, sum,value.var = "Male")
FemalesD            <- reshape2::acast(DT, Year~Cohort, sum,value.var = "Female")

CohortsH            <- sort(unique(DT$Cohort))
MxMn                <- ReduceDimension(MalesD,CohortsH,20,1)
FxMn                <- ReduceDimension(FemalesD,CohortsH,20,1)
RCMales             <- t(apply(MxMn, 1, cumsum))
RCFemales           <- t(apply(FxMn, 1, cumsum))

YearsH              <- sort(unique(DT$Year))
Coh20H              <- CohortsH - CohortsH %% 20
Coh20HN             <- unique(Coh20H)
Hcol                <- gray(seq(from = .9, to = .1,length.out = length(Coh20HN)))
HcolM               <- matrix(Hcol, ncol = length(Coh20HN), nrow = length(YearsH), byrow = TRUE)

# get some coords
bottomsL        <- row(LivCM) + 2008
topsL           <- bottomsL + 1
bottomsP        <- row(ProjCM) + 2008
topsP           <- bottomsP + 1
bottomsH        <- row(RCMales) + 1931
topsH           <- bottomsH + 1
# Years in rows

#dev.new(height = 8, width = 3)
pdf("/home/triffe/git/DISS/latex/Figures/US_DxHist.pdf",height = 8, width = 3)
par(mai = c(.3,.3,.3,.3), xaxs = "i", yaxs = "i", xpd = TRUE)
plot(NULL, type = "n", xlim = c(-3000000, 3000000), ylim = c(1900, 2300), axes = FALSE)
# a year slice:

# draw Lower historical
rect(-cbind(0,RCMales[,-ncol(RCMales)]),bottomsH,-RCMales,topsH, col = HcolM, border = NA)
rect(cbind(0,RCFemales[,-ncol(RCFemales)]),bottomsH,RCFemales,topsH, col = HcolM, border = NA)

# draw Upper projected 
rect(-cbind(0,ProjCM[,-ncol(ProjCM)]),bottomsP,-ProjCM,topsP, col = ProjcolM, border = NA)
rect(cbind(0,ProjCF[,-ncol(ProjCF)]),bottomsP,ProjCF,topsP, col = ProjcolM, border = NA)

# draw Middle Living projected 
rect(-cbind(0,LivCM[,-ncol(LivCM)]),bottomsL,-LivCM,topsL, col = LivcolM, border = NA)
rect(cbind(0,LivCF[,-ncol(LivCF)]),bottomsL,LivCF,topsL, col = LivcolM, border = NA)

# now for dx pyramid on top:
dev.off()







