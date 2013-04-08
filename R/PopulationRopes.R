
# Author: triffe
###############################################################################
source("/home/triffe/git/DISS/R/UtilityFunctions.R")

DT                  <- local(get(load("/home/triffe/git/DISS/Data/RopeData/DTLTUES.Rdata")))
Px                  <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxES.Rdata")))
dxm                 <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmES.Rdata")))
dxf                 <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfES.Rdata")))
DT$Age              <- as.integer(sub(DT$Age, pattern = "\\+", replacement = ""))
NAind               <- which(is.na(DT$Cohort))
DT$Cohort[NAind]    <- DT$Cohort[(NAind - 1)]

# get pop vectors to redistribute:

Males   <- Px$Male[Px$Year == 2009]
Females <- Px$Female[Px$Year == 2009]

# redistribute by future deaths:
MxM     <- ExpectedDx(Males, dxm[, "2009"] )[,111:1]
FxM     <- ExpectedDx(Females, dxf[, "2009"] )[,111:1]
# combine cohorts, summing:
Cohorts2        <- 2008 - 0:110 
MxMn            <- ReduceDimension(MxM,Cohorts2,20,1)
FxMn            <- ReduceDimension(FxM,Cohorts2,20,1)

Coh2            <- Cohorts2 - Cohorts2 %% 20
purples         <- RColorBrewer::brewer.pal(9, "Greens")[3:9]
cols            <- matrix(purples, 
        ncol = length(unique(Coh2)), nrow = 111, byrow=TRUE)
# get cumsums to stack:
MxMat           <- t(apply(MxMn, 1, cumsum))
FxMat           <- t(apply(FxMn, 1, cumsum))
bottoms2         <- row(MxMat) + 2011
tops2            <- bottoms2 + 1




DT$Age              <- as.integer(sub(DT$Age, pattern = "\\+", replacement = ""))
NAind               <- which(is.na(DT$Cohort))
DT$Cohort[NAind]    <- DT$Cohort[(NAind - 1)]


Males <- reshape2::acast(DT, Year~Cohort, sum,value.var = "Male")
Females <- reshape2::acast(DT, Year~Cohort, sum,value.var = "Female")
RCMales <- t(apply(Males,1,cumsum))
RCFemales <- t(apply(Females,1,cumsum))
# Years in rows

#dev.new(height = 8, width = 5)
pdf("/home/triffe/git/DISS/latex/Figures/ES_then_and_now.pdf",height = 8, width = 5)
par(mai = c(.3,.3,.3,.3), xaxs = "i", yaxs = "i", xpd = TRUE)
plot(NULL, type = "n", xlim = c(-200000, 200000), ylim = c(1750, 2120), axes = FALSE)
# a year slice:
# RCMales[10,]

# colors by 10-year cohorts?
Cohorts             <- sort(unique(DT$Cohort))
Years               <- sort(unique(DT$Year))
Coh10               <- Cohorts - Cohorts %% 20
grays               <- gray(seq(from=.9,to=.1,length.out = length(unique(Coh10))))


bottoms             <- row(RCMales) + 1864
tops                <- bottoms + 1
cols                <- matrix(rep(grays, table(Coh10)), 
        ncol = length(Cohorts), nrow = length(Years), byrow=TRUE)
# draw lower rope
rect(-cbind(0,RCMales[,-ncol(RCMales)]),bottoms,-RCMales,tops, border = cols, col = cols, lwd = .2)
rect(cbind(0,RCFemales[,-ncol(RCFemales)]),bottoms,RCFemales,tops, border = cols, col = cols, lwd = .2)

# now for dx pyramid on top:
Cohorts2         <- 2010 + 0:110 
MxMn <- ReduceDimension(MxM,Cohorts2,20,1)
FxMn <- ReduceDimension(FxM,Cohorts2,20,1)

Coh2            <- Cohorts2 - Cohorts2 %% 20
purples         <- rev(RColorBrewer::brewer.pal(9, "Purples")[3:9])
cols            <- matrix(purples, 
                       ncol = length(unique(Coh2)), nrow = 111, byrow=TRUE)

MxMat           <- t(apply(MxMn, 1, cumsum))
FxMat           <- t(apply(FxMn, 1, cumsum))
bottoms         <- row(MxMat) + 2011
tops            <- bottoms + 1
rect(-cbind(0,MxMat[,-ncol(MxMat)]),bottoms,-MxMat,tops, border = cols, col = cols, lwd = .2)
rect(cbind(0,FxMat[,-ncol(FxMat)]),bottoms,FxMat,tops, border = cols, col = cols, lwd = .2)
dev.off()







