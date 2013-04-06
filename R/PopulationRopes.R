
# Author: triffe
###############################################################################
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
# getting some fresh Swedish data. num num
#library(HMDget)
#Px <- HMDget(countries = c("SWE"), wanteditems = c("Population"), 
#        years = list(Population = c(1751:2011)), 
#        drop.tadj = TRUE, format = 0, username = "", password = "")
#colnames(Px) <- c("Year","Age","Female","Male","Total")
#dxm <- HMDget(countries = c("SWE"), wanteditems = c("mltper_1x1"), 
#        years = list(mltper_1x1 = c(1751:2011)), 
#        drop.tadj = TRUE, column = "dx", format = 4, username = "", password = "")
#dxf <- HMDget(countries = c("SWE"), wanteditems = c("fltper_1x1"), 
#        years = list(fltper_1x1 = c(1751:2011)), 
#        drop.tadj = TRUE, column = "dx", format = 4, username = "", password = "")
#save(DT, file = "/home/triffe/git/DISS/Data/SWE/SWEtri.Rdata")
#save(Px, file = "/home/triffe/git/DISS/Data/SWE/Px.Rdata")
#save(dxm, file = "/home/triffe/git/DISS/Data/SWE/dxm.Rdata")
#save(dxf, file = "/home/triffe/git/DISS/Data/SWE/dxf.Rdata")
DT  <- local(get(load("/home/triffe/git/DISS/Data/SWE/SWEtri.Rdata")))
Px  <- local(get(load("/home/triffe/git/DISS/Data/SWE/Px.Rdata")))
dxm <- local(get(load("/home/triffe/git/DISS/Data/SWE/dxm.Rdata")))
dxf <- loca:(get(load("/home/triffe/git/DISS/Data/SWE/dxf.Rdata")))

# get pop vectors to redistribute:
Males   <- Px$Male[Px$Year == 2011]
Females <- Px$Female[Px$Year == 2011]

MxM <- ExpectedDx(Males, dxm[, "2011"] )[,111:1]
FxM <- ExpectedDx(Females, dxf[, "2011"] )[,111:1]

#DT <- read.table("/home/triffe/workspace/CODE_ALL/SWEtri.txt", na.strings = ".",
#        header = FALSE, as.is = TRUE, skip = 3,col.names = c("Year","Age","Cohort","Female","Male","Total"))
DT$Age              <- as.integer(sub(DT$Age, pattern = "\\+", replacement = ""))
NAind               <- which(is.na(DT$Cohort))
DT$Cohort[NAind]    <- DT$Cohort[(NAind - 1)]


Males <- reshape2::acast(DT, Year~Cohort, sum,value.var = "Male")
Females <- reshape2::acast(DT, Year~Cohort, sum,value.var = "Female")
RCMales <- t(apply(Males,1,cumsum))
RCFemales <- t(apply(Females,1,cumsum))
# Years in rows

#dev.new(height = 8, width = 5)
pdf("/home/triffe/git/DISS/latex/Figures/SWE_then_and_now.pdf",height = 8, width = 5)
par(mai = c(.3,.3,.3,.3), xaxs = "i", yaxs = "i", xpd = TRUE)
plot(NULL, type = "n", xlim = c(-60000, 60000), ylim = c(1750, 2120), axes = FALSE)
# a year slice:
# RCMales[10,]

# colors by 10-year cohorts?
Cohorts             <- sort(unique(DT$Cohort))
Years               <- sort(unique(DT$Year))
Coh10               <- Cohorts - Cohorts %% 20
grays               <- gray(seq(from=.9,to=.1,length.out = length(unique(Coh10))))


bottoms             <- row(RCMales) + 1750
tops                <- bottoms + 1
cols                <- matrix(rep(grays, table(Coh10)), 
                            ncol = length(Cohorts), nrow = length(Years), byrow=TRUE)
# draw lower rope
rect(-cbind(0,RCMales[,-ncol(RCMales)]),bottoms,-RCMales,tops, border = cols, col = cols, lwd = .2)
rect(cbind(0,RCFemales[,-ncol(RCFemales)]),bottoms,RCFemales,tops, border = cols, col = cols, lwd = .2)

# now for dx pyramid on top:
Cohorts         <- 2012 + 0:110 

Coh10           <- Cohorts - Cohorts %% 20
purples         <- rev(RColorBrewer::brewer.pal(9, "Purples")[3:9])
cols            <- matrix(rep(purples, table(Coh10)), 
                    ncol = 111, nrow = 111, byrow=TRUE)

MxMat           <- t(apply(MxM, 1, cumsum))
FxMat           <- t(apply(FxM, 1, cumsum))
bottoms         <- row(MxMat) + 2011
tops            <- bottoms + 1
rect(-cbind(0,MxMat[,-ncol(MxMat)]),bottoms,-MxMat,tops, border = cols, col = cols, lwd = .2)
rect(cbind(0,FxMat[,-ncol(FxMat)]),bottoms,FxMat,tops, border = cols, col = cols, lwd = .2)
dev.off()

# -------------------------------------------------------------------------
# oops, data objects were cut down, need another full copy...
# do Spain:
DT  <- local(get(load("/home/triffe/git/DISS/Data/HMD_TLTU/DTLTUES.Rdata")))
Px  <- local(get(load("/home/triffe/git/DISS/Data/HMD_Px/PxES.Rdata")))
dxm <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxmES.Rdata")))
dxf <- local(get(load("/home/triffe/git/DISS/Data/HMD_dx/dxfES.Rdata")))
DT$Age              <- as.integer(sub(DT$Age, pattern = "\\+", replacement = ""))
NAind               <- which(is.na(DT$Cohort))
DT$Cohort[NAind]    <- DT$Cohort[(NAind - 1)]

# get pop vectors to redistribute:
Males   <- Px$Male[Px$Year == 2009]
Females <- Px$Female[Px$Year == 2009]

MxM <- ExpectedDx(Males, dxm[, "2009"] )[,111:1]
FxM <- ExpectedDx(Females, dxf[, "2009"] )[,111:1]


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
Cohorts         <- 2010 + 0:110 

Coh10           <- Cohorts - Cohorts %% 20
purples         <- rev(RColorBrewer::brewer.pal(9, "Purples")[3:9])
cols            <- matrix(rep(purples, table(Coh10)), 
        ncol = 111, nrow = 111, byrow=TRUE)

MxMat           <- t(apply(MxM, 1, cumsum))
FxMat           <- t(apply(FxM, 1, cumsum))
bottoms         <- row(MxMat) + 2010
tops            <- bottoms + 1
rect(-cbind(0,MxMat[,-ncol(MxMat)]),bottoms,-MxMat,tops, border = cols, col = cols, lwd = .2)
rect(cbind(0,FxMat[,-ncol(FxMat)]),bottoms,FxMat,tops, border = cols, col = cols, lwd = .2)
dev.off()







