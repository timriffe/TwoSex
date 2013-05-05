
# let's compare the stability of age and ey pyramids over the whole HMD?
# based on the whole-HMD zipfiles of mtltper_1x1, fltper_1x1 and Population,
# these folders are unzipped and together in a folder called /HMD1x1/

#list.files("/home/triffe/git/DISS/Data/HMD1x1/fltper_1x1")
#list.files("/home/triffe/git/DISS/Data/HMD1x1/mltper_1x1")
#list.files("/home/triffe/git/DISS/Data/HMD1x1/Population")

# split off names of each population
ctries <- unlist(lapply(strsplit(list.files("/home/triffe/git/DISS/Data/HMD1x1/fltper_1x1"),split = "\\."),"[[",1))

Fdx <- do.call(rbind,lapply(ctries, function(ctry){
                    path <- paste0("/home/triffe/git/DISS/Data/HMD1x1/fltper_1x1/",ctry,".fltper_1x1.txt")
                    m <- read.table(path, header = TRUE, skip = 2, as.is = TRUE, na.strings = ".")
                    cbind(CTRY = ctry, m[,c("Year","Age","dx")])
                }))
Mdx <- do.call(rbind,lapply(ctries, function(ctry){
                    path <- paste0("/home/triffe/git/DISS/Data/HMD1x1/mltper_1x1/",ctry,".mltper_1x1.txt")
                    m <- read.table(path, header = TRUE, skip = 2, as.is = TRUE, na.strings = ".")
                    cbind(CTRY = ctry, m[,c("Year","Age","dx")])
                }))
Px <- do.call(rbind,lapply(ctries, function(ctry){
                    path <- paste0("/home/triffe/git/DISS/Data/HMD1x1/Population/",ctry,".Population.txt")
                    m <- read.table(path, header = TRUE, skip = 2, as.is = TRUE, na.strings = ".")
                    cbind(CTRY = ctry, m[,c("Year","Age","Female","Male")])
                }))
Fdx$dx <- unlist(with(Fdx,tapply(dx,paste(CTRY,Year,sep="-"),function(x){
                    x[is.na(x)] <- 0
                    x / sum(x)
                })))
Mdx$dx <- unlist(with(Mdx,tapply(dx,paste(CTRY,Year,sep="-"),function(x){
                            x[is.na(x)] <- 0
                            x / sum(x)
                        })))

rownames(Px)  <-  with(Px,paste(CTRY,Year,Age,sep="-"))
rownames(Fdx) <-  with(Fdx,paste(CTRY,Year,Age,sep="-"))
rownames(Mdx) <-  with(Mdx,paste(CTRY,Year,Age,sep="-"))

Dat <- Fdx
colnames(Dat) <- c("CTRY", "Year", "Age", "fdx")
Dat$mdx <- Mdx$dx
Dat$Pxm <- Px[rownames(Dat),"Male"]
Dat$Pxf <- Px[rownames(Dat),"Female"]

rm(Mdx,Fdx,Px); gc()

# get remaining years counts as well.
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
#Datl <- do.call(rbind,lapply(split(Dat, with(Dat,paste(CTRY,Year,sep="-"))),function(X){
#            X$Pym <- rowSums(ExpectedDx(X$Pxm, X$mdx))
#            X$Pyf <- rowSums(ExpectedDx(X$Pxf, X$fdx))
#            X
#        }))
#rownames(Datl) <- 1:nrow(Datl)
#InnerFun <- compiler::cmpfun(function(X, Njump = .Njump){
#            IDs         <- unique(X$Year)
#            if(length(IDs) > (Njump+1)){
#                DiffCoef   <- matrix(nrow = (length(IDs)-Njump), ncol = 2, dimnames = list(IDs[-(1:Njump)], c("age","ey")))
#                
#                for (i in 1:(length(IDs)-Njump)){
#                    yr1             <- X[X$Year == IDs[i],]
#                    yr2             <- X[X$Year == IDs[i+Njump],]
#                    age1            <- c(yr1$Pxm,yr1$Pxf); age2 <- c(yr2$Pxm,yr2$Pxf)
#                    ey1             <- c(yr1$Pym,yr1$Pyf); ey2 <- c(yr2$Pym,yr2$Pyf)
#                    DiffCoef[i,]   <- c(1 - sum(pmin(age1/sum(age1, na.rm = TRUE),age2/sum(age2, na.rm = TRUE), na.rm = TRUE)), 
#                            1 - sum(pmin(ey1/sum(ey1, na.rm = TRUE),ey2/sum(ey2, na.rm = TRUE), na.rm = TRUE)))
#                }
#                DiffCoef <- DiffCoef[DiffCoef[,1]!=0, ]
#                DiffCoef <- DiffCoef[DiffCoef[,2]!=0, ]
#                DiffCoef <- DiffCoef[!is.na(DiffCoef[,1]),]
#                DiffCoef <- DiffCoef[!is.na(DiffCoef[,2]),]
#                return(DiffCoef)
#            } else {
#                return(c(NA,NA))
#            }
#        })
#hm<- split(Datl, Datl$CTRY)
#
#ResultsAll <- t(sapply(1:50, function(n, .hm){
#            colMeans(do.call(rbind,parallel::mclapply(.hm, InnerFun, Njump = n,mc.cores=2)), na.rm=TRUE)
#        }, .hm = hm))
#
## results for full HMD
#save(ResultsAll, file = "/home/triffe/git/DISS/Data/HMD_Lx/Results/ResultsAll.Rdata")


# now take only years since 1950 and repeat:
#Datl2 <- Datl[Datl$Year >= 1950, ]
#hm2 <- split(Datl2, Datl2$CTRY)
#Results1950p <- t(sapply(1:50, function(n, .hm2){
#            colMeans(do.call(rbind,parallel::mclapply(.hm2, InnerFun, Njump = n,mc.cores=2)), na.rm=TRUE)
#        }, .hm2 = hm2))
#save(Results1950p, file = "/home/triffe/git/DISS/Data/HMD_Lx/Results/Results1950p.Rdata")

ResultsAll   <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/Results/ResultsAll.Rdata")))
Results1950p <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/Results/Results1950p.Rdata")))

pdf("/home/triffe/git/DISS/latex/Figures/PyramidStabilityThetaRatioAll.pdf", height = 5, width = 5)
par(mai = c(.5, .5, .6, .2),xaxs = "i", yaxs = "i")
plot(1:50, ResultsAll[,2]/Results1950p[,1], type = 'l', ylim = c(0,.8), xlim = c(0,51), axes = FALSE,
        col = gray(.2), lwd = 2, xlab = "", ylab = "",
        panel.first = list(rect(0,0,51,.8,col = gray(.95), border=NA),
                abline(h = seq(0,.8,by = .1), col = "white"),
                abline(v = seq(0, 50, by = 5), col = "white"),
                text(0, seq(0, .8, by = .1),seq(0, .8, by = .1), pos = 2, cex = .8, xpd = TRUE),
                text(seq(0, 50, by = 10),0, seq(0, 50, by = 10), pos = 1, cex = .8, xpd = TRUE),
                text(25, -.05, "lag", cex = 1, pos = 1, xpd = TRUE),
                text(0,.87, expression(frac(theta~e[y],theta~age)), cex = 1, xpd = TRUE)))

lines(1:50, Results1950p[,2]/Results1950p[,1], lwd = 2.5, col = gray(.5), lty = 5)

legend(0,.78, lty = c(1,5), col = gray(c(.2,.5)), lwd = c(2,2.5),bty = "n",
        legend = c("All years","Years >= 1950"), xpd = TRUE)
dev.off()


#All <- unique(with(Dat,paste(CTRY,Year,sep="-")))
#All <- as.data.frame(do.call(rbind,strsplit(All,split="-")),stringsAsFactors = FALSE)
#table(All[,1])
#nrow(All) - length(unique(Dat$CTRY))
#
#m50 <- table(All[All[,2]>=1950,1]) - 50
#sum(m50[m50>0])