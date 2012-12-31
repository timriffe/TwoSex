# Author: triffe
###############################################################################

# this function is pretty darn simple. Given: 
# rate matrices 1 and 2
# exposure matrices 1 and 2
# how much of the change in total Births from time 1 to time 2 is due to:
# 1) rates
# 2) exposures
# these will be single output matrices
# the sums of each should add up to the difference in total number of births.
# matrices can be plotted as image()s and time series of the two components can
# be calculated, plotted.

# the idea is to recalculate Fx given the many new mean exposure matrices
# and Ex come from one of the many possible means.
DecomposeBirthsChange <- compiler:::cmpfun(function(Fx1, Fx2, Ex1, Ex2, return.all = FALSE){
    RateComponent       <- (Fx2 - Fx1) * (Ex2 + Ex1) / 2
    ExposureComponent   <- (Fx2 + Fx1) / 2 * (Ex2 - Ex1)
    if (return.all){
        return(list(FxCmat = RateComponent,           # check out plots of these matrices to see how they differ
                ExCmat = ExposureComponent, 
                FxC = sum(RateComponent, na.rm = TRUE),      # these two should add to the change in B
                ExC = sum(ExposureComponent, na.rm = TRUE)))
    } 
    out <- c(sum(RateComponent, na.rm = TRUE), sum(ExposureComponent, na.rm = TRUE))
    names(out) <- c("FxC","ExC")
    out
})

# try out on a given 'p', two consecutive years:
Exyp    <- local(get(load("/home/triffe/git/Dissertation/DISSERTATION/DATA/Exposures/ESexpMeans1x1.Rdata")))
Fxyp    <- local(get(load("/home/triffe/git/Dissertation/DISSERTATION/DATA/Rates.p/ES_Fxyp.Rdata")))
# just loaded in like 1.5 Gb of stuff!
print(object.size(Fxyp),units="Mb");print(object.size(Exyp),units="Mb") # ouch

years  <- 1975:2009
ExCmat <- FxCmat <- matrix(0,ncol = (length(years) - 1), nrow = 199, 
        dimnames = list(dimnames(Exyp[[1]])[[3]], paste(1975:2008, 1976:2009, sep = "-")))

for (n in 1:(length(years) - 1)){
    for (p in 1:199){
        comp.p <- DecomposeBirthsChange(Fxyp[[n]][, , p], Fxyp[[n + 1]][, , p], 
                                        Exyp[[n]][, , p], Exyp[[n + 1]][, , p])
        FxCmat[p, n]    <- comp.p[1]
        ExCmat[p, n]    <- comp.p[2]
    }
}

rm(Exyp);rm(Fxyp);gc()

#image(FxCmat)
#x <- as.numeric(rownames(FxCmat))
#p <- log(x/(1-x))
#plot(p, FxCmat[, 2], type = 'l', ylim = c(-50000,10000))
#lines(p, ExCmat[,2])

#Bchange <- FxCmat + ExCmat
#hist(FxCmat / Bchange)
#hist(ExCmat / Bchange)
#
#which.min(rowSums(FxCmat))
#which.max(rowSums(FxCmat))
#
#which.min(rowSums(ExCmat))
#which.max(rowSums(ExCmat))

save(FxCmat, file = "/home/triffe/git/Dissertation/DISSERTATION/DATA/results/decomposition1/FxCmat.Rdata")
save(ExCmat, file = "/home/triffe/git/Dissertation/DISSERTATION/DATA/results/decomposition1/ExCmat.Rdata")

# ----------------------------------------------------------------------------
# the above doesn't work so well. ummm. how about rescaling Births to total properly
# and decomposing the rate change?
# load("/home/triffe/git/Dissertation/DISSERTATION/DATA/results/decomposition1/FxCmat.Rdata")

