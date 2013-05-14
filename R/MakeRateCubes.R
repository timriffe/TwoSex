# Code from early diagnostics
# Author: triffe
###############################################################################
# a non-trivial, memory-intensive calculation, try not to have to repeat this!
run.this.file <- FALSE
if (run.this.file){
    # 1) load Exposure Cubes (i.e. FAGE~MAGE Stolarsky mean exposures by x (expit of p)- a list of 
    # 3d arrays, very large object!
    E   <- local(get(load("DATA/Exposures/ESexpMeans1x1.Rdata")))
    # 2) load in Bxy
    Bxy <- local(get(load("DATA/ESbirths/ESBxy.Rdata")))
    # 3) create copy of E, to pre-allocate:
    Fxyp <- E
    
    years <- names(Bxy)
    for (i in 1:length(years)){
        out <- apply(E[[i]], 3, function(e.yr.p, .bxy){
                    .bxy / e.yr.p
                }, .bxy = Bxy[[i]])
        dim(out) <- c(111,111,199)
        Fxyp[[i]] <- out
    }
    names(Fxyp) <- 1975:2009
    save(Fxyp, file = "DATA/Rates.p/ES_Fxyp.Rdata")
}

