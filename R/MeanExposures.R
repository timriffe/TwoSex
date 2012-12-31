
# Author: triffe
###############################################################################

# this script simply calculates a given mean exposure. This should be run by hand, as it's non trivial.
# hence the faile brake at the top here: don't want to source it by accident!

# within the script, change N to e.g. 1, 2 or 3, and be sure to change the resultant file name.
# to end in e.g. 1x1, 2x2 or 3x3
run.this.script <- FALSE
if (run.this.script){


source("/home/triffe/git/Dissertation/DISSERTATION/R/MeanFunctions.R")
source("/home/triffe/git/Dissertation/DISSERTATION/R/UtilityFunctions.R")
# SolarskyM() the function of interest
N       <- 3 # change this as necessary

E       <- local(get(load("/home/triffe/git/Dissertation/DISSERTATION/DATA/Exposures/ESexp.Rdata")))
years   <- sort(as.integer(unique(E$Year)))

# these are the different 'p' parameters to iterate over!
delta   <- .005
x       <- seq(delta, 1 - delta, by = delta)
p       <- logit(x) # this way we get more resolution in the steep part of the curve :-)

# preallocate MeanCubeList
Em              <- makeVectorAgeGroups(with(E, Male[Year == years[1]]), ages = 0:110, N = N)
MeanCubeList    <- lapply(years, function(y){
                        array(dim = c(length(Em), length(Em), length(p)), dimnames = list(names(Em), names(Em), x))
                      }
                    )
# note that it's big!
print(object.size(MeanCubeList),units = "Mb")

# begin non-trivial loop
for (i in 1:length(years)){
    Em          <- makeVectorAgeGroups(with(E, Male[Year == years[i]]), ages = 0:110, N = N)
    Ef          <- makeVectorAgeGroups(with(E, Female[Year == years[i]]), ages = 0:110, N = N)
    nages       <- length(Em) 
    Em          <- matrix(Em, nrow = nages, ncol = nages) # males in rows = repeat over columns
    Ef          <- matrix(Ef, nrow = nages, ncol = nages, byrow = TRUE) # repeat over rows
    
    MeanCube    <- array(dim = c(dim(Em), length(p)), dimnames = list(rownames(Em), colnames(Em), x))
    for (z in 1:length(x)){
        MeanCube[, , z]   <- StolarskyM(Em, Ef, p[z])
    }
    MeanCubeList[[i]] <- MeanCube
    cat("\n",years[i],"\n")
}
# did it get bigger or did we preallocate OK?
print(object.size(MeanCubeList), units = "Mb")

# just recall to logit the value x that labels the z levels of the cubes!
#save(MeanCubeList, file = "/home/triffe/git/Dissertation/DISSERTATION/DATA/Exposures/ESexpMeans1x1.Rdata")
#save(MeanCubeList, file = "/home/triffe/git/Dissertation/DISSERTATION/DATA/Exposures/ESexpMeans2x2.Rdata")
save(MeanCubeList, file = "/home/triffe/git/Dissertation/DISSERTATION/DATA/Exposures/ESexpMeans3x3.Rdata")

rm(MeanCubeList); gc()

} # bracket for end of the safety brake 