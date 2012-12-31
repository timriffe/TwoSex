

# Author: triffe
###############################################################################

# Stable marriage problem using the Gale-Shapely algorithm:

males <- c("abe", "bob", "col", "dan", "ed", "fred", "gav", "hal", "ian", "jon")
females <-  c("abi", "bea", "cath", "dee", "eve", "fay", "gay", "hope", "ivy", "jan")

m.prefs <- as.matrix(data.frame(abe = c("abi", "eve", "cath", "ivy", "jan", "dee", "fay", "bea", "hope", "gay"),
bob = c("cath", "hope", "abi", "dee", "eve", "fay", "bea", "jan", "ivy", "gay"),
col = c("hope", "eve", "abi", "dee", "bea", "fay", "ivy", "gay", "cath", "jan"),
dan = c("ivy", "fay", "dee", "gay", "hope", "eve", "jan", "bea", "cath", "abi"),
ed = c("jan", "dee", "bea", "cath", "fay", "eve", "abi", "ivy", "hope", "gay"),
fred = c("bea", "abi", "dee", "gay", "eve", "ivy", "cath", "jan", "hope", "fay"),
gav = c("gay", "eve", "ivy", "bea", "cath", "abi", "dee", "hope", "jan", "fay"),
hal = c("abi", "eve", "hope", "fay", "ivy", "cath", "jan", "bea", "gay", "dee"),
ian = c("hope", "cath", "dee", "gay", "bea", "abi", "fay", "ivy", "jan", "eve"),
jon = c("abi", "fay", "jan", "gay", "eve", "bea", "dee", "cath", "ivy", "hope")))

                
f.prefs <- as.matrix(data.frame(abi = c("bob", "fred", "jon", "gav", "ian", "abe", "dan", "ed", "col", "hal"),
bea = c("bob", "abe(", "col", "fred", "gav", "dan", "ian", "ed", "jon", "hal"),
cath = c("fred", "bob", "ed", "gav", "hal", "col", "ian", "abe", "dan", "jon"),
dee = c("fred", "jon", "col", "abe", "ian", "hal", "gav", "dan", "bob", "ed"),
eve = c("jon", "hal", "fred", "dan", "abe", "gav", "col", "ed", "ian", "bob"),
fay = c("bob", "abe", "ed", "ian", "jon", "dan", "fred", "gav", "col", "hal"),
gay = c("jon", "gav", "hal", "fred", "bob", "abe", "col", "ed", "dan", "ian"),
hope = c("gav", "jon", "bob", "abe", "ian", "dan", "hal", "ed", "col", "fred"),
ivy = c("ian", "col", "hal", "gav", "fred", "bob", "abe", "ed", "jon", "dan"),
jan = c("ed", "hal", "gav", "abe", "bob", "jon", "col", "ian", "fred", "dan")))


MatchMat <- matrix(0, ncol = length(males), nrow = length(females), dimnames = list(females, males))
# 
m.prefs.ind <- matrix(TRUE,ncol = length(males), nrow = length(females), dimnames = list(1:length(females), males))
for (m in 1:length(females)){
    proposers           <- males[colSums(MatchMat) == 0]
    r.ind               <- proposals <- vector(length = length(proposers))
    for (i in 1:length(proposers)) {
        proposals[i]    <- m.prefs[m.prefs.ind[, proposers[i]], proposers[i]][1]
        r.ind[i]        <-  which(m.prefs[, proposers[i]] == proposals[i])  
    }  
    m.prefs.ind[cbind(r.ind,proposers)]     <- FALSE
    MatchMat[cbind(proposals,proposers)]    <- 1
    if (any(rowSums(MatchMat) > 1)){
        deciders        <- females[rowSums(MatchMat) > 1]
        for (i in deciders){
            choice      <- f.prefs[,i][ f.prefs[,i] %in% males[as.logical(MatchMat[i, ])] ][1]
            MatchMat[i, ]       <- 0
            MatchMat[i, choice] <- 1
        }
    }
    if (all(colSums(MatchMat) == 1)){
        break
    }  
}
MatchMat














