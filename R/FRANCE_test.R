# France could have been included but was not for lack of time and because fewer years were available
# these data have not ben fully explored

# Author: triffe
###############################################################################
setwd("/home/triffe/git/DISS/")
path <- "DATA/FRA/BIRTHS/Zipped"
zipped <- list.files(path)
for (f in 1:length(zipped)){
    unzip(file.path(path,zipped[f]), exdir = "DATA/FRA/BIRTHS/UNZIPPED")
}

path <- "DATA/FRA/BIRTHS"

file.names <- paste0("NAIS",1998:2009,".dbf")

library(foreign)
for (i in 1:length(file.names)){
    yr          <- gsub(gsub(file.names[i], pattern = "NAIS", replacement = ""), pattern = ".dbf", replacement = "")
    DAT         <- read.dbf(file.path(path,"UNZIPPED",file.names[i]), as.is = TRUE)
    save(DAT, file = file.path(path, "Rbin", paste0("B", yr, ".Rdata")))
    Bxy         <- with(DAT, table(AGEXACTP,AGEXACTM)) 
    save(Bxy, file = file.path(path, "Rbin", paste0("Bxy", yr, ".Rdata")))
}
# exact age is completed age ; age is age to be completed if year survived. wow use completed for now.
#mat <- reshape2::acast(DAT, AGEXACTP~AGEXACTM)
#with(DAT, table(AGEXACTP,AGEXACTM))
#image(x=17:46+.5,y=17:46+.5,Bxy,xlim=c(17,47),ylim=c(17,48))
#sort(unique(DAT$AGEMERE))
path.rm <- "DATA/FRA/BIRTHS/UNZIPPED"
all.files <- list.files(path.rm)
files.rm <- all.files[!all.files == "varlist_naissances.dbf"]
sapply(files.rm, function(x){
            file.remove(file.path(path.rm,x))
        })
