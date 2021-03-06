
# Author: triffe
###############################################################################
# Also this script was run only once.
# this is not fully reproducible unless you have the huge text files from the CDC and INE
# on disk. This is legacy one-time code. Instead of reproducing the tabulations
# use the data outputs in /DATA/



# those pesky US natality files deserve their own massive script...
# Might have to be creative here...

# ideally, only want male age x female age cross-tables.

# I wonder if sqldf can work from inside a zip file...

library(sqldf)

# put followinb in for loop
zipdir <- tempfile()
dir.create(zipdir)
unzip("/home/triffe/DATA/CDC/BIRTHS/ZIPPED/B1984.zip", exdir=zipdir)
bigfile <- file.path(zipdir, list.files(zipdir)[1])
f <- file(bigfile)
# mother age by father age, inluding NAs?
B.i <- sqldf("select substr(V1, 41, 2) MAGE, substr(V1, 69, 2) FAGE, substr(V1, 35, 1) Sex, sum(substr(V1, 208, 1)) Births 
                            from f group by substr(V1, 41, 2), substr(V1, 69, 2), substr(V1, 35, 1)", 
                    dbname = tempfile())
dim(B.i)
head(B.i)

# clean up file system
close(con)
file.remove(bigfile)
file.remove(zipdir)



widths <- unlist(read.table("/home/triffe/DATA/CDC/BIRTHS/FWF_WIDE/Width2005.txt"))
cnames <- unlist(read.table("/home/triffe/DATA/CDC/BIRTHS/FWF_WIDE/Width2005names.txt",stringsAsFactors = FALSE))

# MAGE, FAGE, SEX
# 1968 no FAGE...
# 1969, 1970, 1971 no weight: multiply tables by 2; just do count of year?
cp <- which(cnames == "SEX")
widths[cp]
sum(widths[1:cp])-widths[cp]+1 # start pos
# count(substr(V1,1,1)) * 2 Births

# 1972-1981:
# FAGE, MAGE, SEX, RECWT
cp <- which(cnames == "RECWT")
widths[cp]
sum(widths[1:cp])-widths[cp]+1 
# 1982-1988
# age_f, age_m, sex, wt

cp <- which(cnames == "age_f")
widths[cp]
sum(widths[1:cp])-widths[cp]+1 
# 1989-2002
# DFAGE, DMAGE, CSEX, RECWT

cp <- which(cnames == "CSEX")
widths[cp]
sum(widths[1:cp])-widths[cp]+1 
# 2003
# FAGERPT or UFAGECOMB; MAGER41, SEX, RECWT
c("UFAGECOMB","MAGER41","SEX","RECWT") %in% cnames
cp <- which(cnames == "RECWT")
widths[cp]
sum(widths[1:cp])-widths[cp]+1 
# 2004-2005
# FAGERPT or UFAGECOMB; MAGER41, SEX, RECWT
c("UFAGECOMB","MAGER","SEX","RECWT") %in% cnames
cp <- which(cnames == "UFAGECOMB")
widths[cp]
sum(widths[1:cp])-widths[cp]+1 
# need to make width files, names for years 2006,7,8,9,10...
# or for speed, just get the positions and widths for the relevant columns

# 2006-2010
# FAGE: 184,2 (99 UNK, as always)
# MAGE: 89,2
# SEX:  436, 1
# no weight (100%), do count of year 15,4

# ---------------------------------------------------------------
# IMPERATIVE CLEANING
# loop to clean beforehand:
#setwd("/home/triffe/DATA/CDC/BIRTHS/ZIPPED")
#new.download <- FALSE
#if (new.download){
#    file.names <- list.files()
#    years      <- 1969:2010
#    for (i in 1:length(file.names)){
#        # i<-5
#        zip.path    <- file.names[i]
#        orig.name   <- as.character(unzip(zip.path, list = TRUE)[1, 1]) # absurd syntax...
#        unzip(zip.path)
#        # basically, keep all ascii chars, replace others with space
#        # system(paste0("perl -i -pe 's/[^[:ascii;]]/ /g' ", orig.name)) 
#        system(paste0("perl -i -pe 's/\\\\/ /g' ", orig.name)) 
#        system(paste0("perl -i -pe 's/\ç/ /g' ", orig.name)) 
#        system(paste0("perl -i -pe 's/,/ /g' ", orig.name)) 
#        new.name    <- paste0("B",years[i],".txt")
#        file.rename(from = orig.name, to = new.name)
#        file.remove(zip.path)
#        gzip(new.name) # gzip, since direct connections can be made
#    }
#}

#setwd("/home/triffe/DATA/CDC/BIRTHS/ZIPPED")
#gz.files<- list.files()
#for (gz in gz.files){
#    gunzip(gz)
#}
# check encodings: note the 5th (1973) isn't ascii...
#txt.files<- list.files()
#encodings <- vector(length = length(txt.files))
#for (txt in txt.files){
#    encodings[txt] <- system(paste0("file -bi ",txt))
#}
## convert, changing ç to c, as far as i know...
system("iconv -f us-ascii -t us-ascii//TRANSLIT -o /DATA/CDC/BIRTHS/B1979.txt")
## now remove backslash (check all)
#system("perl -i -pe 's/\\\\/ /g' B1973.txt" ) 

# CRLF is the key part
#for (i in 1:nrow(postab)){
#    fname <- paste0("B",postab$Year[i],".txt")
#    cat(fname, "\n")
#    system(paste0("file ", fname))
#}



make.tables.from.gz <- FALSE
if (make.tables.from.gz){
    setwd("/home/triffe/DATA/CDC/BIRTHS/ZIPPED")
    library(sqldf)
    library(utils)
    library(R.utils)
    postab      <- read.table("/home/triffe/git/Dissertation/DISSERTATION/R/UStabwidths.txt", 
                      header = TRUE, stringsAsFactors = FALSE)
    file.names  <- list.files()
    
    for (i in 34){
        # i <- 36
        # file.format=list(eol = "\r\n") for years 2003-? 
        # trying "\n" for 2004... (file also moved to different folder
        # complicated sql statement: identifies columns, makes table
       
        cat("starting", postab$Year[i]," ...\n")
        expr    <- with(postab, paste0("select substr(V1,",MAGEstart[i],",",MAGEwidth[i],") MAGE, substr(V1,",
                        FAGEstart[i],",",FAGEwidth[i],") FAGE, substr(V1,",SEXstart[i],",",SEXwidth[i],") Sex, ",
                        wtop[i],"(substr(V1,",Wtstart[i],",",Wtwidth[i],")) Births from f group by MAGE, FAGE, Sex"))
        # establish connection
        f       <- file(file.names[i])
     
        # call to sql, creates data.frame
        Bxy     <- sqldf(expr, dbname = tempfile(), drv = "SQLite", file.format = list(eol = "\n"))
        # clean connection rm(Bxy)
        close(f)
               
        # save
        name.out <- paste0("Bxy", postab$Year[i], ".Rdata")
        save(Bxy, file = file.path("/home/triffe/git/Dissertation/DISSERTATION/DATA",name.out))
        rm(Bxy)
        gc()
    }
}

# there were problems with the newer bigger files- it's faster to read.fwf them than it is to
# figure out SQLite on them.. 2003 onward, file positions are the same. problem apparently due 
# to missing eol at very end...
# ------------------------------------------------------------
library(reshape2)
postab      <- read.table("/home/triffe/git/Dissertation/DISSERTATION/R/UStabwidths.txt", 
                    header = TRUE, stringsAsFactors = FALSE)
            
           
# make Min yri <- 1
for (yri in 2:nrow(postab)){
    porder      <- with(postab[yri,],order(c(FAGEstart, MAGEstart, SEXstart, Wtstart)))
    cols        <- c("FAGE","MAGE","SEX","WT")
    widths      <- with(postab[yri,],c(FAGEwidth, MAGEwidth, SEXwidth, Wtwidth))
    starts      <- with(postab[yri,],c(FAGEstart, MAGEstart, SEXstart, Wtstart))
    blanks      <- diff(c(0,starts[porder])) - c(1,widths[porder][-length(widths)])
    this.year   <- postab$Year[yri]

    # squeeze to ASCII?
    origfile <- paste0("/home/triffe/DATA/CDC/BIRTHS/ZIPPED/B",this.year,".txt.gz")
    saveto   <- paste0("/home/triffe/DATA/CDC/BIRTHS/ZIPPED/B",this.year,"-conv.txt.gz")
    
    # 1) convert unknowns to ASCII.. - to work need to detect original encoding...
    system(paste0("zcat ", origfile, " | iconv -f ISO8859-15 -t US-ASCII//TRANSLIT | gzip >", saveto)) # WORKS! (1969)
    # 2) include only numbers? less necessary now
#    system(paste0("zcat /home/triffe/DATA/CDC/BIRTHS/ZIPPED/B", this.year,
#                    ".txt.gz | sed -e 's/[^0-9]/ /g' | gzip > /home/triffe/DATA/CDC/BIRTHS/ZIPPED/B",  
#                    this.year, "-num.txt.gz"))  
    # 3) parse only needed columns.
    regexp <- paste0("zcat /home/triffe/DATA/CDC/BIRTHS/ZIPPED/B", this.year, "-conv.txt.gz | sed -e 's/^",
                    paste(rep(".", blanks[1]), collapse = ""),                        # initial blanks
                    "\\(", paste(rep(".", widths[porder][1]), collapse = ""),"\\)",   # first column
                            paste(rep(".", blanks[2]), collapse = ""),                        # second blanks
                    "\\(", paste(rep(".", widths[porder][2]), collapse = ""),"\\)",   # second col
                            paste(rep(".", blanks[3]), collapse = ""),                        # third blanks
                    "\\(", paste(rep(".", widths[porder][3]), collapse = ""),"\\)",   # third col
                            paste(rep(".", blanks[4]), collapse = ""),                        # fouth blanks
                    "\\(", paste(rep(".", widths[porder][4]), collapse = ""),"\\)",   # fourth col
                    ".*$/",                                     # closeout
                    paste(rep(" ", blanks[1]), collapse = ""),                        # initial blanks
                    "\\", 1,                     # first column
                    paste(rep(" ", blanks[2]), collapse = ""),                        # second blanks
                    "\\", 2,                    # second col
                    paste(rep(" ", blanks[3]), collapse = ""),                        # third blanks
                    "\\", 3,                    # third col
                    paste(rep(" ", blanks[4]), collapse = ""),                        # fouth blanks
                    "\\", 4,
                    "/' | gzip > /home/triffe/DATA/CDC/BIRTHS/ZIPPED/B",  this.year, "-min.txt.gz")
    cat("\n\n\nStarting",this.year,"....\nRunning sed on file to reduce to needed columns\n")
    system(regexp)
    
    # now read in
    cat("\nReading back into R for reshaping\n")
    DAT <- read.table(paste0("/home/triffe/DATA/CDC/BIRTHS/ZIPPED/B", this.year, "-min.txt.gz"), 
            header = FALSE, as.is = TRUE, stringsAsFactors = FALSE, strip.white = TRUE, col.names = cols[porder])
    if (yri <= 3){
        DAT$WT <- 2
    }
    
    if (yri >= 35){
        DAT$WT <- 1
    }
    # make sure
    if (class(DAT$SEX) == "character"){
        DAT$SEX <- toupper(DAT$SEX)
        DAT$SEX <- as.integer(ifelse(DAT$SEX == "F", 2, 1)) 
    }

    Bxy <- melt(acast(DAT, MAGE~FAGE ~ SEX, sum, value.var = "WT"), value.name="Births")
    colnames(Bxy) <- c("MAGE", "FAGE", "SEX","BIRTHS") 
    file.name <-  paste0("/home/triffe/DATA/CDC/BIRTHS/Bxy/Bxy", this.year,".Rdata") 
    save(Bxy, file = file.name)
    cat("Year", this.year, "done. Saved to:\n",file.name,"\n\n")
}
# add 13 to MAGE in 2003- used MAGER41 variable:
Bxy <- local(get(load("/home/triffe/DATA/CDC/BIRTHS/Bxy/Bxy2003.Rdata")))
Bxy$MAGE <- Bxy$MAGE + 13
save(Bxy, file = "/home/triffe/DATA/CDC/BIRTHS/Bxy/Bxy2003.Rdata")

# read all in
All.paths <- file.path("/home/triffe/DATA/CDC/BIRTHS/Bxy", paste0("Bxy",1969:2010,".Rdata"))
All.Bxy <- lapply(All.paths, function(x){
            local(get(load(x)))
        })
source("R/UtilityFunctions.R")
names(All.Bxy) <- 1969:2010

x <- All.Bxy[["1989"]]
Bxymf <- lapply(All.Bxy, function(x){
              
            Bxym <- reshape2::acast(x[x$SEX == 1, ], MAGE ~ FAGE, sum, value.var = "BIRTHS")
            a99  <- Bxym[,"99"]
            Bxym <- Bxym[, -ncol(Bxym)]
            Bxym <- Bxym + (Bxym / rowSums(Bxym,na.rm=TRUE)) * a99 
            
            Bxyf <- reshape2::acast(x[x$SEX == 2, ], MAGE ~ FAGE, sum, value.var = "BIRTHS")
            a99  <- Bxyf[,"99"]
            Bxyf <- Bxyf[, -ncol(Bxyf)]
            Bxyf <- Bxyf + (Bxyf / rowSums(Bxyf,na.rm=TRUE)) * a99 

            # fix for 1989 mostly
            if ("89" %in% colnames(Bxym) & "89" %in% colnames(Bxyf)){
            if (colSums(Bxym,na.rm=TRUE)["89"] > 100 | colSums(Bxyf,na.rm=TRUE)["89"] > 100){
                a89  <- Bxym[,"89"]
                Bxym <- Bxym[, colnames(Bxym) != "89"]
                Bxym <- Bxym + (Bxym / rowSums(Bxym,na.rm=TRUE)) * a89 
                
                a89  <- Bxyf[,"89"]
                Bxyf <- Bxyf[, colnames(Bxyf) != "89"]
                Bxyf <- Bxyf + (Bxyf / rowSums(Bxyf,na.rm=TRUE)) * a89 
            }}

            list(Bxym = Mna0(t(Bxym)),Bxyf = Mna0(t(Bxyf)))
        })
#BxymfMISSING <- lapply(All.Bxy, function(x){
#            
#            Bxym <- reshape2::acast(x[x$SEX == 1, ], MAGE ~ FAGE, sum, value.var = "BIRTHS")
#            a99m  <- Bxym[,"99"]
#             
#            
#            Bxyf <- reshape2::acast(x[x$SEX == 2, ], MAGE ~ FAGE, sum, value.var = "BIRTHS")
#            a99f  <- Bxyf[,"99"]
#            
#            list(NAm = a99m, NAf = a99f)
#        })
#yearsUS <- 1969:2010
#names(BxymfMISSING) <- yearsUS
#names(Bxymf) <- yearsUS

#USNA <- lapply(as.character(yearsUS),function(yr, .Bxy, .NA){
#                    
#                    NAs <- .NA[[yr]][[1]] + .NA[[yr]][[2]]
#                    Bx <- (colSums(.Bxy[[yr]][[2]]) + colSums(.Bxy[[yr]][[2]]))
#                    NAx <- NAs / Bx
#                    NAt <- sum(NAs) / sum(Bx)
#                    list(NAt,NAx)
#                }, .Bxy = Bxymf, .NA = BxymfMISSING)
#USNAt <- unlist(lapply(USNA,"[[",1))
#plot(yearsUS, USNAt, type = 'l')
#save(USNA, file = "Data/results/USmissings/USNAt.Rdata")
#for (i in 1:length(yearsUS)){
#    plot(USNA[[i]][[2]], type = 'l', main = yearsUS[i])
#    Sys.sleep(1)
#}

Bxymf0_110 <- lapply(Bxymf, function(x){
            Xf <- Xm <- matrix(0, nrow = 111, ncol = 111, dimnames=list(Males = 0:110, Females = 0:110))
            Xm[rownames(x[["Bxym"]]),colnames(x[["Bxym"]])] <- x[["Bxym"]]
            Xf[rownames(x[["Bxyf"]]),colnames(x[["Bxyf"]])] <- x[["Bxyf"]]
            list(Bxym = Xm, Bxyf = Xf)
        })


Bxymf10_65 <- lapply(Bxymf0_110, function(x){   
            xm <- x[["Bxym"]]
            xm[, 11] <- rowSums(xm[, 1:11])
            xm[11, ] <- colSums(xm[1:11, ])
            xm[, 66] <- rowSums(xm[, 66:101])
            xm[66, ] <- colSums(xm[66:101, ])
            
            xf <- x[["Bxyf"]]
            xf[, 11] <- rowSums(xf[, 1:11])
            xf[11, ] <- colSums(xf[1:11, ])
            xf[, 66] <- rowSums(xf[, 66:101])
            xf[66, ] <- colSums(xf[66:101, ])
           
            list(Bxym = xm[11:66, 11:66], Bxyf = xf[11:66, 11:66])
        })

names(Bxymf0_110) <- 1969:2009
names(Bxymf10_65) <- 1969:2009
save(Bxymf0_110 , file="Data/USbirths/USBxymf0_110.Rdata")
save(Bxymf10_65 , file="Data/USbirths/USBxymf10_65.Rdata")


# stick in 0_100 tables and in 10-65 tables (always same dims, named dims)
Bxy0_110 <- lapply(All.Bxy, function(x){
            X <- matrix(0, nrow = 111, ncol = 111, dimnames=list(Males = 0:110, Females = 0:110))
            X[rownames(x),colnames(x)] <- x
            X
        })
Bxy10_65 <- lapply(Bxy0_110, function(x){
            x[, 11] <- rowSums(x[, 1:11])
            x[11, ] <- colSums(x[1:11, ])
            x[, 66] <- rowSums(x[, 66:101])
            x[66, ] <- colSums(x[66:101, ])
            x[11:66, 11:66]
        })
names(Bxy10_65) <-names(Bxy0_110) <- 1969:2010

#save(Bxy0_110, file = "/home/triffe/git/Dissertation/DISSERTATION/DATA/USbirths/USBxy0_110.Rdata")
#save(Bxy10_65, file = "/home/triffe/git/Dissertation/DISSERTATION/DATA/USbirths/USBxy10_65.Rdata")

# making standard tables, including in single list, standardizing ages, distributing missing ages.
#ES <- local(get(load("/home/triffe/git/Dissertation/DISSERTATION/DATA/ESbirths/ESbirths.Rdata")))

# NOW both USBirths and ESBirths are in the same format: 
# lists, with a MAGE x FAGE matrix ages 0-110 in each year

# and finally, make an object, Bxy for spain:
#
#B       <- local(get(load("/home/triffe/git/Dissertation/DISSERTATION/DATA/ESbirths/ESbirths.Rdata")))
#B       <- B[-length(B)] # remove 2010 births for now
#
#Bxy <- lapply(as.character(1975:2009), function(i, B){
#            reshape2:::acast(B[[i]], FAGE ~ MAGE, sum, value.var = "Births")
#        }, B = B)
#names(Bxy) <- 1975:2009
#save(Bxy, file = "/home/triffe/git/Dissertation/DISSERTATION/DATA/ESbirths/ESBxy.Rdata")

BxymfES <- local(get(load("Data/ESbirths/ESBxymf.Rdata")))
BxymfES10_65 <- lapply(BxymfES, function(x){   
            xm <- x[["Bxym"]]
            xm[, 11] <- rowSums(xm[, 1:11])
            xm[11, ] <- colSums(xm[1:11, ])
            xm[, 66] <- rowSums(xm[, 66:101])
            xm[66, ] <- colSums(xm[66:101, ])
            
            xf <- x[["Bxyf"]]
            xf[, 11] <- rowSums(xf[, 1:11])
            xf[11, ] <- colSums(xf[1:11, ])
            xf[, 66] <- rowSums(xf[, 66:101])
            xf[66, ] <- colSums(xf[66:101, ])
            
            list(Bxym = xm[11:66, 11:66], Bxyf = xf[11:66, 11:66])
        })
save(BxymfES10_65, file = "Data/ESbirths/ESBxymf10_65.Rdata")

BxymfES <- local(get(load("Data/ESbirths/ESBxymf.Rdata")))

library(HMDget)
username <- "XXX"
password <- "XXX"
years <- list()
years$mltper_1x1 <- c(1969:2009)
LxmUS <- HMDget(countries = c("USA"), wanteditems = c("mltper_1x1"), 
        years = years, column ="Lx", drop.tadj = TRUE, format = 4, 
        username = username, password = password)
years <- list()
years$fltper_1x1 <- c(1969:2009)
LxfUS <- HMDget(countries = c("USA"), wanteditems = c("fltper_1x1"), 
        years = years, column ="Lx", drop.tadj = TRUE, format = 4, 
        username = username, password = password)

years <- list()
years$mltper_1x1 <- c(1975:2009)
LxmES <- HMDget(countries = c("ESP"), wanteditems = c("mltper_1x1"), 
        years = years, column ="Lx", drop.tadj = TRUE, format = 4, 
        username = username, password = password)
years <- list()
years$fltper_1x1 <- c(1975:2009)
LxfES <- HMDget(countries = c("ESP"), wanteditems = c("fltper_1x1"), 
        years = years, column ="Lx", drop.tadj = TRUE, format = 4, 
        username = username, password = password)
save(LxmUS, file = "Data/HMD_Lx/LxmUS.Rdata")
save(LxfUS, file = "Data/HMD_Lx/LxfUS.Rdata")
save(LxmES, file = "Data/HMD_Lx/LxmES.Rdata")
save(LxfES, file = "Data/HMD_Lx/LxfES.Rdata")

years <- list()
years$mltper_1x1 <- c(1969:2009)
dxmUS <- HMDget(countries = c("USA"), wanteditems = c("mltper_1x1"), 
        years = years, column = "dx", drop.tadj = TRUE, format = 4, 
        username = username, password = password)
years <- list()
years$fltper_1x1 <- c(1969:2009)
dxfUS <- HMDget(countries = c("USA"), wanteditems = c("fltper_1x1"), 
        years = years, column = "dx", drop.tadj = TRUE, format = 4, 
        username = username, password = password)

years <- list()
years$mltper_1x1 <- c(1975:2009)
dxmES <- HMDget(countries = c("ESP"), wanteditems = c("mltper_1x1"), 
        years = years, column = "dx", drop.tadj = TRUE, format = 4, 
        username = username, password = password)
years <- list()
years$fltper_1x1 <- c(1975:2009)
dxfES <- HMDget(countries = c("ESP"), wanteditems = c("fltper_1x1"), 
        years = years, column = "dx", drop.tadj = TRUE, format = 4, 
        username = username, password = password)
save(dxmUS, file = "Data/HMD_dx/dxmUS.Rdata")
save(dxfUS, file = "Data/HMD_dx/dxfUS.Rdata")
save(dxmES, file = "Data/HMD_dx/dxmES.Rdata")
save(dxfES, file = "Data/HMD_dx/dxfES.Rdata")

# get period e0 estimates, all years:
library(HMDget)
years <- list()
years$E0per <- c(1933:2009)
e0perUS <- HMDget(countries = c("USA"), wanteditems = c("E0per"), years = years, 
        column ="", drop.tadj = TRUE, format = 0, username = username, password = password)
years$E0per <- c(1908:2009)
e0perES <- HMDget(countries = c("ESP"), wanteditems = c("E0per"), years = years, 
        column ="", drop.tadj = TRUE, format = 0, username = username, password = password)

colnames(e0perUS) <- colnames(e0perES) <- c("Year","Female","Male","Total")

save(e0perUS, file = "Data/HMD_e0period/e0perUS.Rdata")
save(e0perES, file = "Data/HMD_e0period/e0perES.Rdata")


years <- list()
years$mltper_1x1 <- c(1969:2009)
lxmUS <- HMDget(countries = c("USA"), wanteditems = c("mltper_1x1"), 
        years = years, column = "lx", drop.tadj = TRUE, format = 4, 
        username = username, password = password)
years <- list()
years$fltper_1x1 <- c(1969:2009)
lxfUS <- HMDget(countries = c("USA"), wanteditems = c("fltper_1x1"), 
        years = years, column = "lx", drop.tadj = TRUE, format = 4, 
        username = username, password = password)

years <- list()
years$mltper_1x1 <- c(1975:2009)
lxmES <- HMDget(countries = c("ESP"), wanteditems = c("mltper_1x1"), 
        years = years, column = "lx", drop.tadj = TRUE, format = 4, 
        username = username, password = password)
years <- list()
years$fltper_1x1 <- c(1975:2009)
lxfES <- HMDget(countries = c("ESP"), wanteditems = c("fltper_1x1"), 
        years = years, column = "lx", drop.tadj = TRUE, format = 4, 
        username = username, password = password)
save(lxmUS, file = "Data/HMD_lx/lxmUS.Rdata")
save(lxfUS, file = "Data/HMD_lx/lxfUS.Rdata")
save(lxmES, file = "Data/HMD_lx/lxmES.Rdata")
save(lxfES, file = "Data/HMD_lx/lxfES.Rdata")

years <- list()
years$mltper_1x1 <- c(1969:2009)
muxmUS <- HMDget(countries = c("USA"), wanteditems = c("mltper_1x1"), 
        years = years, column = "mx", drop.tadj = TRUE, format = 4, 
        username = username, password = password)
years <- list()
years$fltper_1x1 <- c(1969:2009)
muxfUS <- HMDget(countries = c("USA"), wanteditems = c("fltper_1x1"), 
        years = years, column = "mx", drop.tadj = TRUE, format = 4, 
        username = username, password = password)
years <- list()
years$mltper_1x1 <- c(1969:2009)
muxmES <- HMDget(countries = c("ESP"), wanteditems = c("mltper_1x1"), 
        years = years, column = "mx", drop.tadj = TRUE, format = 4, 
        username = username, password = password)
years <- list()
years$fltper_1x1 <- c(1969:2009)
muxfES <- HMDget(countries = c("ESP"), wanteditems = c("fltper_1x1"), 
        years = years, column = "mx", drop.tadj = TRUE, format = 4, 
        username = username, password = password)
save(muxmUS, file = "Data/HMD_mux/muxmUS.Rdata")
save(muxfUS, file = "Data/HMD_mux/muxfUS.Rdata")
save(muxmES, file = "Data/HMD_mux/muxmES.Rdata")
save(muxfES, file = "Data/HMD_mux/muxfES.Rdata")

library(HMDget)

head(LxmES)
PxES <- HMDget(countries = c("ESP"), wanteditems = c("Population"), 
        years = list(Population = c(1975:2009)), drop.tadj = TRUE, format = 1, 
        username = username, password = password)
colnames(PxES) <- c("Year","Age","Female","Male","Total")

PxUS <- HMDget(countries = c("USA"), wanteditems = c("Population"), 
        years =  list(Population = c(1969:2009)), drop.tadj = TRUE, format = 1, 
        username = username, password = password)
colnames(PxUS) <- c("Year","Age","Female","Male","Total")
save(PxES, file = "Data/HMD_Px/PxES.Rdata")
save(PxUS, file = "Data/HMD_Px/PxUS.Rdata")


DTLTUES <- HMDget(countries = c("ESP"), wanteditems = c("Deaths_lexis"), 
        years = list(Deaths_lexis = c(1975:2009)), drop.tadj = TRUE, format = 1, 
        username = username, password = password)
colnames(DTLTUES) <- c("Year","Age","Cohort","Female","Male","Total")

DTLTUUS <- HMDget(countries = c("USA"), wanteditems = c("Deaths_lexis"), 
        years =  list(Deaths_lexis = c(1969:2009)), drop.tadj = TRUE, format = 1, 
        username = username, password = password)

colnames(DTLTUUS) <- c("Year","Age","Cohort","Female","Male","Total")


save(DTLTUES, file = "Data/HMD_TLTU/DTLTUES.Rdata")
save(DTLTUUS, file = "Data/HMD_TLTU/DTLTUUS.Rdata")

# --------------------------------------
# add to Ropes folder with full years...

DTLTUES <- HMDget(countries = c("ESP"), wanteditems = c("Deaths_lexis"), 
        years = list(Deaths_lexis = c(1908:2009)), drop.tadj = TRUE, format = 1, 
        username = username, password = password)
colnames(DTLTUES) <- c("Year","Age","Cohort","Female","Male","Total")

DTLTUUS <- HMDget(countries = c("USA"), wanteditems = c("Deaths_lexis"), 
        years =  list(Deaths_lexis = c(1933:2010)), drop.tadj = TRUE, format = 1, 
        username = username, password = password)

colnames(DTLTUUS) <- c("Year","Age","Cohort","Female","Male","Total")


save(DTLTUES, file = "Data/RopeData/DTLTUES.Rdata")
save(DTLTUUS, file = "Data/RopeData/DTLTUUS.Rdata")