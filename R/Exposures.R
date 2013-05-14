# script only run once to get exposures saved


# Author: triffe
###############################################################################
#setwd("/home/triffe/git/DISS/")
## get HMD exposures, see what's available:
#
#library(HMDget)
#years <- list()
#years$Exposures_1x1 <- c(1975:2009)
#ESexp <- HMDget(countries = c("ESP"), wanteditems = c("Exposures_1x1"), years = years, column ="", drop.tadj = TRUE, format = 0, username = "tim.riffe@gmail.com", password = "avery")
#years <- list()
#years$Exposures_1x1 <- c(1969:2009)
#USexp <- HMDget(countries = c("USA"), wanteditems = c("Exposures_1x1"), years = years, column ="", drop.tadj = TRUE, format = 0, username = "tim.riffe@gmail.com", password = "avery")
#
#colnames(USexp) <-colnames(ESexp) <- c("Year","Age","Female","Male","Total")
#
#save(USexp, file = "DATA/Exposures/USexp.Rdata")
#save(ESexp, file = "DATA/Exposures/ESexp.Rdata")


