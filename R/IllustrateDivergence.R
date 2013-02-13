
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
# BxES is 0:110, years 1975:2009
BxUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy10_65.Rdata")))
BxES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65
# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))

# cut ES to same dimentions ages 10:65, square (summing boundary ages)
BxES <- lapply(BxES, function(x){
            x[11, ] <- colSums(x[1:11, ])
            x[, 11] <- rowSums(x[, 1:11])
            x[66, ] <- colSums(x[66:110, ])
            x[, 66] <- rowSums(x[, 66:110])
            x[11:66,11:66]
        })

# get years for each pop
yearsES <- as.integer(names(BxES))
yearsUS <- as.integer(names(BxUS))
ages    <- 10:65
# ASFR sample year 1975:

BxyUS <- BxUS[["1975"]]
BxyES <- BxES[["1975"]]

ExmUS <- with(ExUS, Male[Year == 1975 & Age >= 10 & Age <= 65])
ExfUS <- with(ExUS, Female[Year == 1975 & Age >= 10 & Age <= 65])
ExmES <- with(ExES, Male[Year == 1975 & Age >= 10 & Age <= 65])
ExfES <- with(ExES, Female[Year == 1975 & Age >= 10 & Age <= 65])

ASFRmUS <- rowSums(BxyUS) / ExmUS
ASFRfUS <- colSums(BxyUS) / ExfUS
ASFRmES <- rowSums(BxyES) / ExmES
ASFRfES <- colSums(BxyES) / ExfES

max(c(ASFRmUS,ASFRfUS,ASFRmES,ASFRfES))

par(xaxp = "i", yaxp = "i")
plot(ages, ASFRmUS, type = 'l', ylim = c(0, .21), xlim = c(10,66), axes = FALSE,
        panel.first = list(rect(10,0,66,.21,col = gray(.93), border=NA),
                abline(h = seq(0,.2,by = .02), col = "white")))






# simple series of TFR
TFRm <- vector(length = length())







