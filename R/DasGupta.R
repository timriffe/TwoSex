source("/home/triffe/git/DISS/R/UtilityFunctions.R")


BxUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy10_65.Rdata")))
# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))


# not sure what I'm supposed to see in equation 4.1 everything cancels!
Pxy <- Mxy + Fxy

Pxm <- with(ExUS, Male[Age >= 10 & Age <= 65 & Year == 1970])
Pxf <- with(ExUS, Female[Age >= 10 & Age <= 65 & Year == 1970])        

Bxy <- BxUS[["1970"]]
# eq 5.4
Uxy <- Bxy %row/% rowSums(Bxy)
Vxy <- Bxy %col/% colSums(Bxy)







