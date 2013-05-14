setwd("/home/triffe/git/DISS/")
# this script is what eventually verified for me that there is no hope in empirically finding the 
# best mean function.
# idea: don't compare absolute birth counts but pdfs instead, should remove problem with
# declining tendency, but I'll need to think of whether shape optimization is justifiable

BxUS <- local(get(load("Data/USbirths/USBxy10_65.Rdata")))
# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
ExUS <- local(get(load("Data/Exposures/USexp.Rdata")))


source("R/MeanFunctions.R")

# more helper functions
logit <- function(x){
    log(x/(1-x)) 
}
expit <- function(x){
    exp(x)/(1+exp(x))
}
# experimenting: only one'p' here
min.fun <- compiler::cmpfun(function(p, bm.hat.pdf, bf.hat.pdf, b2.pdf, .Bx2.pdf){
    sqrt(sum((StolarskyM(bm.hat.pdf,bf.hat.pdf, p) - .Bx2.pdf)^2, na.rm=TRUE))
})

years <- as.integer(names(BxUS))
p.hat.vec <- vector(length = length(years) - 2)
for (i in 1:(length(years) - 2)){
    
    # observed births from year 1 and 2
    Bx1  <- BxUS[[i]]
    Bx2  <- BxUS[[i+1]]
    
    # exposures from years 1 and 2
    Exf1 <- ExUS$Female[with(ExUS, Year == years[i] & Age %in% c(10:65))]
    Exm1 <- ExUS$Male[with(ExUS, Year == years[i] & Age %in% c(10:65))]
    Exf2 <- ExUS$Female[with(ExUS, Year == years[i+1] & Age %in% c(10:65))]
    Exm2 <- ExUS$Male[with(ExUS, Year == years[i+1] & Age %in% c(10:65))]
    
    # year 1 male and female rates
    Fxm  <- Bx1 / Exm1
    Fxf  <- t(t(Bx1) / Exf1) # females in columns
    
    # year 2 estimated births using year 1 rates, single sex
    Bxm.hat <- Fxm * Exm2
    Bxf.hat <- t(t(Fxf) * Exf2)
    
    # the corresponding pdfs of male and female estimatd births
    Bxm.hat.pdf <- Bxm.hat / sum(Bxm.hat) 
    Bxf.hat.pdf <- Bxf.hat / sum(Bxf.hat) 
    
    # observed year  pdf
    Bx2.pdf   <- Bx2 / sum(Bx2)
    opt <- optimize(min.fun, bm.hat.pdf = Bxm.hat.pdf, bf.hat.pdf = Bxf.hat.pdf, .Bx2.pdf = Bx2.pdf,
            lower = -200, upper = 200)$minimum  
    p.hat.vec[i] <- opt
}


plot(1970:2009, p.hat.vec, type = 'l')

# that's one 'p' for the whole matrix.
# what if the matrix divides into zones with different 'p.hat' values?
# need to think of how to optimize over blocks within the matrix.
# also need to think of how not to let one or two odd cohorts determine
# the entire fit.

# also try 3x3 fits. see if they're the same.

# 1) easiest: do separate for upper and lower tri. 2 p estimates...
min.fun.ul <- compiler::cmpfun(function(pul, bm.hat.pdf, bf.hat.pdf, b2.pdf, .Bx2.pdf){
    log(
    sqrt(sum((StolarskyM(bm.hat.pdf[lower.tri(bm.hat.pdf)],
                         bf.hat.pdf[lower.tri(bf.hat.pdf)], logit(pul["pl"])) - 
                 .Bx2.pdf[lower.tri(.Bx2.pdf)]) ^ 2, na.rm = TRUE)) +
    sqrt(sum((StolarskyM(bm.hat.pdf[upper.tri(bm.hat.pdf)],
                                    bf.hat.pdf[upper.tri(bf.hat.pdf)], logit(pul["pu"])) - 
                            .Bx2.pdf[upper.tri(.Bx2.pdf)]) ^ 2, na.rm = TRUE)))
})

min.lower <- compiler::cmpfun(function(pl, bm.hat.pdf, bf.hat.pdf, b2.pdf, .Bx2.pdf){
    sum((StolarskyM(bm.hat.pdf[lower.tri(bm.hat.pdf)],
                         bf.hat.pdf[lower.tri(bf.hat.pdf)], logit(pl)) - 
                 .Bx2.pdf[lower.tri(.Bx2.pdf)]) ^ 2, na.rm = TRUE)
        })

p.l.vec <- vector(length = length(years) - 2)
p.hat.mat <- data.frame(pu = vector(length = length(years) - 2),
                        pl = vector(length = length(years) - 2))
for (i in 1:(length(years) - 2)){
    
    # observed births from year 1 and 2
    Bx1  <- BxUS[[i]]
    Bx2  <- BxUS[[i+1]]
    
    # exposures from years 1 and 2
    Exf1 <- ExUS$Female[with(ExUS, Year == years[i] & Age %in% c(10:65))]
    Exm1 <- ExUS$Male[with(ExUS, Year == years[i] & Age %in% c(10:65))]
    Exf2 <- ExUS$Female[with(ExUS, Year == years[i+1] & Age %in% c(10:65))]
    Exm2 <- ExUS$Male[with(ExUS, Year == years[i+1] & Age %in% c(10:65))]
    
    # year 1 male and female rates
    Fxm  <- Bx1 / Exm1
    Fxf  <- t(t(Bx1) / Exf1) # females in columns
    
    # year 2 estimated births using year 1 rates, single sex
    Bxm.hat <- Fxm * Exm2
    Bxf.hat <- t(t(Fxf) * Exf2)
    
    # the corresponding pdfs of male and female estimatd births
    Bxm.hat.pdf <- Bxm.hat / sum(Bxm.hat) 
    Bxf.hat.pdf <- Bxf.hat / sum(Bxf.hat) 
    
    # observed year  pdf
    Bx2.pdf   <- Bx2 / sum(Bx2)
    #starts <- c(pu = .5, pl = .5)
#    opt <- optim(starts, min.fun.ul, 
#            bm.hat.pdf = Bxm.hat.pdf, 
#            bf.hat.pdf = Bxf.hat.pdf, 
#            .Bx2.pdf = Bx2.pdf)
    p.l.vec[i] <- optimize(min.lower, bm.hat.pdf = Bxm.hat.pdf, bf.hat.pdf = Bxf.hat.pdf, .Bx2.pdf = Bx2.pdf,
            lower = 0, upper = 1)$minimum
    #p.hat.mat[i, ] <- opt$par
}

expit(.1)
expit(-4)
logit(-100)
plot(p.hat.mat[, 1])
lines(p.hat.mat[, 2])
plot(p.l.vec)


# single p, but over 3x3
source("R/UtilityFunctions.R")

years <- as.integer(names(BxUS))
p.hat.vec <- vector(length = length(years) - 2)
for (i in 1:(length(years) - 2)){
    
    # observed births from year 1 and 2
    Bx1  <- makeBlockAgeGroups(BxUS[[i]], N = 3)
    Bx2  <- makeBlockAgeGroups(BxUS[[i+1]], N = 3)
    
    # exposures from years 1 and 2
    Exf1 <- makeVectorAgeGroups(ExUS$Female[with(ExUS, Year == years[i] & Age %in% c(10:65))], ages = 10:65, N = 3)
    Exm1 <- makeVectorAgeGroups(ExUS$Male[with(ExUS, Year == years[i] & Age %in% c(10:65))], ages = 10:65, N = 3)
    Exf2 <- makeVectorAgeGroups(ExUS$Female[with(ExUS, Year == years[i+1] & Age %in% c(10:65))], ages = 10:65, N = 3)
    Exm2 <- makeVectorAgeGroups(ExUS$Male[with(ExUS, Year == years[i+1] & Age %in% c(10:65))], ages = 10:65, N = 3)
    
    # year 1 male and female rates
    Fxm  <- Bx1 / Exm1
    Fxf  <- t(t(Bx1) / Exf1) # females in columns
    
    # year 2 estimated births using year 1 rates, single sex
    Bxm.hat <- Fxm * Exm2
    Bxf.hat <- t(t(Fxf) * Exf2)
    
    # the corresponding pdfs of male and female estimatd births
    Bxm.hat.pdf <- Bxm.hat / sum(Bxm.hat) 
    Bxf.hat.pdf <- Bxf.hat / sum(Bxf.hat) 
    
    # observed year  pdf
    Bx2.pdf   <- Bx2 / sum(Bx2)
    opt <- optimize(min.fun, bm.hat.pdf = Bxm.hat.pdf, bf.hat.pdf = Bxf.hat.pdf, .Bx2.pdf = Bx2.pdf,
            lower = -200, upper = 200)$minimum  
    p.hat.vec[i] <- opt
}
plot(1970:2009, p.hat.vec, type = 'l')



########################################################################3
# again, 3x3, but for upper / lower separately.
# make min fun flexible for upper vs lower
min.tri <- compiler::cmpfun(function(.p, bm.hat.pdf, bf.hat.pdf, b2.pdf, .Bx2.pdf, tri.fun = upper.tri){
            sum((StolarskyM(bm.hat.pdf[tri.fun(bm.hat.pdf)],
                                        bf.hat.pdf[tri.fun(bf.hat.pdf)], .p) - 
                                .Bx2.pdf[tri.fun(.Bx2.pdf)]) ^ 2, na.rm = TRUE)
        })

min.tri2 <- compiler::cmpfun(function(.p, bm.hat.pdf, bf.hat.pdf, b2.pdf, .Bx2.pdf, tri.fun = upper.tri){
            sum((StolarskyM(bm.hat.pdf * tri.fun(bm.hat.pdf, diag = TRUE),
                                        bf.hat.pdf * tri.fun(bf.hat.pdf, diag = TRUE), .p) - 
                                .Bx2.pdf *tri.fun(.Bx2.pdf, diag = TRUE)) ^ 2, na.rm = TRUE)
        })

p.l.vec <- vector(length = length(years) - 2)
p.u.vec <- p.l.vec


for (i in 1:(length(years) - 2)){
    
    # observed births from year 1 and 2
    Bx1  <- makeBlockAgeGroups(BxUS[[i]], N = 3)
    Bx2  <- makeBlockAgeGroups(BxUS[[i+1]], N = 3)
    
    # exposures from years 1 and 2
    Exf1 <- makeVectorAgeGroups(ExUS$Female[with(ExUS, Year == years[i] & Age %in% c(10:65))], ages = 10:65, N = 3)
    Exm1 <- makeVectorAgeGroups(ExUS$Male[with(ExUS, Year == years[i] & Age %in% c(10:65))], ages = 10:65, N = 3)
    Exf2 <- makeVectorAgeGroups(ExUS$Female[with(ExUS, Year == years[i+1] & Age %in% c(10:65))], ages = 10:65, N = 3)
    Exm2 <- makeVectorAgeGroups(ExUS$Male[with(ExUS, Year == years[i+1] & Age %in% c(10:65))], ages = 10:65, N = 3)
    
    # year 1 male and female rates
    Fxm  <- Bx1 / Exm1
    Fxf  <- t(t(Bx1) / Exf1) # females in columns
    
    # year 2 estimated births using year 1 rates, single sex
    Bxm.hat <- Fxm * Exm2
    Bxf.hat <- t(t(Fxf) * Exf2)
    
    # the corresponding pdfs of male and female estimatd births
    Bxm.hat.pdf <- Bxm.hat / sum(Bxm.hat) 
    Bxf.hat.pdf <- Bxf.hat / sum(Bxf.hat) 
    
    # observed year  pdf
    Bx2.pdf   <- Bx2 / sum(Bx2)
    #starts <- c(pu = .5, pl = .5)

    p.u.vec[i] <- optimize(min.tri2, bm.hat.pdf = Bxm.hat.pdf, bf.hat.pdf = Bxf.hat.pdf, 
            .Bx2.pdf = Bx2.pdf, tri.fun = upper.tri, lower = -200, upper = 500)$minimum
    p.l.vec[i] <- optimize(min.tri2, bm.hat.pdf = Bxm.hat.pdf, bf.hat.pdf = Bxf.hat.pdf, 
            .Bx2.pdf = Bx2.pdf, tri.fun = lower.tri, lower = -200, upper = 500)$minimum
}

plot(1970:2009, p.u.vec, type = 'l', col = 'red')
lines(1970:2009, p.l.vec, col = "blue")

# now instead of optimizing over the whole pdf, optimize over local pdfs-
# e.g., 3x3 cells of 3x3 blocks? something like that.
wmean(matrix(c(2,3,4,5),2), )
wmeanM2 <- compiler::cmpfun(function(M1, M2, w.frac){
    M1 * w.frac + M2 * w.frac
})

for (i in 1:(length(years) - 2)){
    
    # observed births from year 1 and 2
    Bx1  <- makeBlockAgeGroups(BxUS[[i]], N = 3)
    Bx2  <- makeBlockAgeGroups(BxUS[[i+1]], N = 3)
    
    # exposures from years 1 and 2
    Exf1 <- makeVectorAgeGroups(ExUS$Female[with(ExUS, Year == years[i] & Age %in% c(10:65))], ages = 10:65, N = 3)
    Exm1 <- makeVectorAgeGroups(ExUS$Male[with(ExUS, Year == years[i] & Age %in% c(10:65))], ages = 10:65, N = 3)
    Exf2 <- makeVectorAgeGroups(ExUS$Female[with(ExUS, Year == years[i+1] & Age %in% c(10:65))], ages = 10:65, N = 3)
    Exm2 <- makeVectorAgeGroups(ExUS$Male[with(ExUS, Year == years[i+1] & Age %in% c(10:65))], ages = 10:65, N = 3)
    
    # year 1 male and female rates
    Fxm  <- Bx1 / Exm1
    Fxf  <- t(t(Bx1) / Exf1) # females in columns
    
    # year 2 estimated births using year 1 rates, single sex
    Bxm.hat <- Fxm * Exm2
    Bxf.hat <- t(t(Fxf) * Exf2)
    
    # the corresponding pdfs of male and female estimatd births
    Bxm.hat.pdf <- Bxm.hat / sum(Bxm.hat) 
    Bxf.hat.pdf <- Bxf.hat / sum(Bxf.hat) 
    
    # observed year  pdf
    Bx2.pdf   <- Bx2 / sum(Bx2)
    #starts <- c(pu = .5, pl = .5)
    
    p.u.vec[i] <- optimize(min.tri2, bm.hat.pdf = Bxm.hat.pdf, bf.hat.pdf = Bxf.hat.pdf, 
            .Bx2.pdf = Bx2.pdf, tri.fun = upper.tri, lower = -200, upper = 500)$minimum
    p.l.vec[i] <- optimize(min.tri2, bm.hat.pdf = Bxm.hat.pdf, bf.hat.pdf = Bxf.hat.pdf, 
            .Bx2.pdf = Bx2.pdf, tri.fun = lower.tri, lower = -200, upper = 500)$minimum
}

















