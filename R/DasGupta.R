# preamble:
source("/home/triffe/git/DISS/R/UtilityFunctions.R")
source("/home/triffe/git/DISS/R/MeanFunctions.R")

# BxUS is a list of 56x56 matrices, ages 10-65, males in rows, females in columns
# (1969 - 2010)
# BxES is 0:110, years 1975:2009
BxUS <- local(get(load("/home/triffe/git/DISS/Data/USbirths/USBxy0_110.Rdata")))
BxES <- local(get(load("/home/triffe/git/DISS/Data/ESbirths/ESBxy.Rdata"))) # not cut to 10-65
# exposures, as such, straiht from HMD, all ages 0-110, long form
ExUS <- local(get(load("/home/triffe/git/DISS/Data/Exposures/USexp.Rdata")))
ExES <- local(get(load("/home/triffe/git/DISS/Data/Exposures/ESexp.Rdata")))
LxmUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmUS.Rdata"))) / 1e5
LxfUS <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfUS.Rdata"))) / 1e5
LxmES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxmES.Rdata"))) / 1e5
LxfES <- local(get(load("/home/triffe/git/DISS/Data/HMD_Lx/LxfES.Rdata"))) / 1e5

#--------------------------
# Das Gupta 1972:

# implementing Das Gupta 1978:
pf. <- LxfUS[,"2000"]
pm. <- LxmUS[,"2000"]
Bxy         <- BxUS[["2000"]]

PM <- with(ExUS, Male[Year == 2000])
PF <- with(ExUS, Female[Year == 2000])        

Bxy <- BxUS[["2000"]]
# eq 5.4
Uxy <- Minf0(Mna0(Bxy / rowSums(Bxy)))
Vxy <- Minf0(Mna0(t(t(Bxy) / colSums(Bxy))))


# theta is prop male of births
thetam <- 1.05 / 2.05
thetaf <- 1 / 2.05


mxy <- Minf0(Mna0(Bxy /(Uxy * PM + t(t(Vxy) * PF))))
thetam * pm. * Uxy 
thetaf * pf. * t(Vxy)

Gupta1978Min <- function(r, Uxy, Vxy, lxm, lxf, mxy,.a = .5:110.5){
    (1 - sum(mxy * (thetam * exp(-r * .a) * pm. * Uxy + t(thetaf * exp(-r * .a) * pf. * t(Vxy)))))^2
}
optimize(Gupta1978Min, interval = c(-.02, .02), tol = 1e-15,
        Uxy = Uxy, Vxy = Vxy, lxm = pm., lxf = pf., mxy = mxy)
#mxy <- Mna0(Bxy / (Uxy %row% (Pxm / 1) + Vxy %col% (Pxf / 1)))

unlist(lapply(yearsUS),function())

Gupta1978It <- function()

