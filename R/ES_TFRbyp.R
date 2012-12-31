
# Author: triffe
###############################################################################

Fxyp    <- local(get(load("/home/triffe/git/Dissertation/DISSERTATION/DATA/Rates.p/ES_Fxyp.Rdata")))

TFRp <- do.call(cbind, lapply(Fxyp, function(x){
            apply(x, 3, sum, na.rm = TRUE)
        }))
head(TFRp)
rownames(TFRp) <- seq(.005,.995,.005)

range(TFRp)
plot(NULL, type = "n", ylim = c(1,3), xlim = c(1975,2009))
apply(TFRp, 1, function(y){
            lines(1975:2009,y,col = "#DDDDDD40")
        })
#StolarskyM(20, 40, -2)