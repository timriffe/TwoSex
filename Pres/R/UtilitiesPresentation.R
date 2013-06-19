makeRect <- function(x, y, w, h, ...){
    w2 <- w / 2
    h2 <- h / 2
    rect(x - w2, y - h2, x + w2, y + h2, ...)
}
rotationmat <- compiler::cmpfun(function(theta){
            matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),ncol=2)
        })
rotatevec <- compiler::cmpfun(function(xy,Rmat){
            Rmat %*% c(xy)
        })
# gold here. work on this one
makearc <- compiler::cmpfun(function(x1,y1,x2,y2,bc=.05,nvert=100){
            
            if (x1==x2 & y1 == y2){
                return(cbind(x=rep(x1,nvert),y = rep(y1,nvert)))
            }
            rise    <- y2 - y1 #y2 - y1
            run     <- x2 - x1 #x2-x1
            a       <- .5 * (rise^2 + run^2)^.5 # 1/2 hypotenuse
            b       <- a * bc # how bowed?
            add     <- ifelse(x1 < x2, 0, pi) # direction
            theta   <- add + atan(rise/run)
            Rmat    <- rotationmat(theta)
            pivec   <- seq(from = 0, to = pi, length.out = nvert)
            x       <- a * cos(pivec) + a
            y       <- b * sin(pivec)
            t(apply(cbind(x, y), 1, rotatevec, Rmat = Rmat) + c(x1, y1))   
        })

makearcv <- compiler::cmpfun(function(x1,y1,x2,y2,bc=.05,nvert=100,female = TRUE){
            
            if (x1 == x2 & y1 == y2){
                return(c(rep(x1,nvert),rep(y1,nvert)))
            }
            rise    <- y2 - y1 #y2 - y1
            run     <- x2 - x1 #x2-x1
            a       <- .5 * (rise^2 + run^2)^.5 # 1/2 hypotenuse
            b       <- a * bc # how bowed?
            #add     <- ifelse(female, ifelse(x1<x2,-2*pi,-pi), ifelse(x1<x2,0,pi))
            add     <- ifelse(x1 < x2, 0, pi) # direction
            #add     <- pi * female
            theta   <- add + atan(rise/run) #* ifelse(female,-1,1)
            theta   <- ifelse(female,{-2*pi+theta},theta)
            Rmat    <- rotationmat(theta)
            pivec   <- seq(from = 0, to = pi, length.out = nvert)
            x       <- a * cos(pivec) + a
            y       <- b * sin(pivec)
            c(t(apply(cbind(x, y), 1, rotatevec, Rmat = Rmat) + c(x1, y1)))   
        })
