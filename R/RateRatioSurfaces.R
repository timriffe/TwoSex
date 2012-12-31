# Author: triffe
###############################################################################

RateRatioSurface <- function(B, E, N = 1, min.events = 40, filename, age.range = c(10, 65), zlim = c(-.6, .6)){
    if (missing(filename)){
        stop("\nfirst give the output a name, including the suffix '.pdf'\n")
    }
    source("/home/triffe/git/Dissertation/DISSERTATION/R/UtilityFunctions.R")
    # color ramp function
    colfun          <- grDevices:::colorRampPalette(RColorBrewer:::brewer.pal(9, "RdBu"), space = "Lab")
    # years vectors
    years.char      <- names(B) # just make sure B and E are of same length
    years           <- as.integer(years.char)
    # value breaks
    brks            <- seq(zlim[1], zlim[2], by = .02)
    ages.all        <- 0:110
    ages.range      <- ages.all[ages.all >= min(age.range) & ages.all <= max(age.range)]
    ages.range.N    <- unique(ages.range - ages.range %% N)
    ages5           <- ages.range[ages.range %% 5 == 0]
    # open device 
    pdf(filename)
    for (i in 1:length(years)){
        
        # really, birth matrix just for blanking out cells 
        Bmat                            <- makeBlockAgeGroups(reshape2:::acast(B[[years.char[i]]], FAGE~MAGE, sum, value.var = "Births"), N)
        Em                              <- makeVectorAgeGroups(with(E, Male[Year == years[i]]), ages = 0:110, N = N)
        Ef                              <- makeVectorAgeGroups(with(E, Female[Year == years[i]]), ages = 0:110, N = N)
                      
        R.no.show                       <- Bmat < min.events
        mat2plot                        <- log(Ef %*% (1 / t(Em))) # Identical to log(mFx/fFx)
        dimnames(mat2plot)              <- dimnames(Bmat)
        mat2plot[is.nan(mat2plot)]      <- NA
        mat2plot[is.infinite(mat2plot)] <- NA
        mat2plot[R.no.show ]            <- NA
        mat2plot[mat2plot < zlim[1]]    <- zlim[1]
        mat2plot[mat2plot > zlim[2]]    <- zlim[2]

        fields:::image.plot(x = ages.range.N + N / 2, y = ages.range.N + N / 2, 
                    mat2plot[as.character(ages.range.N), as.character(ages.range.N)], 
                    breaks = brks, 
                    col = rev(colfun(length(brks) - 1)),
                    ylab = "mother age", 
                    xlab = "father age",
                    main = paste0("log(male / female) fertility rates by ", N, "x", N, " age, ",
                            years[i], "\nonly cells based on at least ", min.events, " births shown"),
                    sub = "(identical to log ratio of female to male exposure)",
                    useRaster = TRUE,
                    asp = 1,
                    xlim = age.range,
                    ylim = age.range,
                    zlim = zlim,
                    axes = FALSE
            )
        
        image(x = ages.range.N + N / 2, y = ages.range.N + N / 2, 
                is.na(mat2plot[as.character(ages.range.N), as.character(ages.range.N)]), 
                col = c(NA, gray(.9)), add = TRUE, useRaster = TRUE)
        abline(h = ages5, col = "white")
        abline(v = ages5, col = "white")
        abline(a = 0, b = 1)
        text(min(ages5), ages5, ages5, xpd = TRUE, pos = 2)
        text(ages5, min(ages5), ages5, xpd = TRUE, pos = 1)
    }
    dev.off()
}

# where E is a long df and B is in a list...
B <- local(get(load("/home/triffe/git/Dissertation/DISSERTATION/DATA/ESbirths/ESbirths.Rdata")))
E <- local(get(load("/home/triffe/git/Dissertation/DISSERTATION/DATA/Exposures/ESexp.Rdata")))
B <- B[-length(B)] # remove 2010 births for now

RateRatioSurface(B, E, N = 1, 
        min.events = 40, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateRatioSurfaces/ES_R_1x1.pdf", 
        age.range = c(10,65), zlim = c(-.6, .6))
RateRatioSurface(B, E, N = 2, 
        min.events = 40, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateRatioSurfaces/ES_R_2x2.pdf", 
        age.range = c(10,65), zlim = c(-.6, .6))
RateRatioSurface(B, E, N = 3, 
        min.events = 40, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateRatioSurfaces/ES_R_3x3.pdf", 
        age.range = c(10,65), zlim = c(-.6, .6))

# For US, slight change in data..
###############################################################################

RateRatioSurface <- function(B, E, N = 1, min.events = 40, filename, age.range = c(10, 65), zlim = c(-.6, .6)){
    if (missing(filename)){
        stop("\nfirst give the output a name, including the suffix '.pdf'\n")
    }
    source("/home/triffe/git/Dissertation/DISSERTATION/R/UtilityFunctions.R")
    # color ramp function
    colfun          <- grDevices:::colorRampPalette(RColorBrewer:::brewer.pal(9, "RdBu"), space = "Lab")
    # years vectors
    years.char      <- names(B) # just make sure B and E are of same length
    years           <- as.integer(years.char)
    # value breaks
    brks            <- seq(zlim[1], zlim[2], by = .02)
    ages.all        <- 0:110
    ages.range      <- ages.all[ages.all >= min(age.range) & ages.all <= max(age.range)]
    ages.range.N    <- unique(ages.range - ages.range %% N)
    ages5           <- ages.range[ages.range %% 5 == 0]
    # open device 
    pdf(filename)
    for (i in 1:length(years)){
        
        # really, birth matrix just for blanking out cells 
        Bmat                            <- makeBlockAgeGroups(B[[years.char[i]]], N)
        Em                              <- makeVectorAgeGroups(with(E, Male[Year == years[i]]), ages = 0:110, N = N)
        Ef                              <- makeVectorAgeGroups(with(E, Female[Year == years[i]]), ages = 0:110, N = N)
        R.no.show                       <- Bmat < min.events
        mat2plot                        <- log(Ef %*% (1 / t(Em))) # Identical to log(mFx/fFx)
        dimnames(mat2plot)              <- dimnames(Bmat)
        mat2plot[is.nan(mat2plot)]      <- NA
        mat2plot[is.infinite(mat2plot)] <- NA
        mat2plot[R.no.show ]            <- NA
        mat2plot[mat2plot < zlim[1]]    <- zlim[1]
        mat2plot[mat2plot > zlim[2]]    <- zlim[2]
        
        fields:::image.plot(x = ages.range.N + N / 2, y = ages.range.N + N / 2, 
                mat2plot[as.character(ages.range.N), as.character(ages.range.N)], 
                breaks = brks, 
                col = rev(colfun(length(brks) - 1)),
                ylab = "mother age", 
                xlab = "father age",
                main = paste0("log(male / female) fertility rates by ", N, "x", N, " age, ",
                        years[i], "\nonly cells based on at least ", min.events, " births shown"),
                sub = "(identical to log ratio of female to male exposure)",
                useRaster = TRUE,
                asp = 1,
                xlim = age.range,
                ylim = age.range,
                zlim = zlim,
                axes = FALSE
        )
        
        image(x = ages.range.N + N / 2, y = ages.range.N + N / 2, 
                is.na(mat2plot[as.character(ages.range.N), as.character(ages.range.N)]), 
                col = c(NA, gray(.9)), add = TRUE, useRaster = TRUE)
        abline(h = ages5, col = "white")
        abline(v = ages5, col = "white")
        abline(a = 0, b = 1)
        text(min(ages5), ages5, ages5, xpd = TRUE, pos = 2)
        text(ages5, min(ages5), ages5, xpd = TRUE, pos = 1)
    }
    dev.off()
}

# where E is a long df and B is in a list...
plot(20:40, with(E, Male[Year == 1980])[21:41], type = 'l')
lines(20:40, with(E, Male[Year == 1981])[21:41], col = "red", lty = 2)
lines(20:40, with(E, Male[Year == 1982])[21:41], col = "blue", lty = 2)
abline(v=30)
B <- local(get(load("/home/triffe/git/Dissertation/DISSERTATION/DATA/USbirths/USBxy0_110.Rdata")))
E <- local(get(load("/home/triffe/git/Dissertation/DISSERTATION/DATA/Exposures/USexp.Rdata")))
names(B) <- 1969:2010
B <- B[-length(B)]

RateRatioSurface(B, E, N = 1, 
        min.events = 40, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateRatioSurfaces/US_R_1x1.pdf", 
        age.range = c(10,65), zlim = c(-.6, .6))
RateRatioSurface(B, E, N = 2, 
        min.events = 40, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateRatioSurfaces/US_R_2x2.pdf", 
        age.range = c(10,65), zlim = c(-.6, .6))
RateRatioSurface(B, E, N = 3, 
        min.events = 40, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateRatioSurfaces/US_R_3x3.pdf", 
        age.range = c(10,65), zlim = c(-.6, .6))
