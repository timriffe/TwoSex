
# Author: triffe
###############################################################################

# function needs B and E, good to simply have them available in workspace N <- 3;i<-1
RateSurface <- function(B, E, N = 1, min.events = 40, filename, age.range = c(10,65), sex = "m", zlim = c(0, .035), .log = FALSE){
    if (missing(filename)){
        stop("\nfirst give the output a name, including the suffix '.pdf'\n")
    }
    source("/home/triffe/git/Dissertation/DISSERTATION/R/UtilityFunctions.R")
    # color ramp function
    colfun          <- grDevices:::colorRampPalette(RColorBrewer:::brewer.pal(9, "YlOrRd"), space = "Lab")
    # years vectors
    years.char      <- names(B) # just make sure B and E are of same length
    years           <- as.integer(years.char)
    # value breaks
    brks            <- seq(0, zlim[2], by = .001)
    ages.all        <- 0:110
    ages.range      <- ages.all[ages.all >= min(age.range) & ages.all <= max(age.range)]
    ages.range.N    <- unique(ages.range - ages.range %% N)
    ages5           <- ages.range[ages.range %% 5 == 0]
    # open device 
    pdf(filename)
    for (i in 1:length(years)){
        Bmat                           <- makeBlockAgeGroups(reshape2:::acast(B[[years.char[i]]], FAGE~MAGE, sum, value.var = "Births"), N)
         
        # FAGE in rows = need to transpose for division if sex == "f"
        if (sex == "m"){
            Fx                             <- Bmat / 
                            makeVectorAgeGroups(with(E, Male[Year == years[i]]), ages = 0:110, N = N)
        }
        if (sex == "f"){
            Fx                             <- t(t(Bmat) / 
                            makeVectorAgeGroups(with(E, Female[Year == years[i]]), ages = 0:110, N = N)) 
        }
      
        R.no.show                       <- Bmat < min.events
        mat2plot                        <- Fx
        mat2plot[is.nan(mat2plot)]      <- NA
        mat2plot[is.infinite(mat2plot)] <- NA
        mat2plot[R.no.show ]            <- NA
        mat2plot[mat2plot < 0]          <- 0
        # if logged, forget zlim, use special ticks and brks 
        if (.log){
            ticks <- c(.00001, .00002, .00005, .0001, .0002, .0005, 0.001, .002, .005, .01,.02, .05, .1)
            brks  <- seq(-11.5, -2.5, by = .05)
            fields:::image.plot(x = ages.range.N + N / 2, y = ages.range.N + N / 2, 
                    log(mat2plot[as.character(ages.range.N), as.character(ages.range.N)]), 
                    breaks = brks, 
                    col = colfun(length(brks) - 1),
                    ylab = "mother age", 
                    xlab = "father age",
                    main = paste0(ifelse(sex == "m", "Male", "Female"), "log fertility rates by ",N,"x",N," age, ",
                            years[i],"\nonly cells based on at least ", min.events, " births shown"),
                    useRaster = TRUE,
                    asp = 1,
                    xlim = age.range,
                    ylim = age.range,
                    zlim = range(brks),
                    axes = FALSE,
                    axis.args = list(at = log(ticks), labels = ticks)
            )
        } else {
            mat2plot[mat2plot > zlim[2]]    <- zlim[2]
            fields:::image.plot(x = ages.range.N + N / 2, y = ages.range.N + N / 2, 
                    mat2plot[as.character(ages.range.N), as.character(ages.range.N)], 
                    breaks = brks, 
                    col = colfun(length(brks) - 1),
                    ylab = "mother age", 
                    xlab = "father age",
                    main = paste0(ifelse(sex == "m", "Male", "Female"), " fertility rates by ",N,"x",N," age, ",
                            years[i],"\nonly cells based on at least ", min.events, " births shown"),
                    useRaster = TRUE,
                    asp = 1,
                    xlim = age.range,
                    ylim = age.range,
                    zlim = zlim,
                    axes = FALSE
            #legend.args = list(at = brks[brks %% .1 == 0], labels = brks[brks %% .1 == 0]))
            )
        }
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
RateSurface(B, E, N = 1, 
        min.events = 30, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/ESmFx_1x1.pdf",
        age.range = c(10, 65), sex = "m", zlim = c(0,.04))
RateSurface(B, E, N = 1, 
        min.events = 30, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/ESfFx_1x1.pdf",
        age.range = c(10, 65), sex = "f", zlim = c(0,.04))
RateSurface(B, E, N = 2, 
        min.events = 50, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/ESmFx_2x2.pdf",
        age.range = c(10, 65), sex = "m", zlim = c(0,.06))
RateSurface(B, E, N = 2, 
        min.events = 50, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/ESfFx_2x2.pdf",
        age.range = c(10, 65), sex = "f", zlim = c(0,.06))
RateSurface(B, E, N = 3, 
        min.events = 100, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/ESmFx_3x3.pdf",
        age.range = c(10, 65), sex = "m", zlim = c(0,.08))
RateSurface(B, E, N = 3, 
        min.events = 100, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/ESfFx_3x3.pdf",
        age.range = c(10, 65), sex = "f", zlim = c(0,.08))

# log rate surfaces
RateSurface(B, E, N = 1, 
        min.events = 30, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/ESmlogFx_1x1.pdf",
        age.range = c(10, 65), sex = "m", .log = TRUE)
RateSurface(B, E, N = 1, 
        min.events = 30, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/ESflogFx_1x1.pdf",
        age.range = c(10, 65), sex = "f", .log = TRUE)
RateSurface(B, E, N = 2, 
        min.events = 50, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/ESmlogFx_2x2.pdf",
        age.range = c(10, 65), sex = "m", .log = TRUE)
RateSurface(B, E, N = 2, 
        min.events = 50, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/ESflogFx_2x2.pdf",
        age.range = c(10, 65), sex = "f", .log = TRUE)
RateSurface(B, E, N = 3, 
        min.events = 100, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/ESmlogFx_3x3.pdf",
        age.range = c(10, 65), sex = "m", .log = TRUE)
RateSurface(B, E, N = 3, 
        min.events = 100, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/ESflogFx_3x3.pdf",
        age.range = c(10, 65), sex = "f", .log = TRUE)

# ---------------------------------------------------------------------------------
# try on US data- not sure about dims
RateSurface <- function(B, E, N = 1, min.events = 40, filename, age.range = c(10,65), sex = "m", zlim = c(0, .035), .log = FALSE){
    if (missing(filename)){
        stop("\nfirst give the output a name, including the suffix '.pdf'\n")
    }
    source("/home/triffe/git/Dissertation/DISSERTATION/R/UtilityFunctions.R")
    # color ramp function
    colfun          <- grDevices:::colorRampPalette(RColorBrewer:::brewer.pal(9, "YlOrRd"), space = "Lab")
    # years vectors
    years.char      <- names(B) # just make sure B and E are of same length
    years           <- as.integer(years.char)
    # value breaks
    brks            <- seq(0, zlim[2], by = .001)
    ages.all        <- 0:110
    ages.range      <- ages.all[ages.all >= min(age.range) & ages.all <= max(age.range)]
    ages.range.N    <- unique(ages.range - ages.range %% N)
    ages5           <- ages.range[ages.range %% 5 == 0]
    # open device 
    pdf(filename)
    for (i in 1:length(years)){
        Bmat                           <- makeBlockAgeGroups(B[[years.char[i]]], N)
        
        # FAGE in rows = need to transpose for division if sex == "f"
        if (sex == "m"){
            Fx                             <- Bmat / 
                    makeVectorAgeGroups(with(E, Male[Year == years[i]]), ages = 0:110, N = N)
        }
        if (sex == "f"){
            Fx                             <- t(t(Bmat) / 
                            makeVectorAgeGroups(with(E, Female[Year == years[i]]), ages = 0:110, N = N)) 
        }
        
        R.no.show                       <- Bmat < min.events
        mat2plot                        <- Fx
        mat2plot[is.nan(mat2plot)]      <- NA
        mat2plot[is.infinite(mat2plot)] <- NA
        mat2plot[R.no.show ]            <- NA
        mat2plot[mat2plot < 0]          <- 0
        # if logged, forget zlim, use special ticks and brks 
        if (.log){
            ticks <- c(.00001, .00002, .00005, .0001, .0002, .0005, 0.001, .002, .005, .01,.02, .05, .1)
            brks  <- seq(-11.5, -2.5, by = .05)
            fields:::image.plot(x = ages.range.N + N / 2, y = ages.range.N + N / 2, 
                    log(mat2plot[as.character(ages.range.N), as.character(ages.range.N)]), 
                    breaks = brks, 
                    col = colfun(length(brks) - 1),
                    ylab = "mother age", 
                    xlab = "father age",
                    main = paste0(ifelse(sex == "m", "Male", "Female"), "log fertility rates by ",N,"x",N," age, ",
                            years[i],"\nonly cells based on at least ", min.events, " births shown"),
                    useRaster = TRUE,
                    asp = 1,
                    xlim = age.range,
                    ylim = age.range,
                    zlim = range(brks),
                    axes = FALSE,
                    axis.args = list(at = log(ticks), labels = ticks)
            )
        } else {
            mat2plot[mat2plot > zlim[2]]    <- zlim[2]
            fields:::image.plot(x = ages.range.N + N / 2, y = ages.range.N + N / 2, 
                    mat2plot[as.character(ages.range.N), as.character(ages.range.N)], 
                    breaks = brks, 
                    col = colfun(length(brks) - 1),
                    ylab = "mother age", 
                    xlab = "father age",
                    main = paste0(ifelse(sex == "m", "Male", "Female"), " fertility rates by ",N,"x",N," age, ",
                            years[i],"\nonly cells based on at least ", min.events, " births shown"),
                    useRaster = TRUE,
                    asp = 1,
                    xlim = age.range,
                    ylim = age.range,
                    zlim = zlim,
                    axes = FALSE
            #legend.args = list(at = brks[brks %% .1 == 0], labels = brks[brks %% .1 == 0]))
            )
        }
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

B <- local(get(load("/home/triffe/git/Dissertation/DISSERTATION/DATA/USbirths/USBxy0_110.Rdata")))
E <- local(get(load("/home/triffe/git/Dissertation/DISSERTATION/DATA/Exposures/USexp.Rdata")))
names(B) <- 1969:2010
B <- B[-length(B)]

RateSurface(B, E, N = 1, 
        min.events = 30, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/USmFx_1x1.pdf",
        age.range = c(10, 65), sex = "m", zlim = c(0,.04))
RateSurface(B, E, N = 1, 
        min.events = 30, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/USfFx_1x1.pdf",
        age.range = c(10, 65), sex = "f", zlim = c(0,.04))
RateSurface(B, E, N = 2, 
        min.events = 50, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/USmFx_2x2.pdf",
        age.range = c(10, 65), sex = "m", zlim = c(0,.06))
RateSurface(B, E, N = 2, 
        min.events = 50, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/USfFx_2x2.pdf",
        age.range = c(10, 65), sex = "f", zlim = c(0,.06))
RateSurface(B, E, N = 3, 
        min.events = 100, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/USmFx_3x3.pdf",
        age.range = c(10, 65), sex = "m", zlim = c(0,.08))
RateSurface(B, E, N = 3, 
        min.events = 100, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/USfFx_3x3.pdf",
        age.range = c(10, 65), sex = "f", zlim = c(0,.08))

# log rate surfaces
RateSurface(B, E, N = 1, 
        min.events = 30, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/USmlogFx_1x1.pdf",
        age.range = c(10, 65), sex = "m", .log = TRUE)
RateSurface(B, E, N = 1, 
        min.events = 30, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/USflogFx_1x1.pdf",
        age.range = c(10, 65), sex = "f", .log = TRUE)
RateSurface(B, E, N = 2, 
        min.events = 50, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/USmlogFx_2x2.pdf",
        age.range = c(10, 65), sex = "m", .log = TRUE)
RateSurface(B, E, N = 2, 
        min.events = 50, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/USflogFx_2x2.pdf",
        age.range = c(10, 65), sex = "f", .log = TRUE)
RateSurface(B, E, N = 3, 
        min.events = 100, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/USmlogFx_3x3.pdf",
        age.range = c(10, 65), sex = "m", .log = TRUE)
RateSurface(B, E, N = 3, 
        min.events = 100, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateSurfaces/USflogFx_3x3.pdf",
        age.range = c(10, 65), sex = "f", .log = TRUE)
