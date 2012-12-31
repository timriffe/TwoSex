
# Author: triffe
###############################################################################

B <- local(get(load("/home/triffe/git/Dissertation/DISSERTATION/DATA/ESbirths/ESbirths.Rdata")))
E <- local(get(load("/home/triffe/git/Dissertation/DISSERTATION/DATA/Exposures/ESexp.Rdata")))
B <- B[-length(B)] # remove 2010 births for now
# where E is a long df and B is in a list...

# function needs B and E, good to simply have them available in workspace
RateChangeDiagnostic <- function(B, E, N = 1, min.events = 40, filename, age.range = c(10,65), sex = "m"){
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
    brks            <- seq(-1, 1, by = .01)
    ages.all        <- 0:110
    ages.N          <- unique(ages.all - ages.all %% N)
    age.range.N     <- ages.N[ages.N >= min(age.range) & ages.N <= max(age.range)]
    # open device
    pdf(filename)
    for (i in 1:(length(years) - 1)){
        Bmat1                           <- makeBlockAgeGroups(reshape2:::acast(B[[years.char[i]]], FAGE~MAGE, sum, value.var = "Births"), N)
        Bmat2                           <- makeBlockAgeGroups(reshape2:::acast(B[[years.char[i+1]]], FAGE~MAGE, sum, value.var = "Births"), N)
        
        # FAGE in rows = need to transpose for division if sex == "f"
        if (sex == "m"){
            Fx1                             <- Bmat1 / makeVectorAgeGroups(with(E, Male[Year == years[i]]), ages = 0:110, N = N)
            Fx2                             <- Bmat2 / makeVectorAgeGroups(with(E, Male[Year == years[i+1]]), N = N)
        }
        if (sex == "f"){
            Fx1                             <- t(t(Bmat1) / makeVectorAgeGroups(with(E, Female[Year == years[i]]), ages = 0:110, N = N)) 
            Fx2                             <- t(t(Bmat2) / makeVectorAgeGroups(with(E, Female[Year == years[i+1]]), N = N))
        }
  
        R.no.show                       <- (Bmat1 + Bmat2) < min.events
        chgrate                         <- Fx2 - Fx1
        Rel                             <- chgrate / ((Fx2 + Fx1) / 2)
        mat2plot                        <- Rel
        mat2plot[is.nan(mat2plot)]      <- NA
        mat2plot[is.infinite(mat2plot)] <- NA
        mat2plot[R.no.show ]            <- NA
        mat2plot[mat2plot < -1]         <- -1
        mat2plot[mat2plot > 1]          <- 1
        
        fields:::image.plot(x = age.range.N + N / 2, y = age.range.N + N / 2, 
                mat2plot[as.character(age.range.N), as.character(age.range.N)], 
                breaks = brks, 
                col = rev(colfun(length(brks) - 1)),
                ylab = "mother age", 
                xlab = "father age",
                main = paste0("Relative change in ",ifelse(sex == "m", "male", "female"), " fertility rates by age ",
                        years[i], "-", years[i + 1]),
                useRaster = TRUE,
                asp = 1,
                xlim = age.range,
                ylim = age.range,
                zlim = c(-1,1)
        #legend.args = list(at = brks[brks %% .1 == 0], labels = brks[brks %% .1 == 0]))
        )
        abline(a = 0, b = 1)
    }
    dev.off()
}



RateChangeDiagnostic(B, E, N = 1, 
        min.events = 30, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateChangeSurfaces/ESmR_1x1.pdf",
        age.range = c(10, 65),
        sex = "m")
RateChangeDiagnostic(B, E, N = 1, 
        min.events = 30, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateChangeSurfaces/ESfR_1x1.pdf",
        age.range = c(10, 65),
        sex = "f")
RateChangeDiagnostic(B, E, N = 2, 
        min.events = 50, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateChangeSurfaces/ESmR_2x2.pdf",
        age.range = c(10, 65),
        sex = "m")
RateChangeDiagnostic(B, E, N = 2, 
        min.events = 50, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateChangeSurfaces/ESfR_2x2.pdf",
        age.range = c(10, 65),
        sex = "f")
RateChangeDiagnostic(B, E, N = 3, 
        min.events = 100, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateChangeSurfaces/ESmR_3x3.pdf",
        age.range = c(10, 65),
        sex = "m")
RateChangeDiagnostic(B, E, N = 3, 
        min.events = 100, 
        filename = "/home/triffe/git/Dissertation/DISSERTATION/DiagnosticPlots/RateChangeSurfaces/ESfR_3x3.pdf",
        age.range = c(10, 65),
        sex = "f")

# 1) make rate surfaces
# 2) make relative rate surfaces
# 3) try to decompose changes in rates?
#    a) via B and E for males and females separately?
#    b) via two-sex mean assumptions?? Can decomposition via 2 sex mean assumptions
#       be a viable way to visually determine an optimal mean function?

