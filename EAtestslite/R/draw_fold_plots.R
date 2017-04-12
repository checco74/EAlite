draw_fold_plots <-
function (mydata.full = NULL, myannot.full = NULL, myparameters = init_parameters(),
    myfilter.full = NULL)
{
    plottype = "l"
    library(Hmisc, quietly = T, warn.conflicts = F)
    if (!is.null(mydata.full) && !is.null(myannot.full)) {
        neglog10_P_threshold <- 24
        tiny <- 1e-72
        cat("p-value genome-wide significance level is set to: ")
        gws.level <- max(5e-08, 0.05/sum(!is.na(mydata.full)))
        cat(gws.level, "\n")
        if (myparameters$verbose) cat("setting temporary data storage tables..\n")
        myorder = order(mydata.full)
        mydata <- mydata.full[myorder]
        myannot <- myannot.full[myorder]
        myfilter = NULL
        if ( is.null( myfilter.full ) ) {
            myfilter = array( TRUE, dim = c( length( mydata ), 1 ) )
        } else { myfilter = myfilter.full[myorder,] }
        mymask <- !is.na(mydata) & !is.na(myannot)
        mydata <- mydata[mymask]
        myannot <- myannot[mymask]
        myfilter = myfilter[mymask,]
        myorder <- myorder[mymask]
        neglog10_P_all <- -log10(mydata)
        if (myparameters$verbose) {
            print("draw_qq_plots() neglog10_P_all:")
            print(summary(neglog10_P_all))
        }
        n_all = length(neglog10_P_all)
        histobreaks = seq(
            min(0, floor(min(neglog10_P_all, na.rm = TRUE))),
            ceiling(max(neglog10_P_all, na.rm = TRUE)),
            length.out = 1000
        )
        histotemp_all <- hist(neglog10_P_all, breaks = histobreaks, plot = FALSE)
        cdftemp_all <- cumsum(histotemp_all$counts)
        datafilter <- histotemp_all$mids < -log10(gws.level) | histotemp_all$counts > 1
        datapoints <- histotemp_all$mids[datafilter]
        binconftemp_all <- binconf(cdftemp_all, n_all, method = "exact")
        binconftemp_all <- binconftemp_all[datafilter, ]
        cat("plotting masters..\n")
        plot(datapoints, rep(1, max(10, length(datapoints))),
            type = plottype, lty = 2, lwd = 1.5, pch = 16, cex = 0.5, col = "#777777",
            log = "y", xlab = "Nominal -log10(p)", ylab = "Fold enrichment",
            xlim = c(0, 8), ylim = c(0.1, 10))
        cat("plotting additional points..\n")
        myannotq <- factorize_annotation(myannot, myparameters$mybreaks, myparameters$numof.breaks)
        myannotq_levels <- levels(as.factor(myannotq))
        mycolors <- array(dim = length(myannotq_levels))
        for (k in 1:length(myannotq_levels)) mycolors[k] <-
            mycolorfunction(k, length(myannotq_levels), myparameters$figure.colors)
        for (k in 1:length(myannotq_levels)) {
            cat("plotting ", 'I', k, " -- color=", mycolors[k], "..\n", sep = "")
            neglog10_P <- -log10(mydata[myannotq == myannotq_levels[k]])
            neglog10_P <- neglog10_P[!is.na(neglog10_P)]
            n = length(neglog10_P)
            histotemp <- hist(neglog10_P, breaks = histobreaks, plot = FALSE)
            cdftemp <- cumsum(histotemp$counts)
            binconftemp <- binconf(cdftemp, n, method = "exact")
            binconftemp <- binconftemp[datafilter, ]
            if (length(-log10(1 - binconftemp[, 1])) > 0 && length(datapoints) > 0 && 
                length(-log10(1 - binconftemp[, 1])) == length(datapoints)) {
                foldtemp <- (1 - binconftemp[, 1])/(1 - binconftemp_all[,1] + tiny)
                points(datapoints, foldtemp, type = plottype,
                  lwd = 2, col = mycolors[k], pch = 16)
                if (myparameters$do_draw_cis) {
                  foldtemp_2 <- (1 - binconftemp[, 2])/(1 - binconftemp_all[,2] + tiny)
                  foldtemp_3 <- (1 - binconftemp[, 3])/(1 - binconftemp_all[,3] + tiny)
                  points(datapoints, foldtemp_2, type = plottype, col = mycolors[k], lty = 2)
                  points(datapoints, foldtemp_3, type = plottype, col = mycolors[k], lty = 2)
                }
            }
            else cat("warning: degenerate", 'I', k, "\n")
        }
        cat("writing legends..\n")
        legend("topleft", c("all SNPs", myannotq_levels),
			lwd = c(rep(2, length(myannotq_levels) + 1)), bg = "white",
            col = c("#777777", mycolors[1:length(myannotq_levels)]))
    }
    else {
        return(NULL)
    }
}
