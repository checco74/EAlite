draw_fold_plots <-
function (mydata.full = NULL, myannot.full = NULL, myparameters = init_parameters(),
    myfilter.full = NULL)
{
    plottype = "l"
	defcolor = "#999999"
    if (is.null(mydata.full) || is.null(myannot.full)) {
        return(NULL)
    }
    pd <- get_plot_data(pvals=mydata.full, 
                        annot=myannot.full, 
                        breaks=myparameters$mybreaks, 
                        numof.breaks=myparameters$numof.breaks
                        verbose=myparameters$verbose)
    cat("plotting masters..\n")
    plot(pd$datapoints, rep(1, max(10, length(pd$datapoints))),
        type = plottype, lty = 2, lwd = 1.5, pch = 16, cex = 0.5, col = defcolor,
        log = "y", xlab = "Nominal -log10(p)", ylab = "Fold enrichment",
        xlim = c(0, 8), ylim = c(0.1, 10))
    cat("plotting additional points..\n")
    mycolors <- array(dim = length(pd$folds))
    for (k in 1:length(pd$folds)) mycolors[k] <-
        mycolorfunction(k, length(pd$folds), myparameters$figure.colors)
    for (k in 1:length(pd$folds)) {
        cat("plotting ", 'I', k, " -- color=", mycolors[k], "..\n", sep = "")
        points(pd$datapoints, pd$folds[[k]]$y, type = plottype,
            lwd = 2, col = mycolors[k], pch = 16)
        # draw confidence intervals?
        if (myparameters$do_draw_cis) {
            points(pd$datapoints, pd$folds[[k]]$ci0, type = plottype, col = mycolors[k], lty = 2)
            points(pd$datapoints, pd$folds[[k]]$ci1, type = plottype, col = mycolors[k], lty = 2)
        }
    }
    cat("writing legends..\n")
    annot_names <- unlist(lapply(pd$fold, function(x){x$name}))
    legend("topleft", c("all SNPs", annot_names),
        lwd = c(rep(2, length(annot_names) + 1)), bg = "white",
        col = c(defcolor, mycolors[1:length(annot_names)]))
}


get_plot_data <- 
function(pvals, annot, filter=NULL, breaks='auto', numof.breaks=4, verbose=T) 
{
    # data structure to store plot data
    pd <- NULL
    pd$datapoints <- NULL   # x values
    pd$folds      <- list() # list of y values
    # example for list elements  
    # pd$folds[[i]]$y      <- NULL # y value
    # pd$folds[[i]]$ci0    <- NULL # confidence interval starts
    # pd$folds[[i]]$ci1    <- NULL # confidence interval ends
    # pd$folds[[i]]$name   <- NULL # name of annotation

    library(Hmisc, quietly = T, warn.conflicts = F)
    neglog10_P_threshold <- 24
    tiny <- 1e-72
    cat("p-value genome-wide significance level is set to: ")
    gws.level <- max(5e-08, 0.05/sum(!is.na(pvals)))
    cat(gws.level, "\n")
    if (verbose) cat("setting temporary data storage tables..\n")
    myorder = order(pvals)
    mydata <- pvals[myorder]
    myannot <- annot[myorder]
    myfilter = NULL
    if ( is.null(filter) ) {
        myfilter = array( TRUE, dim = c( length( mydata ), 1 ) )
    } else { myfilter = filter[myorder,] }
    # filter NA values
    mymask <- !is.na(mydata) & !is.na(myannot)
    mydata <- mydata[mymask]
    myannot <- myannot[mymask]
    myfilter = myfilter[mymask,]
    myorder <- myorder[mymask]
    # convert p-values to negative log scale [0,inf]
    neglog10_P_all <- -log10(mydata)
    if (verbose) {
        print("draw_qq_plots() neglog10_P_all:")
        print(summary(neglog10_P_all))
    }
    n_all = length(neglog10_P_all)
    # bin data for plots, so they don't get too heavy
    num_datapoints <- 1000
    histobreaks = seq(
        min(0, floor(min(neglog10_P_all, na.rm = TRUE))),
        ceiling(max(neglog10_P_all, na.rm = TRUE)),
        length.out = num_datapoints
    )
    # make hist for the whole data set (this will span the whole possible range of values)
    histotemp_all <- hist(neglog10_P_all, breaks = histobreaks, plot = FALSE)
    # cumulative sum of histogram bin counts (bin1, bin1+bin2, bin1+bin2+bin3, ...)
    # this is cumulative on negative logarithm of p, i.e. from uninteresting p values to
    # interesting p values
    cdftemp_all <- cumsum(histotemp_all$counts)
    datafilter <- histotemp_all$mids < -log10(gws.level) | histotemp_all$counts > 1
    pd$datapoints <- histotemp_all$mids[datafilter]
    # data struct to store plot data
    binconftemp_all <- binconf(cdftemp_all, n_all, method = "exact")
    binconftemp_all <- binconftemp_all[datafilter, ]
    myannotq <- factorize_annotation(myannot, breaks, numof.breaks)
    myannotq_levels <- levels(as.factor(myannotq))
    for (k in 1:length(myannotq_levels)) {
        fold <- NULL
        fold$name <- as.character(myannotq_levels[k])
        cat("calc ", 'I', k, "..\n", sep = "")
        neglog10_P <- -log10(mydata[myannotq == myannotq_levels[k]])
        neglog10_P <- neglog10_P[!is.na(neglog10_P)]
        n = length(neglog10_P)
        # make hist for the stratum
        histotemp <- hist(neglog10_P, breaks = histobreaks, plot = FALSE)
        # cumulative sum of histogram bin counts (bin1, bin1+bin2, bin1+bin2+bin3, ...)
        # this is cumulative on negative logarithm of p, i.e. from uninteresting p values to
        # interesting p values
        cdftemp <- cumsum(histotemp$counts)
        binconftemp <- binconf(cdftemp, n, method = "exact")
        binconftemp <- binconftemp[datafilter, ]
        # we plot (1 - binconftemp) ratios because we want to emphasize the fold enrichment on
        # the interesting end of the p value spectrum
        if (length(-log10(1 - binconftemp[, 1])) > 0 && length(pd$datapoints) > 0 && 
            length(-log10(1 - binconftemp[, 1])) == length(pd$datapoints)) {
            fold$y <- (1 - binconftemp[, 1])/(1 - binconftemp_all[,1] + tiny)
            # calc confidence intervals
            fold$ci0 <- (1 - binconftemp[, 2])/(1 - binconftemp_all[,2] + tiny)
            fold$ci1 <- (1 - binconftemp[, 3])/(1 - binconftemp_all[,3] + tiny)
        }
        else cat("warning: degenerate", 'I', k, "\n")
        pd$folds[[k]] <- fold
    }
    pd
}
