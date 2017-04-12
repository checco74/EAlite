\name{draw_qq_plots}
\alias{draw_qq_plots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to draw QQ-plots
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
draw_qq_plots(mydata.full = NULL, myannot.full = NULL, myparameters = init_parameters(), myfilter.full = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mydata.full}{
%%     ~~Describe \code{mydata.full} here~~
}
  \item{myannot.full}{
%%     ~~Describe \code{myannot.full} here~~
}
  \item{myparameters}{
%%     ~~Describe \code{myparameters} here~~
}
  \item{myfilter.full}{
%%     ~~Describe \code{myfilter.full} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
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
        if (myparameters$verbose) 
            cat("setting temporary data storage tables..\n")
        myorder = order(mydata.full)
        mydata <- mydata.full[myorder]
        myannot <- myannot.full[myorder]
        myfilter = NULL
        if (is.null(myfilter.full)) {
            myfilter = array(TRUE, dim = c(length(mydata), 1))
        }
        else {
            myfilter = myfilter.full[myorder, ]
        }
        mymask <- !is.na(mydata) & !is.na(myannot)
        mydata <- mydata[mymask]
        myannot <- myannot[mymask]
        myfilter = myfilter[mymask, ]
        myorder <- myorder[mymask]
        neglog10_P_all <- -log10(mydata)
        if (myparameters$verbose) {
            print("draw_qq_plots() neglog10_P_all:")
            print(summary(neglog10_P_all))
        }
        n_all = length(neglog10_P_all)
        histobreaks = seq(min(0, floor(min(neglog10_P_all, na.rm = TRUE))), 
            ceiling(max(neglog10_P_all, na.rm = TRUE)), length.out = 1000)
        histotemp_all <- hist(neglog10_P_all, breaks = histobreaks, 
            plot = FALSE)
        cdftemp_all <- cumsum(histotemp_all$counts)
        datafilter <- histotemp_all$mids < -log10(gws.level) | 
            histotemp_all$counts > 1
        datapoints <- histotemp_all$mids[datafilter]
        binconftemp_all <- binconf(cdftemp_all, n_all, method = "exact")
        binconftemp_all <- binconftemp_all[datafilter, ]
        cat("plotting masters..\n")
        qqplot(-log10(1 - binconftemp_all[, 1]), datapoints, 
            type = plottype, lwd = 2, pch = 16, cex = 0.5, col = "#777777", 
            xlab = "Empirical -log10(p)", ylab = "Nominal -log10(p)", 
            xlim = c(0, 6), ylim = c(0, 8))
        if (myparameters$do_draw_cis) {
            points(-log10(1 - binconftemp_all[, 2]), datapoints, 
                type = plottype, lty = 2)
            points(-log10(1 - binconftemp_all[, 3]), datapoints, 
                type = plottype, lty = 2)
        }
        cat("plotting additional points..\n")
        myannotq <- factorize_annotation(myannot, myparameters$mybreaks, 
            myparameters$numof.breaks)
        myannotq_levels <- levels(as.factor(myannotq))
        mycolors <- array(dim = length(myannotq_levels))
        for (k in 1:length(myannotq_levels)) mycolors[k] <- mycolorfunction(k, 
            length(myannotq_levels), myparameters$figure.colors)
        for (k in 1:length(myannotq_levels)) {
            cat("plotting ", "I", k, " -- color=", mycolors[k], 
                "..\n", sep = "")
            neglog10_P <- -log10(mydata[myannotq == myannotq_levels[k]])
            neglog10_P <- neglog10_P[!is.na(neglog10_P)]
            n = length(neglog10_P)
            histotemp <- hist(neglog10_P, breaks = histobreaks, 
                plot = FALSE)
            cdftemp <- cumsum(histotemp$counts)
            binconftemp <- binconf(cdftemp, n, method = "exact")
            binconftemp <- binconftemp[datafilter, ]
            if (length(-log10(1 - binconftemp[, 1])) > 0 && length(datapoints) > 
                0 && length(-log10(1 - binconftemp[, 1])) == 
                length(datapoints)) {
                points(-log10(1 - binconftemp[, 1]), datapoints, 
                  type = plottype, lwd = 2, col = mycolors[k], 
                  pch = 16)
                if (myparameters$do_draw_cis) {
                  points(-log10(1 - binconftemp[, 2]), datapoints, 
                    type = plottype, col = mycolors[k], lty = 2)
                  points(-log10(1 - binconftemp[, 3]), datapoints, 
                    type = plottype, col = mycolors[k], lty = 2)
                }
            }
            else cat("warning: degenerate", "I", k, "\n")
        }
        cat("writing legends..\n")
        abline(a = 0, b = 1, col = "lightgray", lwd = 2, lty = 2)
        abline(a = -log10(gws.level), b = 0, col = "lightblue", 
            lwd = 2, lty = 3)
        legend("topleft", c("all SNPs", myannotq_levels), lwd = c(rep(2, 
            length(myannotq_levels) + 1)), bg = "white", col = c("#777777", 
            mycolors[1:length(myannotq_levels)]))
        legend("bottomright", c("Expected under null", paste("p = ", 
            format(gws.level))), lty = c(2, 3), lwd = 2, bg = "white", 
            col = c("lightgray", "lightblue"))
    }
    else {
        return(NULL)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
