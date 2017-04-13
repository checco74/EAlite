\name{compute_fdr}
\alias{compute_fdr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to compute FDRs
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compute_fdr(mydata.full = NULL, myannot.full = NULL, myparameters = init_parameters(), myfilter.full = NULL, fdr.file = "")
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
  \item{fdr.file}{
%%     ~~Describe \code{fdr.file} here~~
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
    myfilter.full = NULL, fdr.file = "") 
{
    if (!is.null(mydata) && !is.null(myannot)) {
        myorder = order(mydata.full)
        mydata <- mydata.full[myorder]
        myannot <- myannot.full[myorder]
        myfilter = NULL
        if (is.null(myfilter.full)) {
            myfilter = array(TRUE, dim = c(length(mydata.full), 
                1))
        }
        else {
            myfilter = myfilter.full[myorder, ]
        }
        mymask <- !is.na(mydata) & !is.na(myannot)
        mydata <- mydata[mymask]
        myannot <- myannot[mymask]
        myfilter = myfilter[mymask, ]
        myorder <- myorder[mymask]
        myannotq <- factorize_annotation(myannot, myparameters$mybreaks, 
            myparameters$numof.breaks)
        myannotq_levels <- levels(as.factor(myannotq))
        mycolors <- array(dim = length(myannotq_levels))
        for (k in 1:length(myannotq_levels)) mycolors[k] <- mycolorfunction(k, 
            length(myannotq_levels), myparameters$figure.colors)
        mysummary.fdr <- matrix(NA, nrow = 2 * length(myannotq_levels) + 
            2, ncol = 6)
        colnames(mysummary.fdr) <- c("[=]", "[=]*", "[<]", "[<]*", 
            "[>]", "[>]*")
        rownames(mysummary.fdr) <- c(paste(rep(c("FDR", "GWS"), 
            each = length(myannotq_levels)), "I", 1:length(myannotq_levels), 
            sep = "_"), c("FDR_all", "GWS_all"))
        inames = paste("I", 1:length(myannotq_levels), sep = "_")
        mytable.names <- c("all", paste("FDR", inames, sep = "_"), 
            paste("FDR_lt", inames[-(1:2)], sep = "_"), paste0("FDR_gt", 
                inames[-((length(inames) - 2):length(inames))]))
        mytable <- matrix(NA, nrow = length(mydata), ncol = length(mytable.names))
        colnames(mytable) <- mytable.names
        cat("p-value genome-wide significance level is set to: ")
        gws.level <- max(5e-08, myparameters$myFDR/length(mydata))
        cat(gws.level, "\n")
        mysummary.fdr["GWS_all", ] <- c(0, NA, NA, NA, NA, NA)
        for (i in 1:dim(myfilter)[2]) {
            mytempdata <- mydata[myfilter[, i]]
            gws.entry <- sum(mytempdata < gws.level, na.rm = T)
            mysummary.fdr["GWS_all", 1] <- sum(mysummary.fdr["GWS_all", 
                1], gws.entry, na.rm = T)
        }
        marker.list <- c()
        marker.list.novel <- c()
        temp.fdr <- array(1, dim = dim(myfilter))
        fdr.disc <- logical(length(mydata))
        for (i in 1:dim(myfilter)[2]) {
            temp.fdr[myfilter[, i], i] <- p.adjust(mydata[myfilter[, 
                i]], method = "fdr")
            temp.fdr.disc <- myfilter[, i] & !(myorder \%in\% marker.list) & 
                temp.fdr[, i] < myparameters$myFDR
            temp.fdr.novel.disc <- myfilter[, i] & !(myorder \%in\% 
                marker.list.novel) & temp.fdr[, i] < myparameters$myFDR & 
                mydata >= gws.level
            new.entry <- c(sum(temp.fdr.disc, na.rm = T), sum(temp.fdr.novel.disc, 
                na.rm = T), NA, NA, NA, NA)
            mysummary.fdr["FDR_all", ] <- ifelse(i == 1, 0, mysummary.fdr["FDR_all", 
                ]) + new.entry
            marker.list <- c(marker.list, myorder[temp.fdr.disc & 
                myfilter[, i] & !(myorder \%in\% marker.list)])
            marker.list.novel <- c(marker.list.novel, myorder[temp.fdr.novel.disc & 
                myfilter[, i] & !(myorder \%in\% marker.list.novel)])
            fdr.disc = fdr.disc | temp.fdr.disc
        }
        mytable$all <- rowMeans(temp.fdr, na.rm = T)
        if (myparameters$do_draw_fdr_vs_p_plots) {
            mykey <- c(16)
            mykey.color <- c("#777777")
            mykey.text <- c("all SNPs")
            plot.filter <- runif(length(mydata)) < 1000/length(mydata)
            plot(-log10(mydata[plot.filter]), mytable$all[plot.filter], 
                col = "#777777", pch = 16, xlab = "-Log10(p)", 
                ylab = "FDR", xlim = c(0, 12), type = "b")
        }
        for (k in 1:length(myannotq_levels)) {
            annot_level_name <- paste("I", k, sep = "_")
            cat(" ", "I", k, "[", myannotq, "]..\n")
            temp.fdr <- array(NA, dim = dim(myfilter))
            logtemp_k <- myannotq == myannotq_levels[k]
            cat("   ", sum(logtemp_k, na.rm = T), " [", sum(is.na(logtemp_k)), 
                " dropped because NA] SNPs in quantile..\n", 
                sep = "")
            logtemp_k[is.na(logtemp_k)] <- FALSE
            for (i in 1:dim(myfilter)[2]) {
                logtempfilter_k <- logtemp_k & myfilter[, i]
                temp.fdr[logtempfilter_k, i] <- p.adjust(mydata[logtempfilter_k], 
                  method = "fdr")
            }
            mytable[, paste("FDR", annot_level_name, sep = "_")] <- rowMeans(temp.fdr, 
                na.rm = T)
            if (parameters$do_draw_fdr_vs_p_plots) {
                mykey <- append(mykey, 16)
                mykey.color <- append(mykey.color, mycolors[k])
                mykey.text <- append(mykey.text, annot_level_name)
                points(-log10(mydata[logtemp_k & plot.filter]), 
                  mytable[logtemp_k & plot.filter, paste("FDR", 
                    annot_level_name, sep = "_")], col = mycolors[k], 
                  pch = 16, type = "b")
            }
            marker.filter_k <- FALSE
            for (i in 1:dim(myfilter)[2]) {
                logtempfilter_k <- logtemp_k & myfilter[, i] & 
                  !marker.filter_k
                fdrlogtemp <- logtempfilter_k & temp.fdr[, i] < 
                  myparameters$myFDR
                mysummary.fdr[paste("FDR", annot_level_name, 
                  sep = "_"), "[=]"] <- ifelse(i == 1, 0, mysummary.fdr[paste("FDR", 
                  annot_level_name, sep = "_"), "[=]"]) + sum(fdrlogtemp, 
                  na.rm = T)
                marker.filter_k <- marker.filter_k | fdrlogtemp
            }
            marker.filter_k <- FALSE
            for (i in 1:dim(myfilter)[2]) {
                logtempfilter_k <- logtemp_k & myfilter[, i] & 
                  !marker.filter_k
                gwslogtemp <- logtempfilter_k & mydata < myparameters$myFDR/sum(logtemp_k)
                mysummary.fdr[paste("GWS", annot_level_name, 
                  sep = "_"), "[=]"] <- ifelse(i == 1, 0, mysummary.fdr[paste("GWS", 
                  annot_level_name, sep = "_"), "[=]"]) + sum(gwslogtemp, 
                  na.rm = T)
                marker.filter_k <- marker.filter_k | gwslogtemp
            }
            marker.filter_k <- FALSE
            for (i in 1:dim(myfilter)[2]) {
                logtempfilter_k <- logtemp_k & myfilter[, i] & 
                  !marker.filter_k
                fdrlogtemp <- logtempfilter_k & temp.fdr[, i] < 
                  myparameters$myFDR
                tmpcount <- sum(fdrlogtemp & mydata >= gws.level & 
                  !fdr.disc, na.rm = T)
                mysummary.fdr[paste("FDR", annot_level_name, 
                  sep = "_"), "[=]*"] <- ifelse(i == 1, 0, mysummary.fdr[paste("FDR", 
                  annot_level_name, sep = "_"), "[=]*"]) + tmpcount
                marker.filter_k <- marker.filter_k | (fdrlogtemp & 
                  mydata >= gws.level & !fdr.disc)
            }
            marker.filter_k <- FALSE
            for (i in 1:dim(myfilter)[2]) {
                logtempfilter_k <- logtemp_k & myfilter[, i] & 
                  !marker.filter_k
                tmpcount <- sum(gwslogtemp & mydata >= gws.level, 
                  na.rm = T)
                mysummary.fdr[paste("GWS", annot_level_name, 
                  sep = "_"), "[=]*"] <- ifelse(i == 1, 0, mysummary.fdr[paste("GWS", 
                  annot_level_name, sep = "_"), "[=]*"]) + tmpcount
                marker.filter_k <- marker.filter_k | (gwslogtemp & 
                  mydata >= gws.level)
            }
            if (k > 2) {
                temp.fdr <- array(NA, dim = dim(myfilter))
                logtemp_k <- as.numeric(mydata[, myannotq]) < 
                  k
                logtemp_k[is.na(logtemp_k)] <- FALSE
                for (i in 1:dim(myfilter)[2]) {
                  mytempdata <- mydata[logtemp_k & myfilter[, 
                    i]]
                  temp.fdr[logtemp_k, i] <- p.adjust(mytempdata, 
                    method = "fdr")
                }
                mytable[, paste("FDR_lt", annot_level_name, sep = "_")] <- rowMeans(temp.fdr, 
                  na.rm = T)
                if (parameters$do_draw_fdr_vs_p_plots) {
                  mykey <- append(mykey, 25)
                  mykey.color <- append(mykey.color, mycolors[k])
                  mykey.text <- append(mykey.text, paste("[<]", 
                    annot_level_name))
                  points(-log10(mydata[logtemp_k & plot.filter]), 
                    mytable[logtemp_k & plot.filter, paste("FDR_lt", 
                      annot_level_name, sep = "_")], col = mycolors[k], 
                    pch = 25, type = "b")
                }
                marker.filter_k <- FALSE
                for (i in 1:dim(myfilter)[2]) {
                  logtempfilter_k <- logtemp_k & myfilter[, i] & 
                    !marker.filter_k
                  fdrlogtemp <- logtempfilter_k & temp.fdr[, 
                    i] < myparameters$myFDR
                  mysummary.fdr[paste("FDR", annot_level_name, 
                    sep = "_"), "[<]"] <- ifelse(i == 1, 0, mysummary.fdr[paste("FDR", 
                    annot_level_name, sep = "_"), "[<]"])
                  +sum(fdrlogtemp, na.rm = T)
                  marker.filter_k <- marker.filter_k | fdrlogtemp
                }
                marker.filter_k <- FALSE
                for (i in 1:dim(myfilter)[2]) {
                  logtempfilter_k <- logtemp_k & myfilter[, i] & 
                    !marker.filter_k
                  gwslogtemp <- logtempfilter_k & mydata < myparameters$myFDR/sum(logtemp_k)
                  mysummary.fdr[paste("GWS", annot_level_name, 
                    sep = "_"), "[<]"] <- ifelse(i == 1, 0, mysummary.fdr[paste("GWS", 
                    annot_level_name, sep = "_"), "[<]"])
                  +sum(gwslogtemp, na.rm = T)
                  marker.filter_k <- marker.filter_k | gwslogtemp
                }
                marker.filter_k <- FALSE
                for (i in 1:dim(myfilter)[2]) {
                  logtempfilter_k <- logtemp_k & myfilter[, i] & 
                    !marker.filter_k
                  fdrlogtemp <- logtempfilter_k & temp.fdr[, 
                    i] < myparameters$myFDR
                  tmpcount <- sum(fdrlogtemp & mydata >= gws.level & 
                    !fdr.disc, na.rm = T)
                  mysummary.fdr[paste("FDR", annot_level_name, 
                    sep = "_"), "[<]*"] <- ifelse(i == 1, 0, 
                    mysummary.fdr[paste("FDR", annot_level_name, 
                      sep = "_"), "[<]*"])
                  +tmpcount
                  marker.filter_k <- marker.filter_k | (fdrlogtemp & 
                    mydata >= gws.level & !fdr.disc)
                }
                marker.filter_k <- FALSE
                for (i in 1:dim(myfilter)[2]) {
                  logtempfilter_k <- logtemp_k & myfilter[, i] & 
                    !marker.filter_k
                  gwslogtemp <- logtempfilter_k & mydata < myparameters$myFDR/sum(logtemp_k)
                  tmpcount <- sum(gwslogtemp & mydata >= gws.level, 
                    na.rm = T)
                  mysummary.fdr[paste("GWS", annot_level_name, 
                    sep = "_"), "[<]*"] <- ifelse(i == 1, 0, 
                    mysummary.fdr[paste("GWS", annot_level_name, 
                      sep = "_"), "[<]*"])
                  +tmpcount
                  marker.filter_k <- marker.filter_k | (gwslogtemp & 
                    mydata >= gws.level)
                }
            }
            if (k < length(myannotq_levels) - 1) {
                temp.fdr <- array(NA, dim = dim(myfilter))
                logtemp_k <- as.numeric(mydata[, myannotq]) > 
                  k
                logtemp_k[is.na(logtemp_k)] <- FALSE
                for (i in 1:dim(myfilter)[2]) {
                  mytempdata <- mydata[logtemp_k & myfilter[, 
                    i]]
                  temp.fdr[logtemp_k, i] <- p.adjust(mytempdata[, 
                    pheno], method = "fdr")
                }
                mytable[, paste("FDR_gt", annot_level_name, sep = "_")] <- rowMeans(temp.fdr, 
                  na.rm = T)
                if (parameters$do_draw_fdr_vs_p_plots) {
                  mykey <- append(mykey, 25)
                  mykey.color <- append(mykey.color, mycolors[k])
                  mykey.text <- append(mykey.text, paste("[>]", 
                    annot_level_name))
                  points(-log10(mydata[logtemp_k & plot.filter]), 
                    mytable[logtemp_k & plot.filter, paste("FDR_gt", 
                      annot_level_name, sep = "_")], col = mycolors[k], 
                    pch = 25, type = "b")
                }
                marker.filter_k <- FALSE
                for (i in 1:dim(myfilter)[2]) {
                  logtempfilter_k <- logtemp_k & myfilter[, i] & 
                    !marker.filter_k
                  fdrlogtemp <- logtempfilter_k & temp.fdr[, 
                    i] < myparameters$myFDR
                  mysummary.fdr[paste("FDR", annot_level_name, 
                    sep = "_"), "[>]"] <- ifelse(k == 1, NA, 
                    ifelse(i == 1, 0, mysummary.fdr[paste("FDR", 
                      annot_level_name, sep = "_"), "[>]"]) + 
                      sum(fdrlogtemp, na.rm = T))
                  marker.filter_k <- marker.filter_k | fdrlogtemp
                }
                marker.filter_k <- FALSE
                for (i in 1:dim(myfilter)[2]) {
                  logtempfilter_k <- logtemp_k & myfilter[, i] & 
                    !marker.filter_k
                  gwslogtemp <- logtempfilter_k & mydata[, pheno] < 
                    myparameters$myFDR/sum(logtemp_k)
                  mysummary.fdr[paste("GWS", annot_level_name, 
                    sep = "_"), "[>]"] <- ifelse(k == 1, NA, 
                    ifelse(i == 1, 0, mysummary.fdr[paste("GWS", 
                      annot_level_name, sep = "_"), "[>]"]) + 
                      sum(gwslogtemp, na.rm = T))
                  marker.filter_k <- marker.filter_k | gwslogtemp
                }
                marker.filter_k <- FALSE
                for (i in 1:dim(myfilter)[2]) {
                  logtempfilter_k <- logtemp_k & myfilter[, i] & 
                    !marker.filter_k
                  fdrlogtemp <- logtempfilter_k & temp.fdr[, 
                    i] < myparameters$myFDR
                  tmpcount <- sum(fdrlogtemp & mydata >= gws.level & 
                    !fdr.disc, na.rm = T)
                  mysummary.fdr[paste("FDR", annot_level_name, 
                    sep = "_"), "[>]*"] <- ifelse(k == 1, NA, 
                    ifelse(i == 1, 0, mysummary.fdr[paste("FDR", 
                      annot_level_name, sep = "_"), "[>]*"]) + 
                      tmpcount)
                  marker.filter_k <- marker.filter_k | (fdrlogtemp & 
                    mydata >= gws.level & !fdr.disc)
                }
                marker.filter_k <- FALSE
                for (i in 1:dim(myfilter)[2]) {
                  logtempfilter_k <- logtemp_k & myfilter[, i] & 
                    !marker.filter_k
                  gwslogtemp <- logtempfilter_k & mydata < myparameters$myFDR/sum(logtemp_k)
                  tmpcount <- sum(gwslogtemp & mydata >= gws.level, 
                    na.rm = T)
                  mysummary.fdr[paste("GWS", annot_level_name, 
                    sep = "_"), "[>]*"] <- ifelse(k == 1, NA, 
                    ifelse(i == 1, 0, mysummary.fdr[paste("GWS", 
                      annot_level_name, sep = "_"), "[>]*"]) + 
                      tmpcount)
                  marker.filter_k <- marker.filter_k | (gwslogtemp & 
                    mydata >= gws.level)
                }
            }
        }
        if (parameters$do_draw_fdr_vs_p_plots) {
            legend("topright", mykey.text, col = mykey.color, 
                pch = mykey, bg = "white")
        }
        print(mysummary.fdr)
        if (fdr.file != "") 
            write.table(mytable, file = fdr.file, quote = F, 
                row.names = F)
        return(mytable)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
