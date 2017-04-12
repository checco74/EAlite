\name{perform_hgt}
\alias{perform_hgt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to perform hypergeometric tests
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
perform_hgt(mydata.full = NULL, myannot.full = NULL, myctrlannot.full = NULL, myparameters = init_parameters(), myfilter.full = NULL, P_prctle = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mydata.full}{
%%     ~~Describe \code{mydata.full} here~~
}
  \item{myannot.full}{
%%     ~~Describe \code{myannot.full} here~~
}
  \item{myctrlannot.full}{
%%     ~~Describe \code{myctrlannot.full} here~~
}
  \item{myparameters}{
%%     ~~Describe \code{myparameters} here~~
}
  \item{myfilter.full}{
%%     ~~Describe \code{myfilter.full} here~~
}
  \item{P_prctle}{
%%     ~~Describe \code{P_prctle} here~~
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
function (mydata.full = NULL, myannot.full = NULL, myctrlannot.full = NULL, 
    myparameters = init_parameters(), myfilter.full = NULL, P_prctle = 0.01) 
{
    do_cumulative = FALSE
    if (!is.null(mydata.full) && !is.null(myannot.full)) {
        test.hypothesis = "greater"
        if (myparameters$bidirect) {
            cat("computing two-sided p-values..\n")
            test.hypothesis = "two.sided"
        }
        contingency_table <- function(mytemp, myctrl, mythresh) {
            HGT_pos_signif <- sum((mytemp < mythresh), na.rm = T)
            HGT_pos_nonsignif <- sum((mytemp >= mythresh), na.rm = T)
            HGT_ctrl_signif <- sum((myctrl < mythresh), na.rm = T)
            HGT_ctrl_nonsignif <- sum((myctrl >= mythresh), na.rm = T)
            return(matrix(c(HGT_pos_signif, HGT_pos_nonsignif, 
                HGT_ctrl_signif, HGT_ctrl_nonsignif), nrow = 2))
        }
        mymask <- c(!is.na(mydata.full) & !is.na(myannot.full))
        if (is.null(myctrlannot.full)) 
            myctrlannot.full = array(TRUE, dim = c(length(mydata.full), 
                1))
        if (is.null(myfilter.full)) 
            myfilter.full = array(TRUE, dim = c(length(mydata.full), 
                1))
        myannotq.full <- factorize_annotation(myannot.full, myparameters$mybreaks, 
            myparameters$numof.breaks)
        if (myparameters$verbose) 
            print(head(myannotq.full))
        myannotq_levels <- levels(as.factor(myannotq.full))
        P_thresh = quantile(mydata.full, P_prctle, na.rm = T)
        if (do_cumulative) {
            HGT_test <- array(0, dim = c(length(myannotq_levels), 
                3, 2, 2))
            mypvalues.hgt <- array(NA, dim = c(length(myannotq_levels), 
                3))
            mysummary.hgt <- array(NA, dim = c(length(myannotq_levels), 
                3))
            mymin.hgt <- array(NA, dim = c(length(myannotq_levels), 
                3))
            mymax.hgt <- array(NA, dim = c(length(myannotq_levels), 
                3))
            colnames(mypvalues.hgt) <- c("[=]", "[<]", "[>]")
            colnames(mysummary.hgt) <- c("[=]", "[<]", "[>]")
            colnames(mymin.hgt) <- c("[=]", "[<]", "[>]")
            colnames(mymax.hgt) <- c("[=]", "[<]", "[>]")
        }
        else {
            HGT_test <- array(0, dim = c(length(myannotq_levels), 
                2, 2))
            mypvalues.hgt <- array(NA, dim = c(length(myannotq_levels), 
                1))
            mysummary.hgt <- array(NA, dim = c(length(myannotq_levels), 
                1))
            mymin.hgt <- array(NA, dim = c(length(myannotq_levels), 
                1))
            mymax.hgt <- array(NA, dim = c(length(myannotq_levels), 
                1))
            colnames(mypvalues.hgt) <- c("[=]")
            colnames(mysummary.hgt) <- c("[=]")
            colnames(mymin.hgt) <- c("[=]")
            colnames(mymax.hgt) <- c("[=]")
        }
        rownames(mypvalues.hgt) <- myannotq_levels
        rownames(mysummary.hgt) <- myannotq_levels
        rownames(mymin.hgt) <- myannotq_levels
        rownames(mymax.hgt) <- myannotq_levels
        myfilter = NULL
        if (is.null(myfilter.full)) {
            myfilter = array(TRUE, dim = c(length(mydata.full), 
                1))
        }
        else {
            myfilter = myfilter.full & mymask
        }
        cat("running tests..\n")
        for (i in 1:dim(myfilter)[2]) {
            mydata <- mydata.full[myfilter[, i]]
            myannot <- myannot.full[myfilter[, i]]
            myannotq <- myannotq.full[myfilter[, i]]
            myctrlannot <- myctrlannot.full[myfilter[, i]]
            for (k in 1:length(myannotq_levels)) {
                logtemp_k <- myannotq == myannotq_levels[k]
                logctrl_k <- myctrlannot & !logtemp_k
                if (do_cumulative) {
                  HGT_test[k, 1, , ] <- HGT_test[k, 1, , ] + 
                    contingency_table(mydata[logtemp_k], mydata[logctrl_k], 
                      P_thresh)
                  if (k > 2) {
                    logtemp_k <- as.numeric(myannotq) < k
                    logctrl_k <- myctrlannot & !logtemp_k
                    HGT_test[k, 2, , ] <- HGT_test[k, 2, , ] + 
                      contingency_table(mydata[logtemp_k], mydata[logctrl_k], 
                        P_thresh)
                  }
                  if (k < length(myannotq_levels) - 1) {
                    logtemp_k <- as.numeric(myannotq) > k
                    logctrl_k <- myctrlannot & !logtemp_k
                    HGT_test[k, 3, , ] <- HGT_test[k, 3, , ] + 
                      contingency_table(mydata[logtemp_k], mydata[logctrl_k], 
                        P_thresh)
                  }
                }
                else {
                  HGT_test[k, , ] <- HGT_test[k, , ] + contingency_table(mydata[logtemp_k], 
                    mydata[logctrl_k], P_thresh)
                }
            }
        }
        if (myparameters$verbose) {
            cat("contingency table:\n")
            print(HGT_test)
        }
        for (k in 1:length(myannotq_levels)) {
            myannotq_name <- myannotq_levels[k]
            if (do_cumulative) {
                for (ss in 1:length(colnames(mypvalues.hgt))) {
                  if (all(HGT_test[k, ss, , ] > 0)) {
                    mytest = fisher.test(HGT_test[k, ss, , ], 
                      alternative = test.hypothesis)
                    if (myparameters$verbose) 
                      print(mytest)
                    mysummary.hgt[myannotq_name, ss] <- mytest$estimate
                    mypvalues.hgt[myannotq_name, ss] <- mytest$p.value
                  }
                }
            }
            else {
                if (all(HGT_test[k, , ] > 0)) {
                  mytest = fisher.test(HGT_test[k, , ], alternative = test.hypothesis)
                  if (myparameters$verbose) 
                    print(mytest)
                  mypvalues.hgt[myannotq_name, "[=]"] <- mytest$p.value
                  mysummary.hgt[myannotq_name, "[=]"] <- mytest$estimate
                }
            }
        }
        return(list(pvalues = mypvalues.hgt, summary = mysummary.hgt, 
            min = mymin.hgt, max = mymax.hgt))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
