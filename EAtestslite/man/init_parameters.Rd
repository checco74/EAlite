\name{init_parameters}
\alias{init_parameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to initialize the parameters for the enrichment analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
init_parameters(..., .listed = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{.listed}{
%%     ~~Describe \code{.listed} here~~
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
function (..., .listed = FALSE) 
{
    if (.listed) {
        mylist = (...)
    }
    else mylist = list(...)
    returnlist = list(helpme = FALSE, annotation.files = NULL, 
        bsannotation.file = NULL, covannotation.files = NULL, 
        ctrlannotation.file = "", pheno.files = NULL, bidirect = FALSE, 
        distr = "u", do_compute_fdr = FALSE, do_draw_bs_plots = FALSE, 
        do_draw_fdr_plots = FALSE, do_draw_qq_plots = FALSE, 
        do_draw_cis = FALSE, do_estimate_enrichment = FALSE, 
        do_genomic_correction = FALSE, do_perform_bpt = FALSE, 
        do_perform_hgt = FALSE, expect_rdata = FALSE, interact = FALSE, 
        figure.colors = NULL, N_bins = 100, outprefix = "ea", 
        read.header = FALSE, regdiag = FALSE, verbose = FALSE, 
        filter.names = c(), myP_prctle = c(0.001, 0.01, 0.1, 
            0.5), tiny = 1e-72, testp = 1e-12, neglog10max = 12, 
        N_ytics = 5, myFDR = 0.01, mybreaks = "auto", numof.breaks = 4)
    if (length(mylist) > 0) {
        argcnt = 0
        for (k in 1:length(mylist)) {
            if (is.list(mylist[[k]])) {
                returnlist = init_parameters(mylist[[k]], .listed = TRUE)
            }
            else {
                if (length(names(mylist)) >= k && names(mylist)[k] != 
                  "") {
                  returnlist[[names(mylist)[k]]] = mylist[[k]]
                }
                else if (length(names(returnlist)) - argcnt > 
                  0) 
                  returnlist[[argcnt + 1]] = mylist[[k]]
            }
            argcnt = argcnt + 1
        }
    }
    return(returnlist)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
