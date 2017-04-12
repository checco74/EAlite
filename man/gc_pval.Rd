\name{gc_pval}
\alias{gc_pval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to do genomic correction
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gc_pval(pvalstat, filter = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvalstat}{
%%     ~~Describe \code{pvalstat} here~~
}
  \item{filter}{
%%     ~~Describe \code{filter} here~~
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
function (pvalstat, filter = NULL) 
{
    chisqstat <- qchisq(pvalstat, 1, lower.tail = F)
    if (is.null(filter)) {
        gc.lambda <- median(chisqstat, na.rm = T)/qchisq(0.5, 
            1)
    }
    else {
        gc.lambda = c()
        Nf = dim(filter)[2]
        if (is.null(Nf) || is.na(Nf)) 
            Nf = 1
        for (i in 1:Nf) {
            gc.lambda <- c(gc.lambda, median(chisqstat[filter[, 
                i] & is.finite(chisqstat)], na.rm = T)/qchisq(0.5, 
                1))
        }
        gc.lambda = median(gc.lambda, na.rm = T)
    }
    print(gc.lambda)
    return(pchisq(chisqstat/gc.lambda, 1, lower.tail = F))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
