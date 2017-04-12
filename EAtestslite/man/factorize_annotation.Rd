\name{factorize_annotation}
\alias{factorize_annotation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
make factors from annotations
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
factorize_annotation(annotation, mybreaks = "auto", numof.breaks = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annotation}{
%%     ~~Describe \code{annotation} here~~
}
  \item{mybreaks}{
%%     ~~Describe \code{mybreaks} here~~
}
  \item{numof.breaks}{
%%     ~~Describe \code{numof.breaks} here~~
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
function (annotation, mybreaks = "auto", numof.breaks = 4) 
{
    cat("factorizing annotation..\n", sep = "")
    if (is.factor(annotation) || is.logical(annotation) || mybreaks == 
        "auto") {
        cat("using auto breaks..\n")
        mylevels = levels(as.factor(annotation))
        if (length(mylevels) < 100) {
            return(factor(annotation))
        }
        else mybreaks = "std"
    }
    if (mybreaks == "std") {
        myqdata <- annotation
        mycurrentbreaks <- quantile(myqdata, probs = seq(0, 1, 
            by = 1/numof.breaks), na.rm = T)
        cat("using standard breaks: ", mycurrentbreaks, "\n")
    }
    else if (mybreaks == "unq") {
        myqdata <- unique(annotation)
        mycurrentbreaks <- quantile(myqdata, probs = seq(0, 1, 
            by = 1/numof.breaks), na.rm = T)
        cat("using unique breaks: ", mycurrentbreaks, "\n")
    }
    else {
        cat("breaking on ", mybreaks, "\n")
        mycurrentbreaks <- sort(c(as.numeric(unlist(strsplit(as.character(mybreaks), 
            ","))), min(annotation, na.rm = T), max(annotation, 
            na.rm = T)))
        cat("using custom breaks: ", mycurrentbreaks, "\n")
    }
    if (length(mycurrentbreaks) > length(unique(mycurrentbreaks))) {
        cat("warning: collapsing non-unique breaking points..\n", 
            sep = "")
        mycurrentbreaks = sort(unique(mycurrentbreaks))
    }
    cat("generating factors.. ")
    annotation_quant <- factor(as.character(cut(annotation, mycurrentbreaks, 
        include.lowest = T)))
    cat("done.\n")
    return(annotation_quant)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
