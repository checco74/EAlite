\name{mycolorfunction}
\alias{mycolorfunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to compute colors in QQ-plots
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mycolorfunction(k, n, base.colors = c("#DD1E2F", "#EBB035", "#06A2CB", "#006500"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{base.colors}{
%%     ~~Describe \code{base.colors} here~~
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
function (k, n, base.colors = c("#DD1E2F", "#EBB035", "#06A2CB", 
    "#006500")) 
{
    if (is.null(base.colors)) {
        base.colors <- c("#DD1E2F", "#EBB035", "#06A2CB", "#006500")
    }
    if (length(base.colors) <= 1) {
        base.colors <- c(base.colors, "#DD1E2F", "#EBB035", "#06A2CB", 
            "#006500")
    }
    home.colors <- rgb(t(col2rgb(colors()))/255)
    index = length(home.colors)
    take.doubles <- FALSE
    while (length(base.colors) < n) {
        new.color <- home.colors[index]
        if (!new.color \%in\% base.colors || take.doubles) 
            base.colors <- c(base.colors, new.color)
        index <- index - 1
        if (index == 1) {
            take.doubles = TRUE
            index = length(home.colors)
        }
    }
    return(base.colors[k])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
