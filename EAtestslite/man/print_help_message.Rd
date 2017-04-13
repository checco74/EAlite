\name{print_help_message}
\alias{print_help_message}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to print help for use within an R-script
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print_help_message(parameters = init_parameters())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
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
function (parameters = init_parameters()) 
{
    cat("\n USAGE:\n")
    cat("\n  ealite.Rscript [COMMANDS] [OPTIONS] <input file(s)>\n")
    cat("\n  the script expects all data provided to be coherently sorted;\n")
    cat("  annotation data tracks are specified with '--annot' (see below).\n")
    cat("  the script writes a number of files, all prefixed with <prefix> (see '--o' below).\n")
    cat("\n COMMANDS:\n\n")
    cat("    --help            print help and exit;\n")
    cat("    --annot <annot>   annotation (numeric or word) file; specify as many as desired;\n")
    cat("    --bpt             perform binomial proportion tests [ no tests by default ];\n")
    cat("    --bpts            perform bidirectional binomial proportion tests [ no tests by default ];\n")
    cat("    --cannot <annot>  covariate annotation (numeric or word) file; specify as many as desired;\n")
    cat("    --color <value>   quantile/interval colors in figures; missing ones will be picked semi-randomly;\n")
    cat("    --ctrl <ctrl>     control set file [ default is everything ];\n")
    cat("    --diag            request regression diagnostics [ default ];\n")
    cat("    --distr <value>   expected distribution ('e', 'u') of the phenotype [ default is 'u' (uniform) ];\n")
    cat("    --eest            estimate enrichment and significance [ no estimate by default ];\n")
    cat("    --fdr             compute FDR [ no fdr by default ];\n")
    cat("    --fdrv <value>    threshold for false discovery [ default =", 
        parameters$myFDR, "];\n")
    cat("    --fdrvsp          draw FDR vs p-value plots (forces '--fdr') [ no plots by default ];\n")
    cat("    --filter <file>   logical array to restrict the analyses to: specify as many as desired;\n")
    cat("    --gc              perform genomic correction of p-values [ no gc by default ];\n")
    cat("    --header          if text data files have a header line [ no header expected by default ];\n")
    cat("    --hgt             perform Fisher (hypergeometric) tests [ no tests by default ];\n")
    cat("    --hgts            perform bidirectional Fisher (hypergeometric) tests [ no tests by default ];\n")
    cat("    --i               include annot*cannot interactions in regression [ no interactions by default ]:\n")
    cat("                      non-interaction effects are estimated before including interactions in the model;\n")
    cat("                      interactions are then included singularly for their effects to be estimated;\n")
    cat("    --nbins <value>   number of bins to divide scatter plots (see below) in [ default =", 
        parameters$N_bins, "];\n")
    cat("    --o <prefix>      prefix for output files [ default = '", 
        parameters$outprefix, "' ];\n", sep = "")
    cat("    --prct <value>    p-value percentiles in which to test for significance [ default =", 
        parameters$myP_prctle, "];\n")
    cat("                      only one threshold can be specified per command line option;\n")
    cat("    --qq              draw qq/fold-plots [ no plots by default ];\n")
    cat("    --qqci            draw confidence intervals for qq/fold-plots [ not by default ];\n")
    cat("    --r               expect input files in R data format;\n")
    cat("    --scatter <annot> draw binned scatter plots of squared z-scores VS <annot> for different quantiles of the\n")
    cat("                      annotation of interest (see below); forces --eest option;\n")
    cat("    --strbr <value>   breaking criteria for annotations:\n")
    cat("                      values can be either words 'auto', 'std', 'unq' or comma-separated lists of\n")
    cat("                      break points;\n")
    cat("                        'auto'  use the annotation value as category (categorical annotations),\n")
    cat("                        'std'   make (n) quantiles (see '--strn') using _all_ annotation values,\n")
    cat("                        'unq'   make (n) quantiles removing duplicate annotation values;\n")
    cat("    --strn <value>    n. of annotation quantiles [ default =", 
        parameters$numof.breaks, "];\n")
    cat("                      [note: these are overridden by '--strbr' with explicit break points]\n")
    cat("    --testp <value>   significance tests p-value plots ceiling [ default =", 
        parameters$testp, "];\n")
    cat("    --v               verbose mode on (the program prints lots of stuff on screen);\n")
    cat("\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
