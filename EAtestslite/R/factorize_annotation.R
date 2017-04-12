factorize_annotation <-
function (annotation, mybreaks = "auto", numof.breaks = 4)
{
    cat("factorizing annotation..\n", sep = "")
    if (is.factor( annotation ) || is.logical( annotation ) || mybreaks == "auto") {
        cat("using auto breaks..\n")
        mylevels = levels( as.factor( annotation ) )
        if ( length( mylevels ) < 100 ) {
            return( factor(annotation) )
        } else mybreaks = "std"
    }
    if ( mybreaks == "std" ) {
        myqdata <- annotation
        mycurrentbreaks <- quantile(myqdata, probs = seq(0,
            1, by = 1/numof.breaks), na.rm = T)
        cat("using standard breaks: ", mycurrentbreaks,
            "\n")
    }
    else if ( mybreaks == "unq" ) {
        myqdata <- unique(annotation)
        mycurrentbreaks <- quantile(myqdata, probs = seq(0,
            1, by = 1/numof.breaks), na.rm = T)
        cat("using unique breaks: ", mycurrentbreaks, "\n")
    }
    else {
        cat("breaking on ", mybreaks, "\n")
        mycurrentbreaks <- sort(c(as.numeric(unlist(strsplit(as.character(mybreaks),
            ","))), min(annotation, na.rm = T), max(annotation, na.rm = T)))
        cat("using custom breaks: ", mycurrentbreaks, "\n")
    }
    if (length(mycurrentbreaks) > length(unique(mycurrentbreaks))) {
        cat("warning: collapsing non-unique breaking points..\n", sep = "")
        mycurrentbreaks = sort(unique(mycurrentbreaks))
    }
    cat("generating factors.. ")
    annotation_quant <- factor( as.character(
        cut(annotation, mycurrentbreaks, include.lowest = T)
    ) )
    cat("done.\n")
    return(annotation_quant)
}

