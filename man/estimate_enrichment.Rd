\name{estimate_enrichment}
\alias{estimate_enrichment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to estimate enrichment
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimate_enrichment(mydata.full = NULL, myctrlannot.full = NULL, myregannot.full = NULL, myparameters = init_parameters(), myfilter.full = NULL, mybsannot.full = NULL, diagfile = NULL, dumpfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mydata.full}{
%%     ~~Describe \code{mydata.full} here~~
}
  \item{myctrlannot.full}{
%%     ~~Describe \code{myctrlannot.full} here~~
}
  \item{myregannot.full}{
%%     ~~Describe \code{myregannot.full} here~~
}
  \item{myparameters}{
%%     ~~Describe \code{myparameters} here~~
}
  \item{myfilter.full}{
%%     ~~Describe \code{myfilter.full} here~~
}
  \item{mybsannot.full}{
%%     ~~Describe \code{mybsannot.full} here~~
}
  \item{diagfile}{
%%     ~~Describe \code{diagfile} here~~
}
  \item{dumpfile}{
%%     ~~Describe \code{dumpfile} here~~
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
function (mydata.full = NULL, myctrlannot.full = NULL, myregannot.full = NULL, 
    myparameters = init_parameters(), myfilter.full = NULL, mybsannot.full = NULL, 
    diagfile = NULL, dumpfile = NULL) 
{
    dfcorr = 0.8
    mincount = 5
    do_cumulative = FALSE
    if (!is.null(mydata.full) && !is.null(myregannot.full)) {
        if (is.null(dim(myregannot.full))) 
            myregannot.full = data.frame(defannot = array(myregannot.full, 
                dim = c(length(myregannot.full), 1)))
        if (!is.data.frame(myregannot.full)) 
            myregannot.full = data.frame(myregannot.full)
        cat("removing uninformative covariates..\n")
        myregannot.full = myregannot.full[, apply(!is.na(myregannot.full), 
            2, any), drop = F]
        cat("covariates array dimensions: ")
        print(dim(myregannot.full))
        if (is.null(myctrlannot.full)) 
            myctrlannot.full = array(TRUE, dim = c(length(mydata.full), 
                1))
        mymask <- c(is.finite(mydata.full) & apply(!is.na(myregannot.full), 
            1, all))
        if (sum(mymask) > dim(myregannot.full)[2]) {
            cat(sum(mymask), "samples.\n")
            myregannot_names = names(myregannot.full)
            if (is.null(myregannot_names)) 
                myregannot_names = paste0("V_", 1:dim(myregannot.full)[2])
            for (k in which(myregannot_names == "")) {
                while (myregannot_names[k] == "" || myregannot_names[k] \%in\% 
                  myregannot_names[-k]) myregannot_names[k] = paste0(myregannot_names[k], 
                  "__defannot__", sprintf("\%d", k))
            }
            if (myparameters$do_draw_bs_plots) {
                if (myparameters$verbose) 
                  cat("N. of bins:", myparameters$N_bins, "\n")
                mybsannot_cuts = quantile(mybsannot.full, probs = seq(0, 
                  1, by = 1/myparameters$N_bins), na.rm = T)
                mybsannot.full_bins = cut(mybsannot.full, mybsannot_cuts)
                mybsannot.full_aggregate = as.character(aggregate(mybsannot.full, 
                  list(mybsannot.full_bins), max, na.rm = T)[, 
                  1])
            }
            myannotq_levels = c()
            myannotq.full = factorize_annotation(myregannot.full[, 
                1], myparameters$mybreaks, myparameters$numof.breaks)
            myannotq_levels = levels(myannotq.full)
            Nq = length(myannotq_levels)
            if (length(unique(myregannot.full[, 1])) == Nq) {
                cat("projecting annotation to its factorized version..\n")
                myregannot.full[, 1] = myannotq.full
            }
            if (do_cumulative && Nq > 1) 
                Nq = Nq + 2 * (Nq - 2)
            myfilter = NULL
            if (is.null(myfilter.full)) {
                cat("no filter specified: using all.\n")
                myfilter = array(TRUE, dim = c(length(mydata.full), 
                  1)) & mymask
            }
            else {
                if (myparameters$verbose) {
                  cat("using filters:\n")
                  print(str(myfilter.full))
                  print(str(mymask))
                }
                cat("building filters.. ")
                myfilter = myfilter.full & mymask
                cat("done.\n")
            }
            Nf = dim(myfilter)[2]
            if (myparameters$regdiag) {
                if (myparameters$do_draw_bs_plots) {
                  png(diagfile, height = 250 * Nf, width = 2000)
                  layout(matrix(1:(8 * Nf * (1 + Nq)), nrow = Nf, 
                    ncol = 8 * (1 + Nq), byrow = T))
                }
                else {
                  png(diagfile, height = 250 * Nf, width = 1000)
                  layout(matrix(1:(4 * Nf * (1 + Nq)), nrow = Nf, 
                    ncol = 4 * (1 + Nq), byrow = T))
                  par(mar = c(1, 1, 1, 1))
                }
            }
            cat("creating container tables..\n")
            regression.entries = c("", "(se)", "(t)", "(p)", 
                "(ci_a)", "(ci_b)")
            Nre = length(regression.entries)
            regression.coeff = paste0(rep("intercept", Nre), 
                regression.entries)
            regression.var0 = c()
            regression.var = c()
            cat("processing factors..\n")
            for (k in 1:dim(myregannot.full)[2]) {
                if (myparameters$verbose) {
                  cat("processing factor", k, "..\n")
                  print(summary(myregannot.full[, k]))
                  print(head(myregannot.full[, k]))
                }
                if (is.numeric(myregannot.full[, k]) && length(unique(myregannot.full[, 
                  k])) > 20) {
                  myregannot.full[, k] = (myregannot.full[, k] - 
                    mean(myregannot.full[, k], na.rm = T))/sd(myregannot.full[, 
                    k], na.rm = T)
                }
                else myregannot.full[, k] = factor(myregannot.full[, 
                  k])
                if (is.factor(myregannot.full[, k])) {
                  klevels = levels(myregannot.full[, k])[-1]
                  Nkl = max(length(klevels), 1)
                  regression.coeff = c(regression.coeff, paste0(myregannot_names[k], 
                    "_", rep(klevels, each = Nre), rep(regression.entries, 
                      Nkl)))
                  regression.var = c(regression.var, paste0(myregannot_names[k], 
                    "_", klevels))
                }
                else {
                  regression.coeff = c(regression.coeff, paste0(myregannot_names[k], 
                    regression.entries))
                  regression.var = c(regression.var, myregannot_names[k])
                }
                if (k == 1) {
                  regression.var0 = regression.var
                }
            }
            if (myparameters$interact) {
                for (k in 2:dim(myregannot.full)[2]) {
                  if (is.factor(myregannot.full[, k])) {
                    klevels = levels(myregannot.full[, k])[-1]
                    kbd = rep(paste0(myregannot_names[k], "_", 
                      klevels), each = length(regression.var0))
                    Nkl = max(length(klevels), 1)
                    regression.coeff = c(regression.coeff, paste0(rep(paste(rep(regression.var0, 
                      Nkl), kbd, sep = "*"), each = Nre), rep(regression.entries, 
                      Nkl * length(regression.var0))))
                    regression.var = c(regression.var, paste(rep(regression.var0, 
                      Nkl), rep(paste0(myregannot_names[k], "_", 
                      klevels), each = length(regression.var0)), 
                      sep = "*"))
                  }
                  else {
                    regression.coeff = c(regression.coeff, paste0(rep(paste(regression.var0, 
                      rep(myregannot_names[k], each = length(regression.var0)), 
                      sep = "*"), each = Nre), rep(regression.entries, 
                      length(regression.var0))))
                    regression.var = c(regression.var, paste(regression.var0, 
                      rep(myregannot_names[k], each = length(regression.var0)), 
                      sep = "*"))
                  }
                }
            }
            pLR.names = paste0(myregannot_names[1], "(pLR)")
            if (myparameters$interact && dim(myregannot.full)[2] > 
                1) 
                pLR.names = c(pLR.names, paste0(myregannot_names[1], 
                  "*", myregannot_names[-1], "(pLR)"))
            regression.coeff = c(regression.coeff, pLR.names)
            if (myparameters$verbose) {
                cat("regression coefficients:\n")
                print(regression.coeff)
            }
            N.testq <- length(regression.coeff) + 3
            inames = NULL
            rnames = NULL
            snames = NULL
            if (do_cumulative) {
                snames = c("[=]", "[<]", "[>]")
                rnames = paste0(rep(snames, each = N.testq), 
                  c("", "(F)", "(pF)", regression.coeff))
            }
            else {
                snames = c("[=]")
                rnames = paste0(rep(snames, N.testq), c("", "(F)", 
                  "(pF)", regression.coeff))
            }
            Nm = length(snames)
            inames = paste("I", 1:Nq, sep = "_")
            myenrichment.mean <- array(NA, dim = c(Nm * N.testq, 
                length(inames) + 1))
            myenrichment.median <- array(NA, dim = c(Nm * N.testq, 
                length(inames) + 1))
            myenrichment.max <- array(NA, dim = c(Nm * N.testq, 
                length(inames) + 1))
            myenrichment.min <- array(NA, dim = c(Nm * N.testq, 
                length(inames) + 1))
            myenrichment.sd <- array(NA, dim = c(Nm * N.testq, 
                length(inames) + 1))
            myenrichment.df <- array(NA, dim = c(Nm * N.testq, 
                length(inames) + 1))
            myenrichment.which <- array(NA, dim = c(Nm * N.testq, 
                length(inames) + 1))
            rownames(myenrichment.mean) <- rnames
            colnames(myenrichment.mean) <- c("all", inames)
            rownames(myenrichment.median) <- rnames
            colnames(myenrichment.median) <- c("all", inames)
            rownames(myenrichment.max) <- rnames
            colnames(myenrichment.max) <- c("all", inames)
            rownames(myenrichment.min) <- rnames
            colnames(myenrichment.min) <- c("all", inames)
            rownames(myenrichment.sd) <- rnames
            colnames(myenrichment.sd) <- c("all", inames)
            rownames(myenrichment.df) <- rnames
            colnames(myenrichment.df) <- c("all", inames)
            rownames(myenrichment.which) <- rnames
            colnames(myenrichment.which) <- c("all", inames)
            cat("creating temporary tables..\n")
            myenrichment <- array(NA, dim = c(Nm * N.testq, length(inames) + 
                1, Nf), dimnames = list(entries = rownames(myenrichment.mean), 
                levels = colnames(myenrichment.mean), iterations = NULL))
            myenrichment.numdf <- array(NA, dim = c(Nm * N.testq, 
                length(inames) + 1, Nf), dimnames = list(entries = rownames(myenrichment.mean), 
                levels = colnames(myenrichment.mean), iterations = NULL))
            myenrichment.denomdf <- array(NA, dim = c(Nm * N.testq, 
                length(inames) + 1, Nf), dimnames = list(entries = rownames(myenrichment.mean), 
                levels = colnames(myenrichment.mean), iterations = NULL))
            resnames = c("all")
            if (!is.null(inames)) {
                if (do_cumulative && length(inames) > 2) {
                  resnames = c("all", paste0("[=]", inames), 
                    paste0("[<]", inames[-(1:2)]), paste0("[>]", 
                      inames[-((length(inames) - 2):length(inames))]))
                }
                else resnames = c("all", paste0("[=]", inames))
            }
            mycounts = array(0, dim = c(myparameters$N_bins, 
                1 + Nq), dimnames = list(NULL, resnames))
            myresiduals = array(NA, dim = c(myparameters$N_bins, 
                1 + Nq, Nf), dimnames = list(NULL, resnames, 
                NULL))
            for (i in 1:Nf) {
                mydata <- mydata.full[myfilter[, i]]
                tempz <- qnorm(0.5 * mydata) * (rbinom(length(mydata), 
                  1, 0.5) - 0.5) * 2
                sigma0 = sd(tempz[is.finite(tempz)], na.rm = T)
                if (is.finite(sigma0)) {
                  tempz = tempz/sigma0
                  tempresz = tempz^2
                  N = sum(is.finite(tempresz))
                  cat("enrichment in the complete set:\n")
                  if (myparameters$verbose) 
                    print(myenrichment)
                  myenrichment["[=]", "all", i] <- mean(tempresz, 
                    na.rm = T) - 1
                  myenrichment["[=](F)", "all", i] <- 0
                  myenrichment["[=](pF)", "all", i] <- 1
                  table.entries = paste0(rep("[=]", N.testq - 
                    3), regression.coeff)
                  if (myparameters$verbose) {
                    print("list of table entries:")
                    print(table.entries)
                  }
                  myctrlannot <- myctrlannot.full[myfilter[, 
                    i]]
                  myregannot <- myregannot.full[myfilter[, i], 
                    ]
                  if (is.null(dim(myregannot))) 
                    myregannot = data.frame(defannot = array(myregannot, 
                      dim = c(length(myregannot), 1)))
                  if (!is.data.frame(myregannot)) 
                    myregannot = data.frame(myregannot)
                  if (myparameters$verbose) {
                    cat("response array: ")
                    print(length(tempresz))
                    cat("annotation array: ")
                    print(dim(myregannot))
                  }
                  if (dim(myregannot)[2] > 0 && length(tempresz) > 
                    mincount * dim(myregannot)[2] && N > 1) {
                    cat("full regression model\n")
                    if (myparameters$verbose) {
                      cat("call of regress_enrichment() with arguments:\n")
                      cat("mydata =\n")
                      print(head(tempresz))
                      cat("myregannot =\n")
                      print(head(myregannot))
                      cat("mybreaks = ")
                      print(myparameters$mybreaks)
                      cat("numof.breaks = ")
                      print(myparameters$numof.breaks)
                      cat("interact = ")
                      print(myparameters$interact)
                    }
                    regression.list <- regress_enrichment(tempresz, 
                      myregannot, myparameters)
                    if (myparameters$regdiag) 
                      title(sprintf("set\%d(all)", i))
                    if (myparameters$verbose) {
                      cat("regression list: ")
                      print(str(regression.list))
                      print(rbind(regression.list$reg.names, 
                        table.entries))
                    }
                    if (!is.null(regression.list)) {
                      cat("storing regression stats..\n")
                      myenrichment[table.entries, "all", i] <- regression.list$reg.coeff
                      cat(N, "SNPs analyzed.\n")
                    }
                    regression.redlist = NULL
                    if (myparameters$do_draw_bs_plots) {
                      bsfields = grep(sub("(.RData|.txt)$", "", 
                        basename(myparameters$bsannotation.file)), 
                        names(myregannot))
                      mybsannot <- mybsannot.full[myfilter[, 
                        i]]
                      mybsannot_bins <- mybsannot.full_bins[myfilter[, 
                        i]]
                      myregannotred <- myregannot.full[myfilter[, 
                        i], -c(1, bsfields)]
                      if (all(dim(myregannotred) > 0)) {
                        if (is.null(dim(myregannotred))) 
                          myregannotred = data.frame(array(myregannotred, 
                            dim = c(length(myregannotred), 1)))
                        if (!is.data.frame(myregannotred)) 
                          myregannotred = data.frame(myregannotred)
                        cat("reduced regression model\n")
                        if (myparameters$verbose) {
                          cat("call of regress_enrichment() with arguments:\n")
                          cat("mydata =\n")
                          print(head(tempresz))
                          cat("myregannot =\n")
                          print(head(myregannotred))
                          cat("mybreaks = ")
                          print(myparameters$mybreaks)
                          cat("numof.breaks = ")
                          print(myparameters$numof.breaks)
                          cat("interact = ")
                          print(myparameters$interact)
                        }
                        regression.redlist <- regress_enrichment(tempresz, 
                          myregannotred, myparameters, interact = F)
                      }
                      cat("storing residuals..\n")
                      if (!is.null(regression.redlist)) 
                        tempresz = regression.redlist$residuals
                      mycounts[, "all"] = mycounts[, "all"] + 
                        hist(mybsannot, breaks = mybsannot_cuts, 
                          plot = F)$counts
                      tempaggregate = aggregate(tempresz, list(mybsannot_bins), 
                        sum, na.rm = T)
                      myresiduals[mybsannot.full_aggregate \%in\% 
                        as.character(tempaggregate[, 1]), "all", 
                        i] = tempaggregate[, 2]
                    }
                    myannotq = myannotq.full[myfilter[, i]]
                    cat("annotation strata enrichment (levels=", 
                      Nq, ")\n", sep = "")
                    for (k in 1:Nq) {
                      myannotq_name <- paste("I", k, sep = "_")
                      cat("I", k, "\n")
                      cat("enrichment at level", k, ":\n")
                      logtemp_k <- myannotq == myannotq_levels[k]
                      logtemp_k[!is.finite(logtemp_k)] <- FALSE
                      Nk = sum(is.finite(tempz[logtemp_k]))
                      if (Nk > 1) {
                        if (myparameters$do_draw_bs_plots) {
                          cat("storing residuals..\n")
                          mycounts[, paste0("[=]", myannotq_name)] = mycounts[, 
                            paste0("[=]", myannotq_name)] + hist(mybsannot[logtemp_k], 
                            breaks = mybsannot_cuts, plot = F)$counts
                          tempaggregate = aggregate(tempresz[logtemp_k], 
                            list(mybsannot_bins[logtemp_k]), 
                            sum, na.rm = T)
                          myresiduals[mybsannot.full_aggregate \%in\% 
                            as.character(tempaggregate[, 1]), 
                            paste0("[=]", myannotq_name), i] = tempaggregate[, 
                            2]
                        }
                        myenrichment["[=]", myannotq_name, i] <- mean(tempresz[logtemp_k]^2, 
                          na.rm = T) - 1
                        if (sum(is.finite(tempz) & !logtemp_k & 
                          myctrlannot >= 1) > 1) {
                          zF <- var.test(tempz[is.finite(tempz) & 
                            logtemp_k], tempz[is.finite(tempz) & 
                            !logtemp_k & myctrlannot >= 1], alt = "greater")
                          myenrichment["[=](F)", myannotq_name, 
                            i] <- zF$statistic
                          myenrichment["[=](pF)", myannotq_name, 
                            i] <- zF$p.value
                          myenrichment.denomdf["[=](F)", myannotq_name, 
                            i] <- zF$parameter[2]
                          myenrichment.numdf["[=](F)", myannotq_name, 
                            i] <- zF$parameter[1]
                        }
                      }
                      else cat("no data at this level.\n")
                      if (do_cumulative && k > 2) {
                        cat("enrichment below level", k, ":\n")
                        logtemp_k = FALSE
                        for (kk in 1:(k - 1)) logtemp_k = logtemp_k | 
                          myannotq == myannotq_levels[kk]
                        logtemp_k[!is.finite(logtemp_k)] <- FALSE
                        Nk = sum(is.finite(tempz[logtemp_k]))
                        if (Nk > 1) {
                          if (myparameters$do_draw_bs_plots) {
                            cat("storing residuals..\n")
                            mycounts[, paste0("[<]", myannotq_name)] = mycounts[, 
                              paste0("[<]", myannotq_name)] + 
                              hist(mybsannot[logtemp_k], breaks = mybsannot_cuts, 
                                plot = F)$counts
                            tempaggregate = aggregate(tempresz[logtemp_k], 
                              list(mybsannot_bins[logtemp_k]), 
                              sum, na.rm = T)
                            myresiduals[mybsannot.full_aggregate \%in\% 
                              as.character(tempaggregate[, 1]), 
                              paste0("[<]", myannotq_name), i] = tempaggregate[, 
                              2]
                          }
                          myenrichment["[<]", myannotq_name, 
                            filter.counter] <- mean(tempresz[logtemp_k]^2, 
                            na.rm = T) - 1
                          if (sum(is.finite(tempz) & !logtemp_k & 
                            myctrlannot >= 1) > 1) {
                            zF <- var.test(tempz[is.finite(tempz) & 
                              logtemp_k], tempz[is.finite(tempz) & 
                              !logtemp_k & myctrlannot >= 1], 
                              alt = "greater")
                            myenrichment["[<](F)", myannotq_name, 
                              filter.counter] <- zF$statistic
                            myenrichment["[<](pF)", myannotq_name, 
                              filter.counter] <- zF$p.value
                            myenrichment.denomdf["[<](F)", myannotq_name, 
                              i] <- zF$parameter[2]
                            myenrichment.numdf["[<](F)", myannotq_name, 
                              i] <- zF$parameter[1]
                          }
                        }
                        else cat("no data below this level.\n")
                      }
                      if (do_cumulative && k < Nq - 1) {
                        cat("enrichment above level", k, ":\n")
                        logtemp_k = FALSE
                        for (kk in (k + 1):Nq) logtemp_k = logtemp_k | 
                          myannotq == myannotq_levels[kk]
                        logtemp_k[!is.finite(logtemp_k)] <- FALSE
                        Nk = sum(is.finite(tempz[logtemp_k]))
                        if (Nk > 1) {
                          if (myparameters$do_draw_bs_plots) {
                            cat("storing residuals..\n")
                            mycounts[, paste0("[>]", myannotq_name)] = mycounts[, 
                              paste0("[>]", myannotq_name)] + 
                              hist(mybsannot[logtemp_k], breaks = mybsannot_cuts, 
                                plot = F)$counts
                            tempaggregate = aggregate(tempresz[logtemp_k], 
                              list(mybsannot_bins[logtemp_k]), 
                              sum, na.rm = T)
                            myresiduals[mybsannot.full_aggregate \%in\% 
                              as.character(tempaggregate[, 1]), 
                              paste0("[>]", myannotq_name), i] = tempaggregate[, 
                              2]
                          }
                          myenrichment["[>]", myannotq_name, 
                            filter.counter] <- mean(tempresz[logtemp_k]^2, 
                            na.rm = T) - 1
                          if (sum(is.finite(tempz) & !logtemp_k & 
                            myctrlannot >= 1) > 1) {
                            zF <- var.test(tempz[is.finite(tempz) & 
                              logtemp_k], tempz[is.finite(tempz) & 
                              !logtemp_k & myctrlannot >= 1], 
                              alt = "greater")
                            myenrichment["[>](F)", myannotq_name, 
                              filter.counter] <- zF$statistic
                            myenrichment["[>](pF)", myannotq_name, 
                              filter.counter] <- zF$p.value
                            myenrichment.denomdf["[>](F)", myannotq_name, 
                              i] <- zF$parameter[2]
                            myenrichment.numdf["[>](F)", myannotq_name, 
                              i] <- zF$parameter[1]
                          }
                        }
                        else cat("no data above this level.\n")
                      }
                    }
                  }
                  else cat("insufficient data (", length(tempz), 
                    ")\n", sep = "")
                }
            }
            if (myparameters$verbose) {
                if (!is.null(dumpfile)) 
                  save(myenrichment, file = dumpfile)
                print(myenrichment)
            }
            if (Nf == 1) {
                myenrichment.mean = myenrichment[, , 1]
                myenrichment.median = myenrichment[, , 1]
                myenrichment.max = myenrichment[, , 1]
                myenrichment.min = myenrichment[, , 1]
            }
            else {
                for (m in 1:dim(myenrichment)[1]) {
                  for (n in 1:dim(myenrichment)[2]) {
                    myenrichment.df[m, n] = sum(is.finite(myenrichment[m, 
                      n, ]), na.rm = T)
                    if (myenrichment.df[m, n] > 0) {
                      myenrichment.mean[m, n] = mean(myenrichment[m, 
                        n, ], na.rm = T)
                      myenrichment.median[m, n] = median(myenrichment[m, 
                        n, ], na.rm = T)
                      myenrichment.max[m, n] = max(myenrichment[m, 
                        n, ], na.rm = T)
                      myenrichment.min[m, n] = min(myenrichment[m, 
                        n, ], na.rm = T)
                      myenrichment.sd[m, n] = sd(myenrichment[m, 
                        n, ], na.rm = T)
                      myenrichment.which[m, n] = which.min(abs(myenrichment[m, 
                        n, ] - myenrichment.mean[m, n]))
                    }
                  }
                }
                for (barg in colnames(myenrichment.mean)) {
                  for (marg in snames) {
                    myfstat = paste0(marg, "(F)")
                    myenrichment.mean[paste0(marg, "(pF)"), barg] = pf(myenrichment.mean[myfstat, 
                      barg], myenrichment.numdf[myfstat, barg, 
                      myenrichment.which[myfstat, barg]], myenrichment.denomdf[myfstat, 
                      barg, myenrichment.which[myfstat, barg]], 
                      lower.tail = F)
                    for (tmpdimname in paste0(marg, regression.var)) {
                      if (tmpdimname \%in\% dimnames(myenrichment.mean)[[1]]) {
                        if (myenrichment.df[paste0(tmpdimname, 
                          "(t)"), barg] > 1) {
                          cat("computing approximate p-value ", 
                            paste0(tmpdimname, "(p)"), "..\n", 
                            sep = "")
                          if (myparameters$verbose) {
                            cat("df =", myenrichment.df[paste0(tmpdimname, 
                              "(t)"), barg], "\n")
                            cat("mean =", myenrichment.mean[paste0(tmpdimname, 
                              "(t)"), barg], "\n")
                            cat("stddev =", myenrichment.sd[paste0(tmpdimname, 
                              "(t)"), barg], "\n")
                          }
                          myenrichment.mean[paste0(tmpdimname, 
                            "(p)"), barg] = 2 * pt(abs(sqrt(dfcorr * 
                            myenrichment.df[paste0(tmpdimname, 
                              "(t)"), barg]) * myenrichment.mean[paste0(tmpdimname, 
                            "(t)"), barg]/myenrichment.sd[paste0(tmpdimname, 
                            "(t)"), barg]), myenrichment.df[paste0(tmpdimname, 
                            "(t)"), barg] - 1, lower.tail = FALSE)
                          myenrichment.max[paste0(tmpdimname, 
                            "(p)"), barg] = 2 * pt(abs(sqrt(dfcorr * 
                            myenrichment.df[paste0(tmpdimname, 
                              "(t)"), barg]) * myenrichment.max[paste0(tmpdimname, 
                            "(t)"), barg]), 1, lower.tail = FALSE)
                          myenrichment.min[paste0(tmpdimname, 
                            "(p)"), barg] = 2 * pt(abs(sqrt(dfcorr * 
                            myenrichment.df[paste0(tmpdimname, 
                              "(t)"), barg]) * myenrichment.min[paste0(tmpdimname, 
                            "(t)"), barg]), 1, lower.tail = FALSE)
                        }
                        else {
                          myenrichment.mean[paste0(tmpdimname, 
                            "(p)"), barg] = NA
                          myenrichment.max[paste0(tmpdimname, 
                            "(p)"), barg] = NA
                          myenrichment.min[paste0(tmpdimname, 
                            "(p)"), barg] = NA
                        }
                      }
                    }
                    for (tmpdimname in paste0(marg, pLR.names)) {
                      if (tmpdimname \%in\% dimnames(myenrichment)[[1]]) {
                        mygamma = mean(qchisq(myenrichment[tmpdimname, 
                          barg, ], 1, lower.tail = F), na.rm = T)
                        myenrichment.mean[tmpdimname, barg] = pgamma(mygamma, 
                          shape = Nf/2, scale = 2/Nf, lower.tail = F)
                      }
                    }
                  }
                }
            }
            if (myparameters$regdiag) 
                dev.off()
            myresiduals_mean = rep(NA, dim(myresiduals)[1])
            myresiduals_var = rep(NA, dim(myresiduals)[1])
            myresiduals_array = array(myresiduals[, length(resnames), 
                ], dim = c(myparameters$N_bins, Nf))
            myindex = apply(is.finite(myresiduals_array), 1, 
                any)
            if (sum(myindex) > 1 && myparameters$do_draw_bs_plots) {
                myresiduals_array = array(myresiduals_array[myindex, 
                  ], dim = c(sum(myindex), Nf))
                myresiduals_mean[myindex] = apply(myresiduals_array, 
                  1, sum, na.rm = T)/pmax(mycounts[myindex, length(resnames)], 
                  myparameters$tiny, na.rm = T)
                myresiduals_var[myindex] = apply((myresiduals_array - 
                  myresiduals_mean[myindex])^2, 1, sum, na.rm = T)/pmax(mycounts[myindex, 
                  length(resnames)], myparameters$tiny, na.rm = T)
            }
            return(list(regmean = myenrichment.mean, regmedian = myenrichment.median, 
                regmax = myenrichment.max, regmin = myenrichment.min, 
                resmean = myresiduals_mean, resvar = myresiduals_var))
        }
        else {
            cat("nothing to be done.\n")
        }
    }
    else {
        cat("no data available.\n")
    }
    return(NULL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
