perform_bpt <-
function (mydata.full = NULL, myannot.full = NULL, myctrlannot.full = NULL,
    myparameters = init_parameters(), myfilter.full = NULL, P_prctle = 0.01)
{
    if (!is.null(mydata.full) && !is.null(myannot.full)) {
        mymask <- c( !is.na(mydata.full) & !is.na(myannot.full) )
        if (is.null(myctrlannot.full))
            myctrlannot.full = array(TRUE, dim=c(length(mydata.full), 1))
        if (is.null(myfilter.full))
            myfilter.full = array(TRUE, dim=c(length(mydata.full), 1))
        do_cumulative = FALSE
        myannotq.full <- factorize_annotation(myannot.full, myparameters$mybreaks, myparameters$numof.breaks)
        myannotq_levels <- levels(as.factor(myannotq.full))
        if (any(myannotq.full[mymask] != myannot.full[mymask])) do_cumulative = TRUE
        P_thresh = quantile(mydata.full, P_prctle, na.rm = T)
        if (do_cumulative) {
            BPT_test <- array(dim = c(length(myannotq_levels), 3, dim(myfilter.full)[2]))
            mysummary.bpt <- array(dim = c(length(myannotq_levels), 3))
            mymin.bpt <- array(dim = c(length(myannotq_levels), 3))
            mymax.bpt <- array(dim = c(length(myannotq_levels), 3))
            colnames(mysummary.bpt) <- c("[=]", "[<]", "[>]")
            colnames(mymin.bpt) <- c("[=]", "[<]", "[>]")
            colnames(mymax.bpt) <- c("[=]", "[<]", "[>]")
        }
        else {
            BPT_test <- array(dim = c(length(myannotq_levels), dim(myfilter.full)[2]))
            mysummary.bpt <- array(dim = c(length(myannotq_levels), 1))
            mymin.bpt <- array(dim = c(length(myannotq_levels), 1))
            mymax.bpt <- array(dim = c(length(myannotq_levels), 1))
            colnames(mysummary.bpt) <- c("[=]")
            colnames(mymin.bpt) <- c("[=]")
            colnames(mymax.bpt) <- c("[=]")
        }
        rownames(mysummary.bpt) <- myannotq_levels
        rownames(mymin.bpt) <- myannotq_levels
        rownames(mymax.bpt) <- myannotq_levels
        myfilter = NULL
        if ( is.null( myfilter.full ) ) {
            myfilter = array( TRUE, dim = c( length(mydata.full), 1 ) )
        } else { myfilter = myfilter.full & mymask }
        cat( 'running tests..\n' )
        for ( i in 1 : dim(myfilter)[2] ) {
            mydata <- mydata.full[myfilter[,i]]
            myannot <- myannot.full[myfilter[,i]]
            myannotq <- myannotq.full[myfilter[,i]]
            myctrlannot <- myctrlannot.full[myfilter[,i]]
            if ( length(mydata) > 0
                && length(mydata) == length(myannot)
                && length(mydata) == length(myctrlannot)
            ) {
                BPT_ctrl <- sum(myctrlannot, na.rm = T)
                BPT_ctrl_signif <- sum((myctrlannot & mydata < P_thresh), na.rm = T)
                compute_test_statistic <- function(mytemp) {
                    BPT_positive <- length(mytemp)
                    BPT_positive_signif <- sum((mytemp < P_thresh), na.rm = T)
                    n_a <- BPT_ctrl
                    n_b <- BPT_positive
                    f_a <- BPT_ctrl_signif/(n_a + 0.01)
                    f_b <- BPT_positive_signif/(n_b + 0.01)
                    f_ab <- (f_a * n_a + f_b * n_b)/(n_a + n_b + 0.01)
                    if ( myparameters$verbose ) {
                        cat( 'n_a =', n_a, '\n' )
                        cat( 'n_b =', n_b, '\n' )
                        cat( 'f_a =', f_a, '\n' )
                        cat( 'f_b =', f_b, '\n' )
                        cat( 'f_ab =', f_ab, '\n' )
                    }
                    return((f_a - f_b)/sqrt(f_ab * (1 - f_ab) * (1/(n_a + 0.01) + 1/(n_b + 0.01))))
                }
                for (k in 1:length(myannotq_levels)) {
                    logtemp_k <- myannotq == myannotq_levels[k]
                    if (do_cumulative) {
                      BPT_test[k, 1, i] <- compute_test_statistic(mydata[logtemp_k])
                      if (k > 2) {
                        logtemp_k <- as.numeric(myannotq) < k
                        BPT_test[k, 2, i] <- compute_test_statistic(mydata[logtemp_k])
                      }
                      if (k < length(myannotq_levels) - 1) {
                        logtemp_k <- as.numeric(myannotq) > k
                        BPT_test[k, 3, i] <- compute_test_statistic(mydata[logtemp_k])
                      }
                    } else {
                      BPT_test[k, i] <- compute_test_statistic(mydata[logtemp_k])
                    }
                }
            }
        }
        for (k in 1:length(myannotq_levels)) {
            myannotq_name <- myannotq_levels[k]
            if ( myparameters$bidirect ) {
                if (do_cumulative) {
                    mysummary.bpt[myannotq_name, "[=]"] <- 2*pnorm(-abs(median(BPT_test[k, 1, ])))
                    mysummary.bpt[myannotq_name, "[<]"] <- 2*pnorm(-abs(median(BPT_test[k, 2, ])))
                    mysummary.bpt[myannotq_name, "[>]"] <- 2*pnorm(-abs(median(BPT_test[k, 3, ])))
                    mymin.bpt[myannotq_name, "[=]"] <- 2*pnorm(-abs(min(BPT_test[k, 1, ])))
                    mymin.bpt[myannotq_name, "[<]"] <- 2*pnorm(-abs(min(BPT_test[k, 2, ])))
                    mymin.bpt[myannotq_name, "[>]"] <- 2*pnorm(-abs(min(BPT_test[k, 3, ])))
                    mymax.bpt[myannotq_name, "[=]"] <- 2*pnorm(-abs(max(BPT_test[k, 1, ])))
                    mymax.bpt[myannotq_name, "[<]"] <- 2*pnorm(-abs(max(BPT_test[k, 2, ])))
                    mymax.bpt[myannotq_name, "[>]"] <- 2*pnorm(-abs(max(BPT_test[k, 3, ])))
                } else {
                    mysummary.bpt[myannotq_name, "[=]"] <- 2*pnorm(-abs(median(BPT_test[k, ])))
                    mymin.bpt[myannotq_name, "[=]"] <- 2*pnorm(-abs(min(BPT_test[k, ])))
                    mymax.bpt[myannotq_name, "[=]"] <- 2*pnorm(-abs(max(BPT_test[k, ])))
                }
            } else { 
                if (do_cumulative) {
                    mysummary.bpt[myannotq_name, "[=]"] <- pnorm(median(BPT_test[k, 1, ]))
                    mysummary.bpt[myannotq_name, "[<]"] <- pnorm(median(BPT_test[k, 2, ]))
                    mysummary.bpt[myannotq_name, "[>]"] <- pnorm(median(BPT_test[k, 3, ]))
                    mymin.bpt[myannotq_name, "[=]"] <- pnorm(min(BPT_test[k, 1, ]))
                    mymin.bpt[myannotq_name, "[<]"] <- pnorm(min(BPT_test[k, 2, ]))
                    mymin.bpt[myannotq_name, "[>]"] <- pnorm(min(BPT_test[k, 3, ]))
                    mymax.bpt[myannotq_name, "[=]"] <- pnorm(max(BPT_test[k, 1, ]))
                    mymax.bpt[myannotq_name, "[<]"] <- pnorm(max(BPT_test[k, 2, ]))
                    mymax.bpt[myannotq_name, "[>]"] <- pnorm(max(BPT_test[k, 3, ]))
                } else {
                    mysummary.bpt[myannotq_name, "[=]"] <- pnorm(median(BPT_test[k, ]))
                    mymin.bpt[myannotq_name, "[=]"] <- pnorm(min(BPT_test[k, ]))
                    mymax.bpt[myannotq_name, "[=]"] <- pnorm(max(BPT_test[k, ]))
                }
            }
        }
        return(list( summary = mysummary.bpt, min = mymin.bpt, max = mymax.bpt ))
    }
}
