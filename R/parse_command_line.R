parse_command_line <-
function ( myargs=NULL ) {

	parameters <- init_parameters()

	if ( length( myargs ) > 0 ) {

		opt.vector <- vector()

		index.vector <- 1 : length( myargs )
		cat( "\ncommand line arguments:\n  ", myargs, "\n" )

		sub.index.vector <- index.vector[ myargs == "--help" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$helpme <- TRUE
			opt.vector <- c(
				sub.index.vector
			)
		}
		sub.index.vector <- index.vector[ myargs == "--annot" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$annotation.files <- myargs[ sub.index.vector + 1 ] # array of elements following the flags
			cat( c( "annotations:", parameters$annotation.files ), sep = "\n " )
			parameters$mybreaks = "auto"
			opt.vector <- c(
				opt.vector,
				sub.index.vector,
				sub.index.vector + 1
			)
		}
		sub.index.vector <- index.vector[ myargs == "--bpt" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$do_perform_bpt <- TRUE
			cat( "binomial proportion test requested.\n" )
			opt.vector <- c(
				opt.vector,
				sub.index.vector
			)
		}
		sub.index.vector <- index.vector[ myargs == "--bpts" ]
		if ( length( sub.index.vector ) > 0 ) {
            parameters$bidirect <- TRUE
			parameters$do_perform_bpt <- TRUE
			cat( "bidirectional binomial proportion test requested.\n" )
			opt.vector <- c(
				opt.vector,
				sub.index.vector
			)
		}
		sub.index.vector <- index.vector[ myargs == "--cannot" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$covannotation.files <- myargs[ sub.index.vector + 1 ] # array of elements following the flags
			cat( c( "covariate annotations:", parameters$covannotation.files ), sep = "\n " )
			opt.vector <- c(
				opt.vector,
				sub.index.vector,
				sub.index.vector + 1
			)
		}
		sub.index.vector <- index.vector[ myargs == "--color" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$figure.colors <- myargs[ sub.index.vector + 1 ] # array of elements following the flags
			cat( c( "figure colors:", parameters$figure.colors ), sep = "\n " )
			opt.vector <- c(
				opt.vector,
				sub.index.vector,
				sub.index.vector + 1
			)
		}
		sub.index.vector <- index.vector[ myargs == "--ctrl" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$ctrlannotation.file <- myargs[ sub.index.vector[ length( sub.index.vector ) ] + 1 ]
			cat( parameters$ctrlannotation.file, "set as control annotation.\n" )
			opt.vector <- c(
				opt.vector,
				sub.index.vector,
				sub.index.vector + 1
			)
		}
		sub.index.vector <- index.vector[ myargs == "--dfcorr" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$dfcorr <- as.numeric( myargs[ sub.index.vector[ length( sub.index.vector ) ] + 1 ] )
			cat( c( "correction factor for non-independence set to:", parameters$dfcorr, "\n" ) )
			opt.vector <- c(
				opt.vector,
				sub.index.vector,
				sub.index.vector + 1
			)
		}
		sub.index.vector <- index.vector[ myargs == "--diag" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$regdiag <- TRUE
			cat( "regression diagnostics requested.\n" )
			opt.vector <- c(
				opt.vector,
				sub.index.vector
			)
		}
		sub.index.vector <- index.vector[ myargs == "--distr" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$distr <- myargs[ sub.index.vector[ length( sub.index.vector ) ] + 1 ]
			cat( c( "phenotypes distributions:", parameters$distr ), sep = "\n " )
			opt.vector <- c(
				opt.vector,
				sub.index.vector,
				sub.index.vector + 1
			)
		}
		sub.index.vector <- index.vector[ myargs == "--fdr" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$do_compute_fdr <- TRUE
			cat( "fdr requested.\n" )
			opt.vector <- c(
				opt.vector,
				sub.index.vector
			)
		}
		sub.index.vector <- index.vector[ myargs == "--eest" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$do_estimate_enrichment <- TRUE
			cat( "enrichment estimate requested.\n" )
			opt.vector <- c(
				opt.vector,
				sub.index.vector
			)
		}
		sub.index.vector <- index.vector[ myargs == "--fdrv" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$myFDR <- as.numeric( myargs[ sub.index.vector[ length( sub.index.vector ) ] + 1 ] )
			cat( c( "FDR set to:", parameters$myFDR, "\n" ) )
			opt.vector <- c(
				opt.vector,
				sub.index.vector,
				sub.index.vector + 1
			)
		}
		sub.index.vector <- index.vector[ myargs == "--fdrvsp" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$do_compute_fdr <- TRUE
			parameters$do_draw_fdr_plots <- TRUE
			cat( "FDR_VS_p plots requested.\n" )
			opt.vector <- c(
				opt.vector,
				sub.index.vector
			)
		}
		sub.index.vector <- index.vector[ myargs == "--filter" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$filter.names <- myargs[ sub.index.vector + 1 ] # array of elements following the flags
			cat( c( "filter arrays:", parameters$filter.names ), sep = "\n " )
			opt.vector <- c(
				opt.vector,
				sub.index.vector,
				sub.index.vector + 1
			)
		}
		sub.index.vector <- index.vector[ myargs == "--gc" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$do_genomic_correction <- TRUE
			parameters$genoctrlannotation.file <- myargs[ sub.index.vector[ length( sub.index.vector ) ] + 1 ]
			cat( "genomic control requested based on", parameters$genoctrlannotation.file, "\n" )
			opt.vector <- c(
				opt.vector,
				sub.index.vector,
				sub.index.vector + 1
			)
		}
		sub.index.vector <- index.vector[ myargs == "--header" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$read.header <- TRUE
			cat( "expecting header in input files..\n" )
			opt.vector <- c(
				opt.vector,
				sub.index.vector
			)
		}
		sub.index.vector <- index.vector[ myargs == "--hgt" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$do_perform_hgt <- TRUE
			cat( "Fisher test requested.\n" )
			opt.vector <- c(
				opt.vector,
				sub.index.vector
			)
		}
		sub.index.vector <- index.vector[ myargs == "--hgts" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$bidirect <- TRUE
			parameters$do_perform_hgt <- TRUE
			cat( "bidirectional Fisher test requested.\n" )
			opt.vector <- c(
				opt.vector,
				sub.index.vector
			)
		}
		sub.index.vector <- index.vector[ myargs == "--i" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$interact <- TRUE
			cat( "including interaction terms in regression..\n" )
			opt.vector <- c(
				opt.vector,
				sub.index.vector
			)
		}
		sub.index.vector <- index.vector[ myargs == "--nbins" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$N_bins <- myargs[ sub.index.vector[ length( sub.index.vector ) ] + 1 ]
			cat( c( "number of bins set to:", parameters$N_bins, "\n" ) )
			opt.vector <- c(
				opt.vector,
				sub.index.vector,
				sub.index.vector + 1
			)
		}
		sub.index.vector <- index.vector[ myargs == "--o" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$outprefix <- myargs[ sub.index.vector[ length( sub.index.vector ) ] + 1 ]
			cat( c( "output prefix set to:", parameters$outprefix, "\n" ) )
			opt.vector <- c(
				opt.vector,
				sub.index.vector,
				sub.index.vector + 1
			)
		}
		sub.index.vector <- index.vector[ myargs == "--prct" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$myP_prctle <- as.numeric( myargs[ sub.index.vector + 1 ] ) # array of elements following the flags
			cat( c( "p-value percentiles:", parameters$myP_prctle, "\n" ) )
			opt.vector <- c(
				opt.vector,
				sub.index.vector,
				sub.index.vector + 1
			)
		}
		sub.index.vector <- index.vector[ myargs == "--qq" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$do_draw_qq_plots <- TRUE
			cat( "qq-plots requested.\n" )
			opt.vector <- c(
				opt.vector,
				sub.index.vector
			)
		}
		sub.index.vector <- index.vector[ myargs == "--qqci" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$do_draw_cis <- TRUE
			cat( "confidence intervals requested.\n" )
			opt.vector <- c(
				opt.vector,
				sub.index.vector
			)
		}
		sub.index.vector <- index.vector[ myargs == "--r" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$expect_rdata <- TRUE
			cat( "expecting R data as input..\n" )
			opt.vector <- c(
				opt.vector,
				sub.index.vector
			)
		}
		sub.index.vector <- index.vector[ myargs == "--scatter" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$do_draw_bs_plots <- TRUE
			parameters$do_estimate_enrichment <- TRUE
			parameters$bsannotation.file <- myargs[ sub.index.vector[ length( sub.index.vector ) ] + 1 ]
			cat( "squared z-scores binned scatter plots requested.\n" )
            cat( "  scatter variable: ", parameters$bsannotation.file, "\n" )
			opt.vector <- c(
				opt.vector,
				sub.index.vector,
				sub.index.vector + 1
			)
		}
		sub.index.vector <- index.vector[ myargs == "--strbr" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$mybreaks <- myargs[ sub.index.vector[ length( sub.index.vector ) ] + 1 ]
			cat( c( "annotation strata breaking points:", parameters$mybreaks ), sep = "\n " )
			opt.vector <- c(
				opt.vector,
				sub.index.vector,
				sub.index.vector + 1
			)
		}
		sub.index.vector <- index.vector[ myargs == "--strn" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$numof.breaks <- myargs[ sub.index.vector[ length( sub.index.vector ) ] + 1 ]
			cat( c( " annotation strata quantiles:", parameters$numof.breaks ), sep = "\n " )
			opt.vector <- c(
				opt.vector,
				sub.index.vector,
				sub.index.vector + 1
			)
		}
		sub.index.vector <- index.vector[ myargs == "--testp" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$testp <- as.numeric( myargs[ sub.index.vector[ length( sub.index.vector ) ] + 1 ] )
			cat( c( "test p-value plot ceiling set to:", parameters$testp, "\n" ) )
			if ( parameters$tiny > parameters$testp ) parameters$tiny <- parameters$testp
			parameters$neglog10max <- -log10( parameters$testp )
			if( parameters$bootstrap > 0 ) parameters$neglog10max <- -log10( 1 / parameters$bootstrap )
			opt.vector <- c(
				opt.vector,
				sub.index.vector,
				sub.index.vector + 1
			)
		}
		sub.index.vector <- index.vector[ myargs == "--v" ]
		if ( length( sub.index.vector ) > 0 ) {
			parameters$verbose <- TRUE
			cat( "verbose mode on.\n" )
			opt.vector <- c(
				opt.vector,
				sub.index.vector
			)
		}
		if ( length( opt.vector ) == 0 ) { parameters$pheno.files <- myargs
		} else { parameters$pheno.files <- myargs[ -opt.vector ]
        }

	} else {

		cat( "\nyou may have neglected to provide some necessary input.\n" )
		parameters$helpme <- TRUE

	}

	return( parameters )

}

