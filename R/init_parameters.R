init_parameters <-
function ( ..., .listed=FALSE ) {

    if ( .listed ) {
        mylist = ( ... )
    } else mylist = list( ... )

    returnlist = list(

        helpme = FALSE,

        annotation.files = NULL,
        bsannotation.file = NULL,
        covannotation.files = NULL,
        ctrlannotation.file = "",
        genoctrlannotation.file = "",
        pheno.files = NULL,

        bidirect = FALSE,
        dfcorr = 0.8,
        distr = 'u',
        do_compute_fdr = FALSE,
        do_draw_bs_plots = FALSE,
        do_draw_fdr_plots = FALSE,
        do_draw_qq_plots = FALSE,
        do_draw_cis = FALSE,
        do_estimate_enrichment = FALSE,
        do_genomic_correction = FALSE,
        do_perform_bpt = FALSE,
        do_perform_hgt = FALSE,
        expect_rdata = FALSE,
        interact = FALSE,
        figure.colors = NULL,
        N_bins = 100,
        outprefix = 'ea',
        read.header = FALSE,
        regdiag = FALSE,
        
        verbose = FALSE,

        filter.names = c(),

        # p-value percentiles
        # myP_prctle = c( 0.001 ),
        # myP_prctle = c( 0.001, 0.01, 0.1 ),
        myP_prctle = c( 0.001, 0.01, 0.1, 0.5 ),

        tiny = 1.E-72,  # to prevent NANs of sort
        testp = 1.E-12,  # minimum p-value for significance test plots
        neglog10max = 12,  # negative logarithm of the minimum p-value
        N_ytics = 5,  # number of y-axis tics for significance test plots

        myFDR = 0.01,  # FDR significance level

        mybreaks = 'auto',  # interval break points lists
        numof.breaks = 4  # number of annotation quantiles

    )

    if ( length( mylist ) > 0 ) {
        argcnt = 0
        for ( k in 1 : length( mylist ) ) {
            if ( is.list( mylist[[k]] ) ) {
                returnlist = init_parameters( mylist[[k]], .listed=TRUE )
            } else {
                if ( length( names( mylist ) ) >= k && names( mylist )[k] != '' ) {
                    returnlist[[ names( mylist )[k] ]] = mylist[[k]]
                } else if ( length( names( returnlist ) ) - argcnt > 0 )
                    returnlist[[ argcnt + 1 ]] = mylist[[k]]
            }
            argcnt = argcnt + 1
        }
    }

    return( returnlist )

}

