regress_enrichment <-
function (mydata.full, myregannot.full, ...)
{
    mincount = 5
    myparameters = init_parameters( ... )
    if ( !is.data.frame( myregannot.full ) )
        myregannot.full = data.frame( myregannot.full )
    myregnames = names( myregannot.full )
    cat("regressing on", sprintf( '%d', length( myregnames ) ), "covariates:\n")
    if ( myparameters$verbose ) print( myregnames )
    myname = myregnames[1]
    mycovnames = myregnames[-1]
    if ( myparameters$verbose ) {
      cat( 'regression variables summary:\n' )
      print( length( mydata.full ) ); print( dim( myregannot.full ) );
      print( class( mydata.full ) ); print( class( myregannot.full ) );
      print( typeof( mydata.full ) ); print( typeof( myregannot.full ) );
      print( summary( cbind( mydata.full, myregannot.full ) ) )
    }
    if ( dim( myregannot.full )[2] > 0 && sum( apply( !is.na( mydata.full ) & !is.na( myregannot.full ), 1, all ) ) - 1 > 1 ) {
          mydata = mydata.full[ is.finite( mydata.full ) & apply( !is.na( myregannot.full ), 1, all ) ]
          myregannot = myregannot.full[ is.finite( mydata.full ) & apply( !is.na( myregannot.full ), 1, all ), ]
          if ( is.null( dim( myregannot ) ) )
              myregannot = data.frame( defannot = array( myregannot, dim = c( length( myregannot ), 1 ) ) )
          if ( !is.data.frame( myregannot ) )
              myregannot = data.frame( myregannot )
          myflag = FALSE
          myflagZ = FALSE
          Ncoeff = 1
          Nl = rep( 1, length( myregnames ) )
          for ( k in 1 : length( myregnames ) ) {
            mytmpflag = FALSE
            cat( sprintf( 'checking feasibility for %s..', myregnames[k] ) )
            if ( is.logical( myregannot[,k] ) )
                myregannot[,k] = as.factor( myregannot[,k] )
            if ( is.factor( myregannot[,k] ) ) {
                mytmpflag = length( levels( myregannot[,k] ) ) == 1
                myflag = myflag || mytmpflag
                if ( mytmpflag )
                    cat( ' level degeneracy detected.' )
                if ( !mytmpflag ) {
                    mytmpflag = any( table( myregannot[,k] ) < mincount )
                    myflag = myflag || mytmpflag
                    if ( mytmpflag ) {
                        cat( ' scarcely populated factor level:\n' )
                        print( table( myregannot[,k] ) )
                    }
                }
                if ( !mytmpflag && myparameters$interact && is.factor( myregannot[,1] ) && k > 1 ) {
                    mytmpflag = any( table( interaction(myregannot[,1], myregannot[,k]) ) < mincount )
                    myflagZ = myflagZ || mytmpflag
                    if ( mytmpflag ) {
                        cat( ' scarcely populated interaction factor level:\n' )
                        print( table( interaction(myregannot[,1], myregannot[,k]) ) )
                    }
                }
            }
            if ( is.factor( myregannot[,k] ) )
                Nl[k] = length( levels( myregannot[,k] ) ) - 1
            if ( k > 1 && myparameters$interact ) {
                Ncoeff = Ncoeff + Nl[k] * (Nl[1]+1)
            } else Ncoeff = Ncoeff + Nl[k]
            cat( ifelse( mytmpflag, '\n', ' ok.\n' ) )
          }
          if ( !myflag ) {
              cat("fitting model..\n")
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#               # the following is better suited for array arguments, as opposed to data.frame;                       #
#               # no data.frame attachment is used within the lm call when using array arguments.                     #
#               myannot = myregannot[,1]                                                                              #
#               mycovannot = myregannot[,-1]                                                                          #
#               myannot = array( myregannot[,1], dim=c(length(mydata), 1) )                                           #
#               mycovannot = array( myregannot[,-1], dim=c(length(mydata), length(mycovnames)) )                      #
#                 formula.text <- "mydata~myregannot"                                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
              reg.coeff <- array(NA, dim = c(Ncoeff, 4))
              reg.coeffci <- array(NA, dim = c(Ncoeff, 2))
              pLR = NA
              pLR.names = NA
              formula.text <- paste0( "mydata~", paste0(myregnames, collapse='+') )
              cat('regression: ', formula.text, '\n')
              mydata.lm = lm(as.formula(formula.text), data=myregannot)
              tmpcoefficients = summary(mydata.lm)$coefficients
              tmpcoefficientsci = confint(mydata.lm)
              if ( myparameters$verbose ) {
                  print( mydata.lm$coefficients )
                  print( tmpcoefficients )
                  cat( 'reg.coeff: ' ); print( dim( reg.coeff ) )
                  cat( 'coefficients: ' ); print( length( mydata.lm$coefficients ) );
                  cat( 'summary coefficients: ' ); print( dim( tmpcoefficients ) );
              }
              reg.coeff.names <- c( paste( rep( names(mydata.lm$coefficients), each=dim(tmpcoefficients)[2]+2 ),
                rep( c( colnames(tmpcoefficients), 'CI_low', 'CI_high' ), length(mydata.lm$coefficients) ) ) )
              reg.coeff[which(is.finite(mydata.lm$coefficients)), ] <- tmpcoefficients
              if ( dim(tmpcoefficients)[1] == dim(tmpcoefficientsci)[1] ) {
                  reg.coeffci[which(is.finite(mydata.lm$coefficients)), ] <- tmpcoefficientsci
              } else cat('warning: coefficients dimensions mismatch.\n')
              koffset = sum(is.finite(mydata.lm$coefficients)) + 1
              formula.text <- "mydata~1"
              if ( length( myregannot ) > 1 )
                formula.text <- paste0( "mydata~", paste0(mycovnames, collapse='+') )
              mydata.lm0 = lm(as.formula(formula.text), data=myregannot)
              pLR = anova( mydata.lm0, mydata.lm, test = 'Chisq' )$'Pr(>Chi)'[-1]
              pLR.names <- paste( myname, 'pLR' )
              if ( myparameters$verbose ) { cat( 'pLR: ' ); print( pLR ); }
              if ( myparameters$regdiag ) { cat('plotting diagnostics..\n'); plot(mydata.lm, ask=FALSE) }
              if ( myparameters$interact && dim( myregannot )[2] > 1 ) {
                for ( k in 1 : length( mycovnames ) ) {
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#               # no data.frame attachment is used within the lm call when using array arguments.                     #
#                     formula.text <- "mydata~myregannot+myannot:myregannot[,-1]"                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                    formula.text <- paste0( "mydata~",
                        paste0(c(myregnames, paste(myname, mycovnames[k], sep=':')), collapse='+')
                    )
                    cat('regression with interaction: ', formula.text, '\n')
                    mydata.lmi = lm(as.formula(formula.text), data=myregannot)
                    pLR.names = c( pLR.names, paste0( myname, '*', mycovnames[k], ' pLR' ) )
                    tmpcoefficients = summary(mydata.lmi)$coefficients
                    tmpcoefficientsci = confint(mydata.lmi)
                    if ( myparameters$verbose ) {
                        print( mydata.lmi$coefficients )
                        print( tmpcoefficients )
                        cat( 'reg.coeff: ' ); print( dim( reg.coeff ) )
                        cat( 'coefficients: ' ); print( length( mydata.lmi$coefficients ) );
                        cat( 'summary coefficients: ' ); print( dim( tmpcoefficients ) );
                    }
                    # if the last Nl coefficients (the ones of interest) are defined, store them
                    kint = (length(mydata.lmi$coefficients)-Nl[1]*Nl[k+1]+1) : length(mydata.lmi$coefficients)
                    reg.coeff.names <- c( reg.coeff.names, paste(
                        rep( names(mydata.lmi$coefficients)[ kint ], each=dim(tmpcoefficients)[2]+2 ),
                        rep( c( colnames(tmpcoefficients), 'CI_low', 'CI_high' ), length(kint) ) ) )
                    if ( myflagZ ) {
                        pLR = c( pLR, NA )
                    } else {
                        if ( !any( kint %in% which( is.na(mydata.lmi$coefficients) ) ) ) {
                            for ( h in 0:(Nl[1]*Nl[k+1]-1) ) {
                                reg.coeff[ koffset + h, ] <-
                                    tmpcoefficients[dim(tmpcoefficients)[1] - Nl[1]*Nl[k+1] + h + 1, ]
                                if ( dim(tmpcoefficients)[1] == dim(tmpcoefficientsci)[1] ) {
                                    reg.coeffci[ koffset + h, ] <-
                                        tmpcoefficientsci[dim(tmpcoefficients)[1] - Nl[1]*Nl[k+1] + h + 1, ]
                                } else cat('warning: coefficients dimensions mismatch.')
                            }
                        }
                        pLR = c( pLR, anova( mydata.lm, mydata.lmi, test = 'Chisq' )$'Pr(>Chi)'[-1] )
                    }
                    koffset = koffset + Nl[1]*Nl[k+1]
                }
              }
              if ( myparameters$verbose ) {
                cat( 'returning:\n' )
                print( data.frame( names = c( reg.coeff.names, pLR.names ) ) )
                print( data.frame( coeff = c( t(cbind(reg.coeff, reg.coeffci)), pLR ) ) )
              }
              return( list(
                  reg.names = c( reg.coeff.names, pLR.names ),
                  reg.coeff = c( t(cbind(reg.coeff, reg.coeffci)), pLR ),
                  residuals = mydata.lm$residuals
              ))
          } else {
            cat('warning: factor degeneracy/scarsity detected; analysis dropped.\n')
            return(NULL)
          }
    } else {
      cat('data is broken; analysis dropped\n')
      return(NULL)
    }
}

