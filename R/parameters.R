##' @title Make a list of priors for the catch multiplier to be passed to the SIR function.
##'
##' @param c_mult_1 Catch multiplier for the first catch period; defaults to 1 and used for the entire dataset
##' @param ... Additional catch multipliers as defined using \code{make_prior}. Must be supplied in order of use.
##'
##' @return A named list containing each of the specified priors in a form that
##'   can be used by the SIR function.
make_multiplier_list <- function(c_mult_1 = make_prior(1), ...) {

    add_mult_list <- list(...)
    mult_list <- list(c_mult_1 = c_mult_1)

    if(length(add_mult_list) > 0 ){
        mult_list <- c(mult_list, add_mult_list)
        names(mult_list) <- paste0("c_mult_", 1:length(mult_list))
    }
    mult_list
}



#' Get z from Nk
#'
#' @param Nk Depletion which maximizes yield.
#' @param b branch of lambert W function
#' @param maxiter Number of iterations to solve lambert W
#'
#' @export
#'
#' @examples
#'
#' getz(Nk = 0.6)
getz <- function(Nk = 0.6, b=-1, maxiter=100) {
    z = suppressWarnings(emdbook::lambertW(Nk*log(Nk), b=b, maxiter=maxiter) - log(Nk)) / log(Nk)
    return(z)
}