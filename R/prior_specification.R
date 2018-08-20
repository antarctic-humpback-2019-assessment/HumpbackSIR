##' Generate a prior for use in the Humpback SIR
##'
##' Convenient way to generate a prior that can be used to generate samples for
##' the SIR.
##'
##' @param rfn Function to be used to generate samples from the prior
##'   distribution, such as \code{rnorm} or \code{runif}. Can also be set to a
##'   single number in order to hold that parameter constant.
##' @param par1 First parameter of the prior distribution. For \code{rnorm} this
##'   is the mean, for \code{runif} this is the lower bound.
##' @param par2 Second parameter of the prior distribution. For \code{rnorm}
##'   this is the standard deviation, for \code{runif} this is the upper bound.
##' @param use Boolean, indicates whether this prior is to be used. Defaults to
##'   TRUE.
##' @param label Text label of the distribution. If not set, will try to label
##'   based on the function passed in \code{rfn} Recognizes \code{rnorm}
##'   \code{rlnorm} \code{runif}, and \code{rlunif}. Otherwise will label
##'   "User defined".
##'
##' @return A \code{list} containing the function \code{rfn} for generating a
##'   sample from the prior distribution, a vector \code{pars} containing the
##'   parameters of the distribution, and a boolean \code{use} flag.
##' @export
##'
##' @examples
##' make_prior(rnorm, 0, 1, TRUE)
##' make_prior(runif, 0, 1, TRUE)
##' make_prior(rlunif, 0.01, 0.2, "Log-uniform(0.01, 0.2)")
make_prior <- function(rfn = NA, par1 = NULL, par2 = NULL, use = TRUE, label = NULL) {
  if (is.function(rfn)) {
    fn <- function() rfn(1, par1, par2)
    if (is.null(label)) {
      ## FIXME It would be nice to separate this out to a separate function, but
      ## `substitute` doesn't seem to work right if I'm passing through multiple
      ## functions, even if I specify `env = .GlobalEnv`.
      dist_lab <- switch(deparse(substitute(rfn)),
                        "rnorm"  = "Normal",
                        "rlnorm" = "Log-normal",
                        "runif"  = "Uniform",
                        "rlunif" = "Log-uniform",
                        "User defined")
      label <- paste0(dist_lab, "(", par1, ", ", par2, ")")
    }
  } else if (is.numeric(rfn)) {
    fn <- function() rfn
    par1 <- rfn
    label <- paste0("Constant(", rfn, ")")
  } else {
    fn <- function() NA
    label <- "Not used"
  }
  list(rfn = fn, pars = c(par1, par2), use = use, label = label)
}


##' This is a function to return a random sample from the log-uniform
##' distribution.
##'
##' @title The log-uniform distribution
##'
##' @param n Number of observations.
##' @param min,max Limits of the distribution (must be strictly positive).
##'
##' @return A vector of length \code{n} containing random draws from the
##'   log-uniform distribution.
##' @export
##' @examples
##' rlunif(1, 0.01, 0.2)
rlunif <- function(n, min = 1, max = 2) {
  exp(runif(n, log(min), log(max)))
}
