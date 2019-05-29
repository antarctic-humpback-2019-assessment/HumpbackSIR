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
            ## FIXME It would be nice to separate this out to its own function, but
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


##' Return a random sample from the log-uniform distribution.
##'
##' @title The log-uniform distribution
##'
##' @param n Number of observations.
##' @param min,max Limits of the distribution (must be strictly positive). These
##'   will be the limits of the samples returned.
##'
##' @return A vector of length \code{n} containing random draws from the
##'   log-uniform distribution.
##'
##' Draws independent samples from Uniform(log(min), log(max)) and
##' exponentiates.
##'
##' @export
##' @examples
##' rlunif(1, 0.01, 0.2)
rlunif <- function(n, min = 1, max = 2) {
    exp(runif(n, log(min), log(max)))
}

##' @title Make a list of priors to be passed to the SIR function.
##'
##' @param r_max Population growth rate; defaults to Uniform(0, 0.106).
##' @param K Carrying capacity, defaults to unused, and solution found using
##'   recent observation and sampled \code{r_max}.
##' @param N_obs Prior on a recent abundance estimate. Defaults to Uniform(500,
##'   20,000).
##' @param add_CV Defaults to unused. Additional variability.
##' @param premodern_catch_sample Defaults to unused. Samples between the minimum and maximum premodern catch values.
##' @param z Defaults to constant 2.39. Shape parameter for generalized logistic
##'   population dynamics function.
##' @param q_IA Defaults to unused. Prior on q for indices of abundance. If
##'   \code{use = FALSE}, an analytic solution for q is used.
##' @param q_count Defaults to unused. Prior for q on counts.
##'
##' @return A named list containing each of the specified priors in a form that
##'   can be used by the SIR function.
make_prior_list <- function(r_max = make_prior(runif, 0, 0.106),
                            K = make_prior(use = FALSE),
                            N_obs = make_prior(runif, 500, 20000),
                            add_CV = make_prior(use = FALSE),
                            premodern_catch_sample = make_prior(runif, 0, 1),
                            z = make_prior(2.39),
                            q_IA = make_prior(use = FALSE),
                            q_count = make_prior(use = FALSE)) {
    list(r_max = r_max,
         K = K,
         N_obs = N_obs,
         add_CV = add_CV,
         premodern_catch_sample = premodern_catch_sample,
         z = z,
         q_IA = q_IA,
         q_count = q_count)
}

#' Probability density for correction factor as determined by Best et al 2010
#'
#' @param n number of samples
#' @param crit_95 Struck and loss rate such that there is only a 5% probability of struck and loss rate exceeding this value
#' @param upper Upper limit for struck and loss rate
#'
#' @return
#' @export
#'
#' @examples
rbest <- function(n = 1, crit_95 = 0.139, upper = 0.30){
    dom <- runif(n)
    sapply(dom, function(d) ifelse(d < 0.95, runif(1, 1, 1/(1-crit_95)), runif(1, 1/(1-crit_95), 1/(1-upper))))
}