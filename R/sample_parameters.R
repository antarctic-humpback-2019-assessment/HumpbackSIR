#' Calculate absolute loss for a value of K given current population
#'
#' @param K Pre-expoitation population size in numbers or biomass (depending on
#'   input).
#' @param param_sample List of parameter values that contains named entries for
#'   at least \code{r_max} and \code{z}. Any value of \code{K} will be replaced
#'   by the argument \code{K} above.
#' @param target_N Target current population size.
#' @param data List of data that includes at least named entry \code{catch} of
#'   appropriate length.
#' @param tspan Vector of length two specifying start and end years. Start year
#'   will be pre-exploitation, when population is assumed at carrying capacity.
#'   End year will be the target year.
#'
#' @return Difference between predicted population and observed population.
#' @export
#'
#' @examples
#' K <- 600
#' param_sample <- list(r = 0.08, z = 2.39)
#' target_N <- 750
#' tspan <- c(1990, 2000)
#' data <- list(catch = runif(10, 10, 100))
#' target_K(K, param_sample, target_N, data, tspan)
target_K <- function(K, param_sample,
                     target_N, data,
                     tspan) {
  param_sample$K <- K
  pred_N <- project_population(param_sample = param_sample,
                               data = data,
                               tspan = tspan)
  pred_N$N[length(pred_N$N)] - target_N
}

#' LOGISTIC BISECTION
#'
#' Method of Butterworth and Punt (1995) where the prior distribution of the
#' current absolute abundance $N_{2005}$ and maximum net recruitment rate
#' \code{r_max} are sampled and then used to determine the unique value of the
#' population abundance $N$ in \code{start_Yr} (assumed to correspond to
#' carrying capacity $K$). Requires \code{\link{TARGET.K}} and subsequent
#' dependencies.
#'
#' @param K.low Lower bound for $K$ when preforming the bisection method of Punt
#'   and Butterworth (1995). Default is 1.
#' @param K.high Upper bound for $K$ when preforming the bisection method of
#'   Punt and Butterworth (1995). Default is 500,000.
#' @param r_max The maximum net recruitment rate ($r_{max}$).
#' @param z The parameter that determines the population size where productivity
#'   is maximum (assumed to be 2.39 by the IWC SC).
#' @param num_Yrs The number of projection years. Set as the last year in the
#'   catch or abundance series, whichever is most recent, minus the
#'   \code{start_Yr}.
#' @param start_Yr The first year of the projection (assumed to be the first
#'   year in the catch series).
#' @param target.Pop A sample of the prior on population abundance $N$, in
#'   numbers, set as \code{sample.N.obs} sampled from \code{priors$N.obs}
#' @param catches The time series of catch in numbers or biomass. Currently does
#'   not handle NAs and zeros will have to input a priori for years in which
#'   there were no catches.
#' @param MVP The minimum viable population size in numbers or biomass. Computed
#'   as 4 * \code{\link{num.haplotypes}} to compute minimum viable population
#'   (from Jackson et al., 2006 and IWC, 2007).
#' @param tol The desired accuracy (convergence tolerance) of
#'   \code{\link{stats::uniroot}}.
#'
#' @return A numeric scalar of an estimate of  carrying capacity $K$.
#'
#' @examples
#' LOGISTIC.BISECTION.K(K.low = 1, K.high = 100000, r_max = r_max, z = z,
#'                      num_Yrs = bisection.Yrs, start_Yr = start_Yr,
#'                      target.Pop = target.Pop, catches = catches, MVP = MVP,
#'                      tol = 0.001)
find_K <- function(param_sample, target_N, tspan, data, control) {
  Kmin <- uniroot(target_K,
                  interval = control$K_bisect_lims,
                  tol = control$K_bisect_tol,
                  param_sample = param_sample,
                  target_N = target_N,
                  data = data,
                  tspan = tspan)
  Kmin$root
}

#' Compute analytic estimates of q, the scaling parameter between indices and
#' absolute population size
#'
#' @param rel.Abundance Relative abundance index
#' @param add_CV Coefficient of variation
#' @param Pred_N Predicted population
#' @param start_Yr Initial year
#' @param num.IA Index of abundance
#'
#' @return A numeric estimator for $q$.
#' @export
#'
#' @examples
CALC.ANALYTIC.Q <- function(rel.Abundance, Pred_N, start_Yr,
                            add_CV = 0, num.IA) {
  ## Vector to store the q values
  analytic.Q <- rep(NA, num.IA)

  for (i in 1:num.IA) {
    ## Subseting across each index of abundance
    IA <- Rel.Abundance[Rel.Abundance$Index == i,]
    ## Years for which IAs are available
    IA.yrs <- IA$Year-start_Yr + 1
    ## Computing the value of sigma as in Zerbini et al. 2011
    IA$Sigma <- sqrt(log(1 + IA$CV.IA.obs^2))
    ## Numerator of the analytic q estimator (Zerbini et al., 2011 - eq. (3))
    qNumerator <- sum((log(IA$IA.obs / Pred_N[IA.yrs])) /
                      (IA$Sigma * IA$Sigma + add_CV * add_CV))
    ## Denominator of the analytic q estimator (Zerbini et al., 2011 - eq. (3))
    qDenominator <- sum(1 / (IA$Sigma * IA$Sigma))
    ## Estimate of q
    analytic.Q[i] <- exp(qNumerator / qDenominator)
  }
  analytic.Q
}
