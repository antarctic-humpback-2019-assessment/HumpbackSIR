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
target_K <- function(K, param_sample, data, tspan) {
  param_sample$K <- K
  pred_N <- project_population(param_sample = param_sample,
                               data = data,
                               tspan = tspan)
  tail(pred_N$N, 1) - param_sample$N_obs
}

#' Find K given r_max and N_obs
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
find_K <- function(param_sample, tspan, data, control) {
  Kmin <- uniroot(target_K,
                  interval = control$K_bisect_lims,
                  tol = control$K_bisect_tol,
                  param_sample = param_sample,
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
calc_analytic_q <- function(trajectory, data, data_name) {
  ia_df <- merge(ia_data[[data_name]], trajectory,
                 by = "year", all.x = TRUE, all.y = FALSE)

  q_num <- sum((log(ia_df$obs / ia_df$N)) / (ia_df$sd ^ 2))
  q_den <- sum(ia_df$sd ^ -2)
  q_num / q_den
}

sample_params <- function(priors, data, control) {
  param_sample <- list()
  param_sample$r_max <- priors$r_max$rfn()
  param_sample$z <- priors$z$rfn()

  ## Sample N_obs from prior if necessary
  if (priors$N_obs$use) {
    param_sample$N_obs <- priors$N_obs$rfn()
  }

  ## Find K either using the prior or using N_obs
  if (priors$K$use) {
    param_sample$K <- priors$K$rfn()
  } else {
    param_sample$K <- find_K(param_sample = param_sample,
                             tspan = c(start_year, priors$N_obs$year),
                             data = data,
                             control = control)
  }

  pred_N <- project_population(param_sample = param_sample,
                               data = data)

  ## FIXME Better to specify an S3 type for each parameter, and add methods for
  ## random sampling vs analytic approximations somehow?
  other_pars <- setdiff(names(priors), c("r_max", "K", "z", "N_obs"))
  for (par in other_pars) {
    if (priors[[par]]$use) {
      param_sample[[par]] <- priors[[par]]$rfn()
    } else {
      data_name <- priors[[par]]$data_name
      param_sample[[par]] <- calc_analytic_q(pred_N, data, data_name)
    }
  }
  list(param_sample, pred_N)
}
