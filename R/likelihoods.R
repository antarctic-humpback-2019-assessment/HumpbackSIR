##' Check data frame for required columns
##'
##' Data frames used as data for SIR likelihoods must have columns "year",
##' "obs", and "sd". If standard deviations are not available but coefficients
##' of variation are, the \code{\link{cv_to_sd}} function can be used.
##'
##' @param data List of data passed to a likelihood constructor.
##' @param data_name Name of the entry used by the particular likelihood.
##'
##' @return Nothing; throws an error if required columns are not present.
check_lik_data <- function(data, data_name) {
  if (is.null(data[[data_name]])) {
    stop("Data does not contain \"", data_name, "\"")
  }
  available_names <- names(data[[data_name]])
  required_names <- c("year", "obs", "sd")
  if (!all(required_names %in% available_names)) {
    stop("Data frame must have columns ",
         paste(required_names, collapse = ", "))
  }
  if (!all(complete.cases(data[[data_name]]))) {
    stop("Missing data not allowed.")
  }
  NULL
}

##' Construct log-likelihood for index of abundance
##'
##' @title Index of abundance likelihood function generator
##'
##' @param data_name String specifying the name of the data frame in \code{data}
##'   containing the relevant indices of abundance.
##' @param par_name Name of the parameter to use as \code{q} for this index.
##' @param dfn Density function (e.g. dnorm, dlnorm).
##' @param data List containing at least a data frame named \code{datname} with
##'   columns \code{year}, \code{obs}, and \code{sd}.
##' @param mean_link Mean link function (e.g. \code{log} if \code{dfn} is
##'   \code{dlnorm}). Defaults to \code{identity} for no transform.
##'
##' @return A function that takes arguments \code{trajectory} and
##'   \code{param_sample} and returns the log-likelihood of the observations
##'   given the parameter values.
construct_ia_loglik <- function(data_name, par_name, dfn, data,
                                mean_link = identity) {
  check_lik_data(data, data_name)
  index_data <- data[[data_name]]
  function(trajectory, param_sample) {
    q <- param_sample[[par_name]]
    ## Using merge here ensures that observations and predictions are associated
    ## by year even if they are not sequential.
    ia_df <- merge(index_data, trajectory, by = "year", all.x = TRUE, all.y = FALSE)
    ia_df$pred <- mean_link(q * ia_df$N)
    sum(mapply(dfn, ia_df$obs, ia_df$pred, ia_df$sd, log = TRUE))
  }
}

##' Construct likelihood for observations of absolute abundance
##'
##' @title Absolute abundance likelihood function generator
##'
##' @param data_name String specifying the name of the data frame in \code{data}
##'   containing the relevant indices of abundance.
##' @param dfn Density function (e.g. dnorm, dlnorm).
##' @param data List containing at least a data frame named \code{datname} with
##'   columns \code{year}, \code{obs}, and \code{sd}.
##' @param mean_link Mean link function (e.g. \code{log} if \code{dfn} is
##'   \code{dlnorm}). Defaults to \code{identity} for no transform.
##'
##' @return A function that takes arguments \code{trajectory} and
##'   \code{param_sample}, and returns the log-likelihood of the observed
##'   abundance given the parameter values.
construct_abs_loglik <- function(data_name, dfn, data, mean_link = identity) {
  check_lik_data(data, data_name)
  abs_data <- data[[data_name]]
  function(trajectory, param_sample) {
    ab_df <- merge(abs_data, trajectory,
                   by = "year", all.x = TRUE, all.y = FALSE)
    ab_df$pred <- mean_link(ab_df$N)
    sum(mapply(dfn, ab_df$obs, ab_df$pred, ab_df$sd, log = TRUE))
  }
}

check_growth_data <- function(data, data_name) {
  if (is.null(data[[data_name]])) {
    stop("Data does not contain \"", data_name, "\"")
  }
  available_names <- names(data[[data_name]])
  required_names <- c("start_year", "end_year", "obs", "sd")
  if (!all(required_names %in% available_names)) {
    stop("Data frame must have columns ",
         paste(required_names, collapse = ", "))
  }
  if (!all(complete.cases(data[[data_name]]))) {
    stop("Missing data not allowed.")
  }
  NULL

}

##' Construct likelihood for observations of absolute abundance
##'
##' @title Absolute abundance likelihood function generator
##'
##' @param data_name String specifying the name of the data frame in \code{data}
##'   containing the relevant indices of abundance.
##' @param par_name Name of the parameter to use as \code{q} for this index.
##' @param dfn Density function (e.g. dnorm, dlnorm).
##' @param data List containing at least a data frame named \code{data_name}
##'   with columns \code{start_year}, \code{end_year}, \code{obs}, and
##'   \code{sd}.
##' @param mean_link Mean link function (e.g. \code{log} if \code{dfn} is
##'   \code{dlnorm}). Defaults to \code{identity} for no transformation.
##'
##' @return A function that takes arguments \code{trajectory} and
##'   \code{param_sample}, and returns the log-likelihood of the observed value
##'   given the parameter values.
construct_growth_loglik <- function(data_name, dfn, data, mean_link = identity) {
  check_growth_data(data, data_name)
  growth_data <- data[[data_name]]
  function(trajectory, param_sample) {
    gr_df <- merge(growth_data, trajectory,
                   by.x = "start_year", by.y = "year",
                   all.x = TRUE, all.y = FALSE)
    gr_df <- merge(gr_df, trajectory,
                   by.x = "end_year", by.y = "year",
                   all.x = TRUE, all.y = FALSE)
    gr_df$pred_gr <- calc_growth_rate(gr_df$start_year, gr_df$end_year,
                                      gr_df$N.x, gr_df$N.y)
    gr_df$pred <- mean_link(gr_df$pred_gr)
    sum(mapply(dfn, gr_df$obs, gr_df$pred, gr_df$sd, log = TRUE))
  }

}
##' Gamma density parameterized by mean and standard deviation
##'
##' @title Gamma density parameterized by mean and standard deviation
##' @param x Observation
##' @param mean
##' @param sd
##' @return Density
dgamma_meansd <- function(x, mean, sd) {
  rate <- mean / sd ^ 2
  shape <- mean * rate
  dgamma(x, shape = shape, rate = rate)
}

##' Convert CV to standard deviation for lognormal
##'
##' @title CV to SD
##' @param cv Coefficient of variation
##' @return Numeric standard deviation
cv_to_sd <- function(cv) {
  sqrt(log(1 + cv ^ 2))
}

#' LOG LIKELIHOOD OF ABSOLUTE ABUNDANCE
#'
#' This function computes two estimates of the log-likelihood of the estimated
#' absolute abundance using the equation from Zerbini et al. 2011 (eq. 4) and a
#' lognormal distribution from \code{\link{CALC.LNLIKE}}.
#'
#' @param Obs.N Observed absoluted abundance in numbers as a data.frame
#'   containing year, estimate of absolute abundance, and CV.
#' @param Pred_N Predicted absolute abundance in numbers from
#'   \code{\link{project_population}}.
#' @param start_Yr The first year of the projection (assumed to be the first
#'   year in the catch series).
#' @param add_CV Additional CV to add to variance of lognormal distribution
#'   sampled from \code{priors$add_CV}.
#' @param log Return the log of the likelihood (TRUE/FALSE)
#'
#' @return A list of two numeric scalars of estimates of log-likelihood.
#'
#' @examples
#' Obs.N  <-  data.frame(Year = 2005, Sigma = 5, Obs.N = 1000)
#' Pred_N  <-  1234
#' start_Yr  <-  2005
#' LNLIKE.Ns(Obs.N, Pred_N, start_Yr, add_CV = 0, log=TRUE)
LNLIKE.Ns <- function(Obs.N, Pred_N, start_Yr, add_CV, log = TRUE) {
  loglike.Ns1 <- 0
  loglike.Ns2 <- 0

  ## Years for which Ns are available
  N.yrs <- Obs.N$Year-start_Yr+1
  ## This is the likelihood from Zerbini et al. 2011 (eq. 4)
  loglike.Ns1 <- loglike.Ns1 +
    ((sum(log(Obs.N$Sigma) + log(Obs.N$N.obs) + 0.5 *
          ((((log(Pred_N[N.yrs]) - log(Obs.N$N.obs))^2) /
            (Obs.N$Sigma * Obs.N$Sigma + add_CV * add_CV))))))
  ## This is the log-normal distribution from R (using function dnorm)
  ## FIXME See comments above re: `dlnorm`
  loglike.Ns2 <- loglike.Ns2 + CALC.LNLIKE(Obs.N = Obs.N$N.obs,
                                           Pred_N = (Pred_N[N.yrs]),
                                           CV = sqrt(Obs.N$Sigma * Obs.N$Sigma +
                                                     add_CV * add_CV),
                                           log = log)

  list(loglike.Ns1 = loglike.Ns1, loglike.Ns2 = loglike.Ns2)
}

#' Calculate the log-likelihood of the growth rate
#'
#' Calculates the log-likelihood of the estimated growth rate given the observed
#' growth rate and the standard deviation of the observed growth rate.
#'
#' @param Obs.GR Observed growth rate
#' @param Pred.GR Predicted growth rate
#' @param GR.SD.Obs Standard error of the observed growth rate
#'
#' @return A \code{list} containing \code{loglike.GR1} and \code{loglike.GR2}
#'
#' @examples
#' LNLIKE.GR(0.1, 0.1, 0.1)
LNLIKE.GR <- function(Obs.GR, Pred.GR, GR.SD.Obs) {
  ## TODO Does this need to recalculate and return *both* of these values?
  loglike.GR1 <- 0
  loglike.GR2 <- 0

  ## This is the likelihood from Zerbini et al. 2011 (eq. 6)
  loglike.GR1 <- loglike.GR1 + (((log(GR.SD.Obs) + 0.5 * (((Pred.GR-Obs.GR) / GR.SD.Obs)^2))))

  ## loglike.GR2 <- loglike.GR2 + CALC.LNLIKE(Obs.N = Obs.GR,
  ##                                          Pred_N = Pred.GR,
  ##                                          CV = GR.SD.Obs,
  ##                                          log = FALSE)

  list(loglike.GR1 = loglike.GR1, loglike.GR2 = loglike.GR2)
}

#' Function to calculate the log-likelihood using a lognormal distribution
#'
#' @param Obs.N Time series of observed abundance
#' @param Pred_N Time series of estimated abundance
#' @param CV coefficient of variation
#' @param log whether to export as log-likelihood
#'
#' @return returns a scalar of the likelihood
#'
#' @examples
#' Obs.N <- 2000
#' Pred_N <- 2340
#' CV <- 4
#' CALC.LNLIKE(Obs.N, Pred_N, CV)
CALC.LNLIKE <- function(Obs.N, Pred_N, CV, log = FALSE) {
  sum(dnorm(x = log(Obs.N), mean = log(Pred_N), sd = CV, log = log))
}
