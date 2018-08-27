##' Project population forward from code{N1} for \code{num_years} using
##' Pella-Tomlinson dynamics. Note that this function assumes the initial
##' population is at carrying capacity.
##'
##' @title Project population forward
##'
##' @param param_sample List of parameter values containing at least named
##'   entries for \code{r_max}, \code{K}, and \code{z}.
##' @param data List of data containing at least a data frame \code{catch} with
##'   columns \code{year} and \code{catch}.
##' @param tspan Vector of length two indicating starting and ending years for
##'   projection. If default \code{NULL} is used, the entries of \code{tspan}
##'   will be taken as the first and last rows in \code{data$catch$year}.
##'
##' @return A data frame with columns \code{year} and \code{N}.
##' @export
project_population <- function(param_sample, data, tspan = NULL) {
  if (is.null(tspan)) {
    tspan <- c(head(data$catch$year, 1),
               tail(data$catch$year, 1))
  }
  num_years <- diff(tspan) + 1
  ## TODO Add interface for catch multipliers (struck and loss rates etc. here)
  catch_series <- data$catch$catch
  ## Make sure the catch series starts at the right time
  if (min(data$catch$year[1] != tspan[1])) {
    stop("Catch series and tspan must start the same year")
  }
  ## Make sure the catch series is long enough
  if (length(catch_series) < (num_years - 1)) {
    stop("Catch series must cover all but the final year")
  }
  N <- generalized_logistic(param_sample[["r_max"]],
                            param_sample[["K"]],
                            param_sample[["K"]],
                            param_sample[["z"]],
                            num_years,
                            catch_series)
  data.frame(year = tspan[1]:tspan[2],
             N = N)
}

##' Check whether the population violates the minimum viable population (MVP)
##' constraint provided by the number of haplotypes found. MVP is typically
##' considered to be \code{4 * num_haplotypes} (from Jackson et al. 2006 and
##' IWC, 2007).
##'
##' @title Check for violation of MVP
##'
##' @param trajectory Data frame of population trajectory, as from
##'   \code{\link{project_population}}.
##' @param mvp Minimum viable population
##'
##' @return Logical; TRUE if population trajectory dips below MVP.
check_mvp_violated <- function(trajectory, mvp) {
  any(trajectory$N < mvp)
}

##' Find minimum population and year where it occured.
##'
##' @title Minimum population and year
##'
##' @param trajectory Data frame of population trajectory, as from
##'   \code{\link{project_population}}.
##'
##' @return A data frame with columns \code{year} and \code{N}, corresponding to
##'   the year where the minimum population occured and the population that year
##'   respectively.
minimum_population <- function(trajectory) {
  min_N <- min(trajectory$N)
  min_idx <- which(trajectory$N == min_N)
  trajectory[min_idx, ]
}
