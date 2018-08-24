#' Calculate growth rate based on predicted population
#'
#' \code{calc_growth_rate} computes the predicted growth rate if such
#' information is available from an independent estimate rather than being
#' estimated from data. Growth rate is calculated as:
#'
#' $$r_{t_0 - t_{fin}}^{pred} = \frac{ \sum_{t = t_0} ^{t_{fin - 1}} ln \left(
#' \frac{N_{t+1}^{pred}} { N_t^{pred}} \right) } { t_{fin} - t_0 } = \frac{ ln
#' \left( N_{fin}^{pred} \right) - ln \left( N_{0}^{pred} \right)} { t_{fin} -
#' t_0 }$$
#'
#' where $N^{pred}$ is the model predicted population size, in numbers, at time
#' $t$ or $t+1$ in years, $t_0$ is the start year of the equation (1995 in
#' Zerbini et al. 2011), and $t_{fin}$ is the last year of the equation (1998 in
#' Zerbini et al. 2011).
#'
#' @param years Vector of length 2 indicating the start and end years of the
#'   growth rate observation.
#' @param pred_pop Predicted population in \code{years[1]} and \code{years[2]}.
#'
#' @return A numeric scalar representing predicted growth rate.
#'
#' @examples
#' gr_years <-  c(1995, 1998)
#' gr_pop <- c(1000, 2000)
#' calc_growth_rate(gr_years, gr_pop)
calc_growth_rate <- function(tspan, pred_pop) {
  (log(pred_pop[2]) - log(pred_pop[1])) / diff(tspan)
}

#' Computes the predicted rate of increase for a set of specified years for
#' comparison with trends estimated separately with any of the indices of
#' abundance or count data
#'
#' @param data Count data or relative abundance index to use
#' @param Pred_N Number of individuals predicted
#' @param start_Yr Initial year
#'
#' @return Vector of rates of increase, one per index
#' @export
#'
#' @examples
COMPUTING.ROI <- function(data = data, Pred_N = Pred_N, start_Yr = NULL) {
  ## TODO Can this be combined with the calc_growth_rate function? Or be used
  ## only in post-processing?
  num.indices <- max(data$Index)
  Pred.ROI <- rep(NA, num.indices)

  for (i in 1:num.indices) {
    ## index.ini.year <- (head(subset(data, Index == i)$Year, 1) - start_Yr)
    ## index.final.year <- (tail(subset(data, Index == i)$Year, 1) - start_Yr)
    ## elapsed.years <- index.final.year - index.ini.year
    start_year <- min(data$Year[data$Index == i])
    start_N <- Pred_N$N[Pred_N$year == start_year]
    end_year <- max(data$Year[data$Index == i])
    end_N <- Pred_N$N[Pred_N$year == start_year]

    Pred.ROI[i] <- exp((log(end_N) - log(start_N)) /
                       (end_year - start_year)) - 1
    ## Pred.ROI[i] <- exp((log(Pred_N$N[index.final.year]) -
    ##                     log(Pred_N$N[index.ini.year])) /
    ##                    (elapsed.years)) - 1
  }
  Pred.ROI
}
