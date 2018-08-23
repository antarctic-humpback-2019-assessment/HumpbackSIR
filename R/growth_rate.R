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
calc_growth_rate <- function(years, pred_pop) {
  (log(pred_pop[2]) - log(pred_pop[1])) / diff(years)
}
