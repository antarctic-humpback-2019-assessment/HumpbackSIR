#include <Rcpp.h>
using namespace Rcpp;
//' Generalized logistic (Pella-Tomlinson) population dynamics model
//'

//' This \code{generalized_logistic} returns the population projection using a
//' Pella-Tomlison population dynamics model: $$N_{t+1} =
//' N_{t}+N_{t}*r_{max}*\left[ 1 - \left( \frac{N_{t}}{K} \right) ^z \right] -
//' C_{t}$$ where $N$ is the population size at year $t$ or $t+1$, $r_{max}$ is
//' the maximum net recruitment rate, $K$ is the pre-exploitation population
//' size, $z$ is the parameter that determines the population size where
//' productivity is maximum. For example, a value of 2.39 corresponds to maximum
//' sustainable yield of $0.6K$ and is assumed by the IWC SC, and $C_t$ is the
//' harvest in numbers in year $t$. Population size can be in either numbers or
//' biomass, however, units will have to be the same as the units used for catch,
//' relative abundance, and absolute abundance.
//'
//' @param r_max The maximum net recruitment rate ($r_{max}$).
//' @param K Pre-exploitation population size in numbers or biomass (depending on
//'   input).
//' @param N1 The population size in numbers or biomass at year 1 (generally
//'   assumed to be K).
//' @param z The parameter that determines the population size where productivity
//'   is maximum (assumed to be 2.39 by the IWC SC).
//' @param num_Yrs The number of projection years. Set as the last year in the
//'   catch or abundance series, whichever is most recent, minus the
//'   \code{start_Yr}.
//' @param catches The time series of catch in numbers or biomass. Currently does
//'   not handle NAs and zeros will have to be input a priori for years in which
//'   there were no catches.
//' @param MVP The minimum viable population size in numbers or biomass. Computed
//'   as 4 * \code{\link{num.haplotypes}} to compute minimum viable population
//'   (from Jackson et al., 2006 and IWC, 2007).
//'
//' @return A numeric vector with a population level by year, starting at
//' \code{N1} and running for \code{num_years}.
//'
//' @examples
//' num_Yrs <- 10
//' start_Yr <- 1
//' r_max <- 0.2
//' K <- 1000
//' N1 <- K
//' catches <- round(runif(10, min = 0, max = 150 ), 0)
//' generalized_logistic(r_max, K, N1, z, start_Yr, catch_series)
// [[Rcpp::export]]
NumericVector generalized_logistic(double r_max,
                                   double K,
                                   double N1,
                                   double z,
                                   double num_years,
                                   NumericVector catch_series) {
  NumericVector N(num_years);
  N[0] = N1;                       // The first year in the vector above is N1

  for (int t = 1; t < num_years; t++) {
    N[t] = N[t - 1] +
      r_max * N[t - 1] * (1 - pow(N[t - 1] / K, z)) -
      catch_series[t - 1];
    if(N[t] < 1) {
      N[t] = 1;                    // Make sure the population is positive
    }
  }
  return N;
}

