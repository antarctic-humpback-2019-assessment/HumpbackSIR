#include <Rcpp.h>
using namespace Rcpp;
//' GENERALISED LOGISTIC MODEL
//'
//' \code{GENERALIZED_LOGISTIC} returns the population projection using a
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
//' @param start_Yr The first year of the projection (assumed to be the first
//'   year in the catch series).
//' @param num_Yrs The number of projection years. Set as the last year in the
//'   catch or abundance series, whichever is most recent, minus the
//'   \code{start_Yr}.
//' @param catches The time series of catch in numbers or biomass. Currently does
//'   not handle NAs and zeros will have to input a priori for years in which
//'   there were no catches.
//' @param MVP The minimum viable population size in numbers or biomass. Computed
//'   as 4 * \code{\link{num.haplotypes}} to compute minimum viable population
//'   (from Jackson et al., 2006 and IWC, 2007).
//'
//' @return A list of the minimum population size \code{\link{Min.Pop}}, year of the minimum population size \code{\link{Min.Yr}}, a indicator of wether the minimum population size is below the \code{\link{MVP}}, and the predicted population size \code{Pred.N}.
//'
//' @examples
//' num_Yrs  <-  10
//' start_Yr  <-  1
//' r_max  <-  0.2
//' K  <-  1000
//' N1  <-  K
//' catches  <-  round(runif(10, min = 0, max = 150 ), 0 )
//' MVP  <-  0
//' GENERALIZED_LOGISTIC(r_max, K, N1, z, start_Yr, num_Yrs, catches)
// [[Rcpp::export]]
List GENERALIZED_LOGISTIC(
        double r_max,
        double K,
        double N1,
        double z,
        double start_Yr,
        double num_Yrs,
        NumericVector catches,
        double MVP ) {

    // 1. Setup
    LogicalVector VMVP = false;         // Variable to indicate whether min population is reached
    NumericVector n_hat(num_Yrs);       // Create a vector to hold the model predicted population size
    n_hat[0] = N1;                      // The first year in the vector above is N1

    // 2. Run through population dynamics
    for (int t = 1; t < num_Yrs; t++){
        n_hat[t] = n_hat[t - 1] + r_max * n_hat[t - 1] * (1 - pow(n_hat[t - 1] / K, z) ) - catches[t - 1] ;
        if(n_hat[t] < 1){
            n_hat[t] = 1;               // Make sure the population is positive
        }
    }

    // 3. Summary
    double n_min = min( n_hat );        // Compute Nmin
    IntegerVector y_min ;               //  Compute the year at which Nmin occurred
    for(int t = 0; t < num_Yrs; t++) {
        if (n_hat[t] == n_min) y_min.push_back(t);
    }
    y_min = y_min + start_Yr;
    if (n_min < MVP) { VMVP = true; }   // Determine whether Nmin is below Min Viable Population

    // Compile results
    return List::create(
        Named("Min_Pop") =  n_min ,
        _["Min_Yr"] = y_min ,
        _["Violate_Min_Viable_Pop"] = VMVP,
        _["Pred_N"] = n_hat);
}
