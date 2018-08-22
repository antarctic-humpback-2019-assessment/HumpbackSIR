#include <Rcpp.h>
using namespace Rcpp;
//' Adjusted lognormal likelihood from Zerbini et al. 2011 (eq. 5)
//'
//' \code{dlnorm_zerb} returns the density for the adjusted log normal distribution whose logarith has mean equal to \code{meanlog} and standard deviation equal to \code{sdlog}:
//' $$f\left(x\right) = \frac{1}{\sigma * x} exp \left( \frac{ \left( ln\left( x \right) - \mu \right) ^ 2} {2 \sigma^2} \right)$$
//'
//' @param x vector of quantiles.
//' @param meanlog mean of the distribution on the log scale with a default value of 0.
//' @param sdlog standard deviation of the distribution on the log scale with a default value of 1.
//' @param log logical; if TRUE, probabilities p are given as log(p).
//'
//' @return A vector of densities.
// [[Rcpp::export]]
NumericVector dlnorm_zerb(
        NumericVector x,
        NumericVector meanlog,
        NumericVector sdlog,
        bool return_log = true ) {

    // 1. Setup
    int n_mu = meanlog.length();
    int n_sd = sdlog.length();
    int n_x = x.length();

    NumericVector like(n_x); // Vector of likelihoods

    // 2. Estimate likelihood
    if( (n_mu == 1) & ( n_sd == 1) ){ // 1 mu and 1 sdlog
        for(int i = 0; i < n_x; i++){
            like[i] = 1 / (sdlog[0] * x[i]) * exp(- pow( log( x[i] ) - meanlog[0], 2)/ (2 * sdlog[0] * sdlog[0] ) );
        }

    } else {
        if( (n_mu == n_x) & (n_sd == 1) ){ // Unique mu
            for(int i = 0; i < n_x; i++){
                like[i] = 1 / (sdlog[0] * x[i]) * exp(- pow( log( x[i] ) - meanlog[i], 2)/ (2 * sdlog[0] * sdlog[0] ) );
            }
        }
        else{
            if( (n_mu == n_x) & (n_sd == n_x) ){ // unique mu and sdlog
                for(int i = 0; i < n_x; i++){
                    like[i] = 1 / (sdlog[i] * x[i]) * exp(- pow( log( x[i] ) - meanlog[i], 2)/ (2 * sdlog[i] * sdlog[i] ) );
                }
            }
            else{// Errors
                if((n_sd > 1) & (n_sd != n_x) & (n_mu > 1) & (n_mu != n_x) ){ stop("Error: dimension mismatch between x, meanlog, and sdlog"); }
                if((n_sd > 1) & (n_sd != n_x) ){ stop("Error: dimension mismatch between x and sdlog"); }
                if((n_mu > 1) & (n_mu != n_x) ){ stop("Error: dimension mismatch between x and meanlog"); }
            }
        }
    }

    if(return_log == true){ // log transform
        like = log(like);
    }
    return like;
}
