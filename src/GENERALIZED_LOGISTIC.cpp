#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List GENERALIZED_LOGISTIC(double r_max, double K, double N1, double z, double start_Yr, double num_Yrs, NumericVector catches, double MVP ) {


    int VMVP = 0;     // Variable to indicate whether min population is reached
    NumericVector n_hat(num_Yrs);      // Create a vector to hold the model predicted population size
    n_hat[0] = N1;                     // The first year in the vector above is N1

    for (int t = 1; t < num_Yrs; t++){
        n_hat[t] = n_hat[t - 1] + r_max * n_hat[t - 1] * (1 - pow(n_hat[t - 1] / K, z) ) - catches[t - 1] ;
        if(n_hat[t] < 1){
            n_hat[t] = 1; // Make sure the population is positive
        } 
    }

    double n_min = min( n_hat ); // Compute Nmin
    int y_min = which_min( n_hat ) + start_Yr; //  Compute the year at which Nmin occurred

    if (n_min < MVP) {
        VMVP = 1; // Determine whether Nmin is below Min Viable Population
    }

    // Compile results
    return List::create( Named("Min_Pop") =  n_min , _["Min_Yr"] = y_min , _["Violate_Min_Viable_Pop"] = VMVP, _["Pred_N"] = n_hat);

}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and dev elopment). The R code will be automatically
// run after the compilation.
//

