#' Function to write tables of logistic model parameter and derived quantities for HumpbackSIR similar to Table 5 from Zerbini et al (2011).
#'
#' @param SIR Resample summary from hympbackSIR
#' @param file_name Desired filename to where csv file will be saved. If NULL, will not save.
zerbini_table <- function( SIR, file_name = NULL){

    # Vars of interest
    years <- sort(c( SIR$inputs$target.Yr, SIR$inputs$output.Years))
    vars <- c("r_max", "K", "z", "Nmin", paste0("N", years), "Max_Dep", paste0("status", years))
    vars_latex <- c("$r_{max}$", "$K$", "$z$", "$N_{min}$", paste0("$N_{", years, "}$"), "Max depletion", paste0("Depletion in ", years))
    pop_vars <- c("K", "Nmin", paste0("N", years))
    depletion_vars <- c("Max_Dep", paste0("status", years))

    results <- data.frame(matrix(NA, nrow = length(vars), ncol = 8))
    colnames(results) <- c("Parameter","Mean", "Median", "2.5% CI", "25% CI", "75% CI", "97.5% CI", "Unique")

    x <- SIR$resamples_output[,vars]

    # Get summary statistics
    results[,1] <- vars_latex
    results[,2] <- sapply(x, mean)
    results[,3:7] <- t(sapply(x, quantile, probs= c(0.5, 0.025, 0.25, 0.75, 0.975)))
    results[,8] <- sapply(x, function(x) length(unique(x)))

    # Format things
    results[c(1,3),2:7] <- round(results[c(1,3),2:7], 3)
    results[which(vars %in% depletion_vars),2:7] <- round(results[which(vars %in% depletion_vars),2:7], 3)
    results[which(vars %in% pop_vars),2:7] <- format(round(results[which(vars %in% pop_vars),2:7], 0),big.mark=",",scientific=FALSE)
    results[,8] <- format(round(results[,8], 0),big.mark=",",scientific=FALSE)

    if(!is.null(file_name)){
      write.csv(results, file = paste0(file_name, "_zerbini_table_5.csv"))
    }
    return(results)
}


#' Bayes factor
#'
#' @param SIR SIR Fit model or list of SIR fit models
#' @param prior_probs prior probabilities of each model
#'
#' @return vector of bayes factors
bayes_factor <- function( SIR , prior_probs = NULL){
    # If it is a single SIR, make into a list
    if(class(SIR) == "SIR"){
        stop("Error: only one SIR model provided")
    }

    if(is.null(prior_probs)){
        prior_probs <- rep(1/length(SIR), length(SIR)) # Make uniform prior probs
    }

    # Get average likelihoods
    data_probs <- sapply(SIR, function(x) sum(x$resamples_output$Likelihood)/ length(x$resamples_output$Likelihood) )

    harmonic_mean <- sapply(SIR, function(x) (sum(x$resamples_output$Likelihood ^ -1)/ length(x$resamples_output$Likelihood)) ^-1 ) # Should not be used, facvors parameter rich models

    post_model <- data_probs * prior_probs

    bayes_factor <- post_model/sum(post_model)
    return(bayes_factor)
}




#' Weighted SIR model
#'
#' Function to create a weighted model using bayes factors
#'
#' @param SIR List of SIR fit models
#' @param bayes_factor vector of bayes_factors from \code{\link{bayes_factor}}
#'
#' @return an object of class "SIR"
#' @export
weight_model <- function(SIR, bayes_factor){

    weighted_SIR <- SIR[[1]]

    names_output <- colnames(weighted_SIR$resamples_output)

    # Size of vector
    n_samples <- nrow( weighted_SIR$resamples_trajectories )
    subs_samples <- rmultinom(n = 1, size = n_samples, prob = bayes_factor) # How many values to take from each SIR

    sample_ind <- 1
    for(i in 1:length(SIR)){

        names_resample_sub <- names(SIR[[i]]$resamples_output)
        random_rows <- sample(1:n_samples, size = subs_samples[i], replace = FALSE)

        weighted_SIR$resamples_output[sample_ind:(sample_ind + subs_samples[i] - 1), names_output[names_output %in% names_resample_sub] ] <- SIR[[i]]$resamples_output[random_rows,names_output[names_output %in% names_resample_sub]]
        weighted_SIR$resamples_trajectories[sample_ind:(sample_ind + subs_samples[i] - 1),] <- SIR[[i]]$resamples_trajectories[random_rows, ]
        weighted_SIR$catch_trajectories[sample_ind:(sample_ind + subs_samples[i] - 1),] <- SIR[[i]]$catch_trajectories[random_rows, ]
        sample_ind <- sample_ind + subs_samples[i]
    }

    return(weighted_SIR)
}