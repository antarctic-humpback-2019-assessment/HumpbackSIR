#' Function to write tables of logistic model parameter and derived quantities for HumpbackSIR similar to Table 5 from Zerbini et al (2011).
#'
#' @param SIR Resample summary from hympbackSIR
#' @param file_name Desired filename to where csv file will be saved
zerbini_table <- function( SIR, file_name = "NULL"){

    # Vars of interest
    years <- c( SIR$inputs$target.Yr, SIR$inputs$output.Years)
    vars <- c("r_max", "K", "Nmin", paste0("N", years), "Max_Dep", paste0("status", years))
    vars_latex <- c("$r_{max}$", "$K$", "$N_{min}$", paste0("$N_{", years, "}$"), "Max depletion", paste0("Depletion in ", years))
    pop_vars <- c("K", "Nmin", paste0("N", years))
    depletion_vars <- c("Max_Dep", paste0("status", years))

    results <- data.frame(matrix(NA, nrow = length(vars), ncol = 7))
    colnames(results) <- c("Parameter","Mean", "Median", "2.5% CI", "5% CI", "95% CI", "97.5% CI")

    x <- SIR$resamples_output[,vars]

    # Get summary statistics
    results[,1] <- vars_latex
    results[,2] <- sapply(x, mean)
    results[,3:7] <- t(sapply(x, quantile, probs= c(0.5, 0.025, 0.05, 0.975, 0.95)))

    # Format things
    results[1,2:7] <- round(results[1,2:7], 3)
    results[which(vars %in% depletion_vars),2:7] <- round(results[which(vars %in% depletion_vars),2:7], 3)
    results[which(vars %in% pop_vars),2:7] <- format(round(results[which(vars %in% pop_vars),2:7], 0),big.mark=",",scientific=FALSE)

    write.csv(results, file = paste0(file_name, "_zerbini_table_5.csv"))
    return(results)
}
