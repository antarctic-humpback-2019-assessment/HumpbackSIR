
#' OUTPUT FUNCTION
#'
#' Function that provides a plot of the estimated population trajectory from a SIR outputs model including: median, 95%
#' credible interval, 90% credible interval, catch, and observed absolute abundance.
#'
#' @param SIR A fit SIR model
#' @param object Name of the model object as specified by the user.
#' @param file_name name of a file to identified the files exported by the
#'   function.
#'
#' @return Returns and saves a figure with the population trajectory.
plot_trajectory <- function(SIR, object = "USERDEFINED", file_name = "NULL") {

    # Change name if not supplied
    if(object == "USERDEFINED"){
        object == "trajectory_summary"
    }

    x <- SIR$resamples_trajectories
    abs.abundance <- SIR$inputs$abs.abundance
    catch.data <- SIR$inputs$catch.data
    N_priors <- SIR$inputs$priors_N_obs$pars

    row_names <- c("mean", "median",
                   "2.5%PI", "97.5%PI",
                   "5%PI", "95%PI",
                   "min", "max", "n")

    # Extract values
    output_summary <- matrix(nrow = length(row_names), ncol = dim(x)[2])
    output_summary[1, ] <- sapply(x, mean)
    output_summary[2:6, ] <- sapply(x, quantile, probs= c(0.5, 0.025, 0.975, 0.05, 0.95))
    output_summary[7, ] <- sapply(x, min)
    output_summary[8, ] <- sapply(x, max)
    output_summary[9, ] <- sapply(x, length)

    output_summary <- data.frame(output_summary)
    names(output_summary) <- names(x)
    row.names(output_summary) <- row_names


    # Get 95% CI and range
    Years <- as.numeric(gsub("N_", "", colnames(output_summary)))
    abs.abundance$Upper95 <- qlnorm(0.975, log(abs.abundance$N.obs), abs.abundance$Sigma)
    abs.abundance$Lower95 <- qlnorm(0.025, log(abs.abundance$N.obs), abs.abundance$Sigma)
    ymax <- max(c(max(output_summary[2:6, ]), abs.abundance$N.obs, abs.abundance$Lower95, abs.abundance$Upper95, N_priors))
    ymin <- 0


    # Plot it
    for(i in 1:2){
        if(i == 1){
        filename <- paste0(file_name, "_", object, ".tiff")
        tiff( file = filename , width=169 / 25.4, height = 100 / 25.4, family = "serif", units = "in", res = 300)
        }

        # Plot configuration
        par( mar=c(3, 3 , 0.5 , 0.25) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))
        plot(y = NA, x = NA,
             ylim = c(ymin, ymax),
             xlim = c(min(Years), max(Years)),
             xlab = "Year", ylab = "Number of individuals")

        # Credible interval
        polygon(
            x = c(Years,rev(Years)),
            y = c(output_summary[3, ],rev(output_summary[4, ])),
            col = "Grey80", border = NA) # 95% CI
        polygon( x = c(Years,rev(Years)),
                 y = c(output_summary[5, ], rev(output_summary[6, ])),
                 col = "Grey60", border = NA) # 90% CI

        # Median and catch series
        lines( x = Years, y = output_summary[2, ], lwd = 3) # Median
        lines( x = catch.data$Year, y = catch.data$Catch, lty = 2, lwd = 3, col = 1) # Catch

        # Absolute abundance
        points( x = abs.abundance$Year,
                y = abs.abundance$N.obs,
                col = 1, pch = 16, cex = 2)
        arrows( x0 = abs.abundance$Year,
                y0 = abs.abundance$Lower95,
                x1 = abs.abundance$Year,
                y1 = abs.abundance$Upper95,
                length=0.05, angle=90, code=3, lwd = 3)
        if(i == 1){ dev.off()}
    }
}
