
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
plot_sir <- function(SIR, file_name = "NULL") {

    # Extract SIR objects
    x <- SIR$resamples_trajectories
    abs.abundance <- SIR$inputs$abs.abundance
    rel.abundance <- SIR$inputs$rel.abundance
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

    if(!is.null(rel.abundance)){
        rel.abundance$Upper95 <- qlnorm(0.975, log(rel.abundance$IA.obs), rel.abundance$Sigma)
        rel.abundance$Lower95 <- qlnorm(0.025, log(rel.abundance$IA.obs), rel.abundance$Sigma)
    }


    # Plot trajectory
    for(i in 1:2){
        if(i == 1){
            filename <- paste0(file_name, "_trajectory_summary", ".tiff")
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


    # Plot IOA
    if(!is.null(rel.abundance)){
        IA.yrs <- rel.abundance$Year
        IA.yr.range <- c((min(IA.yrs) - 1):(max(IA.yrs) + 1)) # Range +- 1 of IOA years

        # Predict IOA
        N_hat <- SIR$resamples_trajectories[, paste0("N_", IA.yr.range)] # Estimates of N within IOA years
        q_cols <- grep("q_IA", colnames(SIR$resamples_output)) # Columns of resample Q estimates
        q_est <- SIR$resamples_output[, q_cols]


        IA_pread <- list()
        IA_summary <- list()
        for(i in 1:length(q_cols)){
            # Predict
            IA_pread[[i]] <- N_hat * q_est

            # Summarize
            IA_summary[[i]] <-  matrix(nrow = length(row_names), ncol = dim(IA_pread[[i]])[2])
            IA_summary[[i]][1, ] <- sapply(IA_pread[[i]], mean)
            IA_summary[[i]][2:6, ] <- sapply(IA_pread[[i]], quantile, probs= c(0.5, 0.025, 0.975, 0.05, 0.95))
            IA_summary[[i]][7, ] <- sapply(IA_pread[[i]], min)
            IA_summary[[i]][8, ] <- sapply(IA_pread[[i]], max)
            IA_summary[[i]][9, ] <- sapply(IA_pread[[i]], length)



            IA_summary[[i]] <- data.frame(IA_summary[[i]])
            names(IA_summary[[i]]) <- paste0("IA", i, "_", IA.yr.range)
            row.names(IA_summary[[i]]) <- row_names
        }



        ymax <- max(c(sapply(IA_summary, function(x) max(x[2:6,])), rel.abundance$Lower95, rel.abundance$Upper95))
        ymin <- 0


        for(j in 1:2){
            if(j == 1){
                filename <- paste0(file_name, "_IOA_fit", ".tiff")
                tiff( file = filename , width=169 / 25.4, height = 100 / 25.4, family = "serif", units = "in", res = 300)
            }

            par(mfrow = c(1,length(IA_summary)))

            # Loop through inices
            for(i in 1:length(IA_summary)){
                rel.abundance.sub <- rel.abundance[which(rel.abundance$Index == i),]

                # Plot configuration
                par( mar=c(3, 3 , 0.5 , 0.25) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))
                plot(y = NA, x = NA,
                     ylim = c(ymin, ymax),
                     xlim = c(min(IA.yr.range), max(IA.yr.range)),
                     xlab = "Year", ylab = "Relative abundance")


                # Credible interval
                polygon(
                    x = c(IA.yr.range, rev(IA.yr.range)),
                    y = c(IA_summary[[i]][3, ],rev(IA_summary[[i]][4, ])),
                    col = "Grey80", border = NA) # 95% CI

                polygon( x = c(IA.yr.range, rev(IA.yr.range)),
                         y = c(IA_summary[[i]][5, ], rev(IA_summary[[i]][6, ])),
                         col = "Grey60", border = NA) # 90% CI


                # Median and catch series
                lines( x = IA.yr.range, y = IA_summary[[i]][2, ], lwd = 3) # Median


                # Relative abundance
                points( x = rel.abundance.sub$Year,
                        y = rel.abundance.sub$IA.obs,
                        col = 1, pch = 16, cex = 2)
                arrows( x0 = rel.abundance.sub$Year,
                        y0 = rel.abundance.sub$Lower95,
                        x1 = rel.abundance.sub$Year,
                        y1 = rel.abundance.sub$Upper95,
                        length=0.05, angle=90, code=3, lwd = 3, col = 1)

            }
            if(j == 1){ dev.off()}
        }
    }

}
