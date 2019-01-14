
#' OUTPUT FUNCTION
#'
#' Function that provides a plot of the estimated population trajectory from a SIR outputs model including: median, 95%
#' credible interval, 90% credible interval, catch, and observed absolute abundance.
#'
#' @param SIR A fit SIR model
#' @param file_name name of a file to identified the files exported by the
#'   function. If NULL, does not save.
#' @param Reference A fit SIR model for the reference case
#' @param posterior_pred Logical. If true, includes a posterior predictive distribution of the estimated N
#'
#' @return Returns and saves a figure with the population trajectory.
plot_trajectory <- function(SIR, Reference = NULL, file_name = NULL, posterior_pred = TRUE) {


    # Extract SIR objects
    x <- SIR$resamples_trajectories
    abs.abundance <- SIR$inputs$abs.abundance
    rel.abundance <- SIR$inputs$rel.abundance
    catch.data <- SIR$catch_trajectories
    N_priors <- SIR$inputs$priors_N_obs$pars

    row_names <- c("mean", "median",
                   "2.5%PI", "97.5%PI",
                   "5%PI", "95%PI",
                   "min", "max", "n")

    # Extract N trajectories
    output_summary <- matrix(nrow = length(row_names), ncol = dim(x)[2])
    output_summary[1, ] <- sapply(x, mean)
    output_summary[2:6, ] <- sapply(x, quantile, probs= c(0.5, 0.025, 0.975, 0.05, 0.95))
    output_summary[7, ] <- sapply(x, min)
    output_summary[8, ] <- sapply(x, max)
    output_summary[9, ] <- sapply(x, length)

    output_summary <- data.frame(output_summary)
    names(output_summary) <- names(x)
    row.names(output_summary) <- row_names

    if(!is.null(Reference)){
        ref <- Reference$resamples_trajectories
        reference_summary <- matrix(nrow = length(row_names), ncol = dim(ref)[2])
        reference_summary[1, ] <- sapply(ref, mean)
        reference_summary[2:6, ] <- sapply(ref, quantile, probs= c(0.5, 0.025, 0.975, 0.05, 0.95))
        reference_summary[7, ] <- sapply(ref, min)
        reference_summary[8, ] <- sapply(ref, max)
        reference_summary[9, ] <- sapply(ref, length)
        reference_summary <- data.frame(reference_summary)
        names(reference_summary) <- names(ref)
        row.names(reference_summary) <- row_names
        ref_years <- as.numeric(gsub("N_", "", colnames(reference_summary)))
    }



    if(posterior_pred){
        N_posterior_pred <- matrix(NA, ncol = nrow(abs.abundance), nrow = nrow(x))

        for(i in 1:ncol(N_posterior_pred)){
            N_posterior_pred[,i] <- rlnorm(
                n = nrow(N_posterior_pred),
                meanlog = log( x[,paste0("N_", abs.abundance$Year[i]) ] ),
                sdlog = abs.abundance$Sigma[i] + SIR$resamples_output$add_CV[1]
            )
        }

        N_posterior_pred <- as.data.frame(N_posterior_pred)
        posterior_pred_summary <- matrix(nrow = length(row_names), ncol = dim(N_posterior_pred)[2])
        posterior_pred_summary[1, ] <- sapply(N_posterior_pred, mean)
        posterior_pred_summary[2:6, ] <- sapply(N_posterior_pred, quantile, probs= c(0.5, 0.025, 0.975, 0.05, 0.95))
        posterior_pred_summary[7, ] <- sapply(N_posterior_pred, min)
        posterior_pred_summary[8, ] <- sapply(N_posterior_pred, max)
        posterior_pred_summary[9, ] <- sapply(N_posterior_pred, length)

        colnames(posterior_pred_summary) <- paste0("Posterior_predictive_N_", abs.abundance$Year)
        rownames(posterior_pred_summary) <- row_names
    }

    # Extract catch trajectories
    catch_summary <- matrix(nrow = length(row_names), ncol = dim(catch.data)[2])
    catch_summary[1, ] <- sapply(catch.data, mean)
    catch_summary[2:6, ] <- sapply(catch.data, quantile, probs= c(0.5, 0.025, 0.975, 0.05, 0.95))
    catch_summary[7, ] <- sapply(catch.data, min)
    catch_summary[8, ] <- sapply(catch.data, max)
    catch_summary[9, ] <- sapply(catch.data, length)


    # Get 95% CI and range
    Years <- as.numeric(gsub("N_", "", colnames(output_summary)))
    abs.abundance$Upper95 <- qlnorm(0.975, log(abs.abundance$N.obs), abs.abundance$Sigma)
    abs.abundance$Lower95 <- qlnorm(0.025, log(abs.abundance$N.obs), abs.abundance$Sigma)
    ymax <- max(c(max(output_summary[2:6, ]), abs.abundance$N.obs, abs.abundance$Lower95, abs.abundance$Upper95, N_priors))
    ymin <- 0


    # Plot trajectory
    for(i in 1: (1 + !is.null(file_name))){
        if(i == 2){
            filename <- paste0(file_name, "_trajectory_summary", ".png")
            png( file = filename , width=169 / 25.4, height = 100 / 25.4, family = "serif", units = "in", res = 300)
        }

        # Plot configuration
        par( mar=c(3, 3 , 0.5 , 0.3) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))
        plot(y = NA, x = NA,
             ylim = c(ymin, ymax),
             xlim = c(min(Years), max(Years)),
             xlab = "Year", ylab = "Number of individuals")

        # N Trajectory
        # Credible interval
        polygon(
            x = c(Years,rev(Years)),
            y = c(output_summary[3, ],rev(output_summary[4, ])),
            col = "Grey80", border = NA) # 95% CI
        polygon( x = c(Years,rev(Years)),
                 y = c(output_summary[5, ], rev(output_summary[6, ])),
                 col = "Grey60", border = NA) # 90% CI

        # Median
        lines( x = Years, y = output_summary[2, ], lwd = 3) # Median

        # Reference median
        if(!is.null(Reference)){
            lines( x = ref_years, y = reference_summary[2, ], lwd = 3, lty = 3) # Median
        }

        # Catch trajectory
        # Credible interval
        polygon(
            x = c(Years,rev(Years)),
            y = c(catch_summary[3, ],rev(catch_summary[4, ])),
            col = "Grey80", border = NA) # 95% CI
        polygon( x = c(Years,rev(Years)),
                 y = c(catch_summary[5, ], rev(catch_summary[6, ])),
                 col = "Grey60", border = NA) # 90% CI

        # Median
        lines( x = Years, y = catch_summary[2, ], lty = 2, lwd = 3, col = 1) # Median


        # Absolute abundance
        points( x = abs.abundance$Year,
                y = abs.abundance$N.obs,
                col = 1, pch = 16, cex = 2)
        arrows( x0 = abs.abundance$Year,
                y0 = abs.abundance$Lower95,
                x1 = abs.abundance$Year,
                y1 = abs.abundance$Upper95,
                length=0.05, angle=90, code=3, lwd = 3)

        if(posterior_pred){
            points( x = abs.abundance$Year + 1,
                    y = posterior_pred_summary[2,],
                    col = adjustcolor("Grey80", alpha.f = 1) , pch = 16, cex = 2)
            arrows( x0 = abs.abundance$Year + 1,
                    y0 = as.numeric(posterior_pred_summary[3,]),
                    x1 = abs.abundance$Year + 1,
                    y1 = as.numeric(posterior_pred_summary[4,]),
                    length=0.05, angle=90, code=3, lwd = 3, col = adjustcolor("Grey80", alpha.f = 1))
        }

        if(i == 2){ dev.off()}
    }
}

#' OUTPUT FUNCTION
#'
#' Function that provides a plot of the estimated indices of abundance a SIR  model including: median, 95%
#' credible interval, 90% credible interval, catch, and observed indices of abundance abundance.
#'
#' @param SIR A fit SIR model
#' @param file_name name of a file to identified the files exported by the
#'   function. If NULL, does not save.
#' @param posterior_pred Logical. If true, includes a posterior predictive distribution of the estimated IOA
#' @param ioa_names names of indices of abundance used.
#'
#' @return Returns and saves a figure with the IOA trajectories.
plot_ioa <- function(SIR, file_name = NULL, ioa_names = NULL, posterior_pred = TRUE){

    rel.abundance <- SIR$inputs$rel.abundance
    row_names <- c("mean", "median",
                   "2.5%PI", "97.5%PI",
                   "5%PI", "95%PI",
                   "min", "max", "n")

    if(is.null(rel.abundance)){
        stop("SIR model did not include an IOA")
    }

    rel.abundance$Upper95 <- qlnorm(0.975, log(rel.abundance$IA.obs), rel.abundance$Sigma)
    rel.abundance$Lower95 <- qlnorm(0.025, log(rel.abundance$IA.obs), rel.abundance$Sigma)

    # Predict IOA
    q_cols <- grep("q_IA", colnames(SIR$resamples_output)) # Columns of resample Q estimates
    q_est <- SIR$resamples_output[, q_cols]
    q_est <- as.matrix(q_est, ncol = length(q_cols))


    # Setup objects
    ymax <- c()           # Maximum predicted IOA
    IA.yr.range <- list() # Year range for each IOA
    IA_pread <- list()
    IA_summary <- list()
    IA_posterior_pred <- list()
    IA_posterior_pred_sum <- list()

    # Predict and calculate summary
    for(i in 1:length(q_cols)){

        # Get IOA specifications
        rel.abundance.sub <- rel.abundance[which(rel.abundance$Index == i),]
        IA.yrs <- rel.abundance.sub$Year
        IA.yr.range[[i]] <- c((min(IA.yrs) - 1):(max(IA.yrs) + 1)) # Range +- 1 of IOA years

        # Predict
        N_hat <- SIR$resamples_trajectories[, paste0("N_", IA.yr.range[[i]])] # Estimates of N within IOA years
        IA_pread[[i]] <- N_hat * q_est[,i]

        # Summarize
        IA_summary[[i]] <-  matrix(nrow = length(row_names), ncol = dim(IA_pread[[i]])[2])
        IA_summary[[i]][1, ] <- sapply(IA_pread[[i]], mean)
        IA_summary[[i]][2:6, ] <- sapply(IA_pread[[i]], quantile, probs= c(0.5, 0.025, 0.975, 0.05, 0.95))
        IA_summary[[i]][7, ] <- sapply(IA_pread[[i]], min)
        IA_summary[[i]][8, ] <- sapply(IA_pread[[i]], max)
        IA_summary[[i]][9, ] <- sapply(IA_pread[[i]], length)

        IA_summary[[i]] <- data.frame(IA_summary[[i]])
        names(IA_summary[[i]]) <- paste0("IA", i, "_", IA.yr.range[[i]])
        row.names(IA_summary[[i]]) <- row_names


        # Get posterior predictive
        if(posterior_pred){
            IA_posterior_pred[[i]] <- matrix(NA, nrow = nrow(SIR$resamples_trajectories), ncol = length(IA.yrs))
            IA_posterior_pred_sum[[i]] <-  matrix(nrow = length(row_names), ncol = length(IA.yrs))

            for(j in 1:length(IA.yrs)){
                IA_posterior_pred[[i]][,j] <- rlnorm(
                    n = nrow(IA_posterior_pred[[i]]),
                    meanlog = log( q_est[,rel.abundance.sub$Index[j]] * SIR$resamples_trajectories[, paste0("N_", IA.yrs[j])] ),
                    sdlog = rel.abundance.sub$Sigma[j] + SIR$resamples_output$add_CV[1])
            }

            IA_posterior_pred[[i]] <- data.frame(IA_posterior_pred[[i]])

            # Summarize
            IA_posterior_pred_sum[[i]] <-  matrix(nrow = length(row_names), ncol = dim(IA_posterior_pred[[i]])[2])
            IA_posterior_pred_sum[[i]][1, ] <- sapply(IA_posterior_pred[[i]], mean)
            IA_posterior_pred_sum[[i]][2:6, ] <- sapply(IA_posterior_pred[[i]], quantile, probs= c(0.5, 0.025, 0.975, 0.05, 0.95))
            IA_posterior_pred_sum[[i]][7, ] <- sapply(IA_posterior_pred[[i]], min)
            IA_posterior_pred_sum[[i]][8, ] <- sapply(IA_posterior_pred[[i]], max)
            IA_posterior_pred_sum[[i]][9, ] <- sapply(IA_posterior_pred[[i]], length)

            IA_posterior_pred_sum[[i]] <- data.frame(IA_posterior_pred_sum[[i]])
            names(IA_posterior_pred_sum[[i]]) <- paste0("IA", i, "_", IA.yrs)
            row.names(IA_posterior_pred_sum[[i]]) <- row_names
        }

        # Get plot limits
        ymax[i] <- max(unlist(c(IA_summary[[i]][2:6,], rel.abundance.sub$Lower95, rel.abundance.sub$Upper95))) # Max of posterior
        if(posterior_pred){
            ymax[i] <- max(unlist(c(ymax[i], IA_posterior_pred_sum[[i]][2:6,], rel.abundance.sub$Lower95, rel.abundance.sub$Upper95))) # Max of posterior predictive
        }
    }

    ymin <- 0


    for(j in 1:(1 + !is.null(file_name))){
        if(j == 2){
            filename <- paste0(file_name, "_IOA_fit", ".png")
            png( file = filename , width=169 / 25.4, height = 100 / 25.4, family = "serif", units = "in", res = 300)
        }

        par(mfrow = c(1,length(IA_summary)), mar=c(3, 3 , 0.5 , 0.3) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))

        # Loop through inices
        for(i in 1:length(IA_summary)){
            rel.abundance.sub <- rel.abundance[which(rel.abundance$Index == i),]

            # Plot configuration
            plot(y = NA, x = NA,
                 ylim = c(ymin, ymax[i]),
                 xlim = c(min(IA.yr.range[[i]]), max(IA.yr.range[[i]])),
                 xlab = "Year", ylab = "Relative abundance")


            # Credible interval
            polygon(
                x = c(IA.yr.range[[i]], rev(IA.yr.range[[i]])),
                y = c(IA_summary[[i]][3, ],rev(IA_summary[[i]][4, ])),
                col = "Grey80", border = NA) # 95% CI

            polygon( x = c(IA.yr.range[[i]], rev(IA.yr.range[[i]])),
                     y = c(IA_summary[[i]][5, ], rev(IA_summary[[i]][6, ])),
                     col = "Grey60", border = NA) # 90% CI


            # Median and catch series
            lines( x = IA.yr.range[[i]], y = IA_summary[[i]][2, ], lwd = 3) # Median


            # Relative abundance
            points( x = rel.abundance.sub$Year,
                    y = rel.abundance.sub$IA.obs,
                    col = 1, pch = 16, cex = 2)
            arrows( x0 = rel.abundance.sub$Year,
                    y0 = rel.abundance.sub$Lower95,
                    x1 = rel.abundance.sub$Year,
                    y1 = rel.abundance.sub$Upper95,
                    length=0.05, angle=90, code=3, lwd = 3, col = 1)

            # Posterior predictive
            if(posterior_pred){
                # Mean
                points( x = rel.abundance.sub$Year + 0.25,
                        y = IA_posterior_pred_sum[[i]][1,],
                        col = "Grey60", pch = 16, cex = 2)
                arrows( x0 = rel.abundance.sub$Year + 0.25,
                        y0 = as.numeric(IA_posterior_pred_sum[[i]][3,]),
                        x1 = rel.abundance.sub$Year + 0.25,
                        y1 = as.numeric(IA_posterior_pred_sum[[i]][4,]),
                        length=0.05, angle=90, code=3, lwd = 3, col = "Grey60")

            }

            if(!is.null(ioa_names)){
                legend("topleft", legend = ioa_names[i] ,bty = "n")
            }

        }
        if(j == 2){ dev.off()}
    }
}



#' OUTPUT FUNCTION
#'
#' Function that provides a plot of the estimated posterior densities of parameters from  SIR  model.
#'
#' @param SIR A fit SIR model or list of SIR models. Plots in the order provided.
#' @param file_name name of a file to identified the files exported by the
#'   function. If NULL, does not save.
#' @param multiple_sirs Logical whether or not multiple SIRS are provided as a list.
#' @param lower Vector of lower bounds for x-axis
#' @param upper Vector of upper bounds for x-axis
#' @param priors Default = NULL, wether realized priors are included in SIR
#' @param reference Default = NULL, wether reference case is included in SIR
#'
#' @return Returns and saves a figure with the posterior densities of parameters.
plot_density <- function(SIR, file_name = NULL, multiple_sirs = FALSE, lower = NULL, upper = NULL, prior_list = NULL, inc_reference = FALSE){

    if(multiple_sirs == FALSE){
        sir_list <- list(SIR)
    }
    if(multiple_sirs == TRUE){
        sir_list <- SIR
    }


    posteriors_lwd <- rep(3, length(sir_list))
    posteriors_lty <- c(1, 1:(length(sir_list)-1))
    posteriors_col <- c("grey", rep(1, length(sir_list)-1))

    if(!is.null(prior_list)){
        posteriors_lwd <- c(posteriors_lwd, rep(1, length(prior_list)))
        posteriors_lty <- c(posteriors_lty, c(1, 1:(length(prior_list)-1)))
        posteriors_col <- c(posteriors_col, c("grey", rep(1, length(prior_list)-1)))
        sir_list <- c(sir_list, prior_list)
    }

    # Vars of interest
    years <- sort(unique(c( sapply(sir_list, function(x) x$inputs$target.Yr),
                            sapply(sir_list, function(x) x$inputs$output.Years))))
    vars <- c("r_max", "K", "Nmin", paste0("N", years), "Max_Dep", paste0("status", years))
    vars_latex <- c("$r_{max}$", "$K$", "$N_{min}$", paste0("$N_{", years, "}$"), "Max depletion", paste0("Depletion in ", years))

    # Plot
    for(j in 1:(1 + !is.null(file_name))){
        if(j == 2){
            filename <- paste0(file_name, "_posterior_density", ".png")
            png( file = filename , width=169 / 25.4, height = 100 / 25.4, family = "serif", units = "in", res = 300)
        }

        par(mfrow = c(2,length(vars)/2))
        par( mar=c(3, 3 , 0.5 , 0.3) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))

        # Loop through vars
        for(i in 1:length(vars)){

            # Extract posterio densities
            posterior_dens <- list()
            for(k in 1:length(sir_list)){
                posterior_dens[[k]] <- density(sir_list[[k]]$resamples_output[,vars[i]])
            }

            # Extract prior densities
            # prior_dens <- list()
            # for(k in 1:length(sir_list)){
            #     prior_dens[[k]] <- density(sir_list[[k]]$resamples_output[,vars[i]])
            # }
            # priors_lwd <- rep(3, length(prior_dens))
            # priors_lty <- c(1:length(prior_dens))

            # Get x range
            if(is.null(lower[i])){
                xlow <- quantile(sapply(posterior_dens, "[", "x")$x, probs= c(0.025))
            } else if(is.na(lower[i])){
                xlow <- quantile(sapply(posterior_dens, "[", "x")$x, probs= c(0.025))
            } else{
                xlow <- lower[i]
            }

            if(is.null(upper[i])){
                xup <- quantile(sapply(posterior_dens, "[", "x")$x, probs= c(0.975))
            }
            else if(is.na(upper[i])){
                xup <- quantile(sapply(posterior_dens, "[", "x")$x, probs= c(0.975))
            } else{
                xup <- upper[i]
            }


            # Plot them
            plot(NA,
                 xlim = c(xlow, xup),
                 ylim = range(sapply(posterior_dens, "[", "y")),
                 ylab = "Density", xlab = latex2exp::TeX(vars_latex[i]))
            mapply(lines, posterior_dens, lwd = posteriors_lwd, lty = posteriors_lty, col = posteriors_col)
        }

        if(j == 2){ dev.off()}
    }
}


#' Function to compare posters
#'
#' @param sir_list list of SIR objects
#' @param model_names names of sir objects
#' @param file_name name of a file to identified the files exported by the
#'   function. If NULL, does not save.
#' @param reference_sir Default = NULL, wether reference SIR fit
compare_posteriors <- function(reference_sir = NULL, sir_list, model_names = NULL, file_name = NULL){

    # Extract range of values
    if(!is.null(reference_sir)){
        sir_list <- c(list(reference_sir), sir_list)
    }
    table_results <- list()
    for(i in 1:length(sir_list)){
        table_results[[i]] <- zerbini_table(sir_list[[i]], file_name = NULL)
    }

    # Vars of interest
    years <- sort(unique(c( sapply(sir_list, function(x) x$inputs$target.Yr),
                            sapply(sir_list, function(x) x$inputs$output.Years))))
    vars <- c("r_max", "K", "Nmin", paste0("N", years), "Max_Dep", paste0("status", years))
    vars_latex <- c("$r_{max}$", "$K$", "$N_{min}$", paste0("$N_{", years, "}$"), "Max depletion", paste0("Depletion in ", years))

    for(k in 1:length(vars)){ # Loop through vars
        for(j in 1:(1 + !is.null(file_name))){
            if(j == 2){
                filename <- paste0(file_name, "_posterior_comparison_", vars[k], ".png")
                png( file = filename , width=7, height = 100 / 25.4, family = "serif", units = "in", res = 300)
            }

            par( mar=c(3, 3 , 0.5 , 0.3) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))

            values <- matrix(NA, nrow = nrow(sir_list[[1]]$resamples_output), ncol = length(sir_list))

            for( i in 1:length(sir_list)){
                values[,i] <- sir_list[[i]]$resamples_output[,vars[k]]
            }


            boxplot(values, ylab = TeX(vars_latex[k]), xlab = "Scenario", xaxt = "n", col = "grey", outline = FALSE, cex.axis = 0.75, boxlty = 1, lty = 1)

            # Add x-lab
            if(is.null(model_names)){
                axis(side = 1, at = 1:length(sir_list), labels = as.character(1:length(sir_list)))
            }
            if(!is.null(model_names)){
                axis(side = 1, at = 1:length(sir_list), labels = model_names, cex.axis = 0.75)
            }

            # Add lines for reference
            if(!is.null(reference_sir)){
                ref_box <- boxplot(values[,1], plot = FALSE)

                abline(h = ref_box$stats[1,1], lty = 2, col = "grey30", lwd = 2)
                abline(h = ref_box$stats[3,1], lty = 2, col = "grey30", lwd = 2)
                abline(h = ref_box$stats[5,1], lty = 2, col = "grey30", lwd = 2)
            }

            boxplot(values, ylab = TeX(vars_latex[k]), xlab = "Scenario", xaxt = "n", col = "grey", outline = FALSE, cex.axis = 0.75, add = TRUE)

            if(j == 2){ dev.off()}
        }
    }
}




