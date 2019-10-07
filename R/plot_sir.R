
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
#'
#' @export
plot_trajectory <- function(SIR, Reference = NULL, file_name = NULL, posterior_pred = TRUE, coolors = "#941B0C",     coolors2 = "#104F55") {


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
    output_summary[2:6, ] <- sapply(x, quantile, probs= c(0.5, 0.025, 0.975, 0.25, 0.75))
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
        reference_summary[2:6, ] <- sapply(ref, quantile, probs= c(0.5, 0.025, 0.975, 0.25, 0.75))
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
        posterior_pred_summary[2:6, ] <- sapply(N_posterior_pred, quantile, probs= c(0.5, 0.025, 0.975, 0.25, 0.75))
        posterior_pred_summary[7, ] <- sapply(N_posterior_pred, min)
        posterior_pred_summary[8, ] <- sapply(N_posterior_pred, max)
        posterior_pred_summary[9, ] <- sapply(N_posterior_pred, length)

        colnames(posterior_pred_summary) <- paste0("Posterior_predictive_N_", abs.abundance$Year)
        rownames(posterior_pred_summary) <- row_names
    }

    # Extract catch trajectories
    catch_summary <- matrix(nrow = length(row_names), ncol = dim(catch.data)[2])
    catch_summary[1, ] <- sapply(catch.data, mean)
    catch_summary[2:6, ] <- sapply(catch.data, quantile, probs= c(0.5, 0.025, 0.975, 0.25, 0.75))
    catch_summary[7, ] <- sapply(catch.data, min)
    catch_summary[8, ] <- sapply(catch.data, max)
    catch_summary[9, ] <- sapply(catch.data, length)


    # Get 95% CI and range
    Years <- as.numeric(gsub("N_", "", colnames(output_summary)))
    abs.abundance$Upper95 <- qlnorm(0.975, log(abs.abundance$N.obs), abs.abundance$Sigma)
    abs.abundance$Lower95 <- qlnorm(0.025, log(abs.abundance$N.obs), abs.abundance$Sigma)

    # Set up ranges for plots
    ymax <- max(c(max(output_summary[2:6, ]), abs.abundance$N.obs, abs.abundance$Lower95, abs.abundance$Upper95, N_priors))
    ymin <- -0.2 * ymax

    ymax_catch <- max(catch_summary[2:6, ])
    ymax_catch <- ymax_catch * 2
    ymin_catch <- 0

    # Set up axis
    x_axis <- Years/10
    x_axis <- unique(round(x_axis) * 10)
    x_axis <- x_axis[rep_len(c(TRUE, FALSE), length.out = length(x_axis))]

    y_axis <- seq(0, ymax , length.out = 5)
    y_axis_catch <-  floor(seq(0,  max(catch_summary[2:6, ]) , length.out = 5)/1000) * 1000

    # Plot trajectory
    for(i in 1: (1 + as.numeric(!is.null(file_name)) * 2)){

        # PNG
        if(i == 2){
            filename <- paste0(file_name, "_trajectory_summary", ".png")
            png( file = filename , width=7.5, height = 100 / 25.4, family = "serif", units = "in", res = 300)
        }

        # PDF
        if(i == 3){
            filename <- paste0(file_name, "_trajectory_summary", ".pdf")
            pdf( file = filename , width=7.5, height = 100 / 25.4, family = "serif")
        }


        # Plot configuration
        par( mar=c(3, 3 , 0.5 , 3) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))
        # catch history
        plot(y = NA, x = NA,
             ylim = c(ymin_catch, ymax_catch),
             xlim = c(min(Years), max(Years)),
             xlab = NA, ylab = NA, xaxt = "n",  yaxt = "n")

        # Credible interval
        polygon(
            x = c(Years,rev(Years)),
            y = c(catch_summary[3, ],rev(catch_summary[4, ])),
            col = adjustcolor(coolors, alpha = 0.2), border = NA) # 95% CI
        polygon( x = c(Years,rev(Years)),
                 y = c(catch_summary[5, ], rev(catch_summary[6, ])),
                 col = adjustcolor(coolors, alpha = 0.5), border = NA) # 90% CI

        # Median
        lines( x = Years, y = catch_summary[2, ], lty = 1, lwd = 3, col = coolors) # Median

        axis(side = 4, at = y_axis_catch, col = coolors, col.axis=coolors, font =2)
        mtext(side = 4, "Annual catch (numbers)", line = 1.6, col = coolors, font = 2, adj = 0.4)


        par(new = TRUE)
        # N-trajectory
        plot(y = NA, x = NA,
             ylim = c(ymin, ymax),
             xlim = c(min(Years), max(Years)),
             xlab = NA, ylab = NA, xaxt = "n", col.axis=coolors2, col = coolors2, col.lab= coolors2, font = 2)
        abline(h = 0, col = coolors2)
        mtext(side = 2, "Number of individuals", line = 1.6, font = 2, col=coolors2, adj = 0.6)
        axis(side = 1, x_axis, font = 2)
        mtext(side = 1, "Year", line = 1.6, font = 2)

        # N Trajectory
        # Credible interval
        polygon(
            x = c(Years,rev(Years)),
            y = c(output_summary[3, ],rev(output_summary[4, ])),
            col = adjustcolor(coolors2, alpha = 0.2), border = NA) # 95% CI
        polygon( x = c(Years,rev(Years)),
                 y = c(output_summary[5, ], rev(output_summary[6, ])),
                 col = adjustcolor(coolors2, alpha = 0.5), border = NA) # 90% CI

        # Median
        lines( x = Years, y = output_summary[2, ], lwd = 3, col = coolors2) # Median

        # Reference median
        if(!is.null(Reference)){
            lines( x = ref_years, y = reference_summary[2, ], lwd = 3, lty = 3) # Median
        }

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


        if(i > 1){ dev.off()}
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
#'
#' @export
plot_ioa <- function(SIR, file_name = NULL, ioa_names = NULL, posterior_pred = TRUE, coolors = "#104F55"){

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
        IA_summary[[i]][2:6, ] <- sapply(IA_pread[[i]], quantile, probs= c(0.5, 0.025, 0.975, 0.25, 0.75))
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
            IA_posterior_pred_sum[[i]][2:6, ] <- sapply(IA_posterior_pred[[i]], quantile, probs= c(0.5, 0.025, 0.975, 0.25, 0.75))
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
                col = adjustcolor(coolors, alpha = 0.2), border = NA) # 95% CI

            polygon( x = c(IA.yr.range[[i]], rev(IA.yr.range[[i]])),
                     y = c(IA_summary[[i]][5, ], rev(IA_summary[[i]][6, ])),
                     col = adjustcolor(coolors, alpha = 0.5), border = NA) # 50% CI


            # Median
            lines( x = IA.yr.range[[i]], y = IA_summary[[i]][2, ], col = coolors, lwd = 3) # Median


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
                        y = IA_posterior_pred_sum[[i]][2,],
                        col = "Grey70", pch = 16, cex = 2)
                arrows( x0 = rel.abundance.sub$Year + 0.25,
                        y0 = as.numeric(IA_posterior_pred_sum[[i]][3,]),
                        x1 = rel.abundance.sub$Year + 0.25,
                        y1 = as.numeric(IA_posterior_pred_sum[[i]][4,]),
                        length=0.05, angle=90, code=3, lwd = 3, col = "Grey70")

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
#' @param lower Vector of lower bounds for x-axis
#' @param upper Vector of upper bounds for x-axis
#' @param priors Default = NULL, SIR realized priors or list of SIR realized priors. Plots in the order provided.
#' @param reference Default = NULL, wether reference case is included in SIR
#'
#' @return Returns and saves a figure with the posterior densities of parameters.
#' @export
plot_density <- function(SIR, file_name = NULL, lower = NULL, upper = NULL, priors = NULL, inc_reference = TRUE){

    # Make into list
    if(class(SIR) == "SIR"){
        SIR <- list(SIR)
    }

    # Make into list
    if(class(priors) == "SIR"){
        priors <- list(priors)
    }

    if(inc_reference){
        posteriors_lwd <- rep(3, length(SIR))
        posteriors_lty <- c(1, 1:(length(SIR)-1))
        posteriors_col <- c("grey", rep(1, length(SIR)-1))
    }

    if(inc_reference == FALSE){
        posteriors_lwd <- rep(3, length(SIR))
        posteriors_lty <- c(1:(length(SIR)))
        posteriors_col <- c(rep(1, length(SIR)))
    }

    if(!is.null(priors)){
        if(inc_reference){
            posteriors_lwd <- c(posteriors_lwd, rep(1, length(priors)))
            posteriors_lty <- c(posteriors_lty, c( 1: ifelse(length(priors) > 1, c(1,(length(priors)-1)), 1) ))
            posteriors_col <- c(posteriors_col, c("grey", rep(1, ifelse(length(priors) > 1, (length(priors)-1), 0))))
            SIR <- c(SIR, priors)
        }
        if(inc_reference == FALSE){
            posteriors_lwd <- c(posteriors_lwd, rep(1, length(priors)))
            posteriors_lty <- c(posteriors_lty, c(1:(length(priors))))
            posteriors_col <- c(posteriors_col, c(rep(1, length(priors))))
            SIR <- c(SIR, priors)
        }
    }

    # Vars of interest
    years <- sort(unique(c( sapply(SIR, function(x) x$inputs$target.Yr),
                            sapply(SIR, function(x) x$inputs$output.Years))))
    vars <- c("r_max", "K", "Nmin", paste0("N", years), "Max_Dep", paste0("status", years))
    vars_latex <- c("$r_{max}$", "$K$", "$N_{min}$", paste0("$N_{", years, "}$"), "Max depletion", paste0("Depletion in ", years))

    # Plot
    for(j in 1:(1 + as.numeric(!is.null(file_name)) * 2)){

        # PNG
        if(j == 2){
            filename <- paste0(file_name, "_posterior_density", ".png")
            png( file = filename , width=10, height = 110 / 25.4, family = "serif", units = "in", res = 300)
        }

        # PDF
        if(j == 3){
            filename <- paste0(file_name, "_posterior_density", ".pdf")
            pdf( file = filename , width=10, height = 110 / 25.4, family = "serif")
        }

        par(mfrow = c(2,(length(vars))/2 + 2))
        par( mar=c(3, 0.05 , 0.5 , 0.55) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))

        plot.new()


        # Loop through vars
        for(i in 1:length(vars)){

            # Extract posterio densities
            posterior_dens <- list()
            for(k in 1:length(SIR)){
                posterior_dens[[k]] <- density(SIR[[k]]$resamples_output[,vars[i]])
            }

            # Extract prior densities
            # prior_dens <- list()
            # for(k in 1:length(SIR)){
            #     prior_dens[[k]] <- density(SIR[[k]]$resamples_output[,vars[i]])
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
                 ylim = c(0, range(sapply(posterior_dens, "[", "y"))[2]),
                 ylab = NA, xlab = latex2exp::TeX(vars_latex[i]), yaxt = "n")
            mapply(lines, posterior_dens, lwd = posteriors_lwd, lty = posteriors_lty, col = posteriors_col[1:length(posterior_dens)])

            if(i == (length(vars)/2))  {
                plot.new()
                plot.new()
            }


            if(i %in% c(1, (length(vars))/2 + 1) ) {
                mtext(side = 2, "Density", line = 1, cex= 0.75)
            }


        }

        if(j > 1){ dev.off()}
    }
}



#' Function to compare posters
#'
#' @param SIR list of SIR objects
#' @param model_names names of sir objects
#' @param file_name name of a file to identified the files exported by the
#'   function. If NULL, does not save.
#' @param reference_sir Default = TRUE, reference SIR is in \code{SIR}, should be first if so.
#' @param model_average Default = TRUE, model average SIR is in \code{SIR}, should be last if so.
#' @param years Optional vector of years to compare.
#' @param bayes_factor Optional. Vector of bayesfactors of length SIR
#'
#' @export
compare_posteriors <- function(SIR, model_names = NULL, file_name = NULL, bayes_factor = NULL, reference_sir = TRUE, model_average = TRUE, years = NULL){

    # If it is a single SIR, make into a list
    if(class(SIR) == "SIR"){
        SIR = list(SIR)
    }
    cols <- rep("grey", length(SIR))

    # Extract range of values
    if(reference_sir){
        cols <- c( "#7C9FA2", cols[-1])
    }

    if(model_average){
        cols <- c(cols[-length(cols)] , "#BA6D64")
    }


    library(latex2exp)

    #################################
    # SEPERATE PLOTS
    #################################

    # Vars of interest
    if(is.null(years)){
        years <- sort(unique(c( sapply(SIR, function(x) x$inputs$target.Yr),
                                sapply(SIR, function(x) x$inputs$output.Years))))
    }
    vars <- c("r_max", "K", "Nmin", paste0("N", years), "Max_Dep", paste0("status", years))
    vars_latex <- c("$r_{max}$", "$K$", "$N_{min}$", paste0("$N_{", years, "}$"), "Max depletion", paste0("Depletion in ", years))

    for(k in 1:length(vars)){ # Loop through vars
        for(j in 1:(1 + as.numeric(!is.null(file_name)) * 2)){

            # PNG
            if(j == 2){
                filename <- paste0(file_name, "_posterior_comparison_", vars[k], ".png")
                png( file = filename , width=8, height = 100 / 25.4, family = "serif", units = "in", res = 300)
            }

            # PDF
            if(j == 3){
                filename <- paste0(file_name, "_posterior_comparison_", vars[k], ".pdf")
                pdf( file = filename , width=8, height = 100 / 25.4, family = "serif")
            }

            par( mar=c(5, 3 , 0.5 , 0.3) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))

            values <- matrix(NA, nrow = nrow(SIR[[1]]$resamples_output), ncol = length(SIR))

            for( i in 1:length(SIR)){
                values[,i] <- SIR[[i]]$resamples_output[,vars[k]]
            }


            boxplot(values, ylab = latex2exp::TeX(vars_latex[k]), xlab = NA, xaxt = "n", col = cols, outline = FALSE, cex.axis = 0.75, boxlty = 1, lty = 1)

            # Add x-lab
            if(is.null(model_names)){
                axis(side = 1, at = 1:length(SIR), labels = as.character(1:length(SIR)))
            }
            if(!is.null(model_names)){
                axis(side = 1, at = 1:length(SIR), labels = model_names, cex.axis = 0.75)
            }

            # Add bayes factor
            if(is.null(bayes_factor)){
                mtext(side = 1, text = "Scenario", line = 1.6)
            }
            if(!is.null(bayes_factor)){
                axis(side = 1, at = c(0, 1:length(SIR)), labels = c("BF =", bayes_factor), cex.axis = 0.75, tick = FALSE, line = 1.6)

                mtext(side = 1, text = "Scenario", line = 3.4)
            }

            # Add lines for reference
            if(reference_sir){
                ref_box <- boxplot(values[,1], plot = FALSE)

                abline(h = ref_box$stats[1,1], lty = 2, col = "grey30", lwd = 2)
                abline(h = ref_box$stats[3,1], lty = 2, col = "grey30", lwd = 2)
                abline(h = ref_box$stats[5,1], lty = 2, col = "grey30", lwd = 2)
            }

            boxplot(values, ylab = latex2exp::TeX(vars_latex[k]), xlab = NA, xaxt = "n", col = cols, outline = FALSE, cex.axis = 0.75, add = TRUE)
            # Add means
            mean_vec <- colMeans(values)
            for( i in 1:length(SIR)){
                segments(x0 = .59999 + (i - 1), x1 = 1.39999 + (i - 1), y0 = mean_vec[i], col = 1, lwd = 2, lty = 3)
            }

            if(j > 1){ dev.off()}
        }
    }

    ###########################
    # COMBINED PLOTS
    ###########################

    # Vars of interest
    vars <- c("r_max", "Nmin", paste0("N", years), "K", "Max_Dep", paste0("status", years))
    vars_latex <- c("$r_{max}$", "$N_{min}$", paste0("$N_{", years, "}$"), "$K$", "Max depletion", paste0("Depletion in ", years))

    for(j in 1:(1 + as.numeric(!is.null(file_name)) * 2)){

        # PNG
        if(j == 2){
            filename <- paste0(file_name, "_posterior_comparison", ".png")
            png( file = filename , width=10, height = 8, family = "serif", units = "in", res = 300)
        }

        # PDF
        if(j == 3){
            filename <- paste0(file_name, "_posterior_comparison", ".pdf")
            pdf( file = filename , width=10, height = 8, family = "serif")
        }

        # Set up multiplot
        layout(matrix(c(1:(length(vars) + 2)), (length(vars)/2 + 1), 2, byrow = FALSE), heights = c(rep(1, length(vars)/2), 0.35))
        par( mar=c(0.15, 3.5 , 0.5 , 0.5) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))

        # Loop through vars
        for(k in 1:length(vars)){ # Loop through vars

            values <- matrix(NA, nrow = nrow(SIR[[1]]$resamples_output), ncol = length(SIR))

            for( i in 1:length(SIR)){
                values[,i] <- SIR[[i]]$resamples_output[,vars[k]]
            }

            boxplot(values, ylab = NA, xlab = NA, xaxt = "n", col = cols, outline = FALSE, cex.axis = 1, boxlty = 1, lty = 1)
            mtext(side = 2, text = latex2exp::TeX(vars_latex[k]), line = 1.8)


            # X-Lab
            if(k %in% c(length(vars)/2, length(vars))){
                # Add model name
                if(is.null(model_names)){
                    axis(side = 1, at = 1:length(SIR), labels = as.character(1:length(SIR)))
                }
                if(!is.null(model_names)){
                    axis(side = 1, at = 1:length(SIR), labels = model_names, cex.axis = 0.9, gap.axis = 0)
                }

                # Add bayes factor
                if(is.null(bayes_factor)){
                    mtext(side = 1, text = "Scenario", line = 1.6)
                }
                if(!is.null(bayes_factor)){
                    axis(side = 1, at = c(0, 1:length(SIR)), labels = c("BF =", bayes_factor), cex.axis = 0.9, tick = FALSE, line = 1.6, gap.axis = 0)

                    mtext(side = 1, text = "Scenario", line = 3.4)
                }
            }

            # Add lines for reference
            if(reference_sir){
                ref_box <- boxplot(values[,1], plot = FALSE)

                abline(h = ref_box$stats[1,1], lty = 2, col = "grey30", lwd = 2)
                abline(h = ref_box$stats[3,1], lty = 2, col = "grey30", lwd = 2)
                abline(h = ref_box$stats[5,1], lty = 2, col = "grey30", lwd = 2)
            }

            boxplot(values, ylab = NA, xlab = NA, xaxt = "n", col = cols, outline = FALSE, cex.axis = 1, add = TRUE)
            mtext(side = 2, text = latex2exp::TeX(vars_latex[k]), line = 1.8)

            # Add means
            mean_vec <- colMeans(values)
            for( i in 1:length(SIR)){
                segments(x0 = .59999 + (i - 1), x1 = 1.39999 + (i - 1), y0 = mean_vec[i], col = 1, lwd = 2, lty = 3)
            }

            # Bottom row
            if(k == length(vars)/2){
                plot.new()
            }
        }

        if(j > 1){ dev.off()}
    }
}




