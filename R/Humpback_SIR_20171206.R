# R Script for Bayesian Assessment Model of Humpback Whales - uses the equations
# in Zerbini et al. (2011) Code cleaned up and sent to Andre Punt, John Best and
# Grant Adams on 7 Dec 2017

#' HUMPBACK SIR controls the sampling from the priors, the bisection and
#' likelihoods and the output functions
#'
#' @param file.name name of a file to identified the files exported by the
#'   function
#' @param n.resamples number of resamples to compute the marginal posterior
#'   distributions
#' @param priors List of priors, usually generated using \link{make_prior_list}.
#' @param Klim bounds for K when preforming the bisection method of Punt and
#'   Butterworth (1995). Defined by two elements, the lower and upper bounds.
#'   Default is (1, 500000)
#' @param target.Yr year of the target population estimate for the bisection
#'   method. Default is 2008
#' @param num.haplotypes number of haplotypes to compute minimum viable
#'   population (from Jackson et al., 2006 and IWC, 2007)
#' @param tolerance.for.bisection tolerance value for performing the bisection
#'   method
#' @param output.Yrs year for outputing the predicted abundance estimates.
#'   Default is 2008, but multiple years can be specified. For example, if
#'   outputs for 2005 and 2008 are needed, output.Yrs = c(2005, 2008)
#' @param abs.abundance R object containing year, estimate of absolute
#'   abundance, and CV (see example)
#' @param rel.abundance R object containing years, estimates of relative
#'   abudnance and CVs (see example)
#' @param rel.abundance.key key to speficy if relative abundance data are used
#'   in the likelihood. Default is TRUE
#' @param count.data R object containing years, estimates of counts and effort.
#'   NOT USED
#' @param count.data.key key to speficy in count data are used. Default is
#'   FALSE. NOT USED
#' @param growth.rate.obs observed growth rate (1st element) and standard error
#'   (2nd element) as in Zerbni et al. (2011). If third element is FALSE, the
#'   growth rate is not included in the likelihood
#' @param growth.rate.Yrs Years for which the growth.rate.obs were computed (as
#'   in Zerbini et al., 2011)
#' @param catch.data R object containing the years and catches (see example)
#' @param Threshold threshold for the McCallister et al. (1994) SIR. This is
#'   data-specific. Default is 1e-100.
#' @param Print key to print various pieces of information as the code runs.
#'   Default is 0 (nothing is printed)
#'
#' @return A \code{list} containing posterior samples and metadata
#'
#' TODO: Add the negative binomial likelihood for the count data, which is not
#' currently used even though it is defined in the main function call.
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' HUMPBACK.SIR(file.name = "test.N2005",
#'              n.resamples = 100,
#'              priors = make_prior_list(),
#'              Klim = c(1, 500000),
#'              target.Yr = 2005,
#'              num.haplotypes = 0,
#'              tolerance.for.bisection = 0.0001,
#'              output.Yrs = c(2005, 2006),
#'              abs.abundance = Abs.Abundance.2005,
#'              rel.abundance = Rel.Abundance,
#'              rel.abundance.key = TRUE,
#'              count.data = Count.Data,
#'              count.data.key = FALSE,
#'              growth.rate.obs = c(0.074, 0.033, TRUE),
#'              growth.rate.Yrs = c(1995, 1996, 1997, 1998),
#'              catch.data = Catch.data,
#'              Threshold = 1e-17,
#'              Print = 0)
HUMPBACK.SIR <- function(file.name = "NULL",
                         n.resamples = 1000,
                         priors = make_prior_list(),
                         Klim = c(1, 500000),
                         target.Yr = 2008,
                         num.haplotypes = 66,
                         tolerance.for.bisection = 0.001,
                         output.Yrs = c(2008),
                         abs.abundance = Abs.Abundance,
                         rel.abundance = Rel.Abundance,
                         rel.abundance.key = TRUE,
                         count.data = NULL,
                         count.data.key = FALSE,
                         growth.rate.obs = c(0.074, 0.033, TRUE),
                         growth.rate.Yrs = c(1995, 1996, 1997, 1998),
                         catch.data = Catch.data,
                         Threshold = 1e100,
                         Print = 0) {
  begin.time <- Sys.time()

  ################################
  # Assigning variables
  ################################
  target.Yr <- target.Yr
  ## Use the first year of the projection is set as the first year in the
  ## catch series
  start.Yr <- catch.data$Year[1]
  ## The last year of the projection is set as the last year in the catch or
  ## abundance series, whichever is most recent
  end.Yr <- max(tail(catch.data$Year, 1),
                max(abs.abundance$Year),
                max(rel.abundance$Year))
  ## Setting the target year for the bisection method
  bisection.Yrs <- target.Yr-start.Yr + 1
  ## Setting the years to project
  projection.Yrs <- end.Yr-start.Yr + 1

  ## Assigning the catch data
  catches <- catch.data$Catch
  ## Determining the number of Indices of Abundance available
  num.IA <- max(rel.abundance$Index)
  ## Determining the number of Count Data sets available
  num.Count <- max(count.data$Index)
  ## Computing the value of sigma as in Zerbini et al. 2011
  rel.abundance$Sigma <- sqrt(log(1 + rel.abundance$CV.IA.obs^2))
  ## Computing the value of sigma for the count data as in Zerbini et al. (2011)
  count.data$Sigma <- sqrt(log(1 + count.data$CV.IA.obs^2))
  ## Computing the value of sigma as in Zerbini et al. 2011
  abs.abundance$Sigma <- sqrt(log(1 + abs.abundance$CV.obs^2))
  ## Computing the minimum viable population, if num.haplotypes=0, assumes no MVP
  MVP <- 4 * num.haplotypes

  ## Assigning the value for tolerance in computing K using the uniroot function
  ## (bisection method)
  tol <- tolerance.for.bisection
  ## Start the loop
  i <- 0
  ## Keep track of number of draws
  draw <- 1
  Cumulative.Likelihood <- 0

  #Creating output vectors
  #-------------------------------------
  names <- c("r_max", "K", "sample.N.obs", "add_CV", "Nmin", "YearMin",
             "violate_MVP", paste("N", output.Yrs, sep = ""),
             paste("ROI_IA", unique(rel.abundance$Index), sep = ""),
             paste("q_IA", unique(rel.abundance$Index), sep = ""),
             paste("ROI_Count", unique(count.data$Index), sep = ""),
             paste("q_Count", unique(count.data$Index), sep = ""),
             "NLL.IAs", "NLL.Count", "NLL.N", "NLL.GR", "NLL", "Likelihood",
             "Max_Dep", paste("status", output.Yrs, sep = ""), "draw", "save")


  samples.output <- matrix(0, nrow = 1, ncol = length(names))
  resamples.output <- matrix(0, nrow = 1, ncol = length(names))
  resamples.trajectories <- matrix(NA, nrow = 1, ncol = projection.Yrs)
  final.trajectory <- matrix(NA, nrow = projection.Yrs, ncol = 6)
  Year <- seq(start.Yr, end.Yr, by = 1)

  #Initiating the SIR loop

  while (i < n.resamples) {
    #Sampling from Priors
    #-------------------------------
    save <- FALSE #variable to indicate whether a specific draw is kept

    #Sampling for r_max
    sample.r_max <- 0.2 #setting sample.r_max outside of the bound
    ## FIXME Why is this check necessary; just set the bounds using the prior?
    while (sample.r_max < priors$r_max$pars[1] | sample.r_max > priors$r_max$pars[2]) {
      ## Prior on r_max, keep if within boundaries
      sample.r_max <- priors$r_max$rfn()
    }

    ## Sampling from the N.obs prior
    sample.N.obs <- priors$N_obs$rfn()

    ## Prior on additional CV
    if (priors$add_CV$use) {
      sample.add_CV <- priors$add_CV$rfn()
    } else {
      sample.add_CV <- 0
    }

    ## Sample from prior for `z` (usually constant)
    sample.z <- priors$z$rfn()

    ## Sampling from q priors if q.prior is TRUE; priors on q for indices of
    ## abundance
    if (priors$q_IA$use) {
      q.sample.IA <- replicate(num.IA, priors$q_IA$rfn())
    } else {
      ## FIXME: -9999 is probably not a good sentinel value here; NA?
      q.sample.IA <- rep(-9999, length(unique(rel.abundance$Index)))
    }

    ##priors on q for count data
    if (priors$q_count$use) {
      q.sample.Count <- replicate(num.Count, priors$q_count$rfn())
    } else {
      ## FIXME: Sentinel -9999 again
      q.sample.Count <- rep(-9999, length(unique(count.data$Index)))
    }

    sample.K <- LOGISTIC.BISECTION.K(K.low = Klim[1],
                                     K.high = Klim[2],
                                     r_max = sample.r_max,
                                     z = sample.z,
                                     num.Yrs = bisection.Yrs,
                                     start.Yr = start.Yr,
                                     target.Pop = sample.N.obs,
                                     catches = catches,
                                     MVP = MVP,
                                     tol = tol)

    #Computing the predicted abundances with the samples from the priors
    #----------------------------------------
    Pred.N <- GENERALIZED.LOGISTIC(r_max = sample.r_max,
                                   K = sample.K,
                                   N1 = sample.K,
                                   z = sample.z,
                                   start.Yr = start.Yr,
                                   num.Yrs = projection.Yrs,
                                   catches = catches,
                                   MVP = MVP)
    #Print results if required
    ## FIXME: Print as boolean? Maybe rename `verbose`
    if (Print==1) {
      print(paste("sample.r_max =", sample.r_max,
                  "sample.N.obs =", sample.N.obs,
                  "sample.K =", sample.K,
                  "Pred.N.target =", Pred.N$Pred.N[bisection.Yrs]))}

    #Computing the predicted ROI for the IAs and Count data, if applicable
    #----------------------------------------
    #For IAs
    if (rel.abundance.key) {
      Pred.ROI.IA <- COMPUTING.ROI(data = rel.abundance,
                                   Pred.N = Pred.N$Pred.N,
                                   start.Yr = start.Yr)
    } else {
      Pred.ROI.IA <- rep(0, num.IA)
    }

    #For Count Data
    if (count.data.key) {
      Pred.ROI.Count <- COMPUTING.ROI(data = count.data,
                                      Pred.N = Pred.N$Pred.N,
                                      start.Yr = start.Yr)
    } else {
      Pred.ROI.Count <- rep(0, num.Count)
    }

    #Calculate Analytical Qs if rel.abundance.key is TRUE
    #---------------------------------------------------------
    if (rel.abundance.key) {
      if (!priors$q_IA$use) {
        q.sample.IA <- CALC.ANALYTIC.Q(rel.abundance,
                                       Pred.N$Pred.N,
                                       start.Yr,
                                       sample.add_CV,
                                       num.IA)
      } else {
      q.sample.IA <- q.sample.IA
      }
    }

    #browser()

    ## Calculate Analytical Qs if count.data.key is TRUE
    ## (NOT USED YET - AZerbini, Feb 2013)
    if (rel.abundance.key) {
      if (!priors$q_count$use) {
        q.sample.Count <- CALC.ANALYTIC.Q(count.data,
                                          Pred.N$Pred.N,
                                          start.Yr,
                                          sample.add_CV,
                                          num.Count)
      } else {
      q.sample.Count <- q.sample.Count
      }
    }

    ## FIXME Print -> verbose?
    if (Print == 1) {
      print(paste("q.IAs =", q.sample.IA))
      print(paste("q.Count =", q.sample.Count))
    }

    #Compute the likelihoods
    #--------------------------------
    # (1) relative indices (if rel.abundance.key is TRUE)
    if (rel.abundance.key) {
      lnlike.IAs <- LNLIKE.IAs(rel.abundance,
                               Pred.N$Pred.N,
                               start.Yr,
                               q.sample.IA,
                               sample.add_CV,
                               num.IA,
                               log=TRUE)
    } else {
      lnlike.IAs <- 0
    }
    ## FIXME Print -> verbose
    if (Print==1) {
      print(paste("lnlike.IAs =", lnlike.IAs))
    }

    # (2) count data (if count.data.key is TRUE)
    if (count.data.key) {
      lnlike.Count <- LNLIKE.IAs(count.data,
                                 Pred.N$Pred.N,
                                 start.Yr,
                                 q.sample.Count,
                                 sample.add_CV,
                                 num.Count,
                                 log=TRUE)
    } else {
      lnlike.Count <- 0
    }
    ## FIXME: Print -> verbose (use message?)
    if (Print==1) {
      print(paste("lnlike.Count =", lnlike.Count))
    }

    # (3) absolute abundance
    lnlike.Ns <- LNLIKE.Ns(abs.abundance,
                           Pred.N$Pred.N,
                           start.Yr,
                           sample.add_CV,
                           log=TRUE)
    ## FIXME: Print -> message
    if (Print==1) {
      print(paste("lnlike.Ns =", lnlike.Ns))
    }

    # (4) growth rate if applicable
    if (growth.rate.obs[3]) {
      Pred.GR <- PRED.GROWTH.RATE(growth.rate.Yrs=growth.rate.Yrs,
                                  Pred.N=Pred.N$Pred.N,
                                  start.Yr=start.Yr)
      lnlike.GR <- LNLIKE.GR(Obs.GR=growth.rate.obs[1],
                             Pred.GR=Pred.GR,
                             GR.SD.Obs=growth.rate.obs[2])
    } else {
      lnlike.GR <- 0
    }
    ## FIXME: Print -> message
    if (Print==1) {
      print(paste("lnlike.GR =", lnlike.GR))
    }

    ## These use the likelihoods in Zerbini et al. (2011)
    LL <- lnlike.IAs[[1]] + lnlike.Count[[1]] + lnlike.Ns[[1]] + lnlike.GR[[1]]
    Likelihood <- exp(-LL)
    ## FIXME: Print -> message
    if (Print==1) {
      print(paste("NLL =", LL, "Likelihood =", Likelihood))
    }

    if (Pred.N$Violate.MVP) {
      Likelihood <- 0
      print(paste("MVP violated on draw", draw))
    }

    Cumulative.Likelihood <- Cumulative.Likelihood + Likelihood
    #print(paste("Sample=", i, "Likelihood=", Likelihood))

    if (!Pred.N$Violate.MVP) {
      while (Cumulative.Likelihood > Threshold) {
        print(paste("sample =", i, " & draw =", draw))
        ## FIXME: Print -> message
        if (Print==1) {
          print(paste("draw =", draw,
                      "Likelihood =", Likelihood,
                      "Cumulative =", Cumulative.Likelihood))}
        save <- TRUE
        Cumulative.Likelihood <- Cumulative.Likelihood-Threshold
        resamples.trajectories <- rbind(resamples.trajectories, Pred.N$Pred.N)
        resamples.output <- rbind(resamples.output,
                                  c(sample.r_max,
                                    sample.K,
                                    sample.N.obs,
                                    sample.add_CV,
                                    Pred.N$Min.Pop,
                                    Pred.N$Min.Yr,
                                    Pred.N$Violate.MVP,
                                    c(Pred.N$Pred.N[output.Yrs - start.Yr + 1]),
                                    Pred.ROI.IA,
                                    q.sample.IA,
                                    Pred.ROI.Count,
                                    q.sample.Count,
                                    lnlike.IAs[[1]],
                                    lnlike.Count[[1]],
                                    lnlike.Ns[[1]],
                                    lnlike.GR[[1]],
                                    LL,
                                    Likelihood,
                                    Pred.N$Min.Pop / sample.K,
                                    c(Pred.N$Pred.N[output.Yrs - start.Yr + 1] /
                                      sample.K),
                                    draw,
                                    save))
        i <- i+1
      }
    }

    samples.output <- rbind(samples.output,
                            c(sample.r_max,
                              sample.K,
                              sample.N.obs,
                              sample.add_CV,
                              Pred.N$Min.Pop,
                              Pred.N$Min.Yr,
                              Pred.N$Violate.MVP,
                              c(Pred.N$Pred.N[output.Yrs-start.Yr+1]),
                              Pred.ROI.IA,
                              q.sample.IA,
                              Pred.ROI.Count,
                              q.sample.Count,
                              lnlike.IAs[[1]],
                              lnlike.Count[[1]],
                              lnlike.Ns[[1]],
                              lnlike.GR[[1]],
                              LL,
                              Likelihood,
                              Pred.N$Min.Pop/sample.K,
                              c(Pred.N$Pred.N[output.Yrs-start.Yr+1]/sample.K),
                              draw,
                              save))

    draw <- draw+1
  }

  samples.output <- data.frame(samples.output)
  names(samples.output) <- names
  samples.output <- samples.output[-1, ]
  samples.output.summary <- SUMMARY.SIR(x=samples.output, scenario = file.name)
  ## FIXME Use `paste0` here
  write.csv(samples.output,
            paste(file.name, "_", "samples.output.csv", sep=""))

  resamples.output <- data.frame(resamples.output)
  names(resamples.output) <- names
  resamples.output <- resamples.output[-1, ]
  resamples.output.summary <- SUMMARY.SIR(x = resamples.output,
                                          scenario = file.name)
  ## FIXME Use `paste0` here
  write.csv(resamples.output,
            paste(file.name, "_", "resamples.output.csv", sep=""))

  resamples.trajectories <- data.frame(resamples.trajectories)
  names(resamples.trajectories) <- seq(start.Yr, end.Yr,  1)
  resamples.trajectories <- resamples.trajectories[-1, ]
  ## FIXME Use `paste0` here
  write.csv(resamples.trajectories,
            paste(file.name, "_", "resample.trajectories.csv",  sep=""))

  ## FIXME Use one quantile call to get all?
  final.trajectory[, 1] <- sapply(resamples.trajectories, mean)
  final.trajectory[, 2] <- sapply(resamples.trajectories, median)
  final.trajectory[, 3] <- sapply(resamples.trajectories, quantile,
                                  probs = c(0.025))
  final.trajectory[, 4] <- sapply(resamples.trajectories, quantile,
                                  probs = c(0.975))
  final.trajectory[, 5] <- sapply(resamples.trajectories, quantile,
                                  probs = c(0.05))
  final.trajectory[, 6] <- sapply(resamples.trajectories, quantile,
                                  probs = c(0.95))
  final.trajectory <- data.frame(final.trajectory)
  names(final.trajectory) <- c("mean", "median",
                               "PI.2.5%", "PI.97.5%",
                               "PI.5%", "PI.95%")
  final.trajectory <- data.frame(Year, final.trajectory)

  resamples.per.samples <- dim(samples.output)[1] / dim(resamples.output)[1]

  end.time <- Sys.time()
  print(paste("Time to Compute =", (end.time-begin.time)))

  list(call = call,
       file.name = file.name,
       Date.Time = Sys.time(),
       Time.to.compute.in.minutes = paste((end.time-begin.time) / 60),
       Threshold = Threshold,
       Ratio.Resamples.per.Sample = paste("1 resample",
                                          ":",
                                          resamples.per.samples,
                                          "samples"),
       resamples.output = resamples.output,
       resamples.output.summary = resamples.output.summary$output.table,
       samples.output.summary = samples.output.summary$output.table,
       final.trajectory = final.trajectory,
       inputs = list(draws = draw,
                     n.resamples = n.resamples,
                     prior_r_max = priors$r_max,
                     priors_N.obs = priors$N.obs,
                     target.Yr = target.Yr,
                     MVP = paste("num.haplotypes = ",
                                 num.haplotypes,
                                 "MVP = ",
                                 4 * num.haplotypes),
                     tolerance = tolerance.for.bisection,
                     output.Years = output.Yrs))
}

#' GENERALISED LOGISTIC MODEL
#'
#' \code{GENERALIZED.LOGISTIC} returns the population projection using a
#' Pella-Tomlison population dynamics model: $$N_{t+1} =
#' N_{t}+N_{t}*r_{max}*\left[ 1 - \left( \frac{N_{t}}{K} \right) ^z \right] -
#' C_{t}$$ where $N$ is the population size at year $t$ or $t+1$, $r_{max}$ is
#' the maximum net recruitment rate, $K$ is the pre-exploitation population
#' size, $z$ is the parameter that determines the population size where
#' productivity is maximum. For example, a value of 2.39 corresponds to maximum
#' sustainable yield of $0.6K$ and is assumed by the IWC SC, and $C_t$ is the
#' harvest in numbers in year $t$. Population size can be in either numbers or
#' biomass, however, units will have to be the same as the units used for catch,
#' relative abundance, and absolute abundance.
#'
#' @param r_max The maximum net recruitment rate ($r_{max}$).
#' @param K Pre-exploitation population size in numbers or biomass (depending on
#'   input).
#' @param N1 The population size in numbers or biomass at year 1 (generally
#'   assumed to be K).
#' @param z The parameter that determines the population size where productivity
#'   is maximum (assumed to be 2.39 by the IWC SC).
#' @param start.Yr The first year of the projection (assumed to be the first
#'   year in the catch series).
#' @param num.Yrs The number of projection years. Set as the last year in the
#'   catch or abundance series, whichever is most recent, minus the
#'   \code{start.Yr}.
#' @param catches The time series of catch in numbers or biomass. Currently does
#'   not handle NAs and zeros will have to input a priori for years in which
#'   there were no catches.
#' @param MVP The minimum viable population size in numbers or biomass. Computed
#'   as 4 * \code{\link{num.haplotypes}} to compute minimum viable population
#'   (from Jackson et al., 2006 and IWC, 2007).
#'
#' @return A list of the minimum population size \code{\link{Min.Pop}}, year of the minimum population size \code{\link{Min.Yr}}, a indicator of wether the minimum population size is below the \code{\link{MVP}}, and the predicted population size \code{Pred.N}.
#'
#' @examples
#' num.Yrs  <-  10
#' start.Yr  <-  1
#' r_max  <-  0.2
#' K  <-  1000
#' N1  <-  K
#' catches  <-  round(runif(10, min = 0, max = 150 ), 0 )
#' MVP  <-  0
#' GENERALIZED.LOGISTIC(r_max, K, N1, z, start.Yr, num.Yrs, catches)
GENERALIZED.LOGISTIC <- function(r_max,
                                 K,
                                 N1,
                                 z,
                                 start.Yr,
                                 num.Yrs,
                                 catches,
                                 MVP = 0) {
  ## Variable to indicate whether min population is reached
  Violate.Min.Viable.Pop <- FALSE
  ## Create a vector to hold the model predicted population size
  Pred.N <- rep(NA, num.Yrs)
  ## The first year in the vector above is N1
  Pred.N[1] <- N1

  # Project the population
  for (t in 1:(num.Yrs - 1)) {
    ## FIXME Separate population dynamics function (Rcpp?)
    Pred.N[t + 1] <- max(Pred.N[t] +
                         r_max * Pred.N[t] *
                         (1 - (Pred.N[t] / K)^z) -
                         catches[t],
                         1)
  }

  ## Compute the Nmin
  Min.Pop <- min(Pred.N)
  ## Compute the year at which Nmin occurred
  Min.Yr <- which(Pred.N == Min.Pop) + start.Yr - 1
  if (Min.Pop < MVP) {
    ## Determine whether Nmin is below Min Viable Population
    Violate.Min.Viable.Pop <- TRUE
  }

  list(Min.Pop = Min.Pop,
       Min.Yr = Min.Yr,
       Violate.MVP = Violate.Min.Viable.Pop,
       Pred.N = Pred.N)
}

#' PREDICTED GROWTH RATE
#'
#' \code{PRED.GROWTH.RATE} computes the predicted growth rate if such
#' information is available from an independent estimate rather than being
#' estimated from data. Growth rate is calculated as: $$r_{t_0 - t_{fin}}^{pred}
#' = \frac{ \sum_{t = t_0} ^{t_{fin - 1}} ln \left( \frac{N_{t+1}^{pred}} {
#' N_t^{pred}} \right) } { t_{fin} - t_0 } = \frac{ ln \left( N_{fin}^{pred}
#' \right) - ln \left( N_{0}^{pred} \right)} { t_{fin} - t_0 }$$ where
#' $N^{pred}$ is the model predicted population size, in numbers, at time $t$ or
#' $t+1$ in years, $t_0$ is the start year of the equation (1995 in Zerbini et
#' al. 2011), and $t_{fin}$ is the last year of the equation (1998 in Zerbini et
#' al. 2011).
#'
#' @param growth.rate.Yrs The years to be used for growth rate computation. 1995 - 1998 are used in Zerbini et al. 2011.
#' @param Pred.N Time series of predicted abundance, in numbers, from \code{\link{GENERALIZED.LOGISTIC}}.
#' @param start.Yr The first year of the projection (assumed to be the first year in the catch series).
#'
#' @return A numeric scalar representing predicted growth rate.
#'
#' @examples
#' growth.rate.Yrs  <-  c(1995:1998)
#' Pred.N <- c(1000, 1500, 1500, 2000)
#' start.Yr  <-  1995
#' PRED.GROWTH.RATE(growth.rate.Yrs, Pred.N, start.Yr=start.Yr)
PRED.GROWTH.RATE <- function(growth.rate.Yrs, Pred.N, start.Yr = start.Yr) {
  ## Computing the growth rate years
  GR.Yrs <- growth.rate.Yrs - start.Yr + 1
  Pred.N.GR <- Pred.N[GR.Yrs]

  ## FIXME Just return this line?
  Pred.GR <- (log(Pred.N.GR[length(Pred.N.GR)]) -
              log(Pred.N.GR[1])) / (length(Pred.N.GR) - 1)

  Pred.GR
}

#' Computes the predicted rate of increase for a set of specified years for
#' comparison with trends estimated separately with any of the indices of
#' abundance or count data
#'
#' @param data Count data or relative abundance index to use
#' @param Pred.N Number of individuals predicted
#' @param start.Yr Initial year
#'
#' @return Vector of rates of increase, one per index
#' @export
#'
#' @examples
COMPUTING.ROI <- function(data = data, Pred.N = Pred.N, start.Yr = NULL) {
  num.indices <- max(data$Index)
  Pred.ROI <- rep(NA, num.indices)

  for (i in 1:num.indices) {
    index.ini.year <- (head(subset(data, Index == i)$Year, 1) - start.Yr)
    index.final.year <- (tail(subset(data, Index == i)$Year, 1) - start.Yr)
    elapsed.years <- index.final.year - index.ini.year

    Pred.ROI[i] <- exp((log(Pred.N[index.final.year]) -
                        log(Pred.N[index.ini.year])) /
                       (elapsed.years)) - 1
  }
  Pred.ROI
}

#' Calculate a target K for the bisection method
#'
#' @param r_max The maximum net recruitment rate ($r_{max}$).
#' @param K Pre-expoitation population size in numbers or biomass
#'   (depending on input).
#' @param N1 Population size in numbers or biomass at year 1 (generally
#'   assumed to be K).
#' @param z Generalized logistic shape parameter, determines population
#'   size where productivity is masimum (assumed to be 2.39 by the ISC
#'   SC).
#' @param num.Yrs The number of projection years. Set as the last year
#'   in the catchor abundance series whichever is most recent, minus the
#'   start year.
#' @param start.Yr First year of the projection (assumed to be the first
#'   year in the catch series).
#' @param target.Pop Target population size.
#' @param catches Catch time series. Cannot include NAs,
#' @param MVP Minimum Viable Population Size; `4 * num.haplotypes`
#'
#' @return Vector of differences between predicted population and target
#'   population.
#' @export
#'
#' @examples
#' TARGET.K(r_max, K, N1, z, start.Yr=start.Yr, num.Yrs=bisection.Yrs,
#'          target.Pop=target.Pop, catches=catches, MVP=MVP)
TARGET.K <- function(r_max,
                     K,
                     N1,
                     z,
                     num.Yrs,
                     start.Yr,
                     target.Pop,
                     catches=catches,
                     MVP = 0) {
  Pred.N <- GENERALIZED.LOGISTIC(r_max = r_max,
                                 K = K,
                                 N1 = K,
                                 z = z,
                                 start.Yr = start.Yr,
                                 num.Yrs = num.Yrs,
                                 catches = catches,
                                 MVP = MVP)
  ## FIXME Return this line?
  diff <- Pred.N$Pred.N[num.Yrs] - target.Pop
  #plot(seq(1901, 1900+num.Yrs, 1), Pred.N$Pred.N, ylim=c(0,50000))
  #print(paste("N-Target=", Pred.N$Pred.N[num.Yrs], "Target.Pop=", target.Pop, "K=", K))
  diff
}

#' LOGISTIC BISECTION
#'
#' Method of Butterworth and Punt (1995) where the prior distribution of the
#' current absolute abundance $N_{2005}$ and maximum net recruitment rate
#' \code{r_max} are sampled and then used to determine the unique value of the
#' population abundance $N$ in \code{start.Yr} (assumed to correspond to
#' carrying capacity $K$). Requires \code{\link{TARGET.K}} and subsequent
#' dependencies.
#'
#' @param K.low Lower bound for $K$ when preforming the bisection method of Punt
#'   and Butterworth (1995). Default is 1.
#' @param K.high Upper bound for $K$ when preforming the bisection method of
#'   Punt and Butterworth (1995). Default is 500,000.
#' @param r_max The maximum net recruitment rate ($r_{max}$).
#' @param z The parameter that determines the population size where productivity
#'   is maximum (assumed to be 2.39 by the IWC SC).
#' @param num.Yrs The number of projection years. Set as the last year in the
#'   catch or abundance series, whichever is most recent, minus the
#'   \code{start.Yr}.
#' @param start.Yr The first year of the projection (assumed to be the first
#'   year in the catch series).
#' @param target.Pop A sample of the prior on population abundance $N$, in
#'   numbers, set as \code{sample.N.obs} sampled from \code{priors$N.obs}
#' @param catches The time series of catch in numbers or biomass. Currently does
#'   not handle NAs and zeros will have to input a priori for years in which
#'   there were no catches.
#' @param MVP The minimum viable population size in numbers or biomass. Computed
#'   as 4 * \code{\link{num.haplotypes}} to compute minimum viable population
#'   (from Jackson et al., 2006 and IWC, 2007).
#' @param tol The desired accuracy (convergence tolerance) of
#'   \code{\link{stats::uniroot}}.
#'
#' @return A numeric scalar of an estimate of  carrying capacity $K$.
#'
#' @examples
#' LOGISTIC.BISECTION.K(K.low = 1, K.high = 100000, r_max = r_max, z = z,
#'                      num.Yrs = bisection.Yrs, start.Yr = start.Yr,
#'                      target.Pop = target.Pop, catches = catches, MVP = MVP,
#'                      tol = 0.001)
LOGISTIC.BISECTION.K <- function(K.low,
                                 K.high,
                                 r_max,
                                 z,
                                 num.Yrs,
                                 start.Yr,
                                 target.Pop,
                                 catches,
                                 MVP,
                                 tol = 0.001) {
  Kmin <- uniroot(TARGET.K,
                  tol = tol,
                  c(K.low,
                    K.high),
                  r_max = r_max,
                  z = z,
                  num.Yrs = num.Yrs,
                  start.Yr = start.Yr,
                  target.Pop = target.Pop,
                  catches = catches,
                  MVP = MVP)
  Kmin$root
}

#' Compute analytic estimates of q, the scaling parameter between indices and
#' absolute population size
#'
#' @param rel.Abundance Relative abundance index
#' @param Pred.N Predicted population
#' @param start.Yr Initial year
#' @param add_CV Coefficient of variation
#' @param num.IA Index of abundance
#'
#' @return A numeric estimator for $q$.
#' @export
#'
#' @examples
CALC.ANALYTIC.Q <- function(rel.Abundance, Pred.N, start.Yr,
                            add_CV = 0, num.IA) {
  ## Vector to store the q values
  analytic.Q <- rep(NA, num.IA)

  for (i in 1:num.IA) {
    ## Subseting across each index of abundance
    IA <- subset(rel.Abundance, Index == i)
    ## Years for which IAs are available
    IA.yrs <- IA$Year-start.Yr + 1
    ## Computing the value of sigma as in Zerbini et al. 2011
    IA$Sigma <- sqrt(log(1 + IA$CV.IA.obs^2))
    ## Numerator of the analytic q estimator (Zerbini et al., 2011 - eq. (3))
    qNumerator <- sum((log(IA$IA.obs / Pred.N[IA.yrs])) /
                      (IA$Sigma * IA$Sigma + add_CV * add_CV))
    ## Denominator of the analytic q estimator (Zerbini et al., 2011 - eq. (3))
    qDenominator <- sum(1 / (IA$Sigma * IA$Sigma))
    ## Estimate of q
    analytic.Q[i] <- exp(qNumerator / qDenominator)
  }
  analytic.Q
}

#' Compute the log-likelihood of indices of abundance
#'
#' @param Rel.Abundance Relative abundance
#' @param Pred.N Predicted population size
#' @param start.Yr Initial year
#' @param q.values Scaling parameter
#' @param add_CV Coefficient of variation
#' @param num.IA Number of indices of abundance
#' @param log Boolean, return log likelihood (default TRUE) or
#'   likelihood.
#'
#' @return List of likelihood based on Zerbini et al. (2011) eq. 5 or using `dnorm`
#' @export
#'
#' @examples
LNLIKE.IAs <- function(Rel.Abundance, Pred.N, start.Yr,
                       q.values, add_CV, num.IA, log = TRUE) {
  loglike.IA1 <- 0
  loglike.IA2 <- 0

  for(i in 1:num.IA) {
    ## Subseting across each index of abundance
    IA <- subset(Rel.Abundance, Index == i)
    ## Years for which IAs are available
    IA.yrs <- IA$Year-start.Yr + 1
    ## This is the likelihood from Zerbini et al. 2011 (eq. 5)
    loglike.IA1 <- loglike.IA1 +
      ((sum(log(IA$Sigma) + log(IA$IA.obs) + 0.5 *
            ((((log(q.values[i] * Pred.N[IA.yrs]) - log(IA$IA.obs))^2) /
              (IA$Sigma*IA$Sigma))))))

    ## This is the log-normal distribution from R (using function dnorm)
    ## FIXME Why is this better than using builtin `dlnorm`?
    loglike.IA2 <- loglike.IA2 +
      CALC.LNLIKE(Obs.N = IA$IA.obs,
                  Pred.N = (q.values[i] * Pred.N[IA.yrs]),
                  CV =  sqrt(IA$Sigma * IA$Sigma + add_CV * add_CV),
                  log = log)
  }
  list(loglike.IA1 = loglike.IA1, loglike.IA2 = loglike.IA2)
}

#' LOG LIKELIHOOD OF ABSOLUTE ABUNDANCE
#'
#' This function computes two estimates of the log-likelihood of the estimated
#' absolute abundance using the equation from Zerbini et al. 2011 (eq. 4) and a
#' lognormal distribution from \code{\link{CALC.LNLIKE}}.
#'
#' @param Obs.N Observed absoluted abundance in numbers as a data.frame
#'   containing year, estimate of absolute abundance, and CV.
#' @param Pred.N Predicted absolute abundance in numbers from
#'   \code{\link{GENERALIZED.LOGISTIC}}.
#' @param start.Yr The first year of the projection (assumed to be the first
#'   year in the catch series).
#' @param add_CV Additional CV to add to variance of lognormal distribution
#'   sampled from \code{priors$add_CV}.
#' @param log Return the log of the likelihood (TRUE/FALSE)
#'
#' @return A list of two numeric scalars of estimates of log-likelihood.
#'
#' @examples
#' Obs.N  <-  data.frame(Year = 2005, Sigma = 5, Obs.N = 1000)
#' Pred.N  <-  1234
#' start.Yr  <-  2005
#' LNLIKE.Ns(Obs.N, Pred.N, start.Yr, add_CV = 0, log=TRUE)
LNLIKE.Ns <- function(Obs.N, Pred.N, start.Yr, add_CV, log = TRUE) {
  loglike.Ns1 <- 0
  loglike.Ns2 <- 0

  ## Years for which Ns are available
  N.yrs <- Obs.N$Year-start.Yr+1
  ## This is the likelihood from Zerbini et al. 2011 (eq. 4)
  loglike.Ns1 <- loglike.Ns1 +
    ((sum(log(Obs.N$Sigma) + log(Obs.N$N.obs) + 0.5 *
          ((((log(Pred.N[N.yrs]) - log(Obs.N$N.obs))^2) /
            (Obs.N$Sigma * Obs.N$Sigma + add_CV * add_CV))))))
  ## This is the log-normal distribution from R (using function dnorm)
  ## FIXME See comments above re: `dlnorm`
  loglike.Ns2 <- loglike.Ns2 + CALC.LNLIKE(Obs.N = Obs.N$N.obs,
                                           Pred.N = (Pred.N[N.yrs]),
                                           CV = sqrt(Obs.N$Sigma * Obs.N$Sigma +
                                                     add_CV * add_CV),
                                           log = log)

  list(loglike.Ns1 = loglike.Ns1, loglike.Ns2 = loglike.Ns2)
}

#' Calculate the log-likelihood of the growth rate
#'
#' Calculates the log-likelihood of the estimated growth rate given the observed
#' growth rate and the standard deviation of the observed growth rate.
#'
#' @param Obs.GR Observed growth rate
#' @param Pred.GR Predicted growth rate
#' @param GR.SD.Obs Standard error of the observed growth rate
#'
#' @return A \code{list} containing \code{loglike.GR1} and \code{loglike.GR2}
#'
#' @examples
#' LNLIKE.GR(0.1, 0.1, 0.1)
LNLIKE.GR <- function(Obs.GR, Pred.GR, GR.SD.Obs) {
  ## TODO Does this need to recalculate and return *both* of these values?
  loglike.GR1 <- 0
  loglike.GR2 <- 0

  ## This is the likelihood from Zerbini et al. 2011 (eq. 6)
  loglike.GR1 <- loglike.GR1 + (((log(GR.SD.Obs) + 0.5 * (((Pred.GR-Obs.GR) / GR.SD.Obs)^2))))

  loglike.GR2 <- loglike.GR2 + CALC.LNLIKE(Obs.N = Obs.GR,
                                           Pred.N = Pred.GR,
                                           CV = GR.SD.Obs,
                                           log = FALSE)

  list(loglike.GR1 = loglike.GR1, loglike.GR2 = loglike.GR2)
}

#' Function to calculate the log-likelihood using a lognormal distribution
#'
#' @param Obs.N Time series of observed abundance
#' @param Pred.N Time series of estimated abundance
#' @param CV coefficient of variation
#' @param log whether to export as log-likelihood
#'
#' @return returns a scalar of the likelihood
#'
#' @examples
#' Obs.N <- 2000
#' Pred.N <- 2340
#' CV <- 4
#' CALC.LNLIKE(Obs.N, Pred.N, CV)
CALC.LNLIKE <- function(Obs.N, Pred.N, CV, log = FALSE) {
  sum(dnorm(x = log(Obs.N), mean = log(Pred.N), sd = CV, log = log))
}

#' OUTPUT FUNCTION
#'
#' Function that provides a summary of SIR outputs including: mean, median, 95%
#' credible interval, 90% predicitive interval, max, and sample size.
#'
#' @param x A data.frame of model outputs including: sample.r_max, sample.K,
#'   sample.N.obs, sample.add_CV, Pred.N$Min.Pop, Pred.N$Min.Yr,
#'   Pred.N$Violate.MVP, c(Pred.N$Pred.N[output.Yrs-start.Yr+1]), Pred.ROI.IA,
#'   q.sample.IA, Pred.ROI.Count, q.sample.Count, lnlike.IAs[[1]],
#'   lnlike.Count[[1]], lnlike.Ns[[1]], lnlike.GR[[1]], LL, Likelihood,
#'   Pred.N$Min.Pop/sample.K, c(Pred.N$Pred.N[output.Yrs-start.Yr+1]/sample.K),
#'   draw, save)
#' @param scenario Name of the model run and object as specified by the user.
#'
#' @return Returns a data.frame with summary of SIR outputs
#'
#' @examples
#' x  <-  rnorm(1000, 5, 7)
#' y  <-  rnorm(1000, 6, 9)
#' df <- data.frame(x = x, y = y)
#' SUMMARY.SIR( df , scenario = "example_summary")
SUMMARY.SIR <- function(x, scenario = "USERDEFINED") {
  num.col <- dim(x)[2]
  col.names <- names(x)
  row.names <- c("mean", "median",
                 "2.5%PI", "97.5%PI",
                 "5%PI", "95%PI",
                 "min", "max", "n")

  ## FIXME Only call quantile once?
  output.summary <- matrix(nrow = length(row.names), ncol = num.col)
  output.summary[1, ] <- sapply(x, mean)
  output.summary[2, ] <- sapply(x, median)
  output.summary[3, ] <- sapply(x, quantile, probs=0.025)
  output.summary[4, ] <- sapply(x, quantile, probs=0.975)
  output.summary[5, ] <- sapply(x, quantile, probs=0.05)
  output.summary[6, ] <- sapply(x, quantile, probs=0.95)
  output.summary[7, ] <- sapply(x, min)
  output.summary[8, ] <- sapply(x, max)
  output.summary[9, ] <- sapply(x, length)

  output.table <- data.frame(output.summary)
  names(output.table) <- col.names
  row.names(output.table) <- row.names
  noquote(format(output.table, digits = 3, scientific = FALSE))

  list(scenario=scenario, date=Sys.time(), output.table=output.table)
}

