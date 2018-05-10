# R Script for Bayesian Assessment Model of Humpback Whales - uses the equations in Zerbini et al. (2011)
# Code cleaned up and sent to Andre Punt, John Best and Grant Adams on 7 Dec 2017

#Example call


####################### FUNCTIONS ######################################################
########################################################################################

############ HOUSEKEEPING FUNCTION #############
# HUMPBACK SIR controls the sampling from the priors, the bisection and likelihoods and the output functions
################################################
#' HUMPBACK SIR controls the sampling from the priors, the bisection and likelihoods and the output functions
#'
#' @param file.name name of a file to identified the files exported by the function
#' @param n.samples number of samples for the Rubin SIR: NOT USED
#' @param n.resamples number of resamples to compute the marginal posterior distributions
#' @param prior.K prior for K for future use with the forward method: NOT USED
#' @param prior.r_max prior for r_max. The first element identifies the sampling distribution, the send identifies the lower bound (for uniform distribution) or the mean (for a normal distribution), and the third element corresponds to the upper bound (uniform distribution) or the standard error (normal distribution). The default is Uniform with bounds c(0, 0.12)
#' @param r_max.bound bounds for the r_max prior. Default is c(0, 0.12)
#' @param prior.N.obs prior distribution for a recent abudnance estimate. Elements are equivalent to the prior on r_max
#' @param prior.add.CV prior for additional CV if applicable. It is defined by four elements: (1) the sampling distribution, (2) lower bound or mean, (3) upper bound or SE, (4) boolean variable to specify whether CV add is used or not in the likelihood. Default is a uniform distribution bounded by (0,1) and FALSE (=CV.add not used)
#' @param prior.z prior on the shape parameter. NOT USED, assumed z=2.39 (max productivity at K=0.6)
#' @param q.prior.IA prior on q for indices of abundance. Definition of elements is similar to the prior on CV.add. If the fourth element = FALSE, an analytical solution for q is used (as in Zerbini et al. 2011)
#' @param q.prior.Count similar to q.prior.IA, but for count data. NOT CURRENTLY USED
#' @param Klim bounds for K when preforming the bisection method of Punt and Butterworth (1995). Defined by two elements, the lower and upper bounds. Default is (1, 500000)
#' @param target.Yr year of the target population estimate for the bisection method. Default is 2008
#' @param num.haplotypes number of haplotypes to compute minimum viable population (from Jackson et al., 2006 and IWC, 2007)
#' @param tolerance.for.bisection tolerance value for performing the bisection method
#' @param output.Yrs year for outputing the predicted abundance estimates. Default is 2008, but multiple years can be specified. For example, if outputs for 2005 and 2008 are needed, output.Yrs = c(2005, 2008)
#' @param abs.abundance R object containing year, estimate of absolute abundance, and CV (see example)
#' @param rel.abundance R object containing years, estimates of relative abudnance and CVs (see example)
#' @param rel.abundance.key key to speficy if relative abundance data are used in the likelihood. Default is TRUE
#' @param count.data R object containing years, estimates of counts and effort. NOT USED
#' @param count.data.key key to speficy in count data are used. Default is FALSE. NOT USED
#' @param growth.rate.obs observed growth rate (1st element) and standard error (2nd element) as in Zerbni et al. (2011). If third element is FALSE, the growth rate is not included in the likelihood
#' @param growth.rate.Yrs Years for which the growth.rate.obs were computed (as in Zerbini et al., 2011)
#' @param catch.data R object containing the years and catches (see example)
#' @param Threshold threshold for the McCallister et al. (1994) SIR. This is data-specific. Default is 1e-100.
#' @param Print key to print various pieces of information as the code runs. Default is 0 (nothing is printed)
#'
#' @return A \code{list} containing posterior samples and metadata
#'
#' TO DO:
#' 1. add the negative binomial likelihood for the count data, which is not currently used even though it is defined in the main function call.
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' HUMPBACK.SIR(file.name = "test.N2005",
#'              n.samples = NULL,
#'              n.resamples = 1000,
#'              prior.K = c(NA, NA, NA),
#'              prior.r_max = c("uniform", 0, 0.106),
#'              r_max.bound = c(0, 0.106),
#'              prior.N.obs = c("uniform", 500, 20000),
#'              prior.add.CV = c("uniform", 0, 1, FALSE),
#'              prior.z = c(NA, 2.39, NA),
#'              q.prior.IA = c("uniform", 0, 1, FALSE),
#'              q.prior.Count = c("uniform", 0, 1, FALSE),
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
#' }
HUMPBACK.SIR=function(file.name="NULL", n.samples=1000, n.resamples=1000, prior.K=c(NA, NA, NA), prior.r_max=c("uniform", 0, 0.12), r_max.bound=c(0, 0.12), prior.N.obs=c("uniform", 0, 50000), prior.add.CV=c("uniform", 0, 1, TRUE), prior.z=c(NA, 2.39, NA), q.prior.IA=c("uniform", 0, 1, FALSE), q.prior.Count=c("uniform", 0, 1, FALSE), Klim=c(1, 500000), target.Yr=2008, num.haplotypes=66, tolerance.for.bisection=0.001, output.Yrs=c(2008), abs.abundance=Abs.Abundance, rel.abundance=Rel.Abundance, rel.abundance.key=TRUE, count.data=NULL, count.data.key=FALSE, growth.rate.obs=c(0.074, 0.033, TRUE), growth.rate.Yrs=c(1995, 1996, 1997, 1998), catch.data=Catch.data, Threshold=1e100, Print=0)
{
  require(utils)
  begin.time=Sys.time()

  ################################
  #Assingning variables
  ################################
  n.samples=n.samples #not used as I am not doing the Rubin (1988) SIR anymore.
  target.Yr=target.Yr
  start.Yr=catch.data$Year[1] #the first year of the projection is set as the first year in the catch series
  end.Yr=max(tail(catch.data$Year,1), max(abs.abundance$Year), max(rel.abundance$Year)) #the last year of the projection is set as the last year in the catch or abundance series, whichever is most recent
  bisection.Yrs=target.Yr-start.Yr+1 #setting the target year for the bisection method
  projection.Yrs=end.Yr-start.Yr+1 #setting the years to project

  z=prior.z[2]
  catches=catch.data$Catch #assigning the catch data
  num.IA=max(rel.abundance$Index) #determining the number of Indices of Abundance available
  num.Count=max(count.data$Index) #determining the number of Count Data sets available
  rel.abundance$Sigma=sqrt(log(1+rel.abundance$CV.IA.obs^2)) #computing the value of sigma as in Zerbini et al. 2011
  count.data$Sigma=sqrt(log(1+count.data$CV.IA.obs^2)) #computing the value of sigma for the count data as in Zerbini et al. (2011)
  abs.abundance$Sigma=sqrt(log(1+abs.abundance$CV.obs^2)) #computing the value of sigma as in Zerbini et al. 2011
  MVP=4*num.haplotypes #computing the minimum viable population, if num.haplotypes=0, assumes no MVP

  tol=tolerance.for.bisection #assigning the value for tolerance in computing K using the uniroot function (bisection method)
  i=0 #start the loop
  draw=1 #keep track of number of draws
  Cumulative.Likelihood=0

  #Creating output vectors
  #-------------------------------------
  names=c("r_max", "K", "sample.N.obs", "add.CV", "Nmin", "YearMin", "violate_MVP", paste("N", output.Yrs, sep=""), paste("ROI_IA", unique(rel.abundance$Index), sep=""), paste("q_IA", unique(rel.abundance$Index), sep=""), paste("ROI_Count", unique(count.data$Index), sep=""), paste("q_Count", unique(count.data$Index), sep=""), "NLL.IAs", "NLL.Count", "NLL.N", "NLL.GR", "NLL", "Likelihood", "Max_Dep", paste("status", output.Yrs, sep=""), "draw", "save")

  samples.output=matrix(0, nrow=1, ncol=length(names))
  resamples.output=matrix(0, nrow=1, ncol=length(names))
  resamples.trajectories=matrix(NA, nrow=1, ncol=projection.Yrs)
  final.trajectory=matrix(NA, nrow=projection.Yrs, ncol=6)
  Year=seq(start.Yr, end.Yr, by=1)

  #Initiating the SIR loop

  while(i<n.resamples)
  {

    #Sampling from Priors
    #-------------------------------
    save=FALSE #variable to indicate whether a specific draw is kept

    #Sampling for r_max
    sample.r_max=0.2 #setting sample.r_max outside of the bound
    while(sample.r_max<r_max.bound[1] | sample.r_max>r_max.bound[2])
    {
      sample.r_max=SAMPLE.PRIOR(name=prior.r_max[1], Val.1=prior.r_max[2], Val.2=prior.r_max[3]) #prior on r_max, keep if within boundaries
    }

    #sampling from the N.obs prior
    sample.N.obs=SAMPLE.PRIOR(name=prior.N.obs[1], Val.1=prior.N.obs[2], Val.2=prior.N.obs[3]) #prior on N_obs

    if(prior.add.CV[4]==TRUE) #prior on additional CV
    {
      sample.add.CV=SAMPLE.PRIOR(name=prior.add.CV[1], Val.1=prior.add.CV[2], Val.2=prior.add.CV[3])
    }
    else
    {
      sample.add.CV=0
    }

    #Sampling from q priors if q.prior is TRUE
    if(q.prior.IA[4]==TRUE) #priors on q for indices of abundance
    {
      q.sample.IA=rep(NA, num.IA)

      for (i in 1:num.IA)
      {
        q.sample.IA[i]=SAMPLE.PRIOR(name=q.prior.IA[1], Val.1=q.prior.IA[2], Val.2=q.prior.IA[3])
      }
    }
    else
    {
      q.sample.IA=rep(-9999, length(unique(rel.abundance$Index)))
    }

    if(q.prior.Count[4]==TRUE) #priors on q for count data
    {
      q.sample.Count=rep(NA, 1:num.Count)

      for (i in num.Count)
      {
        q.sample.Count[i]=SAMPLE.PRIOR(name=q.prior.Count[1], Val.1=q.prior.Count[2], Val.2=q.prior.Count[3])
      }
    }
    else
    {
      q.sample.Count=rep(-9999, length(unique(count.data$Index)))
    }

    sample.K=LOGISTIC.BISECTION.K(K.low=Klim[1], K.high=Klim[2], r_max=sample.r_max, z=z, num.Yrs=bisection.Yrs, start.Yr=start.Yr, target.Pop=sample.N.obs, catches=catches, MVP=MVP, tol=tol) #K_prior computed with the bisection/uniroot method

    #Computing the predicted abundances with the samples from the priors
    #----------------------------------------
    Pred.N=GENERALIZED.LOGISTIC(r_max=sample.r_max, K=sample.K, N1=sample.K, z=z, start.Yr=start.Yr, num.Yrs=projection.Yrs, catches=catches, MVP=MVP)
    #Print results if required
    if(Print==1){print(paste("sample.r_max=", sample.r_max, "sample.N.obs=", sample.N.obs, "sample.K=", sample.K, "Pred.N.target=", Pred.N$Pred.N[bisection.Yrs]))}


    #Computing the predicted ROI for the IAs and Count data, if applicable
    #----------------------------------------
    #For IAs
    if(rel.abundance.key==TRUE)
    {
      Pred.ROI.IA=COMPUTING.ROI(data=rel.abundance, Pred.N=Pred.N$Pred.N, start.Yr=start.Yr)
    }
    else
    {
      Pred.ROI.IA=rep(0, num.IA)
    }

    #For Count Data
    if(count.data.key==TRUE)
    {
      Pred.ROI.Count=COMPUTING.ROI(data=count.data, Pred.N=Pred.N$Pred.N, start.Yr=start.Yr)
    }
    else
    {
      Pred.ROI.Count=rep(0, num.Count)
    }



    #Calculate Analytical Qs if rel.abundance.key is TRUE
    #---------------------------------------------------------
    if(rel.abundance.key==TRUE)
      if(q.prior.IA[4]==FALSE)
      {
        q.sample.IA=CALC.ANALYTIC.Q(rel.abundance, Pred.N$Pred.N, start.Yr, sample.add.CV, num.IA)
      }
    else
    {
      q.sample.IA=q.sample.IA
    }

    #browser()

    #Calculate Analytical Qs if count.data.key is TRUE (NOT USED YET - AZerbini, Feb 2013)
    #--------------------------------------------------------
    if(rel.abundance.key==TRUE)
      if(q.prior.Count[4]==FALSE)
      {
        q.sample.Count=CALC.ANALYTIC.Q(count.data, Pred.N$Pred.N, start.Yr, sample.add.CV, num.Count)
      }
    else
    {
      q.sample.Count=q.sample.Count
    }

    if(Print==1){print(paste("q.IAs=", q.sample.IA)); print(paste("q.Count=", q.sample.Count))}


    #Compute the likelihoods
    #--------------------------------
    # (1) relative indices (if rel.abundance.key is TRUE)
    if(rel.abundance.key==TRUE)
    {
      lnlike.IAs=LNLIKE.IAs(rel.abundance, Pred.N$Pred.N, start.Yr, q.sample.IA, sample.add.CV, num.IA, log=TRUE)
    }
    else
    {
      lnlike.IAs=0
    }
    if(Print==1){print(paste("lnlike.IAs=", lnlike.IAs))}

    # (2) count data (if count.data.key is TRUE)
    if(count.data.key==TRUE)
    {
      lnlike.Count=LNLIKE.IAs(count.data, Pred.N$Pred.N, start.Yr, q.sample.Count, sample.add.CV, num.Count, log=TRUE)
    }
    else
    {
      lnlike.Count=0
    }
    if(Print==1){print(paste("lnlike.Count=", lnlike.Count))}

    # (3) absolute abundance
    lnlike.Ns=LNLIKE.Ns(abs.abundance, Pred.N$Pred.N, start.Yr, sample.add.CV, log=TRUE)
    if(Print==1){print(paste("lnlike.Ns=", lnlike.Ns))}

    # (4) growth rate if applicable
    if(growth.rate.obs[3]==TRUE)
    {
      Pred.GR=PRED.GROWTH.RATE(growth.rate.Yrs=growth.rate.Yrs, Pred.N=Pred.N$Pred.N, start.Yr=start.Yr)
      lnlike.GR=LNLIKE.GR(Obs.GR=growth.rate.obs[1], Pred.GR=Pred.GR, GR.SD.Obs=growth.rate.obs[2])
    }
    else
    {
      lnlike.GR=0
    }
    if(Print==1){print(paste("lnlike.GR=", lnlike.GR))}

    LL=lnlike.IAs[[1]]+lnlike.Count[[1]]+lnlike.Ns[[1]]+lnlike.GR[[1]] #these use the likelihoods in Zerbini et al. (2011)
    Likelihood=exp(-LL)
    if(Print==1){print(paste("NLL=", LL, "Likelihood=", Likelihood))}


    if(Pred.N$Violate.MVP==TRUE)
    {
      Likelihood=0
      print(paste("MVP violated on draw", draw))
    }

    Cumulative.Likelihood=Cumulative.Likelihood+Likelihood
    #print(paste("Sample=", i, "Likelihood=", Likelihood))

    if(Pred.N$Violate.MVP==FALSE)
    {
      while(Cumulative.Likelihood>Threshold)
      {
        print(paste("sample=", i, " & draw=", draw))
        if(Print==1){print(paste("draw=", draw, "Likelihood=", Likelihood, "Cumulative=", Cumulative.Likelihood))}
        save=TRUE
        Cumulative.Likelihood=Cumulative.Likelihood-Threshold
        resamples.trajectories=rbind(resamples.trajectories, Pred.N$Pred.N)
        resamples.output=rbind(resamples.output, c(sample.r_max, sample.K, sample.N.obs, sample.add.CV, Pred.N$Min.Pop, Pred.N$Min.Yr, Pred.N$Violate.MVP, c(Pred.N$Pred.N[output.Yrs-start.Yr+1]), Pred.ROI.IA, q.sample.IA, Pred.ROI.Count, q.sample.Count, lnlike.IAs[[1]], lnlike.Count[[1]], lnlike.Ns[[1]], lnlike.GR[[1]], LL, Likelihood, Pred.N$Min.Pop/sample.K, c(Pred.N$Pred.N[output.Yrs-start.Yr+1]/sample.K), draw, save))
        i=i+1
      }
    }

    samples.output=rbind(samples.output, c(sample.r_max, sample.K, sample.N.obs, sample.add.CV, Pred.N$Min.Pop, Pred.N$Min.Yr, Pred.N$Violate.MVP, c(Pred.N$Pred.N[output.Yrs-start.Yr+1]), Pred.ROI.IA, q.sample.IA, Pred.ROI.Count, q.sample.Count, lnlike.IAs[[1]], lnlike.Count[[1]], lnlike.Ns[[1]], lnlike.GR[[1]], LL, Likelihood, Pred.N$Min.Pop/sample.K, c(Pred.N$Pred.N[output.Yrs-start.Yr+1]/sample.K), draw, save))

    draw=draw+1
  }

  samples.output=data.frame(samples.output)
  names(samples.output)=names
  samples.output=samples.output[-1,]
  samples.output.summary=SUMMARY.SIR(x=samples.output, scenario=file.name)
  write.csv(samples.output, paste(file.name,"_","samples.output.csv", sep=""))

  resamples.output=data.frame(resamples.output)
  names(resamples.output)=names
  resamples.output=resamples.output[-1,]
  resamples.output.summary=SUMMARY.SIR(x=resamples.output, scenario=file.name)
  write.csv(resamples.output, paste(file.name,"_","resamples.output.csv", sep=""))

  resamples.trajectories=data.frame(resamples.trajectories)
  names(resamples.trajectories)=seq(start.Yr, end.Yr, 1)
  resamples.trajectories=resamples.trajectories[-1,]
  write.csv(resamples.trajectories, paste(file.name,"_","resample.trajectories.csv", sep=""))

  final.trajectory[,1]=sapply(resamples.trajectories, mean)
  final.trajectory[,2]=sapply(resamples.trajectories, median)
  final.trajectory[,3]=sapply(resamples.trajectories, quantile, probs=c(0.025))
  final.trajectory[,4]=sapply(resamples.trajectories, quantile, probs=c(0.975))
  final.trajectory[,5]=sapply(resamples.trajectories, quantile, probs=c(0.05))
  final.trajectory[,6]=sapply(resamples.trajectories, quantile, probs=c(0.95))
  final.trajectory=data.frame(final.trajectory)
  names(final.trajectory)= c("mean", "median", "PI.2.5%", "PI.97.5%", "PI.5%", "PI.95%")
  final.trajectory=data.frame(Year, final.trajectory)

  resamples.per.samples=dim(samples.output)[1]/dim(resamples.output)[1]

  end.time=Sys.time()
  print(paste("Time to Compute=", (end.time-begin.time)))

  return(list(call=call, file.name=file.name, Date.Time=Sys.time(), Time.to.compute.in.minutes=paste((end.time-begin.time)/60), Threshold=Threshold, Ratio.Resamples.per.Sample=paste("1 resample",":", resamples.per.samples, "samples"),  resamples.output=resamples.output, resamples.output.summary=resamples.output.summary$output.table, samples.output.summary=samples.output.summary$output.table, final.trajectory=final.trajectory, inputs=list(draws=draw, n.resamples=n.resamples, prior.r_max=prior.r_max, prior.N.obs=prior.N.obs, target.Yr=target.Yr, MVP=paste("num.haplotypes=", num.haplotypes, "MVP=", 4*num.haplotypes), tolerance=tolerance.for.bisection, output.Years=output.Yrs)))

}
#### END OF FUNCTION ###################################


#################################################
# FUNCTION TO SAMPLE FROM PRIORS
# This function samples from 4 distributions given two input values (Val.1 and Val.2). If the distribution
# is uniform or log-uniform, then Val.1 and Val.2 correspond to the lower and upper bounds. If normal
# or log-normal, Val.1 corresponds to the mean and Val.2 to the SD.
# -----------------------------------------------
SAMPLE.PRIOR=function(name=NA, Val.1=NA, Val.2=NA)
{
  #uniform prior
  if(name=="uniform"){return(runif(1,as.numeric(Val.1), as.numeric(Val.2)))}
  #log-uniform prior
  else if(name=="log-uniform"){return(exp(runif(1, as.numeric(Val.1), as.numeric(Val.2))))}
  #normal prior
  else if(name=="normal"){return(rnorm(1, as.numeric(Val.1), as.numeric(Val.2)))}
  #log-normal prior
  else if(name=="log-normal"){return(rlnorm(1, as.numeric(Val.1), as.numeric(Val.2)))}
}

# END OF FUNCTION
# -------------------------


###################################################
# GENERALISED LOGISTIC MODEL
# This function does the population projection using a Pella-Tomlison population dynamics model
###################################################
GENERALIZED.LOGISTIC=function(r_max, K, N1, z, start.Yr, num.Yrs, catches, MVP=0)
{
  Violate.Min.Viable.Pop=FALSE #variable to indicate whether min population is reached
  Pred.N=rep(NA, num.Yrs) #create a vector to hold the model predicted population size
  Pred.N[1]=N1 #The first year in the vector above is N1

  #Project the population
  for(t in 1:(num.Yrs-1))
  {
    Pred.N[t+1]=max(Pred.N[t]+r_max*Pred.N[t]*(1-(Pred.N[t]/K)^z)-catches[t],1)
  }

  Min.Pop=min(Pred.N) #Compute the Nmin
  Min.Yr=which(Pred.N==Min.Pop)+start.Yr-1 #Compute the year at which Nmin occurred
  if(Min.Pop<MVP){Violate.Min.Viable.Pop=TRUE} #Determine whether Nmin is below Min Viable Population

  return(list(Min.Pop=Min.Pop, Min.Yr=Min.Yr, Violate.MVP=Violate.Min.Viable.Pop, Pred.N=Pred.N))
}
#END FUNCTION
#---------------------------------

#######################################################
#THiS FUNCTION COMPUTES THE PREDICTED GROWTH RATE IF SUCH INFORMATION IS AVAILABLE
#FROM A INDEPENDENT ESTIMATE (EQUATION 2 in Zerbini et al., 2011) AND WILL NOT BE
#COMPUTED FROM INPUT DATA IN THIS MODEL
#######################################################
PRED.GROWTH.RATE=function(growth.rate.Yrs, Pred.N, start.Yr=start.Yr)
{
  GR.Yrs=growth.rate.Yrs-start.Yr+1 #computing the growth rate years
  Pred.N.GR=Pred.N[GR.Yrs]

  Pred.GR=(log(Pred.N.GR[length(Pred.N.GR)])-log(Pred.N.GR[1]))/(length(Pred.N.GR)-1)

  return(Pred.GR)
}
# END OF FUNCTION
#----------------------------------------

########################################
#THIS FUNCTION COMPUTES THE PREDICTED RATE OF INCREASE FOR A SET OF SPECIFIED
#YEARS FOR COMPARISON WITH TRENDS ESTIMATED SEPARATELY WITH ANY OF THE INDICES OF
#ABUNDANCE OR COUNT DATA
########################################
COMPUTING.ROI=function(data=data, Pred.N=Pred.N, start.Yr=NULL)
{
  num.indices=max(data$Index)
  Pred.ROI=rep(NA, num.indices)

  for(i in 1:num.indices)
  {
    index.ini.year=(head(subset(data, Index==i)$Year,1)-start.Yr)
    index.final.year=(tail(subset(data, Index==i)$Year, 1)-start.Yr)
    elapsed.years=index.final.year-index.ini.year

    Pred.ROI[i]=exp((log(Pred.N[index.final.year])-log(Pred.N[index.ini.year]))/(elapsed.years))-1

  }
  return(Pred.ROI=Pred.ROI)
}
#END OF FUNCTION
#--------------------------


####################################################################
# FUNCTION TO CALCULATE A TARGET K FOR THE BISECTION METHOD
###################################################################

#example call: TARGET.K(r_max, K, N1, z, start.Yr=start.Yr, num.Yrs=bisection.Yrs, target.Pop=target.Pop, catches=catches, MVP=MVP)
#  TARGET.K(r_max=sample.r_max, K, N1, z, start.Yr=start.Yr, num.Yrs=bisection.Yrs, target.Pop=sample.N.obs, catches=catches, MVP=MVP)

TARGET.K = function(r_max, K, N1, z, num.Yrs, start.Yr, target.Pop, catches=catches, MVP=0)
{
  Pred.N=GENERALIZED.LOGISTIC(r_max=r_max, K=K, N1=K, z=z, start.Yr=start.Yr, num.Yrs=num.Yrs, catches=catches, MVP=MVP)
  diff=Pred.N$Pred.N[num.Yrs]-target.Pop
  #plot(seq(1901, 1900+num.Yrs, 1), Pred.N$Pred.N, ylim=c(0,50000))
  #print(paste("N-Target=", Pred.N$Pred.N[num.Yrs], "Target.Pop=", target.Pop, "K=", K))
  return(diff=diff)
}
#END FUNCTION
#---------------------------------


#####################################################################
# THIS FUNCTION DOES THE LOGISTIC BISECTION
#####################################################################
#ex call: LOGISTIC.BISECTION.K(K.low=1, K.high=100000, r_max=r_max, z=z, num.Yrs=bisection.Yrs, start.Yr=start.Yr, target.Pop=target.Pop, catches=catches, MVP=MVP, tol=0.001)
LOGISTIC.BISECTION.K=function(K.low, K.high, r_max, z, num.Yrs, start.Yr, target.Pop, catches, MVP, tol=0.001)
{
  Kmin=uniroot(TARGET.K, tol=tol, c(K.low, K.high), r_max=r_max, z=z, num.Yrs=num.Yrs, start.Yr=start.Yr, target.Pop=target.Pop, catches=catches, MVP=MVP)
  return(Kmin$root)
}
#END OF FUNCTION
#------------------------

######################################################
# THIS FUNCTION COMPUTES THE ANALYTICAL ESTIMATES OF Q
######################################################
CALC.ANALYTIC.Q=function(rel.Abundance, Pred.N, start.Yr, add.CV=0, num.IA)
{
  analytic.Q=rep(NA, num.IA) #vector to store the q values

  for(i in 1:num.IA)
  {
    IA=subset(rel.Abundance, Index==i) #subseting across each index of abundance
    IA.yrs=IA$Year-start.Yr+1 #years for which IAs are available
    IA$Sigma=sqrt(log(1+IA$CV.IA.obs^2)) #computing the value of sigma as in Zerbini et al. 2011
    qNumerator=sum((log(IA$IA.obs/Pred.N[IA.yrs]))/(IA$Sigma*IA$Sigma+add.CV*add.CV)) #numerator of the analytic q estimator (Zerbini et al., 2011 - eq. (3))
    qDenominator=sum(1/(IA$Sigma*IA$Sigma)) #denominator of the analytic q estimator (Zerbini et al., 2011 - eq. (3))
    analytic.Q[i]=exp(qNumerator/qDenominator) #estimate of q
  }
  return(analytic.Q=analytic.Q)
}
# END OF FUNCTION
#-------------------------------------------------------

################################################################################
# THIS FUNCTION COMPUTES THE LN LIKELIHOOD OF THE INDICES OF ABUNDANCE
################################################################################
LNLIKE.IAs=function(Rel.Abundance, Pred.N, start.Yr, q.values, add.CV, num.IA, log=TRUE)
{
  loglike.IA1=0
  loglike.IA2=0

  for(i in 1:num.IA)
  {
    IA=subset(Rel.Abundance, Index==i) #subseting across each index of abundance
    IA.yrs=IA$Year-start.Yr+1 #years for which IAs are available
    loglike.IA1=loglike.IA1+((sum(log(IA$Sigma)+log(IA$IA.obs)+0.5*((((log(q.values[i]*Pred.N[IA.yrs])-log(IA$IA.obs))^2)/(IA$Sigma*IA$Sigma)))))) #this is the likelihood from Zerbini et al. 2011 (eq. 5)

    loglike.IA2=loglike.IA2 + CALC.LNLIKE(Obs.N=IA$IA.obs, Pred.N=(q.values[i]*Pred.N[IA.yrs]), CV= sqrt(IA$Sigma*IA$Sigma + add.CV*add.CV), log=log) #this is the log-normal distribution from R (using function dnorm)

  }
  return(list(loglike.IA1=loglike.IA1, loglike.IA2=loglike.IA2))

}
# END OF FUNCTION
################################################################################

################################################################################
# THIS FUNCTION COMPUTES THE LN LIKELIHOOD OF THE ABSOLUTE ABUNDANCE
################################################################################
LNLIKE.Ns=function(Obs.N, Pred.N, start.Yr, add.CV, log=TRUE)
{
  loglike.Ns1=0
  loglike.Ns2=0

  N.yrs=Obs.N$Year-start.Yr+1 #years for which Ns are available
  loglike.Ns1=loglike.Ns1+((sum(log(Obs.N$Sigma)+log(Obs.N$N.obs)+0.5*((((log(Pred.N[N.yrs])-log(Obs.N$N.obs))^2)/(Obs.N$Sigma*Obs.N$Sigma+add.CV*add.CV)))))) #this is the likelihood from Zerbini et al. 2011 (eq. 4)
  loglike.Ns2=loglike.Ns2 + CALC.LNLIKE(Obs.N=Obs.N$N.obs, Pred.N=(Pred.N[N.yrs]), CV= sqrt(Obs.N$Sigma*Obs.N$Sigma + add.CV*add.CV), log=log) #this is the log-normal distribution from R (using function dnorm)

  return(list(loglike.Ns1=loglike.Ns1, loglike.Ns2=loglike.Ns2))

}

# END OF FUNCTION
################################################################################

################################################################################
# THIS FUNCTION COMPUTES THE LN LIKELIHOOD OF THE GROWTH RATE
################################################################################
LNLIKE.GR=function(Obs.GR, Pred.GR, GR.SD.Obs)
{
  loglike.GR1=0
  loglike.GR2=0

  loglike.GR1=loglike.GR1+(((log(GR.SD.Obs)+0.5*(((Pred.GR-Obs.GR)/GR.SD.Obs)^2)))) #this is the likelihood from Zerbini et al. 2011 (eq. 6)

  loglike.GR2=loglike.GR2 + CALC.LNLIKE(Obs.N=Obs.GR, Pred.N=Pred.GR, CV=GR.SD.Obs, log=FALSE)

  return(list(loglike.GR1=loglike.GR1, loglike.GR2=loglike.GR2))

}
# END OF FUNCTION
################################################################################

################################################################################
# THIS FUNCTION COMPUTES THE LOG LIKELIHOOD
################################################################################
CALC.LNLIKE=function(Obs.N, Pred.N, CV, log=F)
{
  return(sum(dnorm(x=log(Obs.N), mean=log(Pred.N), sd=CV, log=F)))
}
# END OF FUNCTION
################################################################################

################################################################################
# OUTPUT FUNCTION
# FUNCTION TO PROVIDE AN SUMMARY OF THE SIR.HUMPBACK OUTPUT
################################################################################
SUMMARY.SIR=function(x, scenario="USERDEFINED")
  # x corresponds to the output resample table from the HUMPBACK.SIR function
{
  num.col=dim(x)[2]
  col.names=names(x)
  row.names=c("mean", "median", "2.5%PI", "97.5%PI", "5%PI", "95%PI", "min", "max", "n")

  output.summary=matrix(nrow=length(row.names), ncol=num.col)
  output.summary[1, ]=sapply(x, mean)
  output.summary[2, ]=sapply(x, median)
  output.summary[3, ]=sapply(x, quantile, probs=0.025)
  output.summary[4, ]=sapply(x, quantile, probs=0.975)
  output.summary[5, ]=sapply(x, quantile, probs=0.05)
  output.summary[6, ]=sapply(x, quantile, probs=0.95)
  output.summary[7, ]=sapply(x, min)
  output.summary[8, ]=sapply(x, max)
  output.summary[9, ]=sapply(x, length)

  output.table=data.frame(output.summary)
  names(output.table)=col.names
  row.names(output.table)=row.names
  noquote(format(output.table, digits=3, scientific=FALSE))

  return(list(scenario=scenario, date=Sys.time(), output.table=output.table))

}
#################################
# END OF FUNCTION
#################################



