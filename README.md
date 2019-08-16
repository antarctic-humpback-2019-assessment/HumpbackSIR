# HumpbackSIR
This is a package to implement a deterministic generalized logistic model in a Bayesian framework using a Sampling-Importance-Resampling (SIR) algorithm as implemented by McAllister et al. (1994). The population dynamics are modeled as follows:

$N_{t+1}=N_t+N_t*r_{max}  \left[1 - \left(N_t/K \right)^z \right] - C_t *SLR_{p(t)}$

where Nt is the estimated population abundance in year t, K is the estimated population carrying capacity, z is the assumed shape parameter corresponding to the percentage of K at which maximum sustainable yield production is achieved, rmax is the estimated maximum population growth rate, Ct is the annual catch, and SLRp(t) is a correction factor for the period of years that includes year t to account for whales that were struck and lost. The population was assumed to be at equilibrium (carrying capacity) in 1830, prior to the onset of historical whaling.
The estimable parameters of this model are K, rmax and , where  determines the true catch for the pre-modern era given uncertainty in landings numbers, i.e.
C_t= C_(t,min) +θ*(C_(t,max)-C_(t,min) )     (eq. 2)

where Ct,min is the minimum estimate of catch in year t, Ct,max is the maximum estimated catch in year t. The parameter is K is not assigned a prior. Rather, abundance was projected using a “backwards” approach [71], which avoids explicitly defining a prior for K by instead assigning a prior to a recent abundance, Nrecent and back-calculating the abundance trajectory. 
The based priors for the parameters of the model wereare defined below. :
Priors for rmax : were bounded at 11.8%/year to prevent unrealistic rates of population growth [72].
Ncurrent:???
θ ~ U[0,1]
SLRP(t): ??? 
In the population model, abundance was projected using a “backwards” approach [71], which avoids explicitly defining a prior for K by assigning a prior to a recent abundance and back-calculating the abundance trajectory. To account for uncertainty in pre-modern catches, catch was estimated between the minimum and maximum estimated catches for each year as:

C_t= C_(t,min)+θ*(C_(t,max)-C_(t,min) )     (eq. 2)

where Ct,min is the minimum estimate of catch in year t, Ct,max is the maximum estimated catch in year t, and θ is a parameter uniformly distributed between 0 and 1.
In the population model, abundance was projected using a “backwards” approach [71], which avoids explicitly defining a prior for K by assigning a prior to a recent abundance and back-calculating the abundance trajectory. Priors are therefore defined for a recent estimate of absolute abundance, rmax, SLRp(t), and θ. Priors for rmax were bounded at 11.8%/year to prevent unrealistic rates of population growth [72]. Model lLikelihoods were constructed for the absolute abundance and relative indices of relative abundance data assuming log-normal distributions (eq. [4] and [5], p.134 in [31]). The catchability coefficients for the iIndices of relative abundance were analytically integrated out to produce a marginal likelihoods, related to population abundance by analytically estimating a catchability coefficient q assuming a U[-, ] uniform prior on the log-catchabilityscale for each index (eq. [3], p. 134 in [31], [73]). 

