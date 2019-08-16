# HumpbackSIR
This is a package to implement a deterministic generalized logistic model in a Bayesian framework using a Sampling-Importance-Resampling (SIR) algorithm as implemented by McAllister et al. (1994). The population dynamics are modeled as follows:

N_(t+1)= N_t + N_t * r_max * [1 - (N_t / K) ^z] - C_t * SLR_p(t)

where Nt is the estimated population abundance in year t, K is the estimated population carrying capacity, z is the assumed shape parameter corresponding to the percentage of K at which maximum sustainable yield production is achieved, r_max is the estimated maximum population growth rate, C_t is the annual catch, and SLR_p(t) is a correction factor for the period of years that includes year t to account for whales that were struck and lost. The estimable parameters of this model are K, rmax and θ, where θ determines the true catch for the pre-modern era given uncertainty in landings numbers:

C_t = C_(t,min) + θ * (C_(t,max) - C_(t,min) )

where C_(t,min) is the minimum estimate of catch in year t, C_(t,max) is the maximum estimated catch in year t. The parameter is K is not assigned a prior. Rather, abundance was projected using a “backwards” approach, which avoids explicitly defining a prior for K by instead assigning a prior to a recent abundance, Nrecent and back-calculating the abundance trajectory.  Priors are therefore defined for a recent estimate of absolute abundance, rmax, SLRp(t), and θ. Model likelihoods were constructed for the absolute abundance and relative indices of relative abundance data assuming log-normal distributions. Catchability coefficients for the indices of relative abundance are analytically integrated out to produce a marginal likelihoods.

# References
McAllister, M. K., Pikitch, E. K., Punt, A. E., Hilborn, R. 1994 A Bayesian approach to stock assessment and harvest decisions using the sampling/importance resampling algorithm. Canadian Journal of Fisheries and Aquatic Sciences. 12, 2673-2687. 

Zerbini, A. N., Ward, E., Engel, M., Andriolo, A., Kinas, P. G. 2011 A Bayesian assessment of the conservation status of humpback whales (Megaptera novaeangliae) in the western South Atlantic Ocean (Breeding Stock A). J. Cetacean Res. Manage. (special issue 3). 131-144. 