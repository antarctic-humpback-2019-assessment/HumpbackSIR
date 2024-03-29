% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{GENERALIZED_LOGISTIC}
\alias{GENERALIZED_LOGISTIC}
\title{GENERALISED LOGISTIC MODEL}
\usage{
GENERALIZED_LOGISTIC(r_max, K, N1, z, start_yr, num_Yrs, catches, MVP)
}
\arguments{
\item{r_max}{The maximum net recruitment rate ($r_{max}$).}

\item{K}{Pre-exploitation population size in numbers or biomass (depending on
input).}

\item{N1}{The population size in numbers or biomass at year 1 (generally
assumed to be K).}

\item{z}{The parameter that determines the population size where productivity
is maximum (assumed to be 2.39 by the IWC SC).}

\item{start_yr}{The first year of the projection (assumed to be the first
year in the catch series).}

\item{num_Yrs}{The number of projection years. Set as the last year in the
catch or abundance series, whichever is most recent, minus the
\code{start_yr}.}

\item{catches}{The time series of catch in numbers or biomass. Currently does
not handle NAs and zeros will have to input a priori for years in which
there were no catches.}

\item{MVP}{The minimum viable population size in numbers or biomass. Computed
as 4 * \code{\link{num.haplotypes}} to compute minimum viable population
(from Jackson et al., 2006 and IWC, 2007).}
}
\value{
A list of the minimum population size \code{\link{Min.Pop}}, year of the minimum population size \code{\link{Min.Yr}}, a indicator of wether the minimum population size is below the \code{\link{MVP}}, and the predicted population size \code{Pred.N}.
}
\description{
\code{GENERALIZED_LOGISTIC} returns the population projection using a
Pella-Tomlison population dynamics model: $$N_{t+1} =
N_{t}+N_{t}*r_{max}*\left[ 1 - \left( \frac{N_{t}}{K} \right) ^z \right] -
C_{t}$$ where $N$ is the population size at year $t$ or $t+1$, $r_{max}$ is
the maximum net recruitment rate, $K$ is the pre-exploitation population
size, $z$ is the parameter that determines the population size where
productivity is maximum. For example, a value of 2.39 corresponds to maximum
sustainable yield of $0.6K$ and is assumed by the IWC SC, and $C_t$ is the
harvest in numbers in year $t$. Population size can be in either numbers or
biomass, however, units will have to be the same as the units used for catch,
relative abundance, and absolute abundance.
}
\examples{
num_Yrs  <-  10
start_yr  <-  1
r_max  <-  0.2
K  <-  1000
N1  <-  K
catches  <-  round(runif(10, min = 0, max = 150 ), 0 )
MVP  <-  0
GENERALIZED_LOGISTIC(r_max, K, N1, z, start_yr, num_Yrs, catches)
}
