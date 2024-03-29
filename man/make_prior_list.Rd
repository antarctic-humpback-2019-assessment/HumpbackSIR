% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior_specification.R
\name{make_prior_list}
\alias{make_prior_list}
\title{Make a list of priors to be passed to the SIR function.}
\usage{
make_prior_list(r_max = make_prior(runif, 0, 0.118), K = make_prior(use
  = FALSE), N_obs = make_prior(runif, 500, 40000),
  add_CV = make_prior(use = FALSE),
  premodern_catch_sample = make_prior(runif, 0, 1),
  z = make_prior(2.39), q_IA = make_prior(use = FALSE),
  q_count = make_prior(use = FALSE))
}
\arguments{
\item{r_max}{Population growth rate; defaults to Uniform(0, 0.106).}

\item{K}{Carrying capacity, defaults to unused, and solution found using
recent observation and sampled \code{r_max}.}

\item{N_obs}{Prior on a recent abundance estimate. Defaults to Uniform(500,
20,000).}

\item{add_CV}{Defaults to unused. Additional variability.}

\item{premodern_catch_sample}{Defaults to unused. Samples between the minimum and maximum premodern catch values.}

\item{z}{Defaults to constant 2.39. Shape parameter for generalized logistic
population dynamics function.}

\item{q_IA}{Defaults to unused. Prior on q for indices of abundance. If
\code{use = FALSE}, an analytic solution for q is used.}

\item{q_count}{Defaults to unused. Prior for q on counts.}
}
\value{
A named list containing each of the specified priors in a form that
  can be used by the SIR function.
}
\description{
Make a list of priors to be passed to the SIR function.
}
