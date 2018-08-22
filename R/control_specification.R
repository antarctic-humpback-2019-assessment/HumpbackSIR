##' Set fitting parameters for SIR with reasonable defaults.

##'
##' @title Set parameters for SIR
##'
##' @param K_bisect_lims Bounds for K when performing the bisection method of
##'   Pund and Butterworth (1995). Two-element vector denoting upper and lower
##'   bounds. Default is \code{c(1, 5e5)}.
##' @param K_bisect_tol Tolerance passed to \code{uniroot} when using the
##'   bisection method to find K.
##' @param threshold Threshold for McCallister et al. (1994) SIR. Data specific.
##'   Defaults to 1e-17.
##' @param progress_bar Logical; display progress bar?
##' @param verbose Integer, defaults to 0; See Details for specifics. Note that
##'   more intermediate output will slow down sampling.
##'
##' Use of \code{progress_bar} and \code{verbose} together is not allowed. If
##' both are set, only \code{progress_bar} will be shown. The \code{verbose}
##' argument takes an integer, controlling how much is output. All output uses
##' \code{message}, so works with functions like \code{supressMessages}. These
##' levels are
##'
##' \itemize{
##'   \item \code{0} no output;
##'   \item \code{1} sample count and MVP violations;
##'   \item \code{2} above plus likelihood totals;
##'   \item \code{3} above plus likelihood components;
##'   \item \code{4} above plus sampled parameter values.
##' intermediate
##'
##' }
##'
##' @return A list to be passed to the \code{control} argument of \code{HUMPBACK.SIR}.
sir_control <- function(K_bisect_lims = c(1, 5e5),
                        K_bisect_tol = 1e-4,
                        threshold = 1e-17,
                        progress_bar = FALSE,
                        verbose = 0) {
  if (progress_bar && verbose) {
    warning("Progress bar supercedes verbose output.")
    verbose <- 0

  }
  list(K_bisect_lims = K_bisect_lims,
       K_bisect_tol = K_bisect_tol,
       threshold = threshold,
       progress_bar = progress_bar,
       verbose = verbose)
}
