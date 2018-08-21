##' Set parameters for SIR
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
##' @param verbose Integer, defaults to 0; print intermediate values during
##'   sampling? Higher values give more output. Note that printing more
##'   intermediate values will slow down computation.
##'
##' @return A list to be passed to the \code{control} argument of \code{HUMPBACK.SIR}.
sir_control <- function(K_bisect_lims = c(1, 5e5),
                        K_bisect_tol = 1e-4,
                        threshold = 1e-17,
                        progress_bar = FALSE,
                        verbose = 0) {
  if (progress_bar) {
    stop("Progress bar not implemented.")
  }
  if (progress_bar && verbose) {
    warn("Progress bar may not be visible if verbose output used.")
  }
  list(K_bisect_lims = K_bisect_lims,
       K_bisect_tol = K_bisect_tol,
       threshold = threshold,
       progress_bar = progress_bar,
       verbose = verbose)
}
