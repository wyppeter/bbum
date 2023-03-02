#' BBUM FDR correction of p values
#'
#' Fits the BBUM model on the dataset, and transforms raw p values into values
#'   corrected for false discovery rate (FDR) of null **and** secondary signal
#'   according to the BBUM model, as multiple testing correction. Optionally, it
#'   automatically detects extreme value outliers among the background set, and
#'   resolves correction issues by trimming the outliers.
#'
#' @param pvals A vector of all numerical p values, including both signal and
#'   background sets.
#' @inheritParams BBUM_fit
#' @param add_starts List of named vectors for additional starts of fitting
#'   algorithm beyond the default set.
#' @param only_start Whether the algorithm should only use the given starts
#'   (\code{add_starts}) to fit.
#' @param auto_outliers Toggle automatic outlier trimming.
#' @param rtrimmax Maximum fraction of data points allowed to be outliers in the
#'   background set of data (to be trimmed).
#' @param atrimmax Maximum absolute number of data points allowed to be outliers
#'   in the background set of data (to be trimmed).
#' @param two_tailed Toggle the "two-tailed" case of BBUM correction, if the
#'   background assumption is weak and bona fide hits in the background class
#'   are relevant. See Details. Default behavior is off.
#' @param quiet Suppress printed messages and warnings.
#'
#' @return A named list with the following items:
#' * \code{pvals}: Vector of input p values.
#' * \code{pBBUMs}: Vector of p values corrected for FDR by BBUM modeling.
#' * \code{estim}: A named list of fitted parameter values.
#' * \code{LL}: Value of the maximized log-likelihood.
#' * \code{convergence}: Convergence code from \code{optim}.
#' * \code{outlier_trim}: Number of outliers trimmed in the background set.
#' * \code{r.passed}: Boolean for whether the fitted \code{r} value was under
#'   the threshold for flagging outliers.
#'
#' @details \code{pBBUM} represents the expected overall FDR level if the cutoff
#'   were set at that particular p value. This is similar to the interpretation
#'   of p values corrected through the typical \code{p.adjust(method = "fdr")}.
#' @details \code{pBBUM} values are designed for the signal set p values only,
#'   Values for the background set are given but not valid as significance
#'   testing adjustment, and so should **not** be used to call any hits. They
#'   are provided primarily to compare the equivalent transformation against the
#'   signal set to assess the adjustment strategy. The background set should
#'   **not** be considered for hits.
#' @details \code{BBUM_corr} functions best with p values filtered for poor
#'   quality data points in prior. Such points tend to have high p values and
#'   may disrupt the uniform null distribution.
#' @details Default starts for BBUM fitting are implemented. If additional
#'   starts should included, or only custom starts should be considered, make
#'   use of \code{add_starts} and/or \code{only_start} arguments.
#' @details If more than one start achieved the identical likelihood, a random
#'   start is chosen among them.
#' @details Automatic outlier detection relies on the model fitting a value of
#'   \code{r > 1}. Such a result suggests that a stronger signal (presumably
#'   outliers) exists in the background set than in the signal set, which
#'   violates the assumptions of the model. This is a *conservative* strategy.
#'   The ideal way to deal with outliers is to identify and handle them before
#'   any statistical analyses. For benchmarking of the trimming strategy, see
#'   Wang & Bartel, 2022.
#' @details Adding too many starts or allowing too much outlier trimming can
#'   increase computation time.
#' @details  If the background assumption is weak, such that a small number
#'   of bona fide hits are anticipated and relevant to the hypothesis at
#'   hand among the data points designated "background class", the FDR could be
#'   made to include the background class. This is akin to a two-tailed test
#'   (despite a one-tailed assumption to begin with). This would allow the
#'   generation of genuine FDR-corrected p values for the background class
#'   points as well. Toggle this using the \code{two_tailed} value.
#' @details Due to the asymptotic behavior of the function when any
#'   p values = 0, any p values < \code{.Machine$double.xmin*10} would be
#'   constrained to \code{.Machine$double.xmin*10}.
#'
#'
#' @examples
#' BBUM_corr(
#'   pvals         = c(0.501, 0.203, 0.109, 0.071, 0.019, 0.031, 0.001,
#'                     0.000021, 0.00010, 0.03910,
#'                     0.0001,
#'                     0.11, 0.27, 0.36, 0.43, 0.50, 0.61, 0.77, 0.87, 0.91,
#'                     0.13, 0.21, 0.38, 0.42, 0.52, 0.60, 0.73, 0.81, 0.97),
#'   signal_set    = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE,
#'                     TRUE,  TRUE,  TRUE,
#'                     FALSE,
#'                     FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
#'                     TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE),
#'   add_starts = list(c(lambda = 0.9, a = 0.6, theta = 0.1, r = 0.1)),
#'   limits = list(a = c(0.1,0.7))
#' )
#'
#' @export
BBUM_corr = function(
  pvals, signal_set,
  add_starts = list(), only_start = FALSE,
  limits = list(),
  pBBUM.alpha = 0.05,
  auto_outliers = TRUE, rthres = 1,
  rtrimmax = 0.05, atrimmax = 10,
  two_tailed = FALSE,
  quiet = FALSE) {

  # Check inputs and initiate fit ----
  if(length(pvals) != length(signal_set)){
    stop("Length of pvals does not match length of signal_set!")
    }

  # Get the p value vectors
  dt_bg_set     = pvals[!signal_set & !is.na(pvals)] %>% sort()
  dt_signal_set = pvals[ signal_set & !is.na(pvals)] %>% sort()

  # Get list of starts to try
  if(!only_start){
    # Use multi-start
    start.tries = c(
      # Concatenate lists of different starts
      BBUM_startgen(limits = limits, pBBUM.alpha = pBBUM.alpha),
      add_starts
    )
  } else if (length(add_starts) > 0){
    # Only start at given starts
    start.tries = add_starts
  } else {
    stop("No starts were given to fit!")
  }
  start.tries.input = start.tries %>%
    lapply(., BBUM_params_toReal, limits = limits, rcap = FALSE, pBBUM.alpha = pBBUM.alpha)
    # turn into logit-transformed values
  start.tries.input.cap = start.tries %>%
    lapply(., BBUM_params_toReal, limits = limits, rcap = TRUE,  pBBUM.alpha = pBBUM.alpha)
    # turn into logit-transformed values

  # Fitting ----

  # First pass with no trimming
  outlier_trim = 0
  best.fit.i = BBUM_fit(
    dt_signal_set = dt_signal_set,
    dt_bg_set = dt_bg_set,
    starts = start.tries.input,
    limits = limits,
    rcap = FALSE,
    pBBUM.alpha = pBBUM.alpha,
    outlier_trim = outlier_trim,
    rthres = rthres
  )

  r.pass.i = unname(best.fit.i$r.pass)

  if(!auto_outliers){

    # Auto trim DISABLED: just check for issues if we are not in quiet mode

    print_additional_outliers = FALSE

    if(!quiet && !r.pass.i) {
      warning(">> r may converge at >= 1: Possible violation of distribution assumptions. r constrained to < 1. Make sure to check fit model; results may not be accurate if r ~ 1.")
    }

    # Re-fit, now cap theta
    best.fit = BBUM_fit(
      dt_signal_set = dt_signal_set,
      dt_bg_set = dt_bg_set,
      starts = start.tries.input.cap,
      limits = limits,
      rcap = TRUE,
      pBBUM.alpha = pBBUM.alpha,
      outlier_trim = outlier_trim,
      rthres = rthres
    )

  } else {

    # Auto trim ENABLED
    # Iterative step-up outlier trimming from bg set

    trim_max = min(floor(length(dt_bg_set) * rtrimmax), atrimmax)
    print_additional_outliers = FALSE

    while(!r.pass.i && outlier_trim < trim_max) {

      # If r does not pass the test, start trimming from the bg set
      print_additional_outliers = TRUE
      outlier_trim = outlier_trim + 1  # step up!

      if(!quiet){ print(paste0(">> Trying to trim outliers x ", outlier_trim, "...")) }

      best.fit.i = BBUM_fit(
        dt_signal_set = dt_signal_set,
        dt_bg_set = dt_bg_set,
        starts = start.tries.input,
        limits = limits,
        rcap = FALSE,
        pBBUM.alpha = pBBUM.alpha,
        outlier_trim = outlier_trim,
        rthres = rthres
      )

      r.pass.i = unname(best.fit.i$r.pass)  # new pass boolean
    }

    # We are done. Either it now passes, or we reached the limit.
    if(!r.pass.i){
      # Fail
      outlier_trim = 0
      print_additional_outliers = FALSE
      if(!quiet) {
        warning(">> Failed to trim outliers to resolve the possible violation of distribution assumptions. r may converge at >= 1. r constrained to < 1. Make sure to check fit model; results may not be accurate if r ~ 1.")
      }
    }

    # Re-fit, now cap theta
    best.fit = BBUM_fit(
      dt_signal_set = dt_signal_set,
      dt_bg_set = dt_bg_set,
      starts = start.tries.input.cap,
      limits = limits,
      rcap = TRUE,
      pBBUM.alpha = pBBUM.alpha,
      outlier_trim = outlier_trim,
      rthres = rthres
    )
  }

  # Load fit results ----
  # Log-likelihood
  LL = unname(best.fit$LL)
  # Coefs
  vals = best.fit$estim

  val.l  = unname(vals$lambda)
  val.a  = unname(vals$a)
  val.th = unname(vals$theta)
  val.r  = unname(vals$r)

  if(!quiet){
    print(">> Fit results:")
    print(paste0(
      "lambda = ",  signif(val.l, 3),
      "; a = ",     signif(val.a, 3),
      "; theta = ", signif(val.th, 3),
      "; r = ",     signif(val.r, 3)
      ))


  }

  # BBUM transform and multiple testing correction (FDR) ----
  #####
  if (two_tailed) {
    dataratio = length(dt_signal_set)/length(dt_bg_set)
    pBBUMs = BBUM_FDR(q = pvals,
                      lambda = val.l, a = val.a, theta = val.th, r = val.r,
                      dtratio = dataratio
    )
  } else {
    pBBUMs = BBUM_FDR(q = pvals,
                      lambda = val.l, a = val.a, theta = val.th, r = val.r
    )
  }

  # Additional outlier notice
  if(!quiet && print_additional_outliers){
    print(paste0(
      ">> Additional outliers trimmed from background class (# = ",
      outlier_trim,
      ")"))
  }

  # Return ----
  list(
    pvals = pvals, pBBUMs = pBBUMs,
    estim = vals, LL = LL,
    convergence = best.fit$convergence,
    outlier_trim = outlier_trim,
    r.passed = r.pass.i
  )
}
