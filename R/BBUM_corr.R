#' BBUM FDR correction of p values
#'
#' Fits the BBUM model on the data, and transforms raw p values into values
#'   corrected for false discovery rate (FDR) of null **and** secondary effects
#'   according to the BBUM model, as multiple testing correction. Optionally, it
#'   automatically detects extreme value outliers among the background set, and
#'   resolves correction issues by trimming the outliers.
#'
#' @param pvals A vector of all numerical p values, including both signal and
#'   background sets.
#' @inheritParams BBUM_fit
#' @param add_starts List of named vectors for additional starts of fitting
#'   algorithm.
#' @param only_start Whether the algorithm should only use the given starts to fit.
#' @param auto_outliers Toggle automatic outlier trimming.
#' @param rtrimmax Maximum fraction of data points allowed to be outliers in the
#'   background set of data (to be trimmed).
#' @param atrimmax Maximum absolute number of data points allowed to be outliers
#'  in the background set of data (to be trimmed).
#' @param quiet Suppress printed messages and warnings.
#'
#' @return A named list with the following items:
#' * \code{pvals}: Vector of input p values.
#' * \code{pBBUMs}: Vector of p values corrected for FDR by BBUM modeling.
#' * \code{BBUM.l}: Fitted \code{lambda} parameter.
#' * \code{BBUM.a}: Fitted \code{a} parameter.
#' * \code{BBUM.th}: Fitted \code{theta} parameter.
#' * \code{BBUM.r}: Fitted \code{r} parameter.
#' * \code{outlier_trim}: Number of outliers trimmed in the background set.
#' * \code{r.passed}: Boolean for whether the fitted \code{r} value was under
#'   the threshold for flagging outliers.
#' * \code{BBUM.LL}: Value of the maximized log-likelihood.
#' * \code{fitdist_obj}: Fit output object of \code{optim}.
#'
#' @details \code{pBBUM} represents the expected overall FDR level if the cutoff
#'   were set at that particular p value. This is similar to the interpretation
#'   of p values corrected through the typical \code{p.adjust(method = "fdr")}.
#' @details Default starts for BBUM fitting are implemented. If additional
#'   starts should included, or only custom starts should be considered, make
#'   use of \code{add_starts} and/or \code{only_start} arguments.
#' @details If more than one start achieved the identical maximum likelihood,
#'   A random start is chosen among them.
#' @details Automatic outlier detection relies on the model fitting a value of
#'   \code{r > 1}. For benchmarking of the trimming strategy, see
#'   Wang & Bartel, 2021.
#' @details Adding too many starts or allowing too much outlier trimming can
#'   increase computation time.
#'
#' @examples
#' BBUM_corr(
#'   pvals        = c(0.701, 0.503, 0.109, 0.0071, 0.019, 0.031, 0.0014,
#'                     0.0003, 0.0001, 0.049, 0.0001,
#'                     0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
#'                     0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
#'   signal_set    = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE,
#'                     TRUE, TRUE, TRUE, FALSE,
#'                     FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
#'                     TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
#'   add_starts = list(c(lambda = 0.9, a = 0.6, theta = 0.1, r = 0.1)),
#'   limits = list(a = c(0.50,0.75))
#' )
#'
#' @export
BBUM_corr = function(
  pvals, signal_set,
  add_starts = list(), only_start = FALSE,
  limits = list(),
  auto_outliers = TRUE, rthres = 1,
  rtrimmax = 0.05, atrimmax = 10,
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
      BBUM_startgen(limits = limits),
      add_starts
    )
  } else if (length(add_starts) > 0){
    # Only start at given starts
    start.tries = add_starts
  } else {
    stop("No starts were given to fit!")
  }
  start.tries.input = start.tries %>%
    lapply(., BBUM_params_toReal, limits = limits, rcap = F)
    # turn into logit-transformed values
  start.tries.input.cap = start.tries %>%
    lapply(., BBUM_params_toReal, limits = limits, rcap = T)
    # turn into logit-transformed values

  # Fitting ----
  # Iterative step-up outlier trimming from bg set
  # First pass: no trimming
  outlier_trim = 0
  best.fit = BBUM_fit(
    dt_signal_set = dt_signal_set,
    dt_bg_set = dt_bg_set,
    starts = start.tries.input,
    limits = limits,
    rcap = F,
    outlier_trim = outlier_trim,
    rthres = rthres
  )
  r.pass.i = unname(best.fit$r.pass)
  # Check for issues if auto trim is disabled
  if(!quiet && !auto_outliers && !r.pass.i){
    warning(">> r >= 1: Probable violation of distribution assumptions. Results are not accurate.")
    # Re-fit, now cap theta
    best.fit = BBUM_fit(
      dt_signal_set = dt_signal_set,
      dt_bg_set = dt_bg_set,
      starts = start.tries.input.cap,
      limits = limits,
      rcap = T,
      outlier_trim = outlier_trim,
      rthres = rthres
    )
  }
  # Auto trim section
  trim_max = min(floor(length(dt_bg_set) * rtrimmax), atrimmax)
  print_additional_outliers = F
  if(auto_outliers && !r.pass.i){
    # Activate trimming
    # If r does not pass the test, start trimming from the bg set
    while(!r.pass.i && outlier_trim <= trim_max) {
      outlier_trim = outlier_trim + 1
      best.fit.i = BBUM_fit(
        dt_signal_set = dt_signal_set,
        dt_bg_set = dt_bg_set,
        starts = start.tries.input,
        limits = limits,
        rcap = F,
        outlier_trim = outlier_trim,
        rthres = rthres
      )
      r.pass.i = unname(best.fit.i$r.pass)
    }

    # If we failed to fix it with the max number of trims, we give up and go back to the first pass
    # Otherwise, adopt the one that worked
    if(r.pass.i && outlier_trim != trim_max){
      # Success
      best.fit = best.fit.i
      print_additional_outliers = T
    } else {
      outlier_trim = 0
      if(!r.pass.i && !quiet){
        # Fail
        warning(">> Failed to trim outliers to resolve the violation of distribution assumptions. r >= 1. Results are not accurate.")
        # Re-fit, now cap theta
        best.fit = BBUM_fit(
          dt_signal_set = dt_signal_set,
          dt_bg_set = dt_bg_set,
          starts = start.tries.input.cap,
          limits = limits,
          rcap = T,
          outlier_trim = outlier_trim,
          rthres = rthres
        )
      }
    }
  }

  # Load fit results ----
  # Log-likelihood
  LL = unname(best.fit$LL)
  # Coefs
  val.l  = unname(best.fit$val.l)
  val.a  = unname(best.fit$val.a)
  val.th = unname(best.fit$val.th)
  val.r  = unname(best.fit$val.r)

  if(!quiet){
    print(">> Fit results:")
    print(paste0(
      "lambda = ",  signif(val.l, 4),
      ", a = ",     signif(val.a, 4),
      ", theta = ", signif(val.th, 4),
      ", r = ",     signif(val.r, 4)
    ))
  }

  # BBUM transform and multiple testing correction (FDR) ----
  pBBUMs = BBUM_FDR(pvals, val.l, val.a, val.th, val.r)

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
    BBUM.l  = val.l,  BBUM.a = val.a,
    BBUM.th = val.th, BBUM.r = val.r,
    outlier_trim = outlier_trim, r.passed = r.pass.i,
    BBUM.LL = LL,
    fitdist_obj = best.fit$fit.obj
  )
}
