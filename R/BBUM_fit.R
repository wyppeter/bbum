#' BBUM statistical modeling
#'
#' Fitting the BBUM model on data containing a set with primary effects (signal
#'   set) and a set without (background set), using maximum likelihood estimation
#'   (MLE) with the Nelder-Mead algorithm. It chooses the best solution among
#'   all starts provided.
#'
#' @details Either use \code{dt_signal_set} and \code{dt_bg_set} to input data
#'   separately, or use \code{dt_all} and \code{signal_set} to input data.
#'   When both pairs are defined, \code{dt_all} and \code{signal_set} are
#'   ignored.
#' @details If more than one start achieved the identical maximum likelihood,
#'   A random start is chosen among them.
#' @details \code{rcap} is used internally to decide on the default limits for
#'   \code{r}.
#' @details A failed \code{r.pass} code is not
#'   triggered if \code{lambda} is too big for a reliable fitting of \code{a} to
#'   begin with.
#'
#' @param dt_signal_set,dt_bg_set Vectors of numerical p values, belonging to
#'   the signal set and background set respectively.
#' @param dt_all A vector of all numerical p values, including both signal and
#'   background sets.
#' @param signal_set A vector of booleans signifying which values among
#'   \code{dt_all} are signal set data points. Used only in conjunction with
#'   \code{dt_all}.
#' @param starts A list of named vectors of starts for the four BBUM parameters.
#' @inheritParams BBUM_loglik
#' @param outlier_trim Number of strongest points among the background class to be
#'   trimmed as outliers.
#' @param rthres Threshold value of \code{r} parameter to trigger a failed
#'   \code{r.pass} value.
#'
#' @return A named list with the following items:
#' * \code{estim}: A named list of fitted parameter values.
#' * \code{LL}: Value of the maximized log-likelihood.
#' * \code{convergence}: Convergence code from \code{optim}.
#' * \code{outlier_trim}: The input value of the \code{outlier_trim} argument.
#' * \code{r.passed}: Boolean for whether the fitted \code{r} value was under
#'   the threshold for flagging outliers.
#'
#' @examples
#' BBUM_fit(
#'   dt_signal_set = c(0.000021, 0.00010, 0.03910, 0.031, 0.001,
#'                     0.13, 0.21, 0.38, 0.42, 0.52, 0.60, 0.73, 0.81, 0.97),
#'   dt_bg_set     = c(0.501, 0.203, 0.109, 0.071, 0.019,
#'                     0.11, 0.27, 0.36, 0.43, 0.50, 0.61, 0.77, 0.87, 0.91),
#'   starts = list(c(lambda = 0.9, a = 0.6, theta = 0.1, r = 0.1))
#' )
#' BBUM_fit(
#'   dt_all        = c(0.501, 0.203, 0.109, 0.071, 0.019, 0.031, 0.001,
#'                     0.000021, 0.00010, 0.03910,
#'                     0.0001,
#'                     0.11, 0.27, 0.36, 0.43, 0.50, 0.61, 0.77, 0.87, 0.91,
#'                     0.13, 0.21, 0.38, 0.42, 0.52, 0.60, 0.73, 0.81, 0.97),
#'   signal_set    = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE,
#'                     TRUE,  TRUE,  TRUE,
#'                     FALSE,
#'                     FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
#'                     TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE),
#'   starts = list(c(lambda = 0.9, a = 0.6, theta = 0.1, r = 0.1)),
#'   outlier_trim = 1
#' )
#'
#' @export
BBUM_fit = function(
  dt_signal_set = NULL,
  dt_bg_set     = NULL,
  dt_all        = NULL,
  signal_set    = NULL,
  starts,
  limits = list(),
  rcap = TRUE,
  outlier_trim = 0,
  rthres = 1
) {

  # Get data input ----
  if(is.null(dt_signal_set) && is.null(dt_bg_set)){
    if(is.null(dt_all) && is.null(signal_set)){
      stop("No inputs for BBUM_fit!")
    } else if(length(dt_all) != length(signal_set)) {
      stop("Different lengths for dt_all and signal_set!")
    } else{
      dt_signal_set = dt_all[ signal_set & !is.na(dt_all)] %>% sort()
      dt_bg_set     = dt_all[!signal_set & !is.na(dt_all)] %>% sort()
    }
  }
  # Check input
  if(length(dt_bg_set) < 10 || length(dt_signal_set) < 10){
    stop("Too few data points for reliable model fitting!")
  }

  # If lambda is too big, r may be unreliable (as a is unreliable)
  # In this case, no need to check for r being over 1.
  lambda.thres = 1 - 0.1/(length(dt_signal_set)+length(dt_bg_set))

  # For every start, conduct fitting ----
  all.fit.tries = starts %>%
    lapply(., function(start.i) {
      tryCatch({

        # Optimization
        BBUM_corr.i = stats::optim(
          par = start.i,
          fn = BBUM_loglik,
          posSet = dt_signal_set,
          negSet = dt_bg_set[(outlier_trim+1):length(dt_bg_set)],  # trim any outliers
          limits = limits,
          rcap = rcap,
          method = "Nelder-Mead",
          control = list(maxit = 500,
                         fnscale = -1)  # maximization
        )

        # Extract values and calculate estimated confidence intervals
        par.raw.i = BBUM_corr.i$par

        params.estim.i = BBUM_params_recoverLin(
          par.raw.i,
          limits = limits, rcap = rcap)

        list(LL     = BBUM_corr.i$value,
             conv   = BBUM_corr.i$convergence,
             val.l  = unname(params.estim.i["lambda"]),
             val.a  = unname(params.estim.i["a"]),
             val.th = unname(params.estim.i["theta"]),
             val.r  = unname(params.estim.i["r"])
             )
      },
      error = function(m){
        warning(m)
        list(LL  = -Inf,
             conv = 999,
             val.l  = NA_real_,
             val.a  = NA_real_,
             val.th = NA_real_,
             val.r  = NA_real_
             )
      },
      warning = function(m){
        warning(m)
        list(LL  = -Inf,
             conv = 999,
             val.l  = NA_real_,
             val.a  = NA_real_,
             val.th = NA_real_,
             val.r  = NA_real_
             )
      }
      )
    }) %>%
    data.table::rbindlist(.)  # arrange all results into a table

  # Select best fit ----
  # If more than one identical LL, take one randomly
  succ.fits = all.fit.tries %>%
    dplyr::filter(conv == 0)  # successful convergence
  if(length(succ.fits) == 0){

    # Nothing that worked
    stop("All attempts to fit failed!")

  } else {

    # Load results and r check ----
    best.fit = succ.fits %>%
      dplyr::filter(LL == max(LL)) %>%
      dplyr::sample_n(1) %>%
      dplyr::mutate(outlier_trim = outlier_trim,
                    r.pass = val.r < rthres &
                      val.l < lambda.thres)

    list(
      estim = list( lambda = best.fit$val.l,
                    a      = best.fit$val.a,
                    theta  = best.fit$val.th,
                    r      = best.fit$val.r    ),
      LL = best.fit$LL,
      convergence = best.fit$conv,
      outlier_trim = outlier_trim,
      r.pass = best.fit$r.pass
    )
  }
}
