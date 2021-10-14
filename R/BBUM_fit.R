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
#' * \code{LL}: Value of the maximized log-likelihood.
#' * \code{val.l}: Fitted \code{lambda} parameter.
#' * \code{val.a}: Fitted \code{a} parameter.
#' * \code{val.th}: Fitted \code{theta} parameter.
#' * \code{val.r}: Fitted \code{r} parameter.
#' * \code{fit.obj}: Fit output object of \code{optim}.
#'
#' @examples
#' BBUM_fit(
#'   dt_signal_set = c(0.031, 0.0014, 0.0003, 0.0001, 0.049,
#'                     0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
#'   dt_bg_set     = c(0.701, 0.503, 0.109, 0.0071, 0.019,
#'                     0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
#'   starts = list(c(lambda = 0.9, a = 0.1, theta = 0.1, r = 0.1))
#' )
#' BBUM_fit(
#'   dt_all        = c(0.701, 0.503, 0.109, 0.0071, 0.019, 0.031, 0.0014,
#'                     0.0003, 0.0001, 0.049, 0.0001,
#'                     0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
#'                     0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
#'   signal_set    = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE,
#'                     TRUE, TRUE, TRUE, FALSE,
#'                     FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
#'                     TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
#'   starts = list(c(lambda = 0.9, a = 0.1, theta = 0.1, r = 0.1)),
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
        LL.i = BBUM_corr.i$value
        params.estim.i = BBUM_params_recoverLin(BBUM_corr.i$par, limits = limits, rcap = rcap)
        val.l.i  = params.estim.i["lambda"]
        val.a.i  = params.estim.i["a"]
        val.th.i = params.estim.i["theta"]
        val.r.i  = params.estim.i["r"]
        list(LL  = LL.i,
             val.l  = val.l.i,
             val.a  = val.a.i,
             val.th = val.th.i,
             val.r  = val.r.i,
             fit.obj = I(BBUM_corr.i))
      },
      error = function(m){
        warning(m)
        list(LL  = -Inf,
             val.l  = NA_real_,
             val.a  = NA_real_,
             val.th = NA_real_,
             val.r  = NA_real_,
             fit.obj = NA)
      },
      warning = function(m){
        warning(m)
        list(LL  = -Inf,
             val.l  = NA_real_,
             val.a  = NA_real_,
             val.th = NA_real_,
             val.r  = NA_real_,
             fit.obj = NA)
      }
      )
    }) %>%
    data.table::rbindlist(.)

  # Select best fit ----
  # If more than one identical LL, take one randomly
  best.fit = all.fit.tries %>%
    dplyr::filter(LL == max(LL)) %>%
    dplyr::sample_n(1)

  # Load results and r check ----
  best.fit %>%
    dplyr::mutate(outlier_trim = outlier_trim,
                  r.pass = val.r < rthres & val.l < lambda.thres)
}
