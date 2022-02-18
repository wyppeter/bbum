#' Expected sensitivity and specificity by the BBUM model
#'
#' \code{BBUM_expperf} computes the expected performance in terms of sensitivity
#'   and specificity values at quantiles \code{q} (raw p-values) according to
#'   the given parameters of the BBUM distribution.
#'
#' This is primarily used for plotting in \code{BBUM_plot(option = "confusion")}.
#'
#' @details Sensitivity is the fraction of positives that are correctly called
#'   as positives.
#' @details Specificity is the fraction of negatives that are correctly called
#'   as negatives.
#'
#' @param q Vector of quantiles.
#' @inheritParams BBUM_FDR
#' @inheritParams BBUM_DEcorr
#'
#' @details  If the background assumption is weak, such that a small number
#'   of bona fide hits are anticipated and relevant to the hypothesis at
#'   hand among the data points designated "background class", the FDR could be
#'   made to include the background class. This is akin to a two-tailed test
#'   (despite a one-tailed assumption to begin with). This would allow the
#'   generation of genuine FDR-corrected p values for the background class
#'   points as well. Toggle this using the \code{dtratio} value.
#'
#' @return A named list of two vectors of the same length as \code{q},
#'   under names \code{sensitivity} and \code{specificity}.
#'
#' @examples
#' # Default
#' BBUM_expperf(q = c(0.001, 0.007, 0.19, 0.5, 0.99),
#'   lambda = 0.65, a = 0.1, theta = 0.02, r = 0.07)
#' # "Two-tailed"
#' BBUM_expperf(q = c(0.001, 0.007, 0.19, 0.5, 0.99),
#'   lambda = 0.65, a = 0.1, theta = 0.02, r = 0.07, dtratio = 1.13)
#'
#' @export
BBUM_expperf = function(q,
                        lambda, a, theta, r,
                        dtratio = Inf) {
  frac_signaldata = 1/(1+dtratio^-1)
  sensi = dplyr::if_else(
    q <= 0, 0,
    q^(a*r)*frac_signaldata
  )
  speci = dplyr::if_else(
    q <= 0, 1,
    1 - (lambda * q + (1-lambda) * q^a)*frac_signaldata
  )
  list(
    sensitivity = sensi,
    specificity = speci
  )
}
