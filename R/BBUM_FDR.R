#' False discovery rate by the BBUM model
#'
#' \code{BBUM_FDR} computes the false discovery rate (FDR) value at quantiles
#'   \code{q} according to the given parameters of the BBUM distribution.
#'
#' @details This modifed FDR is defined to include primary signal and exclude
#'   both secondary signal and null.
#'
#' @param q Vector of quantiles.
#' @inheritParams BBUM_distribution
#' @param dtratio If using the "two-tailed" case of BBUM correction, provide
#'   the ratio of number of data points in the signal class over that in the
#'   background class. Leave as \code{Inf} for the "one-tailed" case (the
#'   default). See Details.
#'
#' @details  If the background assumption is weak, such that a small number
#'   of bona fide hits are anticipated and relevant to the hypothesis at
#'   hand among the data points designated "background class", the FDR could be
#'   made to include the background class. This is akin to a two-tailed test
#'   (despite a one-tailed assumption to begin with). This would allow the
#'   generation of genuine FDR-corrected p values for the background class
#'   points as well. Toggle this using the \code{dtratio} value.
#'
#' @return A vector of FDR values at each value of \code{q}.
#'
#' @examples
#' # Default
#' BBUM_FDR(q = c(0.001, 0.007, 0.19, 0.5, 0.99),
#'   lambda = 0.65, a = 0.1, theta = 0.02, r = 0.07)
#' # "Two-tailed"
#' BBUM_FDR(q = c(0.001, 0.007, 0.19, 0.5, 0.99),
#'   lambda = 0.65, a = 0.1, theta = 0.02, r = 0.07, dtratio = 1.13)
#'
#' @export
BBUM_FDR = function(q, lambda, a, theta, r, dtratio = Inf) {
  if (dtratio == Inf) {
    dplyr::if_else(
      q <= 0, 0,
      1 - ( (theta * q^(a*r)) / pbbum(q, lambda, a, theta, r) )
    )
  } else {
    frac_signaldata = 1/(1+dtratio^-1)
    dplyr::if_else(
      q <= 0, 0,
      1 - (
        ( theta * q^(a*r))*frac_signaldata /
          (pbbum(q, lambda, a, theta, r)*(frac_signaldata) +
             pbbum(q, lambda, a, 0, r)*(1-frac_signaldata) )
        )
    )
  }
}
