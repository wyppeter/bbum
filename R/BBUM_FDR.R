#' False discovery rate by the BBUM model
#'
#' \code{BBUM_FDR} computes the false discovery rate (FDR) value at quantiles
#'   \code{q} according to the given parameters of the BBUM distribution.
#'
#' @param q Vector of quantiles.
#' @inheritParams BBUM_distribution
#'
#' @return A vector of FDR values at each value of \code{q}.
#'
#' @examples
#' BBUM_FDR(q = c(0.001, 0.007, 0.19, 0.5, 0.99),
#'   lambda = 0.65, a = 0.1, theta = 0.02, r = 0.07)
#'
#' @export
BBUM_FDR = function(q, lambda, a, theta, r) {
  dplyr::if_else(
    q <= 0, 0,
    1 - ( (theta * q^(a*r)) / pbbum(q, lambda, a, theta, r) )
  )
}
