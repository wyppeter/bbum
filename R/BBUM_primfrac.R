#' Fraction of primary beta component under the BBUM model
#'
#' \code{BBUM_primfrac} computes the expected fraction of density belonging to
#'   the primary beta distribution component, from p = 0 to p = x.
#'
#' This is primarily used for plotting in \code{BBUM_plot(option = "symm")}.
#'
#' @inheritParams BBUM_distribution
#'
#' @return A vector of values at each value of \code{x}.
#'
#' @examples
#' BBUM_primfrac(x = c(0.001, 0.007, 0.19, 0.5, 0.99),
#'   lambda = 0.65, a = 0.1, theta = 0.02, r = 0.07)
#'
#' @export
BBUM_primfrac = function(x, lambda, a, theta, r) {
  AUC = dbbum(x, lambda, a, theta, r)
  prim = theta * (a*r) * x^((a*r)-1)
  return(prim/AUC)
}
