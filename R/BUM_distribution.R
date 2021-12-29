#' Beta-uniform mixture (BUM) distribution
#'
#' Density and distribution function, and random generation for a BUM
#'   distribution. The quantile function cannot be simply expressed and is thus
#'   not available.
#'
#' @name BUM_distribution
#'
#' @param x,q Vector of quantiles.
#' @param lambda Vector of BUM parameter \code{lambda}. *lambda* is the
#'   fraction of null (uniform distribution density) over all density.
#' @param a Vector of BUM parameter \code{a}, which corresponds to the *a*
#'   shape parameter of the beta distribution component. It describes the
#'   steepness of the beta distribution.
#'
#' @return \code{dbum} gives the density (PDF), and \code{pbum} gives the
#'   distribution function (CDF).
#'
#' @details The \code{b} parameter of the typical beta distribution function is
#'   constrained to be \code{1} here to limit it to a one-sided, monotonic case,
#'   peaking at 0, when \code{0 < a < 1}.
#'
#' @examples
#' dbum(x = 0.096, lambda = 0.65, a = 0.1)
#' pbum(q = 0.096, lambda = 0.65, a = 0.1)
#'
#' dbum(x = c(0.013, 0.04, 0.93, 0.8), lambda = 0.73, a = 0.02)
#' pbum(q = c(0.013, 0.04, 0.93, 0.8), lambda = 0.73, a = 0.02)
#'
NULL

#' @rdname BUM_distribution
#' @export
dbum = function(x, lambda, a) { lambda     + (1 - lambda) * a * x^(a-1) }

#' @rdname BUM_distribution
#' @export
pbum = function(q, lambda, a) { lambda * q + (1 - lambda)     * q^a     }

