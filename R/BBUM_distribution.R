#' Bi-beta-uniform mixture (BBUM) distribution
#'
#' Density, distribution function, and random generation for the BBUM
#'   distribution, with four parameters, \code{lambda}, \code{a}, \code{theta},
#'   and \code{r}. The quantile function cannot be simply expressed and is thus
#'   not available.
#'
#' @name BBUM_distribution
#'
#' @param x,q Vector of quantiles.
#' @param n Number of observations.
#' @param lambda Vector of BBUM parameter \code{lambda}. *lambda* is the
#'   fraction of null (uniform distribution density) over all density except the
#'   primary beta distribution, i.e. null plus secondary beta.
#' @param a Vector of BBUM parameter \code{a}, which corresponds to the *a*
#'   shape parameter of the secondary beta distribution component. It describes
#'   the steepness of the second beta distribution.
#' @param theta Vector of BBUM parameter \code{theta}. *theta* is the fraction
#'   of primary beta distribution density over all density.
#' @param r Vector of BBUM parameter \code{r}, which is the ratio of the *a*
#'   shape parameter of the primary beta distribution over that of the secondary
#'   beta distribution. In other words, \code{a*r} is the shape parameter
#'   for the primary beta distribution component.
#'
#' @return \code{dbbum} gives the density (PDF), and \code{pbbum} gives the
#'   distribution function (CDF). \code{rbbum} generates random values from the
#'   distribution, and \code{rbbum.ID} generates random values while identifying
#'   the part of BBUM they belong to, as a \code{data.frame} of two columns with
#'   length \code{n}. \code{pvalue} corresponds to the generated values, while
#'   \code{cate} is \code{0} for null, \code{1} for primary beta, and \code{2}
#'   for secondary beta.
#'
#' @details The \code{b} parameter of the typical beta distribution function is
#'   constrained to be \code{1} here to limit it to a one-sided, monotonic case,
#'   peaking at 0, when \code{0 < a < 1}.
#'
#' @examples
#' dbbum(x = 0.096, lambda = 0.65, a = 0.1, theta = 0.02, r = 0.07)
#' pbbum(q = 0.096, lambda = 0.65, a = 0.1, theta = 0.02, r = 0.07)
#'
#' dbbum(x = c(0.013, 0.04, 0.93, 0.8), lambda = 0.73, a = 0.02, theta = 0.11, r = 0.003)
#' pbbum(q = c(0.013, 0.04, 0.93, 0.8), lambda = 0.73, a = 0.02, theta = 0.11, r = 0.003)
#'
#' rbbum(n = 100, lambda = 0.65, a = 0.1, theta = 0.02, r = 0.07)
#'
#' rbbum.ID(n = 100, lambda = 0.65, a = 0.1, theta = 0.02, r = 0.07)
#'
NULL

#' @rdname BBUM_distribution
#' @export
dbbum = function(x, lambda, a, theta, r) {
  theta                    * (a*r) * x^((a*r)-1) +  # primary
    (1-theta) * lambda     * 1                   +  # null
    (1-theta) * (1-lambda) * a * x^(a-1)            # secondary
}

#' @rdname BBUM_distribution
#' @export
pbbum = function(q, lambda, a, theta, r) {
  theta                    * q^(a*r) +  # primary
    (1-theta) * lambda     * q       +  # null
    (1-theta) * (1-lambda) * q^a        # secondary
}

#' @rdname BBUM_distribution
#' @export
rbbum = function(n, lambda, a, theta, r) {
  n_prim = stats::rbinom(1, size = n,        prob = theta)
  n_null = stats::rbinom(1, size = n-n_prim, prob = lambda)
  n_seco = n - n_prim - n_null

  pts_null =         stats::runif(n_null, min = 0, max = 1)
  pts_seco = qbeta_a(stats::runif(n_seco, min = 0, max = 1), a)
  pts_prim = qbeta_a(stats::runif(n_prim, min = 0, max = 1), a*r)

  sample(c(pts_null, pts_seco, pts_prim))
}

#' @rdname BBUM_distribution
#' @export
rbbum.ID = function(n, lambda, a, theta, r) {
  n_prim = stats::rbinom(1, size = n,        prob = theta)
  n_null = stats::rbinom(1, size = n-n_prim, prob = lambda)
  n_seco = n - n_prim - n_null

  pts_null =         stats::runif(n_null, min = 0, max = 1)
  pts_seco = qbeta_a(stats::runif(n_seco, min = 0, max = 1), a)
  pts_prim = qbeta_a(stats::runif(n_prim, min = 0, max = 1), a*r)

  ordered.p = data.frame(
    pvalue = c(    pts_null,     pts_seco,     pts_prim),
    cate   = c(rep(0,n_null),rep(2,n_seco),rep(1,n_prim))
  )

  shuffled_order = sample(1:n)
  ordered.p[shuffled_order, ]
}

#' Quantile function of beta distribution with b = 1
#'
#' Quantile function for the special case of the beta distribution, where b = 1.
#'
#' @param p Vector of probabilities.
#' @param a Vector of a / shape1.
#'
#' @return The quantile function.
#'
#' @examples
#' qbeta_a(p = 0.013, a = 0.1)
#'
#' @export
qbeta_a = function(p, a) { p^(1/a) }
