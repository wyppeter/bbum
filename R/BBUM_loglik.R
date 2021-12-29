#' Log-likelihood for the BBUM model
#'
#' \code{BBUM_loglik} computes the total log-likelihood of the given BBUM
#'   parameters, given data points with and without the primary beta component.
#'   This is the maximization target for fitting the BBUM model to data (MLE).
#'
#' @details p values lower than \code{.Machine$double.xmin*10} are constrained
#'   to that limit value before evaluation to avoid zero and machine limit
#'   issues.
#'
#' @param params Named vector of BBUM parameters.
#' @param posSet Vector of values following a BBUM distribution including the
#'  primary beta component ("signal set" or "sample set").
#' @param negSet Vector of values following a BBUM distribution without the
#'  primary beta component, i.e. a BUM distribution with the secondary beta
#'  distribution ("background set").
#' @param limits Named list of custom limits for specific paramters. Parameters
#'   not mentioned would be default values.
#' @param rcap Whether the parameter r should have a stringent upper bound in
#'   this instance (for smart toggling of outlier detection).
#' @param pBBUM.alpha Cutoff level of BBUM-FDR-adjusted p values for
#'   significance testing. Only used here to generate appropriate default
#'   limits.
#'
#' @return The value of the total log-likelihood. The logarithm used is the
#'   natural logarithm.
BBUM_loglik = function(params, posSet, negSet, limits = list(), rcap, pBBUM.alpha) {

  # Remove NAs
  posSet = posSet[!is.na(posSet)]
  negSet = negSet[!is.na(negSet)]

  # Get params
  params.trans = BBUM_params_recoverLin(params,
                                        limits = limits,
                                        rcap = rcap,
                                        pBBUM.alpha = pBBUM.alpha)  # identify any non-default limits defined
  lambda = params.trans["lambda"]
  a      = params.trans["a"]
  theta  = params.trans["theta"]
  r      = params.trans["r"]

  # Resolve ~0 -> NaNs/Infs near machine limit
  posSet[posSet < .Machine$double.xmin*10] = .Machine$double.xmin*10
  negSet[negSet < .Machine$double.xmin*10] = .Machine$double.xmin*10

  # Calculate log-likelihoods
  LL.total = sum(
    log(dbbum(posSet, lambda, a, theta, r))
  ) +
    sum(
      log(dbbum(negSet, lambda, a, 0,     1))
    )

  return(LL.total)
}
