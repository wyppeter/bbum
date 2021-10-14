#' Log-likelihood for the BBUM model
#'
#' \code{BBUM_loglik} computes the total log-likelihood of the given BBUM
#'   parameters, given data points with and without the primary beta component.
#'   This is the maximization target for fitting the BBUM model to data.
#'
#' @param params Named vector of BBUM parameters.
#' @param posSet Vector of values following a BBUM distribution including the
#'  primary beta component.
#' @param negSet Vector of values following a BBUM distribution without the
#'  primary beta component, i.e. a BUM distribution with the secondary beta
#'  distribution.
#' @param limits Named list of custom limits for specific paramters. Parameters
#'   not mentioned would be default values.
#' @param rcap Whether the parameter r should have a stringent upper bound in
#'   this instance.
#'
#' @return The value of the total log-likelihood. The logarithm used is the
#'   natural logarithm.
BBUM_loglik = function(params, posSet, negSet, limits = list(), rcap) {

  # Remove NAs
  posSet = posSet[!is.na(posSet)]
  negSet = negSet[!is.na(negSet)]

  # Get params
  params.trans = BBUM_params_recoverLin(params, limits = limits, rcap = rcap)  # identify any non-default limits defined
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
      log(dbbum(negSet, lambda, a, 0,     0))
    )

  return(LL.total)
}
