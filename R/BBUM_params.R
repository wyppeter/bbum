#' Transform fitting values to actual values
#'
#' For internal use of the package.
#'
#' @inheritParams BBUM_loglik
#'
BBUM_params_recoverLin = function(params, limits = list(), rcap, pBBUM.alpha) {
  # Retrieve custom lims
  parlims = utils::modifyList(BBUM_default_lim(rcap = rcap, pBBUM.alpha = pBBUM.alpha), limits)
  return(c(
    lambda = inv.logit( unname(params["lambda"]), left = parlims[["lambda"]][1], right = parlims[["lambda"]][2]),
    a      = inv.logit( unname(params["a"])     , left = parlims[["a"]][1],      right = parlims[["a"]][2]),
    theta  = inv.logit( unname(params["theta"]) , left = parlims[["theta"]][1],  right = parlims[["theta"]][2]),
    r      = inv.logit( unname(params["r"])     , left = parlims[["r"]][1],      right = parlims[["r"]][2])
  ))
}

#' Transform actual values to fitting values
#'
#' For internal use of the package.
#'
#' @inheritParams BBUM_loglik
#'
BBUM_params_toReal     = function(params, limits = list(), rcap, pBBUM.alpha) {
  # Retrieve custom lims
  parlims = utils::modifyList(BBUM_default_lim(rcap = rcap, pBBUM.alpha = pBBUM.alpha), limits)
  return(c(
    lambda = logit( unname(params["lambda"]), left = parlims[["lambda"]][1], right = parlims[["lambda"]][2]),
    a      = logit( unname(params["a"])     , left = parlims[["a"]][1],      right = parlims[["a"]][2]),
    theta  = logit( unname(params["theta"]) , left = parlims[["theta"]][1],  right = parlims[["theta"]][2]),
    r      = logit( unname(params["r"])     , left = parlims[["r"]][1],      right = parlims[["r"]][2])
  ))
}
