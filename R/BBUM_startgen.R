# Guess generator for the fitting algorithm ----

#' Rescaling guess values
#'
#' Transformer of guess values defined in \code{[0, 1]} to whatever the limits
#'   are.
#'
#' @param props Named list of parameter guess as proportions defined in
#'   \code{[0, 1]}.
#' @inheritParams BBUM_loglik
#'
BBUM_limapply = function(props, limits = list(), pBBUM.alpha) {
  # Retrieve custom lims
  parlims = utils::modifyList(BBUM_default_lim(rcap = FALSE, pBBUM.alpha = pBBUM.alpha), limits)
  # Rescale
  # If right > 1, limit to 1
  c(
    lambda = lin_rescale(props[["lambda"]], parlims[["lambda"]][1], min(parlims[["lambda"]][2],1)),
    a      = lin_rescale(props[["a"]],      parlims[["a"]][1],      min(parlims[["a"]][2]     ,1)),
    theta  = lin_rescale(props[["theta"]],  parlims[["theta"]][1],  min(parlims[["theta"]][2] ,1)),
    r      = lin_rescale(props[["r"]],      parlims[["r"]][1],      min(parlims[["r"]][2]     ,1))   # always capped
  )
}

#' Generating guess values
#'
#' Generates properly rescaled start values for BBUM fitting, including 6 static
#'   configurations by default (\code{4*3/2}).
#'
#' @inheritParams BBUM_loglik
#'
BBUM_startgen = function(limits = list(), pBBUM.alpha) {
  c(0.1, 0.9) %>%  # fixed start proportions
    purrr::map(function(x) {  # combinations of relative magnitudes
      list(
        BBUM_limapply(limits = limits, pBBUM.alpha = pBBUM.alpha, props = c(
          lambda = 1-x,
          a      = 1-x,
          theta  = x,
          r      = x)
        ),
        BBUM_limapply(limits = limits, pBBUM.alpha = pBBUM.alpha, props = c(
          lambda = 1-x,
          a      = x,
          theta  = 1-x,
          r      = x)
        ),
        BBUM_limapply(limits = limits, pBBUM.alpha = pBBUM.alpha, props = c(
          lambda = 1-x,
          a      = x,
          theta  = x,
          r      = 1-x)
        )

      )
    }) %>%  # this gives a list of lists of combos
    unlist(recursive = F, use.names = F)  # flatten one level; a list of all combos
}
