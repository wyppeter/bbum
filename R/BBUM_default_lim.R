#' Default static limits for BBUM fitting
#'
#' For internal use of the package.
#'
#' @inheritParams BBUM_loglik
#'
BBUM_default_lim = function(rcap, pBBUM.alpha) {
  if(rcap){
    list(
      lambda = c( 0, 1 ),
      a      = c( 0, 1 ),
      theta  = c( 0, 1 - pBBUM.alpha * 2 ),
      r      = c( 0, 1 )
    )
  } else {
    list(
      lambda = c( 0, 1 ),
      a      = c( 0, 1 ),
      theta  = c( 0, 1 - pBBUM.alpha * 2 ),
      r      = c( 0, 10)
    )
  }
}
