#' Expected hits count under the BBUM model
#'
#' \code{BBUM_expectedhits} computes the expected number of true hits, according
#'   to the given BBUM distribution and the pBBUM value cutoff (\code{alpha})
#'   selected.
#'
#' @param params Named list of BBUM parameters.
#' @param alpha The cutoff selected for BBUM-FDR-corrected p values.
#' @param n Total number of considered data points. This should only include the
#'   signal set.
#'
#' @details \code{params} should be named as \code{BBUM.l}, \code{BBUM.a},
#'   \code{BBUM.th}, \code{BBUM.a}.
#'
#' @return The numeric value of the expected hits count, as a float.
#'
#' @examples
#' BBUM_expectedhits(
#'   params = list(BBUM.l = 0.65, BBUM.a = 0.1, BBUM.th = 0.02, BBUM.r = 0.07),
#'   alpha = 0.05,
#'   n = 432
#'   )
#'
#' @export
BBUM_expectedhits = function(params, alpha, n) {
  unname((params$BBUM.th * alpha^(params$BBUM.a * params$BBUM.r))*n)
}
