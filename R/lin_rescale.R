#' Linear rescaling
#'
#' Rescale numeric values from 0 to 1 to a specified interval.
#'
#' @details Values smaller than 0 or larger than 1 are rescaled to values
#'   linearly extrapolated from the specified interval.
#'
#' @param x A numeric vector.
#' @param left,right The lower and upper boundaries for rescaled interval.
#'
#' @return A vector with the same length as \code{x}.
#'
#' @examples
#' lin_rescale(0.4, left = 100, right = 200)
#' lin_rescale(c(-2,0.3,0.6,10), left = -10, right = 10)
#'
#' @export
lin_rescale = function(x, left, right) {
  left + x * (right-left)
}
