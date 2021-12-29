#' Logit function
#'
#' \code{logit} computes values transformed by the logit or log-odds function.
#'   \code{inv.logit} is the inverse of \code{logit}. The logarithm base used
#'   can be defined (default is natural logarithm). The left and right bounds
#'   of the logit function can also be defined, with the default domain being
#'   \code{(0, 1)}.
#'
#' @param x A numeric vector.
#' @param base The base used for the logarithm.
#' @param left,right The lower and upper domain boundaries of the corresponding
#'   logit function. Defaults are 0 and 1.
#' @param suppress.nan Whether to output \code{NaN} or \code{-Inf}/\code{Inf}
#'   when \code{x} is out of bounds. If \code{TRUE}, an out-of-bounds warning is
#'   issued while suppressing the \code{NaN} into either \code{Inf}
#'   (\code{x > right}) or \code{-Inf} (\code{x < left}).
#'
#' @return A vector with the same length as \code{x}.
#'
#' @examples
#' logit(0.7)
#' logit(6.2, left = 0, right = 10)
#' logit(-0.07, suppress.nan = FALSE)
#' logit(0.01, base = 10)
#'
#' inv.logit(2.31)
#' inv.logit(5.1, left = -1, right = 1)
#' inv.logit(2, base = 2)
#'
#' @export
logit = function(x, base = exp(1), left = 0, right = 1, suppress.nan = TRUE) {
  # Map (left, right) to (-Inf, Inf)
  # If x is out of bounds, returns -Inf or Inf instead of NaN (suppress.nan == TRUE)
  ifelse(x > right & suppress.nan, {warning("logit() out of bounds"); Inf},
         ifelse(x < left & suppress.nan, {warning("logit() out of bounds"); -Inf},
                log((x-left)/(right-x))
         )
  )
}

#' @rdname logit
#' @export
inv.logit = function(x, base = exp(1), left = 0, right = 1) {
  # Map (-Inf, Inf) to (left, right)
  if(right <= left) { stop("Incorrect bounds in inv.logit() call.") }
  e = exp(x)
  r = ifelse(is.infinite(e), 1, (e/(1+e)))
  r * (right-left) + left
}
