#' Variance matrix of multinomial distribution
#'
#' Estimates the variance-covariance matrix of a multinomial distribution.
#'
#' @param p A vector of prevalences which sums to 1, or an integer vector of
#'   counts if \code{n} is missing.
#' @param n The size of the sample used to find \code{p}.
#'
#' @return a square covariance matrix
#' @export
VarMultinom <- function(p, n) {
  if (missing(n)) {
    # treat p as counts
    n <- sum(p)
    p <- p / n
  }
  (diag(p) - p %o% p) / n
}
