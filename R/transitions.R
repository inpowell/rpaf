#' Simple to/from transition matrix
#'
#' Generates a simple transition matrix where members with a particular risk
#' factor level are distributed amongst the other levels, holding all other risk
#' factors constant.
#'
#' Suppose that \code{df} is a data frame of risk factors, with the risk factor
#' of interest (RFOI) in \code{df[[column]]}. We are interested in the
#' counterfactual scenario when those taking RFOI level \code{from} are replaced
#' with RFOI level \code{to}.
#'
#' This function allows specification of multiple \code{to} levels. This can be
#' done by specifying \code{to} as a character vector.
#'
#' By default, the current \code{from} population will be evenly dispersed among
#' the \code{to} population(s). This can be changed by specifying the
#' \code{prop} argument, which specifies the proportions of the \code{from}
#' population that move to the equivalent \code{to} populations. \code{prop}
#' should be a non-negative numeric vector of the same length as \code{to} that
#' sums to at most 1. If it sums to less than one, then the counterfactual
#' scenario will include the remainder in the original \code{from} population.
#'
#' @section Warning:
#'
#'   If \code{df} contains any extraneous columns, then the dispersion among
#'   groups may not be accurate and this function has undefined behaviour.
#'
#' @param df a data frame of risk factors; see Warning
#' @param column the name of the risk factor of interest
#' @param from the source population level in \code{column}
#' @param to the target population level(s) in \code{column}
#' @param prop counterfactual dispersion proportions; see Details
#'
#' @return a square transition matrix
#' @export
#'
#' @family transition functions
TransitionSimple <- function(df, column, from, to, prop = NULL) {
  if (any(duplicated(df)))
    stop("Input data frame cannot have duplicate rows.")

  if (is.null(prop))
    prop <- rep(1/length(to), length(to))

  remain <- 1 - sum(prop)
  if (remain < 0)
    stop("Proportions add up to more than 1")

  n <- nrow(df)
  i <- seq_len(n)
  f <- df[, names(df) != column]

  # these two splits should have the same order, we hope
  i_spl <- split(i, f)
  col_spl <- split(df[[column]], f)

  submats <- lapply(col_spl, function(col) {
    col <- as.character(col)
    # from index
    fi <- which(col == from)
    # to indices
    ti <- match(to, col, nomatch = 0)

    # transition submatrix
    A <- diag(nrow = length(col))
    A[fi, fi] <- remain
    A[ti, fi] <- prop
    A
  })

  # final transition matrix
  tmat <- matrix(0, nrow = n, ncol = n)
  for (K in seq_along(i_spl)) {
    ii <- i_spl[[K]]
    A <- submats[[K]]
    tmat[ii, ii] <- A
  }
  tmat
}
