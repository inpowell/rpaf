#' Value matching in data frames
#'
#' Find row indices of data frame that match a reference
#'
#' @param x the data frame to be matched
#' @param table reference rows
#' @param by columns of interest in each data frame
#' @param nomatch,incomparables see \code{\link{match}}
#'
#' @return an integer vector of length \code{nrow(x)}
match.data.frame <- function(x, table, by = intersect(names(x), names(table)),
                             nomatch = NA_integer_, incomparables = NULL) {
  # select relevant columns only -- also ensures same order
  x <- subset(x, select = by)
  table <- subset(table, select = by)

  # build interactions
  x <- interaction(x)
  table <- interaction(table)

  # assertion to test that ordering is the same
  stopifnot(all.equal.character(levels(x), levels(table)))

  # final match
  match(as.numeric(x), as.numeric(table), nomatch, incomparables)
}

#' Convert normal-based confidence intervals into standard errors
#'
#' This function returns the standard errors that would give a
#' normally-based confidence interval of given width with the
#' specified level.
#'
#' @param ci a 2-row matrix of confidence intervals, with lower limit
#'   in the first row and upper limit in the second row
#' @param trans a function that transforms the confidence intervals
#'   into normal space.
#' @param level the level of the confidence intervals supplied
#'
#' @return a vector of standard errors
#' @export
#'
#' @examples inst/examples/eg_risk_relative.R
ci2se <- function(ci, trans = identity, level = 0.95) {
  trans <- match.fun(trans)
  ci <- trans(ci)
  dist <- apply(ci, 2, diff)
  a <- (1 - level) / 2
  return(dist / (2 * qnorm(a, lower.tail = FALSE)))
}

#' Rearrange vector to matrix in order by factors
#'
#' This has undefined behaviour when
#' \code{anyDuplicated(interaction(rowf,colf))}, which is not tested.
#'
#' @param data vector of matrix entries
#' @param rowf,colf factors for row and column respectively of the
#'   same length as data
#'
#' @return a matrix
vec2mat <- function(data, rowf, colf) {
  rowf <- as.factor(rowf)
  colf <- as.factor(colf)

  arr <- array(dim = c(nlevels(rowf), nlevels(colf)),
               dimnames = list(levels(rowf), levels(colf)))
  arr[] <- data[order(colf, rowf)]
  arr
}

#' Convert a matrix into a 3d array in order by factors
#'
#' This has undefined behaviour when
#' \code{anyDuplicated(interaction(rowf,colf))}, which is not tested.
#'
#' We want \code{rowf} for dim 1, \code{colf} for dim 2, and the
#' columns of \code{data} for dim 3.
#'
#' @param data the matrix
#' @param rowf,colf factors of length \code{nrow(data)}
#'
#' @return
mat2arr <- function(data, rowf, colf) {
  rowf <- as.factor(rowf)
  colf <- as.factor(colf)

  nlay <- ncol(data)
  nn <- colnames(data)

  arr <- array(dim = c(nlevels(rowf), nlevels(colf), nlay),
               dimnames = list(levels(rowf), levels(colf), nn))

  rowf <- rep(rowf, nlay)
  colf <- rep(colf, nlay)

  arr[] <- data[order(col(data), colf, rowf)]
  arr
}
