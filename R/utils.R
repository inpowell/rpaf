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
