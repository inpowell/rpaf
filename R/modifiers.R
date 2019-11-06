#' Simple to/from dataframe modification function
#'
#' Returns a function that changes all \code{from} levels of a column in a data
#' frame to the \code{to} level, and returns the modified data frame. If
#' \code{from} is not specified, then the entire column is replaced.
#'
#' @param factor the name of the factor column to be changed
#' @param to the target level of the factor
#' @param from,... the source levels of the factor (defaulting to all levels if
#'   unspecified)
#'
#' @seealso \code{\link{paf}} and \code{\link{paf_ci}}, where this function can
#'   be used
#'
#' @return a function that takes a data frame and outputs a modified data frame
#' @export
ModifyCategorical <- function(factor, to, from = NULL, ...) {
  from <- c(from, ...)
  if (!is.null(from))
    function(data) {
      data[data[[factor]] %in% from, factor] <- to
      data
    }
  else
    function(data) {
      # data[[factor]] <- factor(to, levels = levels(data[[factor]]))
      # the following handles missing data better:
      from <- setdiff(levels(data[[factor]]), to)
      data[data[[factor]] %in% from, factor] <- to
      data
    }
}
