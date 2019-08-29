#' Risk function information
#'
#' A standardised object containing information regarding risk
#' functions needed to calculate PAFs.
#'
#' A risk object is a list that contains the following elements:
#'
#' \describe{
#'
#' \item{\code{riskfn(df, ...)}}{A function which takes a population
#' as described in \code{df} and optional further arguments, and
#' outputs the event risk for each member in the population.}
#'
#' \item{\code{dtransvar(df, ...)}}{A function that returns the
#' derivatives of risks with respect to any underlying random
#' parameters. The output should be a matrix with element in row
#' \eqn{i} and column \eqn{j} \deqn{J_{i,j} = \frac{\partial r(x_i;
#' \beta)}{\partial \beta_j}} where \eqn{x_i} corresponds to row
#' \eqn{i} of the input dataframe.}
#'
#' \item{\code{var(d)}}{A function that for any vector \code{d},
#' returns the variance of the dot product \eqn{d^\top \beta}.}
#'
#' \item{\code{source_df}}{(optional) The source population used for
#' risk function calculation.}
#'
#' }
#'
#' @name risk-object
#' @family risk functions
NULL

#' Risk information for binary GLMs
#'
#' This function extracts information on the risk from generalised linear models
#' with a binary predictor. This can be used for logistic regression, or other
#' models with binary dependent variables, such as the probit model.
#'
#' @param model a fitted GLM
#'
#' @return a \code{\link{risk-object}} list
#' @export
#'
#' @importFrom stats delete.response model.frame model.matrix predict.glm terms vcov
#' @family risk functions
risk_glm <- function(model) {
  Terms <- delete.response(terms(model))
  mu.eta <- model$family$mu.eta

  list(
    riskfn = function(df, ...) {
      predict.glm(model, newdata = df, type = "response")
    },
    dtransvar = function(df, ...) {
      X <- model.matrix(Terms, data = df, xlev = model$xlevels, contrasts.arg = model$contrasts)
      X * mu.eta(drop(X %*% coef(model)))
    },
    var = function(d) drop(t(d) %*% vcov(model) %*% d),
    source_df = model.frame(model)
  )
}
