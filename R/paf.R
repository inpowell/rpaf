#' PAF estimates
#'
#' Calculates a point estimate (optionally with confidence interval) for a PAF
#' given information about the risk and a counterfactual distribution. External
#' prevalence estimates can be used in the calculation of the PAF.
#'
#' If \code{df} is missing and the original population (used for risk
#' estimation) can be accessed in \code{risk}, then the PAF will be calculated
#' from the original population.
#'
#' The counterfactual population can be specified in two ways: by modifying the
#' data itself, or by altering the prevalence estimates (if using). These two
#' methods should not be used together.
#'
#' First, the data frame is modified by calling \code{mod_df}, which should be a
#' function that accepts a single input data frame and outputs a similar data
#' frame. The output data frame should have the same column structure and number
#' of rows as the original data frame.
#'
#' The second method involves changing external prevalence estimates. This is
#' done using a transition matrix (\code{mod_prev}) where each column sums to 1.
#' Then, the counterfactual prevalence estimate is \code{new_prev <- mod_prev
#' \%*\% prevalence}. By construction of \code{mod_prev}, we note that
#' \code{sum(new_prev) == sum(old_prev)}.
#'
#' @param risk a \code{\link{risk-object}}
#' @param df a data frame describing the new population
#' @param mod_df a function to modify \code{df}; see Details
#' @param prevalence a vector of external prevalence estimates corresponding to
#'   rows in \code{df}
#' @param var_prev variance-covariance matrix for the prevalence estimates
#' @param mod_prev a transition matrix to modify \code{df}; see Details
#' @param method method to calculate confidence intervals (currently the delta
#'   method alone is supported)
#' @param level the confidence level required, default 0.95
#' @param ... further arguments passed to \code{risk$riskfn} and
#'   \code{risk$dtransvar}
#'
#' @return the PAF estimate as a length-one numeric vector
#' @export
#' @name paf
NULL

#' @rdname paf
paf <- function(risk, df, mod_df = identity, prevalence = 1, mod_prev = NULL, ...) {
  if (mod_df != identity && !is.null(mod_prev))
    warning("The behaviour of this function is uncertain for tandem counterfactual methods.")

  if (missing(df))
    df <- risk$source_df

  # build counterfactual scenarios
  if (is.null(mod.prev)) {
    mprev <- prevalence # should be 1, if the warning above is heeded
    mdf <- mod_df(df)
  } else {
    mprev <- drop(mod_prev %*% prevalence)
    mdf <- df
  }

  # estimate risk and incidence
  r <- risk$riskfn(df, ...)
  mr <- risk$riskfn(mdf, ...)
  inc <- sum(r * prevalence)
  minc <- sum(mr * mprev)

  return(1 - minc / inc)
}

#' @rdname paf
paf_confidence <- function(risk, df, mod_df = identity, prevalence = 1, var_prev = 0, mod_prev = NULL,
                           method = c("delta"), level = 0.95, ...) {
  if (mod_df != identity && !is.null(mod_prev))
    warning("The behaviour of this function is uncertain for tandem counterfactual methods.")

  if (missing(df))
    df <- risk$source_df

  # build counterfactual scenarios
  if (is.null(mod.prev)) {
    mprev <- prevalence # should be 1, if the warning above is heeded
    mdf <- mod_df(df)
    mod.prev <- 1 # hack so later calculations work
  } else {
    mprev <- drop(mod_prev %*% prevalence)
    mdf <- df
  }

  # estimate risk and incidence
  r <- risk$riskfn(df, ...)
  mr <- risk$riskfn(mdf, ...)
  inc <- sum(r * prevalence)
  minc <- sum(mr * mprev)

  # calculate gradients
  der_r <- risk$dtransvar(df, ...)
  der_mr <- risk$dtransvar(df, ...)

  # ... with respect to coefficients (beta) and probabilities (p):
  der_paf_b <- colSums(der_mr * mprev) / minc - colSums(der_r * prev) / inc
  der_paf_p <- t(mod_prev) %*% mr / minc - r / inc

  # transformed PAF estimate and variance
  paf_est <- log(minc) - log(minc)
  paf_var <- drop(der_paf_b %*% risk$vcov %*% t(der_paf_b) +
                    der_paf_p %*% var_prev %*% t(der_paf_p))

  # confidence interval calculation
  a <- (1 - level) / 2
  a <- c(a, 1 - a)

  paf_est <- c(paf_est, paf_est + sqrt(paf_var) * qnorm(a))
  strucutre(paf_est, names = c("PAF", "upr", "lwr"))
}
