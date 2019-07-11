#' PAF estimates in a specific population
#'
#' Calculates a point estimate (\code{paf}) or a point estimate with confidence
#' interval (\code{paf_ci}) for a PAF given an explicit population, information
#' about the risk and a counterfactual distribution.
#'
#' If \code{df} is missing and the original population (used for risk
#' estimation) can be accessed in \code{risk}, then the PAF will be calculated
#' from the original population.
#'
#' In this function, the counterfactual population is defined by modifying
#' entries in the population data frame. This is done through the user-specified
#' function \code{mod_df}, which should accept a single input data frame and
#' outputs a similar data frame. The output data frame should have the same
#' column structure and number of rows as the original data frame.
#'
#' @seealso \code{\link{paf_ext}} for PAF calculations with external prevalences
#'
#' @param risk a \code{\link{risk-object}}
#' @param df a data frame describing the new population
#' @param mod_df a function to modify \code{df}; see Details
#' @param method method to calculate confidence intervals (currently the delta
#'   method alone is supported)
#' @param level the confidence level required, default 0.95
#' @param ... further arguments passed to \code{risk$riskfn} and
#'   \code{risk$dtransvar}
#'
#' @return the PAF estimate as a length-one numeric vector
#' @name paf
NULL

#' @rdname paf
#' @export
paf <- function(risk, df, mod_df = identity, ...) {
  if (missing(df))
    df <- risk$source_df

  # build counterfactual scenario
  mdf <- mod_df(df)

  # estimate risk and incidence
  r <- risk$riskfn(df, ...)
  r_ <- risk$riskfn(mdf, ...)
  I <- sum(r)
  I_ <- sum(r_)

  return(1 - I_ / I)
}

#' @rdname paf
#' @export
paf_ci <- function(risk, df, mod_df = identity, ...,
                           method = c("delta"), level = 0.95) {
  if (missing(df))
    df <- risk$source_df

  # build counterfactual scenarios
  df_ <- mod_df(df)

  # estimate risk and incidence
  r <- risk$riskfn(df, ...)
  r_ <- risk$riskfn(df_, ...)
  I <- sum(r)
  I_ <- sum(r_)

  # calculate gradients
  der_r <- risk$dtransvar(df, ...)
  der_r_ <- risk$dtransvar(df_, ...)

  # ... with respect to coefficients:
  der_paf_b <- colSums(der_r_) / I_ - colSums(der_r) / I

  # transformed PAF estimate and variance
  PAF <- log(I_) - log(I)
  PAFs2 <- drop(t(der_paf_b) %*% risk$vcov %*% der_paf_b)

  # confidence interval calculation
  a <- (1 - level) / 2
  a <- c(1 - a, a) # reverse order because -expm1 is a negative transformation

  PAF <- c(PAF, PAF + sqrt(PAFs2) * qnorm(a))
  PAF <- -expm1(PAF)

  structure(PAF, names = c("PAF", "lwr", "upr"))
}

#' PAF estimates with external prevalences
#'
#' Calculates a point estimate (\code{paf_ext}) or a point estimate with
#' confidence interval (\code{paf_ext_ci}) for a PAF given an population and
#' prevalence estimates, information about the risk and a counterfactual
#' distribution.
#'
#' If \code{df} is missing and the original population (used for risk
#' estimation) can be accessed in \code{risk}, then the PAF will be calculated
#' from the original population.
#'
#' In this function, the counterfactual population is defined using a transition
#' matrix \eqn{A} = \code{mod_prev} such that, if \eqn{p} = \code{prevalence} is
#' a vector of prevalences, then the counterfactual prevalence is \eqn{Ap}. Note
#' that the columns of \code{mod_prev} should all sum to 1, and that
#' \code{prevalence} is reweighted to sum to one.
#'
#' @seealso \code{\link{paf}} for PAF calculations with explicit populations;
#'   \code{\link{TransitionSimple}} and related functions to generate transition
#'   matrices.
#'
#' @inheritParams paf
#' @param prevalence a vector of external prevalence estimates corresponding to
#'   rows in \code{df}
#' @param var_prev variance-covariance matrix for the prevalence estimates
#' @param mod_prev a transition matrix to modify \code{df}; see Details
#'
#' @name paf_ext
NULL

#' @rdname paf_ext
#' @export
paf_ext <- function(risk, df, prevalence, mod_prev, ...) {
  r <- risk$riskfn(df, ...)
  p <- prevalence / sum(prevalence)
  p_ <- drop(mod_prev %*% p)

  I <- sum(r * p)
  I_ <- sum(r * p_)

  return(1 - I_ / I)
}

#' @rdname paf_ext
#' @export
paf_ext_ci <- function(risk, df, prevalence, mod_prev, var_prev, ...,
                               method = c("delta"), level = 0.95) {
  # risk and prevalence -- underscore denotes modification
  r <- risk$riskfn(df, ...)
  p <- prevalence
  p_ <- drop(mod_prev %*% p)

  # incidence
  I <- sum(r * p)
  I_ <- sum(r * p_)

  # point estimate
  PAF <- log(I_) - log(I)

  # derivatives/gradient/Jacobian of risk wrt coefficients
  der_r <- risk$dtransvar(df, ...)

  # PAF dertivatives wrt coefficients (b) and prevalences (p)
  der_paf_b <- colSums(der_r * p_) / I_ - colSums(der_r * p) / I
  der_paf_p <- t(mod_prev) %*% r / I_ - r / I

  # Variance estimate
  PAFs2 <- drop(t(der_paf_b) %*% risk$vcov %*% der_paf_b +
                  t(der_paf_p) %*% var_prev %*% der_paf_p)

  # confidence interval
  a <- (1 - level) / 2
  a <- c(1 - a, a) # reversed due to transformation
  PAF <- c(PAF, PAF + sqrt(PAFs2) * qnorm(a))
  PAF <- -expm1(PAF)

  structure(PAF, names = c("PAF", "lwr", "upr"))
}
