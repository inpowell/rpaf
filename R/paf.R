#' PAF estimates in a specific population
#'
#' Calculates a point estimate (\code{paf}) or a point estimate with
#' confidence interval (\code{paf_ci}) for a PAF given an explicit
#' population, information about the risk and a counterfactual
#' distribution.
#'
#' If \code{data} is missing and the original population (used for risk
#' estimation) can be accessed in \code{risk}, then the PAF will be
#' calculated from the original population.
#'
#' In this function, the counterfactual population is defined by
#' modifying entries in the population data frame. This is done
#' through the user-specified function \code{mod_data}, which should
#' accept a single input data frame and outputs a similar data frame.
#' The output data frame should have the same column structure and
#' number of rows as the original data frame.
#'
#' @seealso \code{\link{paf_ext}} for PAF calculations with external
#'   prevalences
#'
#' @param risk a \code{\link{risk-object}}
#' @param data a data frame describing the new population
#' @param mod_data a function to modify \code{data}; see Details
#' @param se_fit if \code{TRUE}, then return the estimated standard
#'   error in transformed space
#' @param method method to calculate confidence intervals (currently
#'   the delta method alone is supported)
#' @param level the confidence level required, default 0.95
#' @param ... further arguments passed to \code{risk$riskfn} and
#'   \code{risk$dtransvar}
#' @param na.action how NAs are treated, as per
#'   \code{\link{na.action}} using the terms from
#'   \code{\link{risk-object}}.
#'
#' @return the PAF estimate as a length-one numeric vector
#' @name paf
NULL

#' @rdname paf
#' @export
paf <- function(risk, data, mod_data = identity, ..., na.action = options("na.action")) {
  if (missing(data))
    data <- risk$source_data

  if (!is.null(risk$terms) && missing(na.action)) {
    data <- model.frame(risk$terms, data = data)
  } else if (!is.null(risk$terms) && !missing(na.action)) {
    data <- model.frame(risk$terms, data = data, na.action = na.action)
  } else {
    data <- match.fun(na.action)(data)
  }

  # build counterfactual scenario
  mdata <- mod_data(data)

  # estimate risk and incidence
  r <- risk$riskfn(data, ...)
  r_ <- risk$riskfn(mdata, ...)
  I <- sum(r)
  I_ <- sum(r_)

  return(1 - I_ / I)
}

#' @rdname paf
#' @export
paf_ci <- function(risk, data, mod_data = identity, ..., se_fit = TRUE,
                           method = c("delta"), level = 0.95, na.action = options("na.action")) {
  if (missing(data))
    data <- risk$source_data

  if (!is.null(risk$terms) && missing(na.action)) {
    data <- model.frame(risk$terms, data = data)
  } else if (!is.null(risk$terms) && !missing(na.action)) {
    data <- model.frame(risk$terms, data = data, na.action = na.action)
  } else {
    data <- match.fun(na.action)(data)
  }

  # build counterfactual scenarios
  data_ <- mod_data(data)

  # estimate risk and incidence
  r <- risk$riskfn(data, ...)
  r_ <- risk$riskfn(data_, ...)
  I <- sum(r)
  I_ <- sum(r_)

  # calculate gradients
  der_r <- risk$dtransvar(data, ...)
  der_r_ <- risk$dtransvar(data_, ...)

  # ... with respect to coefficients:
  der_paf_b <- colSums(der_r_) / I_ - colSums(der_r) / I

  # transformed PAF estimate and variance
  PAF <- log(I_) - log(I)
  PAFs2 <- risk$var(der_paf_b)

  # confidence interval calculation
  a <- (1 - level) / 2
  a <- c(1 - a, a) # reverse order because -expm1 is a negative transformation

  PAF <- c(PAF, PAF + sqrt(PAFs2) * qnorm(a))
  PAF <- -expm1(PAF)

  ret <- structure(PAF, names = c("PAF", "lwr", "upr"))
  if (se_fit)
    ret <- c(ret, "se.trans" = sqrt(PAFs2))

  ret
}

#' PAF estimates with external prevalences
#'
#' Calculates a point estimate (\code{paf_ext}) or a point estimate with
#' confidence interval (\code{paf_ext_ci}) for a PAF given an population and
#' prevalence estimates, information about the risk and a counterfactual
#' distribution.
#'
#' If \code{data} is missing and the original population (used for risk
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
#'   rows in \code{data}
#' @param var_prev variance-covariance matrix for the prevalence estimates
#' @param mod_prev a transition matrix to modify \code{data}; see Details
#'
#' @name paf_ext
NULL

#' @rdname paf_ext
#' @export
paf_ext <- function(risk, data, prevalence, mod_prev, ..., na.action = options("na.action")) {
  if (!is.null(risk$terms) && missing(na.action)) {
    data <- model.frame(risk$terms, data = data)
  } else if (!is.null(risk$terms) && !missing(na.action)) {
    data <- model.frame(risk$terms, data = data, na.action = na.action)
  } else {
    data <- match.fun(na.action)(data)
  }

  r <- risk$riskfn(data, ...)
  p <- prevalence / sum(prevalence)
  p_ <- drop(mod_prev %*% p)

  I <- sum(r * p)
  I_ <- sum(r * p_)

  return(1 - I_ / I)
}

#' @rdname paf_ext
#' @export
paf_ext_ci <- function(risk, data, prevalence, mod_prev, var_prev, ..., se_fit = TRUE,
                               method = c("delta"), level = 0.95, na.action = options("na.action")) {
  if (!is.null(risk$terms) && missing(na.action)) {
    data <- model.frame(risk$terms, data = data)
  } else if (!is.null(risk$terms) && !missing(na.action)) {
    data <- model.frame(risk$terms, data = data, na.action = na.action)
  } else {
    data <- match.fun(na.action)(data)
  }

  # risk and prevalence -- underscore denotes modification
  r <- risk$riskfn(data, ...)
  p <- prevalence
  p_ <- drop(mod_prev %*% p)

  # incidence
  I <- sum(r * p)
  I_ <- sum(r * p_)

  # point estimate
  PAF <- log(I_) - log(I)

  # derivatives/gradient/Jacobian of risk wrt coefficients
  der_r <- risk$dtransvar(data, ...)

  # PAF dertivatives wrt coefficients (b) and prevalences (p)
  der_paf_b <- colSums(der_r * p_) / I_ - colSums(der_r * p) / I
  der_paf_p <- t(mod_prev) %*% r / I_ - r / I

  # Variance estimate
  PAFs2 <- drop(risk$var(der_paf_b) +
                  t(der_paf_p) %*% var_prev %*% der_paf_p)

  # confidence interval
  a <- (1 - level) / 2
  a <- c(1 - a, a) # reversed due to transformation
  PAF <- c(PAF, PAF + sqrt(PAFs2) * qnorm(a))
  PAF <- -expm1(PAF)

  ret <- structure(PAF, names = c("PAF", "lwr", "upr"))

  if (se_fit)
    ret <- c(ret, "se.trans" = sqrt(PAFs2))

  ret
}
