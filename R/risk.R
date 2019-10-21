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
#' \item{\code{terms}}{(optional) A \code{formula} or \code{terms}
#' object from the underlying model that can be used for missing data
#' handling. If this is not specified, then the PAF methods will use
#' only complete cases.}
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
    terms = Terms,
    source_df = model.frame(model)
  )
}

#' Risk information for location-scale survival models
#'
#' This function extracts information on risk from survival models based on a
#' location-scale distribution.
#'
#' @param model a \code{survreg} object
#'
#' @return a \code{\link{risk-object}} list
#' @export
#'
#' @importFrom stats delete.response model.frame model.matrix vcov
#' @importFrom survival strata survreg.distributions untangle.specials
#'
#' @family risk functions
risk_survreg <- function(model) {
  # See https://github.com/therneau/survival/blob/master/R/predict.survreg.R as
  # the inspiration for much of this function

  # What do we need?
  #  - survival quantile: time, transformation, linear predictor, stratum scale
  #  - gradients: time, transformation, linear predictor, scale, stratum index

  Terms <- model$terms
  if (!inherits(Terms, "terms"))
    stop("Inavlid terms component of model")

  strata <- attr(Terms, "specials")$strata
  Terms <- delete.response(Terms)
  coef <- model$coefficients
  nvar <- length(model$coefficients)
  fixedscale <- (nvar == ncol(model$var))

  #
  # Grab the distribution -- copied with minor modification for Therneau
  #
  if (is.character(model$dist)) dd <- survreg.distributions[[model$dist]]
  else dd <- model$dist
  if (is.null(dd$trans)) {
    trans <- function(x) x  # identity transformation
  } else {
    trans <- dd$trans
    dtrans <- dd$dtrans
  }
  if (!is.null(dd$dist)) dd <- survreg.distributions[[dd$dist]]

  # this block defines get_strata <- function(df) {...} to get indices
  if (length(strata) && !fixedscale) {
    #
    # We need to reconstruct the original "strata" variable
    #
    mf <- model.frame(model)
    strat.unt <- untangle.specials(Terms, 'strata', 1)
    if (length(strat.unt$vars)==1) strata.keep <- mf[[strat.unt$vars]]
    else strata.keep <- strata(mf[,strat.unt$vars], shortlabel=TRUE)
    strata <- as.numeric(strata.keep)
    nstrata <- max(strata)

    get_strata <- function(df) {
      if (length(strat.unt$vars)==1) newstrat <- df[[strat.unt$vars]]
      else newstrat <- strata(df[,strat.unt$vars], shortlabel=TRUE)
      match(newstrat, levels(strata.keep))
    }
  } else {
    nstrata <- 1
    get_strata <- function(df) rep(1L, nrow(df))
  }

  list(
    riskfn = function(df, time, ...) {
      X <- model.matrix(Terms, data = df, xlev = model$xlevels, contrasts.arg = model$contrasts)
      lp <- drop(X %*% coef)
      gt <- trans(time)
      sigma <- model$scale[get_strata(df)]
      return(dd$density((gt - lp) / sigma)[,1]) # first column corresponds to CDF
    },
    dtransvar = function(df, time, ...) {
      X <- model.matrix(Terms, data = df, xlev = model$xlevels, contrasts.arg = model$contrasts)
      lp <- drop(X %*% coef)
      gt <- trans(time)
      strata_i <- get_strata(df)
      sigma <- model$scale[strata_i]

      qmean <- (gt - lp) / sigma
      density <- dd$density(qmean)[,3] # third column corresponds to PDF

      J <- -density / sigma * X # X is a matrix, but multiplication goes by columns

      if (!fixedscale) {
        dlsigma <- matrix(NA, nrow = nrow(J), ncol = nstrata)
        dlsigma[] <- ifelse(col(dlsigma) == strata_i, -qmean * density, 0)
        J <- cbind(J, dlsigma)
      }
      return(J)
    },
    hazardfn = function(df, time, ...) {
      X <- model.matrix(Terms, data = df, xlev = model$xlevels, contrasts.arg = model$contrasts)
      lp <- drop(X %*% coef)
      gt <- trans(time)
      dgt <- dtrans(time)
      sigma <- model$scale[get_strata(df)]

      FF <- dd$density((gt - lp) / sigma)

      return(dgt / sigma * FF[,3] / FF[,2])
    },
    dhazardfn = function(df, time, ...) {
      X <- model.matrix(Terms, data = df, xlev = model$xlevels, contrasts.arg = model$contrasts)
      lp <- drop(X %*% coef)
      gt <- trans(time)
      dgt <- dtrans(time)
      sigma <- model$scale[get_strata(df)]

      FF <- dd$density((gt - lp) / sigma)

      return(X * (-dgt/sigma^2) * FF[,3] / FF[,2] * (FF[,4] + FF[,3] / FF[,2]))
    },
    terms = Terms,
    var = function(d) drop(t(d) %*% vcov(model) %*% d)
  )
}

#' Risk information for constant piecewise hazards model
#'
#' This function returns a \code{\link{risk-object}} based on a
#' parametric piecewise constant hazards model. The risk calculation
#' is performed largely by an underlying risk object from
#' \code{\link{risk_survreg}}, with some convenience added.
#'
#' This function creates time-varying covariates from \code{data}, and
#' for that reason, the column names \code{.period}, \code{.ID},
#' \code{.tstart}, \code{.tstop} and \code{.event} should *not* be
#' used. These column names can be used in \code{formula} however,
#' and any arguments specified in \code{...} passed on to
#' \code{\link[survival]{survreg}}, and are taken to mean \describe{
#'
#' \item{\code{.period}}{the time period varying as a factor}
#'
#' \item{\code{.ID}}{a unique identifier for each input row}
#'
#' \item{\code{.tstart}, \code{.tstop}}{start and end times of each
#' observation period, often the endpoints of the time period (if
#' censoring occurs)}
#'
#' \item{\code{.event}}{an event indicator}
#'
#' }
#'
#' Note that the response in \code{formula} should be specified
#' without any reference to these intermediate variables.
#'
#' @param formula the model formula, which may include \code{.period}
#'   (see Details)
#' @param data the data frame from which observations are taken
#' @param breaks a vector of time points demarcating episodes of
#'   constant hazard, *including* beginning and end points
#' @param ... further arguments to be passed to
#'   \code{\link[survival]{survreg}}
#'
#' @return a \code{\link{risk-object}} list
#' @export
risk_pch <- function(formula, data, breaks, ...) {
  newvars <- c('.tstart', '.tstop', '.event', '.period', '.ID')
  Terms <- terms(formula, data = data)
  ttf <- attr(Terms, "factors")
  ttn <- dimnames(ttf)[[1]]
  dropx <- which(as.logical(colSums(ttf[ttn %in% newvars, , drop = FALSE])))
  s_tt <- drop.terms(Terms, dropx = dropx,
                     keep.response = TRUE)

  dt <- diff(breaks)
  cuts <- breaks[-c(1, lengths(breaks))]
  period_labels <- paste0("(", breaks[-length(breaks)], ",", breaks[-1], "]")
  period_base <- factor(1:length(dt), levels = 1:length(dt),
                        labels = period_labels)

  longdata <- survSplit(s_tt, data = data, cut = cuts, zero = breaks[1],
                        start = '.tstart', end = '.tstop', event = '.event',
                        episode = '.period', id = '.ID')
  longdata$.period <- period_base[longdata$.period] # applies levels

  lform <- update(formula, Surv(.tstop - .tstart, .event) ~ .)
  model <- survreg(lform, data = longdata, ..., dist = "exponential")
  risk <- risk_survreg(model)

  lengthen <- function(df) {
    df$.ID <- seq_len(nrow(df))
    merge(df, data.frame(".period" = period_base), by = NULL)
  }

  list(
    riskfn = function(df, time, ...) {
      ldf <- lengthen(df)
      # ignore periods outside of time range
      keep <- (breaks[-length(breaks)] <= time)[ldf$.period]
      ldf <- subset(ldf, keep)

      # calculate start times, end times
      tstart <- breaks[-length(breaks)][ldf$.period]
      tstop <- breaks[-1][ldf$.period]
      t <- pmin(time, tstop) - tstart

      r <- risk$riskfn(ldf, t)
      1 - tapply(1 - r, ldf$.ID, FUN = prod)
    },
    dtransvar = function(df, time, ...) {
      ldf <- lengthen(df)
      # ignore periods outside of time range
      keep <- (breaks[-length(breaks)] <= time)[ldf$.period]
      ldf <- subset(ldf, keep)

      # calculate start times, end times
      tstart <- breaks[-length(breaks)][ldf$.period]
      tstop <- breaks[-1][ldf$.period]
      t <- pmin(time, tstop) - tstart

      r <- risk$riskfn(ldf, t)
      J <- risk$dtransvar(ldf, t)

      P <- tapply(1 - r, ldf$.ID, FUN = prod)
      S <- (1/(1 - r)) * J # summand matrix

      # sums complete -- rows correspond to people, cols to regression
      # coefficients
      tS <- tapply(S, list(ldf$.ID[row(S)], col(J)), FUN = sum)

      apply(tS, 2, "*", P)
    },
    hazardfn = function(df, time, ...) {
      # get all hazards up to time
      ldf <- lengthen(df)
      l <- risk$hazardfn(ldf, time, ...)
      vec2mat(l, ldf$.ID, ldf$.period)
    },
    dhazardfn = function(df, time, ...) {
      ldf <- lengthen(df)
      l <- risk$dhazardfn(ldf, time, ...)
      mat2arr(l, ldf$.ID, ldf$.period)
    },
    terms = delete.response(s_tt),
    var = function(d) drop(t(d) %*% vcov(model) %*% d)
  )
}

#' Risk information from relative risks and their standard errors
#'
#' This function synthesises relative risk estimates and standard
#' errors from a specified set of factors. The covariance matrix is
#' estimated conservatively according to the \code{eta} parameter:
#' \eqn{\eta = 1} corresponds to the "worst-case" variance of the PAF,
#' while \eqn{\eta = 0} corresponds to uncorrelated log-relative
#' risks.
#'
#' @param RR a named list of relative risk estimates; see Examples
#' @param se.trans a named list of standard error estimates for log-relative risks
#' @param eta number between 0 and 1, ranging from least to most conservative
#'
#' @return a \code{\link{risk-object}} list
#' @export
#'
#' @family risk functions
#'
#' @example inst/examples/eg_risk_relative.R
risk_relative <- function(RR, se.trans, eta = 1.0) {
  if (eta > 1.0 || eta < 0.0)
    stop("Parameter eta must be between 0 and 1.")

  if (names(RR) != names(se.trans) ||
      !all.equal(lapply(RR, names), lapply(se.trans, names)))
    stop("Lists RR and se.trans must have same name order.")

  # calculate total relative risks using outer product (note:
  # reference is upper left entry)
  orr <- Reduce(outer, RR)
  names(dimnames(orr)) <- names(RR)

  # convert to data frame and remove risk column, which should always
  # appear at the end
  ref <- as.data.frame.table(orr)
  rr <- ref[[length(RR) + 1]]
  ref <- subset(ref, select = -(length(RR) + 1))

  # Model matrix equivalent
  Xref <- lapply(ref, function(fact) {
    mat <- matrix(0, nrow = length(fact), ncol = nlevels(fact))
    mat[col(mat) == as.numeric(fact)] <- 1
    mat
  })
  Xref <- do.call(cbind, Xref)

  # Jacobian reference
  Jref <- Xref * rr

  # Standard error reference
  SEref <- do.call(c, se.trans)

  list(
    riskfn = function(df, ...) {
      rr[match.data.frame(df, ref)]
    },
    dtransvar = function(df, ...) {
      Jref[match.data.frame(df, ref), ]
    },
    var = function(d) {
      s <- sign(d)
      R <- eta * s %o% s
      diag(R) <- 1
      V <- (SEref %o% SEref) * R
      drop(t(d) %*% V %*% d)
    }
  )
}
