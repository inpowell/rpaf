#' Risk information for constant piecewise hazards model with
#' competing risks
#'
#' This function returns a \code{\link{risk-object}} based on two
#' parametric piecewise constant hazards models, one each for the
#' primary and secondary events.
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
#' @family risk functions
#'
#' @return a \code{\link{risk-object}} list
#' @export
risk_cr <- function(resp1, resp2, predictors, data, breaks, ...) {
  # get Surv objects from formulas
  mr1 <- model.response(model.frame(resp1, data, na.action = na.pass))
  mr2 <- model.response(model.frame(resp2, data, na.action = na.pass))
  # convert to matrices
  mmr1 <- as.matrix(mr1)
  mmr2 <- as.matrix(mr2)
  # censor secondary event by primary
  mmr2[,1] <- pmin(mmr1[,1], mmr2[,1])
  mmr2[as.logical(mmr1[,2]), 2] <- 0

  f1 <- Surv(mmr1[,1], mmr1[,2]) ~ .
  f2 <- Surv(mmr2[,1], mmr2[,2]) ~ .

  # new model.frames with data
  data1 <- model.frame(f1, data, na.action = na.pass)
  data2 <- model.frame(f2, data, na.action = na.pass)

  risk1 <- risk_pch(update(f1, predictors), data1, breaks)
  risk2 <- risk_pch(update(f2, predictors), data2, breaks)

  list(
    # SCOPING WARNING: conflict in `data`
    riskfn = function(data, time, ...) {
      # hazards lambda
      l1 <- risk1$hazardfn(data, time)
      l2 <- risk2$hazardfn(data, time)

      diffs <- diff(pmin(time, breaks))
      S1 <- diffs * t(l1)
      S1 <- apply(S1, 2, cumsum)
      S1 <- t(exp(-S1))
      S2 <- diffs * t(l2)
      S2 <- apply(S2, 2, cumsum)
      S2 <- t(exp(-S2))

      S1_lag <- cbind(1, S1[,-ncol(S1)])
      S2_lag <- cbind(1, S2[,-ncol(S2)])

      # final risk
      rowSums(l1 / (l1 + l2) * (S1_lag * S2_lag - S1 * S2))
    },
    # SCOPING WARNING: conflict in `data`
    dtransvar = function(data, time, ...) {
      l1 <- risk1$hazardfn(data, time)
      l2 <- risk2$hazardfn(data, time)

      diffs <- diff(pmin(time, breaks))
      S1 <- diffs * t(l1)
      S1 <- apply(S1, 2, cumsum)
      S1 <- t(exp(-S1))
      S2 <- diffs * t(l2)
      S2 <- apply(S2, 2, cumsum)
      S2 <- t(exp(-S2))

      S1_lag <- cbind(1, S1[,-ncol(S1)])
      S2_lag <- cbind(1, S2[,-ncol(S2)])

      # hazard gradients
      # d1: ID; d2: time; d3: coef
      dl1 <- risk1$dhazardfn(data, time)
      dl2 <- risk2$dhazardfn(data, time)

      # survival gradients
      dS1 <- aperm(dl1, c(2, 1, 3)) * diffs # d1: time; d2: ID; d3: coef
      dS1 <- apply(dS1, 2:3, cumsum)
      dS1 <- aperm(dS1, c(2, 1, 3)) * as.vector(S1)
      dS2 <- aperm(dl2, c(2, 1, 3)) * diffs
      dS2 <- apply(dS2, 2:3, cumsum)
      dS2 <- aperm(dS2, c(2, 1, 3)) * as.vector(S2)

      dS1_lag <- dS2_lag <- array(dim = dim(dS1))
      dS1_lag[,1,] <- dS2_lag[,1,] <- 0
      dS1_lag[, 2:(dim(S1)[2]), ] <- dS1[, 1:(dim(S1)[2] - 1), ]
      dS2_lag[, 2:(dim(S2)[2]), ] <- dS2[, 1:(dim(S2)[2] - 1), ]

      dr1 <- dl1 * as.vector(l2 / (l1 + l2)^2) * as.vector(S1_lag * S2_lag - S1 * S2) +
        (dS1_lag * as.vector(S2_lag) - dS1 * as.vector(S2)) * as.vector(l1 / (l1 + l2))
      dr1 <- apply(dr1, c(1, 3), sum) # d1: ID; d2: coef

      dr2 <- dl2 * as.vector(-l1 / (l1 + l2)^2) * as.vector(S1_lag * S2_lag - S1 * S2) +
        (dS2_lag * as.vector(S1_lag) - dS2 * as.vector(S1)) * as.vector(l1 / (l1 + l2))
      dr2 <- apply(dr2, c(1, 3), sum) # d1: ID; d2: coef

      cbind(dr1, dr2)
    },
    terms = risk1$terms,
    var = function(d) {
      p <- length(d) %/% 2L
      risk1$var(d[1:p]) + risk2$var(d[-(1:p)])
    }
  )
}
