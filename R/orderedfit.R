#' Ordered Logistic/Probit Model Estimation with Weighted Likelihood
#'
#' Estimates an ordinal regression model using a weighted likelihood approach.
#' This function is adapted with minor modifications from Dong, Nakaya, and Brunsdon (2020),
#' who proposed it as a local estimation procedure in their geographically weighted ordinal regression model.
#'
#' @param y An ordered factor vector of responses.
#' @param X A matrix or data frame of independent variables (without intercept).
#' @param link The link function to use. Either \code{"probit"} (default) or \code{"logit"}.
#' @param Gweights A numeric vector of weights for each observation (e.g., kernel weights from spatial models).
#' @param id.fixed.vars Integer vector specifying which coefficients to hold fixed during optimization. Use \code{NA} to indicate no fixed coefficients.
#' @param start.values Initial values for the parameters to be optimized.
#'
#' @return A list containing:
#' \item{LogL}{Log-likelihood at the optimum.}
#' \item{coefficients}{Estimated coefficients including thresholds.}
#' \item{std.err}{Standard errors of the coefficients.}
#' \item{t.value}{t-statistics of the coefficients.}
#' \item{Betas}{Estimated slope coefficients.}
#' \item{zeta}{Estimated threshold values (not the reparameterized form).}
#' \item{V.g}{Variance-covariance matrix for the slope coefficients.}
#'
#' @details
#' This function is adapted with only minor modifications from Dong, G., Nakaya, T., & Brunsdon, C. (2020).
#' It retains their original formulation of the likelihood, gradient, and Hessian for efficient estimation of an ordered response model.
#'
#' The following notes are quoted directly from the original authors:
#' \emph{
#' " First, write the likelihood function of a standard ordered model and its gradient and hessian for fast optimisation.
#' We assume there are five ordered responses as in our SWB data, but extension to data with fewer or more categories is trivial.
#' ALSO, we DO NOT have an intercept in the model. Then, we maximise the log-likelihood function and output estimated parameters."
#' }
#'
#' @references
#' Dong, G., Nakaya, T., & Brunsdon, C. (2020).
#' Geographically weighted regression models for ordinal categorical response variables: An application to geo-referenced life satisfaction data.
#' \emph{Computers, Environment and Urban Systems}, 80, 101428.
#'
#' @importFrom stats dnorm dlogis pnorm
#' @importFrom maxLik maxLik
#' @importFrom MASS ginv
#'
#' @export
orderedfit <- function(y, X, link = "probit", Gweights, id.fixed.vars, start.values) {
  # param: model parameters--[betas;gamm1;gamm2;gamm3;gamm4]
  # y: responses--[1,2,3,4,5]
  # X: independent variables--ideally they are of data.frame objects
  # link: link function--it could a probit link or a logit link
  # later I might use a simple formula object to extract these arguments

  ## a lot of preliminary quantities that are used in evaluations of likelihood and gradient
  X <- X
  ind.names <- colnames(X)
  X <- as.matrix(X)

  # the number of observations and indenpendent variables
  n <- nrow(X)
  p <- ncol(X)

  lev.y <- levels(y)
  ly <- nlevels(y) - 1
  y <- unclass(y)


  # the log-likelihood function
  loglordered <- function(param) {
    # betas come first then threshold parameters (gamm1, gamm2, gamm3, gamm4)
    # note to satisfy the increasing order of threshold parameters, we use a reparameterisation procedure
    betas <- param[1:p]
    theta <- param[p + seq_len(ly)]

    # threshold values augumented by gamm0 = -Inf and gamm5 = Inf
    # using an exponential reparameterisation
    gamm <- c(-Inf, cumsum(c(theta[1], exp(theta[-1]))), Inf)

    # the linear predictor effect
    eta <- drop(X %*% betas)

    # calculate the probabilities corresponding the actual response
    # note pfun is the chosen link function, which we shall see later. pfun calculate the cumulative probabilities and
    # dfun calculate the densities
    z1 <- pmin(100, gamm[y + 1] - eta)
    z2 <- pmax(-100, gamm[y] - eta)

    if (link == "probit") {
      pr <- pnorm(z1) - pnorm(z2)
      p1 <- dnorm(z1)
      p2 <- dnorm(z2)
    } else {
      pr <- plogis(z1) - plogis(z2)
      p1 <- dlogis(z1)
      p2 <- dlogis(z2)
    }


    # the log-likelihood values
    # log.l <- -sum(Gweights * log(pr))
    if (all(pr > 0)) log.l <- sum(Gweights * log(pr)) else log.l <- -Inf

    ### calculate the gradient of the log-likelihood function

    # First calculate the Jacobian matrix due to our use of reparameterisation in optimisation
    k <- length(theta)
    exp.theta <- exp(theta)
    jacobian <- matrix(0, k, k)
    jacobian[, 1] <- rep(1, k)
    for (i in 2:k) jacobian[i:k, i] <- exp.theta[i]

    # the gradient with respect to betas
    # densities using dfun


    grad.betas <- t(X) %*% (Gweights * (p1 - p2) / pr)

    # the gradient with respect to gamms
    # we have 4 threshold parameters to estimate

    y.temp <- matrix(0, n, ly)
    ind.Y1 <- col(y.temp) == y
    ind.Y2 <- col(y.temp) == (y - 1)
    XX <- ind.Y1 * p1 - ind.Y2 * p2
    grad.gamm <- -t(XX) %*% (Gweights / pr)
    grad.gamm <- t(grad.gamm) %*% jacobian

    # combine them
    gradient.logl <- -c(grad.betas, grad.gamm)
    # assign a gradient attribute
    attr(log.l, "gradient") <- gradient.logl
    return(log.l)
  }

  #### Now optimise the log-likelihood function
  # using an unified maximisation approach provided in maxLik R package
  # check if we have global coefficients
  if (all(is.na(id.fixed.vars))) {
    # optimisation without global coefficients
    res.opt <- maxLik(logLik = loglordered, start = start.values, method = "BFGS")
  } else {
    # optimisation with global coefficients constraints
    res.opt <- maxLik(
      logLik = loglordered, start = start.values, method = "BFGS",
      fixed = id.fixed.vars
    )
  }

  ### save results---coefficients, standard errors, t-values, log-l
  LogL <- res.opt$maximum
  Betas <- res.opt$estimate[1:p]
  theta <- res.opt$estimate[p + seq_len(ly)]
  zeta <- cumsum(c(theta[1], exp(theta[-1])))
  names(Betas) <- ind.names
  names(zeta) <- paste(lev.y[-length(lev.y)], lev.y[-1], sep = "|")
  coefficients <- c(Betas, zeta)
  # get inference for zetas rather than thetas
  # final jacobian matrix
  vc <- ginv(-res.opt$hessian)
  k <- length(theta)
  J <- matrix(0, k, k)
  J[, 1] <- rep(1, k)
  for (i in 2:k) J[i:k, i] <- exp(theta)[i]
  A <- diag(p + k)
  A[p + seq_len(k), p + seq_len(k)] <- J
  # the final variance-covariance matrix
  V <- A %*% vc %*% t(A)
  std.err <- drop(sqrt(diag(V)))
  t.value <- drop(coefficients / std.err)

  ### as we need to calcuate the willingness-to-pay measures we need to save the covariance between income and pollution
  V.g <- V[1:p, 1:p]
  if (!all(is.na(id.fixed.vars))) {
    temp.ind <- seq_len(p)[-id.fixed.vars]
    V.g <- V.g[temp.ind, temp.ind]
  }
  res <- list(
    LogL = LogL, coefficients = coefficients, std.err = std.err,
    t.value = t.value, Betas = Betas, zeta = zeta, V.g = V.g
  )
  return(res)
}
