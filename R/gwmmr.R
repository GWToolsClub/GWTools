#' Geographically Weighted Multivariate Multiple Regression (GWMMMR)
#'
#' This function fits a Geographically Weighted Multivariate Multiple Regression (GWMMR) model.
#' It provides local estimates of regression coefficients, standard errors, fitted values,
#' and conducts local multivariate hypothesis testing using Pillai's trace.
#'
#' @param formula A formula object specifying the model structure (e.g., \code{cbind(y1, y2) ~ x1 + x2 + x3}).
#' @param cordxy A matrix of spatial coordinates, with columns representing X and Y.
#' @param distmethod A character string indicating the distance calculation method to use.
#' Supported options are \code{"euclidean"} (default), \code{"greatcircle"} (great-circle distance), and \code{"manhattan"} (city-block distance).
#' @param h The bandwidth used in the kernel function.
#' @param kernel A character string specifying the kernel function used to compute spatial weights.
#' Options are \code{"global"}, \code{"gaussian"}, \code{"exponential"}, \code{"bisquare"}, \code{"tricube"}, and \code{"boxcar"}.
#' @param adaptive Logical. If \code{TRUE} (default), uses adaptive bandwidth. If \code{FALSE}, uses fixed bandwidth.
#' @param data A data frame containing all variables used in the model.
#'
#' @return A list of class \code{"gwmmr_model"} containing:
#' \itemize{
#'    \item \code{formula}: Model formula.
#'    \item \code{test}: A data.frame with local coefficients, standard errors, fitted values, prediction errors, and multivariate test statistics.
#'    \item \code{SSCP_E}: Residual sum of squares and cross-products matrix.
#'    \item \code{sigmahat}: Estimated residual covariance matrix.
#'    \item \code{epp}: Trace of the local hat matrix (S), representing the effective number of parameters across all response variables.
#'    \item \code{epp2}: Trace of S'S, representing degrees of freedom used in estimating residual covariance.
#'    \item \code{enp}: the effective number of parameters per response variable.
#'    \item \code{enp2}: A bias-adjusted effective degrees of freedom used for AICc and local multivariate testing.
#'    \item \code{df2}: Degrees of freedom used for standard error and AICc calculation.
#'    \item \code{AIC}, \code{AICc}: Model selection criteria values.
#'    \item \code{call}: The matched call.
#' }
#'
#' @references
#' Chen, V. Y.-J., Yang, T.-C., & Jian, H.-L. (2022).
#' Geographically Weighted Regression Modeling for Multiple Outcomes.
#' *Annals of the American Association of Geographers*, 112(1), 1–18.
#'
#' @examples
#' set.seed(123)
#'
#' n <- 300
#' x <- runif(n, 0, 100)
#' y <- runif(n, 0, 100)
#'
#' X1 <- rnorm(n, 50, 10)
#' X2 <- rnorm(n, 30, 5)
#' X3 <- rnorm(n, mean = 10, sd = 2)
#'
#' y1 <- 1.2 * X1 - 0.5 * X2 + 1.5 * X3 + rnorm(n, 0, 2)
#' y2 <- 0.8 * X1 + 0.3 * X2 + 1.2 * X3 + rnorm(n, 0, 2)
#'
#' data <- data.frame(x = x, y = y, X1 = X1, X2 = X2, X3 = X3, y1 = y1, y2 = y2)
#'
#' formula <- cbind(y1, y2) ~ X1 + X2 + X3
#' coords <- cbind(data$x, data$y)
#'
#' gwmmr_result <- gwmmr(
#'   formula = formula,
#'   cordxy = coords,
#'   h = 200,
#'   kernel = "bisquare",
#'   adaptive = TRUE,
#'   data = data
#' )
#' gwmmr_result
#' @export
gwmmr <- function(formula, cordxy, distmethod = "euclidean",
                  h, kernel, adaptive = TRUE, data) {
  start <- proc.time()

  if (!is.matrix(cordxy)) {
    stop("cordxy must be a matrix.")
  }

  # ----- kernel function ----
  kernel_weights <- function(d, h, nr, kernel, adaptive) {
    if (adaptive == TRUE) {
      rh <- round(h)
      sd <- sort(d)
      bw <- sd[rh]
    } else {
      bw <- h
    }

    # w <- rep(0, length(d))
    if (kernel == "global") {
      w <- rep(1, nr)
    } else if (kernel == "gaussian") { # /*gaussian kernel*/
      w <- exp((-0.5) * ((d / bw)**2))
    } else if (kernel == "exponential") {
      w <- exp(-d / bw) # /*exponential kernel*/
    } else if (kernel == "bisquare") {
      w <- (1 - (d / bw)^2)^2
      index <- which(d > bw)
      w[index] <- 0 # /*bisquare nearest neighbor*/
    } else if (kernel == "tricube") {
      w <- (1 - abs(d / bw)^3)^3
      index <- which(d > bw)
      w[index] <- 0
    } else if (kernel == "boxcar") {
      w <- ifelse(d <= bw, 1, 0)
    } else {
      stop("Unsupported kernel type. Choose from
           'global', 'gaussian', 'exponential',
           'bisquare', 'tricube' or 'boxcar'.")
    }
    w <- as.numeric(w)
    return(w)
  }

  # ---- Prepare data ----
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)
  ox <- model.matrix(formula, data = data)[, -1, drop = FALSE]
  n <- nrow(data)
  nr <- nrow(ox)
  X <- as.matrix(cbind(1, ox))
  nv <- ncol(X)
  ny <- nrow(y)
  nc <- ncol(y)
  Y <- array(y, dim = c(ny, nc))
  beta <- array(0, dim = c(nv * nc, nr))
  strr <- array(0, dim = c(nr, nv))
  pred <- array(0, dim = c(ny, nc))
  leverage <- array(0, dim = c(nr, 1))
  leverage2 <- array(0, dim = c(nr, 1))

  for (i in 1:n) {
    d <- gwdist(cordxy, method = distmethod, target_idx = i)
    w <- kernel_weights(d, h, nr, kernel, adaptive)

    newx <- X * sqrt(w) # n × p
    newy <- Y * sqrt(w) # n × m
    xpx <- t(newx) %*% newx # p × p
    xpxi <- MASS::ginv(xpx) # inverse crossproducts
    b <- xpxi %*% t(newx) %*% newy # p × m
    beta[, i] <- b
    yhat <- X[i, ] %*% b # 1 × m
    pred[i, ] <- yhat
    hatmatrix <- newx %*% xpxi %*% t(newx) # n × n
    leverage[i, ] <- diag(hatmatrix)[i] # S[i]

    C <- xpxi %*% t(newx * sqrt(w)) # p x n
    varb <- C %*% t(C) # p x p
    strr[i, ] <- sqrt(diag(varb)) # standard error of the coefficients
    leverage2[i, ] <- (X[i, ] %*% C) %*% t(X[i, ] %*% C) # (S'S)[i]
  }

  error <- Y - pred # n x m
  SSCP_E <- t(error) %*% error

  epp <- sum(leverage) * nc # trace(S)*y個數 (effective number of parameters)
  epp2 <- sum(leverage2) # trace(S'S)
  df2 <- epp2 * nc

  sigmahat <- SSCP_E / (nr - epp)
  enp <- epp / nc # trace(S)
  enp2 <- 2 * enp - epp2 # 2*trace(S) - trace(S'S)

  AIC <- nr * (log(det(SSCP_E / nr))) + nr * nc * (log(2 * 3.141593) + 1) + 2 * enp * nc + (nc * (nc + 1))
  AICc <- nr * (log(det(SSCP_E / nr))) + nr * nc * log(2 * 3.141593) + (nr * nc * (nr + enp)) / (nr - enp - nc - 1)

  # manova
  L <- diag(nv)
  Pillai <- array(0, dim = c(nr, nv))
  Fvalue_Pillai <- array(0, dim = c(nr, nv))
  P_value <- array(0, dim = c(nr, nv))
  for (i in 1:nv) {
    L_i <- array(L[i, ], dim = c(1, nv))
    M <- diag(nc)
    manova <- manova_gwmmr(
      formula = formula, cordxy = cordxy,
      distmethod = distmethod, h = h, kernel = kernel,
      adaptive = adaptive, data = data,
      L = L_i, M = M, v = enp2
    )

    Pillai[, i] <- manova$test$Pillai
    Fvalue_Pillai[, i] <- manova$test$Fvalue_Pillai
    P_value[, i] <- manova$test$P_value
  }

  betase <- array(0, dim = c(nr, nv * nc))
  bname <- array(0, dim = c(1, nv * nc))
  sname <- array(0, dim = c(1, nv * nc))

  for (k in 1:nc) {
    bname[, (k * (ncol(ox) + 1) - ncol(ox)):(k * (ncol(ox) + 1))] <- stringr::str_trim(stringr::str_c("estx", as.character(0:ncol(ox)), "m", as.character(k)))
    sname[, (k * (ncol(ox) + 1) - ncol(ox)):(k * (ncol(ox) + 1))] <- stringr::str_trim(stringr::str_c("se_estx", as.character(0:ncol(ox)), "m", as.character(k)))
    betase[, (k * (ncol(ox) + 1) - ncol(ox)):(k * (ncol(ox) + 1))] <- sqrt(as.numeric(diag(sigmahat)[k])) * strr
  }
  # test <- data.frame(t(exp(beta)), betase, pred, error, Pillai, Fvalue_Pillai, P_value, X, Y, cordxy)  # beta 做 exp 轉換
  test <- data.frame(t(beta), betase, pred, error, Pillai,
                     Fvalue_Pillai, P_value, X, Y, cordxy)

  # name
  pre <- stringr::str_trim(stringr::str_c("pred", as.character(1:nc)))
  Pillainame <- stringr::str_trim(stringr::str_c("Pillai", as.character(0:ncol(ox))))
  Fvalue_Pillainame <- stringr::str_trim(stringr::str_c("Fvalue_Pillai", as.character(0:ncol(ox))))
  P_valuename <- stringr::str_trim(stringr::str_c("P_value", as.character(0:ncol(ox))))
  xxname <- stringr::str_trim(stringr::str_c("x", as.character(0:ncol(ox))))
  yyname <- stringr::str_trim(stringr::str_c("y", as.character(1:nc)))
  errorname <- stringr::str_trim(stringr::str_c("error", as.character(1:nc)))
  cordname <- c("cordx", "cordy")
  colnames(test) <- c(bname, sname, pre, errorname, Pillainame,
                      Fvalue_Pillainame, P_valuename, xxname, yyname, cordname)

  time <- proc.time()-start
  print(time)

  out <- list(
    formula = formula, test = test, SSCP_E = SSCP_E, sigmahat = sigmahat,
    epp = epp, epp2 = epp2, enp = enp, enp2 = enp2,
    df2 = df2, AIC = AIC, AICc = AICc,
    call = match.call()
  )
  class(out) <- "gwmmr_model"
  return(out)
}
