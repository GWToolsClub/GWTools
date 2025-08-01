#' Simplified Estimation for Geographically Weighted Multivariate Regression (GWMMR)
#'
#' This function performs a lightweight estimation for Geographically Weighted Multivariate Multiple Regression (GWMMMR) models.
#' It returns local regression coefficients, fitted values, and a Frobenius-norm-based variability statistic for each coefficient,
#' without computing AIC, or hypothesis testing. This function is designed for internal reuse, such as in permutation
#' or Monte Carlo procedures.
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
#' @return A list containing:
#' \itemize{
#'   \item \code{beta}: A matrix of local regression coefficients.
#'   \item \code{pred}: A matrix of fitted values for each observation and response.
#'   \item \code{residual}: A matrix of residuals, defined as observed minus fitted values.
#'   \item \code{stat}: A vector of Frobenius norms calculated from the local covariance matrices of selected pairs of coefficients. Used as a summary statistic in Monte Carlo testing.
#' }
#'
#' @details
#' This function omits computationally intensive components such as leverage matrix, AIC, and local hypothesis testing.
#' It is particularly useful for scenarios requiring repeated model fitting, such as resampling-based significance testing.
#'
#' @references
#' Chen, V. Y.-J., Yang, T.-C., & Jian, H.-L. (2022).
#' Geographically Weighted Regression Modeling for Multiple Outcomes.
#' *Annals of the American Association of Geographers*, 112(1), 1–18.
#'
#' @importFrom MASS ginv
#' @importFrom stats cov model.frame model.response model.matrix
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
#' gwmmr_simple_result <- gwmmr_simple(
#'   formula = formula,
#'   cordxy = coords,
#'   h = 200,
#'   kernel = "bisquare",
#'   adaptive = TRUE,
#'   data = data
#' )
#' gwmmr_simple_result
#' @export
gwmmr_simple <- function(formula, cordxy, distmethod = "euclidean",
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
  pred <- array(0, dim = c(ny, nc))

  for (i in 1:n) {
    d <- gwdist(cordxy, method = distmethod, target_idx = i)
    w <- kernel_weights(d, h, nr, kernel, adaptive)

    newx <- X * sqrt(w) # n × p
    newy <- Y * sqrt(w) # n × m
    xpx <- t(newx) %*% newx # p × p
    xpxi <- ginv(xpx) # inverse crossproducts
    b <- xpxi %*% t(newx) %*% newy # p × m
    beta[, i] <- b
    yhat <- X[i, ] %*% b # 1 × m
    pred[i, ] <- yhat
  }
  beta <- t(beta)
  residual <- Y - pred
  stat <- array(0, dim = c(nv, 1))

  for (k in 1:nv) {
    covb <- cov(beta[, c(k, k + nv)])
    stat[k, ] <- norm(covb, type = "F")
  }

  output <- list(beta = beta, pred = pred, residual = residual, stat = stat)
  return(output)
}
