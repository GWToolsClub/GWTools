#' Bandwidth Selection via AICc for Geographically Weighted Multivariate Multiple Regression (GWMMMR)
#'
#' This function selects the optimal bandwidth for a Geographically Weighted Multivariate Multiple Regression (GWMMMR) model
#' by minimizing the multivariate AICc using golden section search algorithm.
#'
#' @param formula A formula object specifying the model structure (e.g., \code{cbind(y1, y2) ~ x1 + x2 + x3}).
#' @param cordxy A matrix of spatial coordinates, with columns representing X and Y.
#' @param distmethod A character string indicating the distance calculation method to use.
#' Supported options are \code{"euclidean"} (default), \code{"greatcircle"} (great-circle distance), and \code{"manhattan"} (city-block distance).
#' @param kernel A character string specifying the kernel function used to compute spatial weights.
#' Options are \code{"global"}, \code{"gaussian"}, \code{"exponential"}, \code{"bisquare"}, \code{"tricube"}, and \code{"boxcar"}.
#' @param adaptive Logical. If \code{TRUE} (default), uses adaptive bandwidth. If \code{FALSE}, uses fixed bandwidth.
#' @param data A data frame containing all variables used in the model.
#' @param bw_min Minimum bandwidth value for the AICc-based search.
#' If \code{NA} (default), it is automatically set to the 10th percentile of non-zero distances (for fixed bandwidth) or 10% of the number of observations (for adaptive bandwidth).
#' @param bw_max Maximum bandwidth value for the AICc-based search.
#' If \code{NA} (default), it is automatically set to the 90th percentile of non-zero distances (for fixed bandwidth) or 90% of the number of observations (for adaptive bandwidth).
#'
#' @return A numeric value representing the optimal bandwidth that minimizes the multivariate AICc.
#'
#' @details
#' This function selects the optimal bandwidth for the Geographically Weighted Multivariate Multiple Regression (GWMMR) model
#' by minimizing a multivariate AICc, as proposed in Chen et al. (2022). The GWMMR model extends traditional GWR
#' to handle multiple correlated outcomes, allowing for spatially varying associations across multiple response variables.
#'
#' @references
#' Chen, V. Y.-J., Yang, T.-C., & Jian, H.-L. (2022).
#' Geographically Weighted Regression Modeling for Multiple Outcomes.
#' *Annals of the American Association of Geographers*, 112(1), 1–18.
#'
#' @importFrom stats model.matrix model.frame model.response
#' @importFrom MASS ginv
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
#' bw <- bw_gwmmr_aic(
#'   formula = formula,
#'   cordxy = coords,
#'   kernel = "bisquare",
#'   adaptive = TRUE,
#'   data = data
#' )
#' @export
bw_gwmmr_aic <- function(formula, cordxy, distmethod = "euclidean",
                         kernel, adaptive = TRUE, data,
                         bw_min = NA, bw_max = NA) {
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

  # ---- distance matrix ----
  dist_mat <- gwdist(cordxy, method = distmethod)

  all_dist <- sort(dist_mat)[which(sort(dist_mat) != 0)]
  all_dist <- all_dist[!duplicated(all_dist)]

  if (adaptive == TRUE) {
    auto_bw_min <- floor(0.1 * nrow(data))
    auto_bw_max <- ceiling(0.9 * nrow(data))
  } else {
    auto_bw_min <- all_dist[max(1, floor(length(all_dist) * 0.1))]
    auto_bw_max <- all_dist[ceiling(length(all_dist) * 0.9)]
  }

  if (!is.na(bw_min) && !is.na(bw_max)) {
    if (bw_min > bw_max) {
      stop("Error: bw_min cannot be greater than bw_max")
    }
  } else if (!is.na(bw_min) && is.na(bw_max)) {
    bw_max <- auto_bw_max
    if (bw_min > bw_max) {
      stop("Error: The computed bw_max is smaller than the provided bw_min.")
    }
  } else if (is.na(bw_min) && !is.na(bw_max)) {
    bw_min <- auto_bw_min
    if (bw_min > bw_max) {
      stop("Error: The computed bw_min is larger than the provided bw_max.")
    }
  } else {
    bw_min <- auto_bw_min
    bw_max <- auto_bw_max
  }

  # ---- Prepare data ----
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)
  ox <- model.matrix(formula, data = data)[, -1, drop = FALSE]
  nr <- nrow(ox)
  X <- as.matrix(cbind(rep(1, nr), ox))
  nv <- ncol(X)
  ny <- nrow(y)
  nc <- ncol(y)
  Y <- array(y, dim = c(ny, nc))
  beta <- array(0, dim = c(nv * nc, nr))
  # strr <- array(0, dim = c(nr, nv))
  pred <- array(0, dim = c(ny, nc))
  leverage <- array(0, dim = c(nr, 1))
  # leverage2 <- array(0, dim = c(nr, 1))
  error <- array(y, dim = c(ny, nc))

  aich <- function(h) {
    for (i in 1:nr) {
      d <- gwdist(cordxy, method = distmethod, target_idx = i)
      w <- kernel_weights(d, h, nr, kernel, adaptive)

      newx <- X * sqrt(w)
      newy <- Y * sqrt(w)
      xpx <- t(newx) %*% newx
      xpxi <- ginv(xpx) # inverse crossproducts
      b <- xpxi %*% (t(newx) %*% newy)
      beta[, i] <- b
      yhat <- X[i, ] %*% b
      pred[i, ] <- yhat
      hatmatrix <- newx %*% xpxi %*% t(newx)
      leverage[i, ] <- diag(hatmatrix)[i] # S
    }

    error <- Y - pred
    SSCP_E <- t(error) %*% error
    epp <- sum(leverage) * nc # trace(S)*y個數 (effective number of psrameters)
    # sigmahat=SSCP.E/(nr-epp)

    sigmahat <- SSCP_E / nr
    k <- epp / nc
    AIC <- nr * log(det(sigmahat)) + nr * nc * (log(2 * 3.14159) + 1) + 2 * k * nc + (nc * (nc + 1))
    AICc <- nr * log(det(sigmahat)) + nr * nc * log(2 * 3.14159) + (nr * nc * (nr + k)) / (nr - k - nc - 1)
    print(sprintf("bandwidth: %.4f, AICc value: %.4f", h, AICc))
    return(AICc)
  }

  ## Golden Section Search##
  eps <- 1
  r <- (sqrt(5) - 1) / 2
  a0 <- bw_min
  b0 <- bw_max
  x1 <- r * a0 + (1 - r) * b0
  x2 <- (1 - r) * a0 + r * b0
  fx1 <- aich(x1)
  # print (fx1)
  fx2 <- aich(x2)
  # print (fx2)
  it <- 0
  while (b0 - a0 > eps) {
    it <- it + 1
    if (fx1 < fx2) {
      b0 <- x2
      x2 <- x1
      x1 <- r * a0 + (1 - r) * b0
      fx2 <- fx1
      fx1 <- aich(x1)
    } else {
      a0 <- x1
      x1 <- x2
      x2 <- (1 - r) * a0 + r * b0
      fx1 <- fx2
      fx2 <- aich(x2)
    }
    #print(cbind( it, a0 ,fx1 ,b0 ,fx2))
  }

  h <- ifelse(fx1 <= fx2, a0, b0)
  # print(it)
  bw_time <- proc.time() - start
  # print(bw_time)
  return(h)
}
