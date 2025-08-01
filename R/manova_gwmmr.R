#' Local MANOVA Test for Geographically Weighted Multivariate Multiple Regression (GWMMMR)
#'
#' This function performs local multivariate hypothesis testing for each spatial location in a
#' Geographically Weighted Multivariate Multiple Regression (GWMMMR) model. It tests user-specified
#' linear hypotheses of the form \code{H0: LBM = C} using Pillai's trace statistic, and returns local
#' test results including Pillai's trace, F-statistics, and p-values at each location.
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
#' @param L A numeric matrix specifying linear constraints on the coefficient matrix \code{B}. Each row defines a linear combination to be tested.
#' @param M A numeric matrix specifying linear combinations among response variables.
#' @param v A numeric scalar representing the effective number of parameters. Use the value of \code{enp2} returned from \code{\link{gwmmr}}.
#'
#' @return A list with class \code{"gwmmr_manova"} containing:
#' \itemize{
#'    \item \code{test}: A data frame with columns Pillai, Fvalue_Pillai, and P_value.
#'    \item \code{E}: Within-group variation matrix.
#'    \item \code{H}: Between-group variation matrix.
#'    \item \code{df1}, \code{df2}: Local degrees of freedom used in the F-test approximation of Pillai's trace.
#'    \item \code{call}: The matched call.
#' }
#'
#' @details
#' This function implements a localized MANOVA framework for GWMMR models, based on the multivariate general linear hypothesis (MGLH)
#' formulation proposed in Chen et al. (2022). At each spatial location, the function fits a local multivariate regression model and
#' tests the hypothesis \code{H0: L B M = C}, where \code{B} is the coefficient matrix, and \code{C} is typically assumed to be a zero matrix.
#' The Pillai trace statistic is computed and transformed into an approximate F distribution under multivariate normality assumptions.
#'
#'
#' @references
#' Chen, V. Y.-J., Yang, T.-C., & Jian, H.-L. (2022).
#' Geographically Weighted Regression Modeling for Multiple Outcomes.
#' *Annals of the American Association of Geographers*, 112(1), 1–18.
#'
#' @importFrom stats pf model.frame model.response model.matrix
#' @importFrom Matrix rankMatrix
#' @importFrom MASS ginv
#' @importFrom matlib tr
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
#'
#' manova_result <- manova_gwmmr(
#'   formula = formula,
#'   cordxy = coords,
#'   h = 200,
#'   kernel = "bisquare",
#'   adaptive = TRUE,
#'   data = data,
#'   L = array(c(0,1,0,0),dim=c(1,4)),
#'   M = diag(2),
#'   v = gwmmr_result$enp
#' )
#' manova_result
#' @export
manova_gwmmr <- function(formula, cordxy, distmethod = "euclidean",
                         h, kernel, adaptive = TRUE, data,
                         L, M, v) {
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
  X <- as.matrix(cbind(rep(1, nr), ox))
  nv <- ncol(X)
  ny <- nrow(y)
  nc <- ncol(y)
  Y <- array(y, dim = c(ny, nc))
  F_value <- array(0, dim = c(nr, 1))
  stat <- array(0, dim = c(nr, 1))
  P_value <- array(0, dim = c(nr, 1))

  for (i in 1:n) {
    d <- gwdist(cordxy, method = distmethod, target_idx = i)
    w <- kernel_weights(d, h, nr, kernel, adaptive)

    newx <- X * sqrt(w)
    newy <- Y * sqrt(w)
    xpx <- t(newx) %*% newx
    xpxi <- ginv(xpx) # inverse crossproducts
    b <- xpxi %*% (t(newx) %*% newy) # p x m

    # test
    E <- t(M) %*% (t(newy) %*% newy - (t(b) %*% (t(newx) %*% newx) %*% b)) %*% M # 組內變異 b x b

    H <- t(M) %*% t(L %*% b) %*% ginv(L %*% xpxi %*% t(L)) %*% (L %*% b) %*% M # 組間變異 b x b
    HEi <- H %*% ginv(E)
    LAMDA <- eigen(HEi)$value
    Pillai <- tr(H %*% ginv(H + E))

    stat[i, ] <- Pillai

    # F-value
    a <- nrow(L)
    q <- ncol(L)
    p <- nrow(M)
    b <- ncol(M)
    r <- rankMatrix(X)[1]

    s <- min(a, b)
    mm <- 0.5 * (abs(a - b) - 1)
    nn <- 0.5 * (v - b - 1)

    PF <- ((2 * nn + s + 1) / (2 * mm + s + 1)) * (Pillai / (s - Pillai)) # Pillai F_value
    F_value[i, ] <- PF
    P_value[i, ] <- 1 - pf(PF, df1 = s * (2 * mm + s + 1), df2 = s * (2 * nn + s + 1)) # Pillai P_value
  }
  df1 <- rep(s * (2 * mm + s + 1), n)
  df2 <- rep(s * (2 * nn + s + 1), n)
  test <- data.frame(stat, F_value, P_value)

  # name
  statname <- "Pillai"
  fname <- "Fvalue_Pillai"
  P_valuename <- "P_value"
  colnames(test) <- c(statname, fname, P_valuename)

  output <- list(test = test, E = E, H = H,
                 df1 = s * (2 * mm + s + 1), df2 = s * (2 * nn + s + 1),
                 call = match.call())
  class(output) <- "gwmmr_manova"
  return(output)
}

