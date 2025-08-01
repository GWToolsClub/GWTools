#' Geographically Weighted Quantile Regression (GWQR)
#'
#' This function fits a Geographically Weighted Quantile Regression (GWQR) model
#' using kernel-weighted local quantile regression at each spatial location.
#' It supports both local constant and local linear estimation, and parallel computing for faster execution.
#'
#' @param formula A formula object specifying the model structure (response ~ predictors).
#' @param tau The quantile level to estimate. Default is 0.5 (median regression).
#' @param cordxy A matrix of spatial coordinates, with columns representing X and Y.
#' @param distmethod A character string indicating the distance calculation method to use.
#' Supported options are \code{"euclidean"} (default), \code{"greatcircle"} (great-circle distance), and \code{"manhattan"} (city-block distance).
#' @param h A numeric value specifying the bandwidth. This can represent a fixed distance or a number of nearest neighbors (if \code{adaptive = TRUE}).
#' @param kernel A character string specifying the kernel function used to compute spatial weights.
#' Options are \code{"global"}, \code{"gaussian"}, \code{"exponential"}, \code{"bisquare"}, \code{"tricube"}, and \code{"boxcar"}.
#' @param adaptive Logical. If \code{TRUE} (default), uses adaptive bandwidth. If \code{FALSE}, uses fixed bandwidth.
#' @param data A data frame containing all variables used in the model.
#' @param locallin Logical. If \code{TRUE}, includes local linear terms (interactions with spatial coordinates).
#' @param parallel Logical. If \code{TRUE} (default), enables parallel computation for estimation.
#' @param core Integer. The number of CPU cores to use when \code{parallel = TRUE}.
#' If \code{NULL}, the number of cores is set to the total number of detected cores minus 2.
#'
#' @return
#' A list of class \code{"gwqr_model"} containing:
#' \itemize{
#'   \item \code{formula}: The model formula.
#'   \item \code{beta}: Matrix of local coefficient estimates.
#'   \item \code{std}: Matrix of standard errors for each coefficient.
#'   \item \code{tvalue}: Matrix of local t-statistics.
#'   \item \code{muhat}: Vector of fitted values.
#'   \item \code{error}: Vector of residuals.
#'   \item \code{wsum}: The overall check loss.
#'   \item \code{call}: The matched function call.
#' }
#'
#' @details
#' This function supports two estimation types:
#'
#' **1. Local Constant (default)**:
#' Estimates a separate set of coefficients at each location, assuming constant relationships within the local neighborhood.
#'
#' **2. Local Linear (set `locallin = TRUE`)**:
#' Optionally includes interactions between predictors and spatial coordinates.
#' This allows the model to better capture spatial gradients and improve robustness near boundaries or in sparse areas.
#'
#' @importFrom quantreg rq
#' @importFrom parallel detectCores makeCluster clusterExport clusterEvalQ clusterSplit clusterApply stopCluster
#'
#' @references
#' Chen, V. Y.-J., Deng, W.-S., Yang, T.-C., & Matthews, S. A. (2012).
#' Geographically weighted quantile regression (GWQR): An application to U.S. mortality data.
#' *Geographical Analysis*, 44(2), 134–150.
#'
#' @examples
#' \dontrun{
#' data(boston, package = "spData")
#' boston <- boston.c
#' cordxy <- cbind(boston$LON, boston$LAT)
#'
#' formula <- CMEDV ~ LSTAT + RM
#'
#' bw <- bw_gwqr_cv(
#'   formula = formula, tau = 0.5,
#'   cordxy = cordxy, distmethod = "euclidean",
#'   kernel = "gaussian", adaptive = TRUE,
#'   data = boston, locallin = FALSE,
#'   parallel = TRUE, core = NULL)
#'
#' result <- gwqr(
#'   formula = formula,
#'   tau = 0.5,
#'   cordxy = cordxy,
#'   distmethod = "euclidean",
#'   h = bw,
#'   kernel = "gaussian",
#'   adaptive = TRUE,
#'   data = boston,
#'    locallin = FALSE,
#'    parallel = TRUE,
#'    core = NULL
#' )
#' result
#' }
#' @export
gwqr <- function(formula, tau = 0.5, cordxy, distmethod = "euclidean",
                 h, kernel, adaptive = TRUE, data, locallin = FALSE,
                 parallel = TRUE, core = NULL) {
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
      w <- as.numeric(w)
    } else if (kernel == "exponential") {
      w <- exp(-d / bw)
      w <- as.numeric(w) # /*exponential kernel*/
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
  mat <- model.frame(formula, data = data)
  y <- mat[, 1]
  x <- as.matrix(model.matrix(formula, data = data)) # designed x matrix including intercept (1's)
  ox <- x[, -1] # designed x matrix excluding intercept (1's)
  nvar <- ncol(x)
  nr <- length(y)
  cordx <- cordxy[, 1]
  cordy <- cordxy[, 2]


  # ---- estimate function ----
  estimate <- function(times) {
    beta <- c()
    muhat <- c()
    std <- c()
    tvalue <- c()
    for (i in times) {
      d <- gwdist(cordxy, method = distmethod, target_idx = i)
      w <- kernel_weights(d, h, nr, kernel, adaptive)

      regdata <- data.frame(ox)

      # for local linear estimation
      if (locallin == TRUE) {
        xstart <- cordx[i]
        ystart <- cordy[i]

        xd <- cordx - xstart
        yd <- cordy - ystart

        xdz <- xd * ox
        ydz <- yd * ox

        allx <- cbind(ox, xd, yd, xdz, ydz)

        regdata <- data.frame(allx)
      }
      fit <- quantreg::rq( y ~ ., tau = tau, weights = w,
                           method = "fn", data = regdata)
      fit1 <- quantreg::summary.rq(fit, se = "nid")
      bmat <- coef(fit1)[1:nvar, 1]
      semat <- coef(fit1)[1:nvar, 2]
      tmat <- coef(fit1)[1:nvar, 3]

      std <- rbind(std, semat)
      tvalue <- rbind(tvalue, tmat)
      yhat <- x[i, ] %*% as.matrix(bmat[1:nvar])
      beta <- rbind(beta, bmat)
      muhat <- rbind(muhat, yhat)
      # print(muhat)
    }
    # colnames(beta)[ncol(beta)]<- c("lambda")
    output <- list(beta = beta, std = std, muhat = muhat, tvalue = tvalue)
    return(output)
  }

  if (parallel == TRUE) {
    if (is.null(core)) {
      core <- parallel::detectCores() - 2
    }

    cl <- parallel::makeCluster(core)
    parallel::clusterExport(cl, c(
      "nr", "nvar", "tau", "x", "formula", "data", "y", "cordxy", "kernel", "h",
      "kernel_weights", "gwdist"
    ), envir = environment())
    parallel::clusterEvalQ(cl, c(library(quantreg), library(MASS)))
    sp <- parallel::clusterSplit(cl, 1:nr)

    res <- parallel::clusterApply(cl = cl, sp, fun = estimate)
    # str(res)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  } else {
    res <- lapply(split(1:nr, 1:nr), estimate)
  }

  beta <- c()
  muhat <- c()
  std <- c()
  tvalue <- c()
  for (l in 1:length(res)) {
    beta <- rbind(beta, res[[l]]$beta)
    std <- rbind(std, res[[l]]$std)
    tvalue <- rbind(tvalue, res[[l]]$tvalue)
    muhat <- rbind(muhat, res[[l]]$muhat)
  }

  error <- y - muhat
  # loss function
  wsum <- sum((error * (tau - ifelse(error < 0, 1, 0))))
  out <- list(
    formula = formula, beta = beta, std = std, muhat = muhat,
    tvalue = tvalue, error = error, wsum = wsum,
    call = match.call()
  )
  class(out) <- "gwqr_model"
  return(out)
}
