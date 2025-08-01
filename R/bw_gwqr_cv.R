#' Cross-Validation Bandwidth Selection for Geographically Weighted Quantile Regression (GWQR)
#'
#' This function selects the optimal bandwidth for a Geographically Weighted Quantile Regression (GWQR)
#' by minimizing the cross-validation (CV) score using a golden section search algorithm.
#'
#' @param formula A formula object specifying the model structure (response ~ predictors).
#' @param tau The quantile level to estimate. Default is 0.5 (median regression).
#' @param cordxy A matrix of spatial coordinates, with columns representing X and Y.
#' @param distmethod A character string indicating the distance calculation method to use.
#' Supported options are \code{"euclidean"} (default), \code{"greatcircle"} (great-circle distance), and \code{"manhattan"} (city-block distance).
#' @param kernel A character string specifying the kernel function used to compute spatial weights.
#' Options are \code{"global"}, \code{"gaussian"}, \code{"exponential"}, \code{"bisquare"}, \code{"tricube"}, and \code{"boxcar"}.
#' @param adaptive Logical. If \code{TRUE} (default), uses adaptive bandwidth. If \code{FALSE}, uses fixed bandwidth.
#' @param data A data frame containing all variables used in the model.
#' @param locallin Logical. If \code{TRUE}, includes local linear terms (interactions with spatial coordinates).
#' @param bw_min Minimum bandwidth value for the cross-validation search.
#' If \code{NA} (default), it is automatically set to the 10th percentile of non-zero distances (for fixed bandwidth) or 10% of the number of observations (for adaptive bandwidth).
#' @param bw_max Maximum bandwidth value for the cross-validation search.
#' If \code{NA} (default), it is automatically set to the 90th percentile of non-zero distances (for fixed bandwidth) or 90% of the number of observations (for adaptive bandwidth).
#' @param parallel Logical. If \code{TRUE} (default), enables parallel computation during the cross-validation procedure to improve speed.
#' @param core Integer. The number of CPU cores to use when \code{parallel = TRUE}.
#' If \code{NULL}, the number of cores is set to the total number of detected cores minus 2.
#'
#' @return A numeric value representing the optimal bandwidth.
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
#' @references
#' Chen, V. Y.-J., Deng, W.-S., Yang, T.-C., & Matthews, S. A. (2012).
#' Geographically weighted quantile regression (GWQR): An application to U.S. mortality data.
#' *Geographical Analysis*, 44(2), 134–150.
#'
#' @importFrom stats model.frame model.matrix terms coef
#' @importFrom quantreg rq
#' @importFrom parallel detectCores makeCluster clusterExport clusterEvalQ clusterSplit clusterApply stopCluster
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
#' bw
#' }
#' @export
bw_gwqr_cv <- function(formula, tau = 0.5, cordxy, distmethod = "euclidean",
                       kernel, adaptive = TRUE, data, locallin = FALSE,
                       bw_min = NA, bw_max = NA, parallel = TRUE, core = NULL) {
  start <- proc.time()

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
    if (bw_min > bw_max)
      stop("Error: bw_min cannot be greater than bw_max")
  } else if (!is.na(bw_min) && is.na(bw_max)) {
    bw_max <- auto_bw_max
    if (bw_min > bw_max)
      stop("Error: The computed bw_max is smaller than the provided bw_min.")
  } else if (is.na(bw_min) && !is.na(bw_max)) {
    bw_min <- auto_bw_min
    if (bw_min > bw_max)
      stop("Error: The computed bw_min is larger than the provided bw_max.")
  } else {
    bw_min <- auto_bw_min
    bw_max <- auto_bw_max
  }

  # ---- CV function ----
  cvh <- function(h) {
    # ---- Prepare data ----
    mat <- model.frame(formula, data = data)
    y <- mat[, 1]
    x <- as.matrix(model.matrix(formula, data = data))
    if (attr(terms(formula), "intercept")) x <- x[, -1] else x <- x
    nr <- length(y)
    xx <- cbind(rep(1, nr), x)
    nvar <- ncol(xx)
    cordx <- cordxy[, 1]
    cordy <- cordxy[, 2]

    # ---- estimate function ----
    estimate <- function(times) {
      beta <- c()
      muhat <- c()
      for (i in times) {
        d <- gwdist(cordxy, method = distmethod, target_idx = i)
        # leave-one-out
        cid <- which(d != 0)
        dd <- d[cid]

        w <- kernel_weights(dd, h, nr - 1, kernel, adaptive)

        regdata <- data.frame(data[cid, ], w)
        # for local linear estimation
        if (locallin == TRUE) {
          xstart <- cordx[i]
          ystart <- cordy[i]

          xd <- cordx - xstart
          yd <- cordy - ystart

          xdz <- xd * x
          ydz <- yd * x

          allx <- cbind(x, xd, yd, xdz, ydz)
          callx <- allx[cid, ]

          regdata <- data.frame(callx)
        }

        fit <- rq(formula, tau = tau, weights = w,
                  method = "fn", data = regdata)
        fit1 <- summary(fit,se="nid")
        bmat <- coef(fit1)[1:nvar,1]
        #bmat <- fit$coefficients[1:nvar]
        beta <- rbind(beta, bmat)

        yhat <- xx[i, ] %*% as.matrix(bmat[1:nvar])
        muhat <- rbind(muhat, yhat)
      }
      output <- list(beta = beta, muhat = muhat)
      return(output)
    }

    # ---- parallel ----
    if (parallel == TRUE) {
      if (is.null(core)) {
        core <- detectCores() - 2
      }

      cl <- makeCluster(core)
      clusterExport(cl, c(
        "y", "nr", "cordxy", "cordx", "cordy",
        "formula", "xx", "x", "nvar", "kernel", "tau",
        "kernel_weights", "gwdist", "data", "bw_min",
        "bw_max", "h", "adaptive", "locallin"
      ),
      envir = environment()
      )
      clusterEvalQ(cl, c(library(quantreg), library(MASS)))
      sp <- clusterSplit(cl, 1:nr)

      res <- clusterApply(cl = cl, sp, fun = estimate)
      # str(res)
      on.exit(stopCluster(cl), add = TRUE)
    } else {
      res <- lapply(split(1:nr, 1:nr), estimate)
    }

    beta <- c()
    muhat <- c()
    for (l in 1:length(res)) {
      beta <- rbind(beta, res[[l]]$beta)
      muhat <- rbind(muhat, res[[l]]$muhat)
    }

    error <- y - muhat
    wsum <- sum((error * (tau - ifelse(error < 0, 1, 0))))
    print(sprintf("bandwidth: %.4f, wsum value: %.4f", h, wsum))
    out <- list(y = y, x = x, beta = beta, muhat = muhat,
                error = error, wsum = wsum)
    return(out)
  }

  # /*Golden Section Search*/
  eps <- 1
  r <- (sqrt(5) - 1) / 2
  a0 <- bw_min
  b0 <- bw_max
  x1 <- r * a0 + (1 - r) * b0
  x2 <- (1 - r) * a0 + r * b0
  fx1 <- cvh(x1)$wsum
  # print(cbind(fx1))
  fx2 <- cvh(x2)$wsum

  it <- 0
  while ((b0 - a0) > eps) {
    it <- it + 1
    if (fx1 < fx2) {
      b0 <- x2
      x2 <- x1
      fx2 <- fx1
      x1 <- r * a0 + (1 - r) * b0
      fx1 <- cvh(x1)$wsum
    } else {
      a0 <- x1
      x1 <- x2
      fx1 <- fx2
      x2 <- (1 - r) * a0 + r * b0
      fx2 <- cvh(x2)$wsum
    }
  }
  h <- ifelse(fx1 <= fx2, a0, b0) #* 將fx1小於fx2的b0改成a0;
  # print(paste(time1,Sys.time()))
  # print(paste(a0,b0))
  # print(cvh(h)$wsum)
  bw.time <- proc.time() - start
  # print(bw.time)
  names(h) <- c("lambda")
  return(h)
}
