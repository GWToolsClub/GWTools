#' Cross-Validation Bandwidth Selection for GWOLR
#'
#' This function selects the optimal bandwidth for a Geographically Weighted Ordinal Logistic Regression (GWOLR) model using cross-validation,
#' based on the method proposed by Dong, Nakaya, and Brunsdon (2020). This implementation also supports parallel computing to enhance computational efficiency.
#'
#' @param formula A model formula specifying the ordinal response and predictor variables.
#' @param cordxy A matrix of spatial coordinates, with columns representing X and Y.
#' @param distmethod A character string indicating the distance calculation method to use.
#' Supported options are \code{"euclidean"} (default), \code{"greatcircle"} (great-circle distance), and \code{"manhattan"} (city-block distance).
#' @param kernel A character string specifying the kernel function used to compute spatial weights.
#' Options are \code{"global"}, \code{"gaussian"}, \code{"exponential"}, \code{"bisquare"}, \code{"tricube"}, and \code{"boxcar"}.
#' @param adaptive Logical. If \code{TRUE} (default), uses adaptive bandwidth. If \code{FALSE}, uses fixed bandwidth.
#' @param data A data frame containing all variables used in the model.
#' @param link A character string specifying the link function for the ordinal logistic regression. Options are \code{"probit"} (default) and \code{"logit"}.
#' @param fixed_vars Optional. A character vector of predictor variable names to be treated as fixed (i.e., non-spatially varying) in the GWOLR model.
#' Default is \code{NULL}, in which case all predictors are treated as spatially varying.
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
#' This function builds upon the implementation by Dong, Nakaya, and Brunsdon (2020), which uses a golden-section search to minimize leave-one-out cross-validation (CV) error.
#' It introduces several practical improvements and extensions over the original version:
#'
#' \itemize{
#'    \item \strong{Parallel Computing}: Supports parallel computation using the \code{parallel} package to accelerate the CV process, especially for large datasets. The number of cores can be set via \code{core}.
#'    \item \strong{Automatic Distance Matrix Calculation}: Automatically computes the distance matrix using the provided coordinates and selected \code{distmethod}.
#'    \item \strong{Flexible Distance Metrics}: Adds support for \code{"manhattan"} and \code{"greatcircle"} distances, in addition to the original Euclidean metric.
#'    \item \strong{Flexible Bandwidth Ranges}: Uses a more flexible default bandwidth search range, and allows user-specified \code{bw_min} and \code{bw_max} values.
#'    \item \strong{Expanded Kernel Options}: Supports \code{"global"}, \code{"gaussian"}, \code{"exponential"}, \code{"bisquare"}, \code{"tricube"}, and \code{"boxcar"} kernels.
#' }
#'
#' These enhancements make the function more flexible, robust, and suitable for modern spatial modeling workflows in R.
#'
#' @references
#' Dong, G., Nakaya, T., & Brunsdon, C. (2020). Geographically weighted regression models for ordinal categorical response variables: An application to geo-referenced life satisfaction data. \emph{Computers, Environment and Urban Systems}, 80, 101428.
#'
#' @importFrom stats model.frame model.response model.matrix terms lm coefficients runif pnorm plogis
#' @importFrom MASS polr
#' @importFrom maxLik maxLik
#' @importFrom parallel detectCores makeCluster clusterExport clusterEvalQ clusterSplit clusterApply stopCluster
#'
#' @examples
#' # This example is adapted from the demo script provided in the appendix of Dong et al. (2020).
#' \dontrun{
#' #if (require("gstat")){
#'   # Simulate spatial coordinates
#'   cordx<- seq(100, 199, 10)
#'   cordy <- seq(200, 299, 10)
#'   coords <- expand.grid(cordx, cordy)
#'   cordxy <- as.matrix(coords)
#'   n.total <- nrow(coords) # 100
#'
#'   g.dummy <- gstat(formula = z ~ 1, locations = ~Var1 + Var2, dummy = TRUE,
#'                     beta = 0, model = vgm(psill = 1, model = "Exp", range = 50), nmax = 20)
#'   sim <- predict(g.dummy, newdata = coords, nsim = 2)
#'
#'   # Generate spatially varying coefficients
#'   coe_x1 <- sim$sim1 - mean(sim$sim1)
#'   coe_x2 <- sim$sim2 - mean(sim$sim2)
#'   coe_mat <- cbind(coe_x1, coe_x2)
#'
#'   # Generate independent variables
#'   X <- matrix(runif(2 * n.total, min = 1, max = 10), n.total)
#'   colnames(X) <- c("x1", "x2")
#'
#'   # Generate latent variable and ordinal response
#'   eta <- rowSums(X * coe_mat) + rlogis(n.total)
#'   cut_points <- as.numeric(quantile(eta))
#'   y <- as.integer(cut(eta, breaks = cut_points, include.lowest = TRUE))
#'
#'   # Assemble data
#'   data <- data.frame(x1 = X[, 1], x2 = X[, 2], y = factor(y))
#'   formula <- y ~ x1 + x2
#'
#'   bw <- bw_gwolr_cv(formula = formula, cordxy = cordxy, kernel = "bisquare",
#'                     adaptive = FALSE, data = data,
#'                     link = "logit", fixed_vars = NULL,
#'                     parallel = TRUE, core = 4)
#' }
#' }
#' @export
bw_gwolr_cv <- function(formula, cordxy, distmethod = "euclidean", kernel,
                        adaptive = TRUE, data,
                        link = "probit", fixed_vars = NULL,
                        bw_min = NA, bw_max = NA,
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

    #w <- rep(0, length(d))
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
    if (bw_min > bw_max) stop("Error: bw_min cannot
                              be greater than bw_max")
  } else if (!is.na(bw_min) && is.na(bw_max)) {
    bw_max <- auto_bw_max
    if (bw_min > bw_max) stop("Error: The computed bw_max is
                              smaller than the provided bw_min.")
  } else if (is.na(bw_min) && !is.na(bw_max)) {
    bw_min <- auto_bw_min
    if (bw_min > bw_max) stop("Error: The computed bw_min is
                              larger than the provided bw_max.")
  } else {
    bw_min <- auto_bw_min
    bw_max <- auto_bw_max
  }

  # ---- cv function ----
  cvh <- function(h) {
    # ---- Prepare data ----
    mat <- model.frame(formula, data = data)
    y <- model.response(mat)
    x <- as.matrix(model.matrix(formula, data = data))
    if (attr(terms(formula), "intercept")) x <- x[, -1] else x <- x
    nr <- length(y)
    # xx <- cbind(rep(1, nr), x)
    # nvar <- ncol(x)

    # ---- start value ----
    global_model <- try(suppressWarnings(polr(formula, data = data,
                                              method = link)), silent = TRUE)
    if (inherits(global_model, "try-error")) {
      # using simple regression coefficients
      ols_model <- lm(as.numeric(y) ~ x - 1)
      start_betas <- drop(coefficients(ols_model))
      start_zeta <- runif(nlevels(y) - 1)
      start_zeta <- cumsum(start_zeta)
      start_values <- c(start_betas, start_zeta)
    } else {
      start_betas <- drop(global_model$coefficients)
      start_zeta <- drop(global_model$zeta)
      start_theta <- drop(c(start_zeta[1], log(diff(start_zeta))))
      start_values <- drop(c(start_betas, start_theta))
    }

    if (is.null(fixed_vars)) {
      index_fixed <- NA
    } else {
      index_fixed <- match(fixed_vars, colnames(x))
    }

    estimate <- function(times) {
      beta <- c()
      theta <- c()
      eta <- c()
      cv <- c()
      for (i in times) {
        d <- gwdist(cordxy, method = distmethod, target_idx = i)
        cid <- which(d != 0)
        dd <- d[cid]

        w <- kernel_weights(dd, h, nr - 1, kernel, adaptive)
        fit <- try(orderedfit(y = y[cid], X = x[cid, ], link = link,
                              Gweights = w, id.fixed.vars = index_fixed,
                              start.values = start_values), silent = TRUE)

        if (!inherits(fit, "try-error")) {
          bmat <- fit$Betas
          zeta <- fit$zeta
          eta <- x[i, ] %*% as.matrix(bmat)

          if (link == "probit") {
            cum <- pnorm(zeta - eta)
          } else {
            cum <- plogis(zeta - eta)
          }

          prob <- c(cum[1], diff(cum), 1 - max(cum))
          cv <- c(cv, 1 - prob[y[i]])
        } else {
          cv <- Inf
          break
        }
      }
      output <- list(beta = beta, theta = theta, eta = eta, cv = cv)
      return(output)
    }

    # ---- parallel ----
    if (parallel == TRUE) {
      if (is.null(core)) {
        core <- detectCores() - 2
      }

      cl <- makeCluster(core)
      clusterExport(cl, c(
        "y", "nr", "cordxy", "formula", "x", "kernel", "link", "dist_mat",
        "kernel_weights", "gwdist", "data", "bw_min", "bw_max", "h", "adaptive",
        "orderedfit"
      ), envir = environment())
      clusterEvalQ(cl, c(library(MASS), library(maxLik)))
      sp <- clusterSplit(cl, 1:nr)

      res <- clusterApply(cl = cl, sp, fun = estimate)
      # str(res)
      on.exit(stopCluster(cl), add = TRUE)
    } else {
      res <- lapply(split(1:nr, 1:nr), estimate)
    }

    beta <- c()
    theta <- c()
    eta <- c()
    cv <- c()
    for (l in seq_along(res)) {
      beta <- rbind(beta, res[[l]]$beta)
      theta <- rbind(theta, res[[l]]$theta)
      eta <- rbind(eta, res[[l]]$eta)
      cv <- c(cv, res[[l]]$cv)
    }

    if (any(is.infinite(cv))) {
      print(sprintf("bandwidth: %.4f, CV score: Inf", h))
      return(list(cv = Inf))
    }

    cv_score <- sum(cv^2)
    print(sprintf("bandwidth: %.4f, CV score: %.4f", h, cv_score))
    out <- list(y = y, x = x, beta = beta, theta = theta,
                eta = eta, cv_score = cv_score)
    return(out)
  }

  # /*Golden Section Search*/
  eps <- 1
  r <- (sqrt(5) - 1) / 2
  a0 <- bw_min
  b0 <- bw_max
  x1 <- r * a0 + (1 - r) * b0
  x2 <- (1 - r) * a0 + r * b0
  fx1 <- cvh(x1)$cv
  fx2 <- cvh(x2)$cv

  it <- 0
  while ((b0 - a0) > eps) {
    it <- it + 1
    if (fx1 < fx2) {
      b0 <- x2
      x2 <- x1
      fx2 <- fx1
      x1 <- r * a0 + (1 - r) * b0
      fx1 <- cvh(x1)$cv
    } else {
      a0 <- x1
      x1 <- x2
      fx1 <- fx2
      x2 <- (1 - r) * a0 + r * b0
      fx2 <- cvh(x2)$cv
    }
  }
  h <- ifelse(fx1 <= fx2, a0, b0)
  #bw_time <- proc.time() - start
  #print(bw_time)
  return(h)
}
