#' Geographically Weighted Ordinal Logistic Regression (GWOLR)
#'
#' Estimates a Geographically Weighted Ordinal Logistic Regression (GWOLR) model.
#' The implementation is based on Dong, Nakaya, and Brunsdon (2020), and has been extended to support parallel computing for large datasets.
#'
#' @param formula A model formula specifying the ordinal response and predictor variables.
#' @param cordxy A matrix of spatial coordinates, with columns representing X and Y.
#' @param distmethod A character string indicating the distance calculation method to use.
#' Supported options are \code{"euclidean"} (default), \code{"greatcircle"} (great-circle distance), and \code{"manhattan"} (city-block distance).
#' @param h Bandwidth value for the spatial weighting kernel.
#' @param kernel A character string specifying the kernel function used to compute spatial weights.
#' Options are \code{"global"}, \code{"gaussian"}, \code{"exponential"}, \code{"bisquare"}, \code{"tricube"}, and \code{"boxcar"}.
#' @param adaptive Logical. If \code{TRUE} (default), uses adaptive bandwidth. If \code{FALSE}, uses fixed bandwidth.
#' @param data A data frame containing all variables used in the model.
#' @param link A character string specifying the link function for the ordinal logistic regression. Options are \code{"probit"} (default) and \code{"logit"}.
#' @param fixed_vars Optional. A character vector of predictor variable names to be treated as fixed (i.e., non-spatially varying) in the GWOLR model.
#' Default is \code{NULL}, in which case all predictors are treated as spatially varying.
#' @param parallel Logical. If \code{TRUE} (default), enables parallel computation during the estimation procedure to improve speed.
#' @param core Integer. The number of CPU cores to use when \code{parallel = TRUE}.
#' If \code{NULL}, the number of cores is set to the total number of detected cores minus 2.
#'
#' @return An object of class \code{"gwolr_model"}.
#' The default print method displays the call, model formula, and a summary of the local coefficients.
#' Additional components are accessible via list indexing:
#'
#' \itemize{
#'    \item{\code{call}: The original function call.}
#'    \item{\code{formula}: The model formula.}
#'    \item{\code{beta}: A matrix of locally estimated regression coefficients (one row per location).}
#'    \item{\code{std}: A matrix of standard errors for the local coefficients.}
#'    \item{\code{tvalue}: A matrix of pseudo t-values.}
#'    \item{\code{cov}: A list of local covariance matrices.}
#' }
#'
#' @details
#' This function builds upon the original implementation by Dong, Nakaya, and Brunsdon (2020) for fitting Geographically Weighted Ordinal Logistic Regression (GWOLR) models.
#' It introduces several practical improvements and extensions over the original version:
#'
#' \itemize{
#'    \item \strong{Parallel Computing}: Supports parallel computation using the \code{parallel} package to accelerate the estimation process, especially for large datasets. The number of cores can be set via \code{core}.
#'    \item \strong{Automatic Distance Matrix Calculation}: Automatically computes the distance matrix using the provided coordinates and selected \code{distmethod}.
#'    \item \strong{Flexible Distance Metrics}: Adds support for \code{"manhattan"} and \code{"greatcircle"} distances, in addition to the original Euclidean metric.
#'    \item \strong{Expanded Kernel Options}: Supports \code{"global"}, \code{"gaussian"}, \code{"exponential"}, \code{"bisquare"}, \code{"tricube"}, and \code{"boxcar"} kernels.
#' }
#'
#' These enhancements make the function more flexible, robust, and suitable for modern spatial modeling workflows in R.
#'
#' @references
#' Dong, G., Nakaya, T., & Brunsdon, C. (2020). Geographically weighted regression models for ordinal categorical response variables: An application to geo-referenced life satisfaction data. \emph{Computers, Environment and Urban Systems}, 80, 101428.
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
#'                    beta = 0, model = vgm(psill = 1, model = "Exp",
#'                    range = 50), nmax = 20)
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
#'
#'   est <- gwolr(formula = formula, cordxy = cordxy, h = bw,
#'                 kernel = "bisquare", adaptive = FALSE, data = data,
#'                 link = "logit", fixed_vars = NULL, parallel = TRUE, core = 4)
#'
#'   print(est)
#' }
#' }
#' @export
gwolr <- function(formula, cordxy, distmethod = "euclidean", h,
                  kernel, adaptive = TRUE, data,
                  link = "probit", fixed_vars = NULL,
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

  # ---- Prepare data ----
  mat <- model.frame(formula, data = data)
  y <- model.response(mat)
  x <- as.matrix(model.matrix(formula, data = data))
  if (attr(terms(formula), "intercept")) x <- x[, -1] else x <- x
  nr <- length(y)
  # xx <- cbind(rep(1, nr), x)
  # p <- ncol(x)
  # q <- nlevels(y) - 1

  # ---- start value ----
  global_model <- try(suppressWarnings(MASS::polr(formula, data = data,
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

  # ----- fixed regression coefficients ----
  if (is.null(fixed_vars)) {
    index_fixed <- NA
  } else {
    index_fixed <- match(fixed_vars, colnames(x))
  }

  estimate <- function(times) {
    beta <- c()
    std <- c()
    tvalue <- c()
    cov <- list()
    for (i in times) {
      d <- gwdist(cordxy, method = distmethod, target_idx = i)
      w <- kernel_weights(d, h, nr, kernel, adaptive)

      fit <- try(orderedfit(y, x, link, Gweights = w,
                            id.fixed.vars = index_fixed,
                            start.values = start_values), silent = TRUE)

      if (!inherits(fit, "try-error")) {
        bmat <- fit$coefficients
        semat <- fit$std.err
        tmat <- fit$t.value
        covmat <- fit$V.g

        beta <- rbind(beta, bmat)
        std <- rbind(std, semat)
        tvalue <- rbind(tvalue, tmat)
        cov <- append(cov, list(covmat))
      } else {
        cat("The", i, "-th observation without estimates", "\n")
      }
    }
    output <- list(beta = beta, std = std, tvalue = tvalue, cov = cov)
  }

  # ---- parallel ----
  if (parallel) {
    if (is.null(core)) {
      core <- parallel::detectCores() - 2
    }

    cl <- parallel::makeCluster(core)
    parallel::clusterExport(cl, c(
      "y", "nr", "cordxy", "formula", "x", "kernel",
      "link", "kernel_weights", "gwdist", "data", "h",
      "adaptive", "orderedfit"
    ), envir = environment())
    parallel::clusterEvalQ(cl, c(library(MASS), library(maxLik)))
    sp <- parallel::clusterSplit(cl, 1:nr)

    res <- parallel::clusterApply(cl = cl, sp, fun = estimate)
    # str(res)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  } else {
    res <- lapply(split(1:nr, 1:nr), estimate)
  }


  beta <- c()
  std <- c()
  tvalue <- c()
  cov <- list()
  for (l in seq_along(res)) {
    beta <- rbind(beta, res[[l]]$beta)
    std <- rbind(std, res[[l]]$std)
    tvalue <- rbind(tvalue, res[[l]]$tvalue)
    cov <- c(cov, res[[l]]$cov)
  }

  #time <- proc.time() - start
  #print(time)
  out <- list(
    formula = formula, beta = beta, std = std, tvalue = tvalue, cov = cov,
    call = match.call()
  )
  class(out) <- "gwolr_model"
  return(out)
}
