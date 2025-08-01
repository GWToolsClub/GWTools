#' Monte Carlo Test fot Spatial Variability for GWMMR Coefficients
#'
#' This function performs a Monte Carlo significance test for spatial variability of local coefficients in the Geographically Weighted Multivariate Regression (GWMMR) Model.
#' It compares observed coefficient variability (measured using the Frobenius norm) with a distribution obtained by randomly permuting spatial coordinates.
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
#' @param n_permutations Number of Monte Carlo permutations. Default is 99.
#'
#' @return A list with one element:
#' \itemize{
#'    \item \code{p_values}: A vector of p-values for each coefficient.
#' }
#'
#' @details
#' The test statistic is the Frobenius norm of the covariance matrix of local coefficients (for each variable).
#' By randomly permuting spatial coordinates, the spatial structure is destroyed and a null distribution is generated.
#' The p-value for each coefficient is computed as the proportion of permuted norms greater than the observed norm.
#'
#' @references
#' Chen, V. Y.-J., Yang, T.-C., & Jian, H.-L. (2022).
#' Geographically Weighted Regression Modeling for Multiple Outcomes.
#' *Annals of the American Association of Geographers*, 112(1), 1–18.
#'
#' @examples
#' \dontrun{
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
#' monte_result <- montecarlo_gwmmr(
#'   formula = formula,
#'   cordxy = coords,
#'   h = 200,
#'   kernel = "bisquare",
#'   adaptive = TRUE,
#'   data = data
#' )
#' monte_result
#' }
#' @export
montecarlo_gwmmr <- function(formula, cordxy, distmethod = "euclidean",
                             h, kernel, adaptive = TRUE,
                             data, n_permutations = 99){
  start <- proc.time()

  if (!is.matrix(cordxy)) {
    stop("cordxy must be a matrix.")
  }

  mf <- model.frame(formula, data = data)
  ox <- model.matrix(formula, data = data)[, -1, drop = FALSE]
  param_names <- c("(Intercept)", colnames(ox))
  nv <- ncol(as.matrix(ox))+1

  frob_stat_perm <- array(0, dim = c(nv, n_permutations))
  p_values <- array(0, dim = c(nv, 1))

  # ------------ Step1: GWMMR variability -------------
  gwmmr_result <- gwmmr_simple(formula = formula, cordxy = cordxy,
                             distmethod = distmethod, h = h,
                             kernel = kernel, adaptive = adaptive,
                             data = data)

  frob_stat_obs <- gwmmr_result$stat

  # ------------ Step2~4: Monte Carlo Test ------------
  pb <- progress::progress_bar$new(total = n_permutations)

  for (i in 1:n_permutations) {
    idx <- sample(1:nrow(data), replace = FALSE)
    permuted_coords <- cordxy[idx, ]

    gwmmr_permuted <- gwmmr_simple(formula = formula, cordxy = permuted_coords,
                                   distmethod = distmethod, h = h,
                                   kernel = kernel, adaptive = adaptive,
                                   data = data)
    frob_stat_perm[, i] <- gwmmr_permuted$stat
    pb$tick()
  }
  for (j in 1:nv) {
    r <- length(which(frob_stat_perm[j, ] > frob_stat_obs[j, ]))
    p_values[j, ] <- r / (n_permutations+1)
  }

  monte_time <- proc.time() - start
  print(monte_time)

  p_values_named <- data.frame(
    variable = param_names,
    p_value = as.vector(p_values)
  )

  out <- list(
    p_values = p_values_named
  )

  return(out)
}
