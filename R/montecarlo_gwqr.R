#' Monte Carlo Test for Spatial Variability in Geographically Weighted Quantile Regression (GWQR)
#'
#' This function performs a Monte Carlo significance test for spatial variability of local coefficients in the Geographically Weighted Quantile Regression (GWQR) Model.
#' This implementation also supports parallel computing to enhance computational efficiency.
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
#' @param n_permutations Integer. Number of Monte Carlo permutations to perform. Default is 99.
#' @param parallel Logical. If \code{TRUE} (default), enables parallel computation.
#' @param core Integer. The number of CPU cores to use when \code{parallel = TRUE}.
#' If \code{NULL}, the number of cores is set to the total number of detected cores minus 2.
#'
#' @return
#' A named list with one element: \code{p_values}, a named vector of p-values for each coefficient,
#' indicating whether its spatial variation is statistically significant.
#'
#' @details
#' The test statistic is the spatial variance of each coefficient.
#' By randomly permuting spatial coordinates, the spatial structure is destroyed and a null distribution is generated.
#' The p-value is the proportion of permuted variances exceeding the observed variance.
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
#' bw<- bw_gwqr_cv(
#'   formula = formula, tau = 0.5,
#'   cordxy = cordxy, distmethod = "euclidean",
#'   kernel = "gaussian", adaptive = TRUE,
#'   data = boston, locallin = FALSE,
#'   parallel = TRUE, core = NULL)
#'
#' result <- montecarlo_gwqr(
#'   formula = formula,
#'   tau = 0.5,
#'   cordxy = cordxy,
#'   distmethod = "euclidean",
#'   h = bw,
#'   kernel = "gaussian",
#'   adaptive = TRUE,
#'   data = boston,
#'   locallin = FALSE,
#'   parallel = TRUE,
#'   core = 3)
#' result
#' }
#' @export
montecarlo_gwqr <- function(formula, tau = 0.5, cordxy, distmethod ="euclidean",
                            h, kernel, adaptive = TRUE, data, locallin = FALSE,
                            n_permutations = 99, parallel = TRUE, core = NULL){
  start_time <- proc.time()

  if (!is.matrix(cordxy)) {
    stop("cordxy must be a matrix.")
  }

  # ------------ Step1: GWQR variability -------------
  gwqr_result <- gwqr(formula = formula, tau = tau,
                      cordxy = cordxy, distmethod = distmethod,
                      h = h, kernel = kernel, adaptive = adaptive,
                      data = data, locallin = locallin,
                      parallel = parallel, core = core)
  beta_j_obs <- gwqr_result$beta
  nu_j_obs <- apply(beta_j_obs, 2, var)

  var_names <- colnames(beta_j_obs)

  nu_j_permuted <- matrix(NA, n_permutations, length(nu_j_obs))

  # ------------ Step2~4: Monte Carlo Test ------------
  pb <- progress::progress_bar$new(total = n_permutations)

  for (i in 1:n_permutations) {
    idx <- sample(1:nrow(data), replace = FALSE)
    permuted_coords <- cordxy[idx, ]

    gwqr_permuted <- gwqr(formula = formula, tau = tau,
                          cordxy = permuted_coords, distmethod = distmethod,
                          h = h, kernel = kernel, adaptive = adaptive,
                          data = data, locallin = locallin,
                          parallel = parallel, core = core)

    beta_j_perm <- gwqr_permuted$beta
    nu_j_permuted[i, ] <- apply(beta_j_perm, 2, var)

    pb$tick()
  }

  # Step 5 : Compute p-value
  p_values <- setNames(sapply(1:length(nu_j_obs), function(j) {
    (sum(nu_j_permuted[ ,j] >= nu_j_obs[j])) / n_permutations
  }), var_names)

  monte_time <- proc.time()-start_time
  print(monte_time)

  out <- list(
    p_values = p_values
  )

  return(out)
}
