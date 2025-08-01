#' Monte Carlo Test for Spatial Variablility in GWGLM
#'
#' This function performs a Monte Carlo test to examine whether the local coefficients
#' in a Geographically Weighted Generalized Linear Model (GWGLM) vary significantly across different locations.
#'
#' The function first computes the observed variance of each local coefficient across space.
#' Then, it repeatedly permutes the spatial coordinates and refits the model to simulate a null distribution
#' of coefficient variances under spatial randomness.
#' For each coefficient, a \code{p}-value is calculated by comparing the observed variance
#' to the distribution of permuted variances.
#'
#' @param formula A formula object specifying the model structure (response ~ predictors).
#' @param family A character string specifying the distribution family. Supported options are \code{"gaussian"}, \code{"poisson"}, and \code{"binomial"}.
#' @param cordxy A matrix of spatial coordinates, with columns representing X and Y.
#' @param distmethod A character string indicating the distance calculation method to use.
#' Supported options are \code{"euclidean"} (default), \code{"greatcircle"} (great-circle distance), and \code{"manhattan"} (city-block distance).
#' @param h A numeric value specifying the bandwidth. This can represent a fixed distance or a number of nearest neighbors (if \code{adaptive = TRUE}).
#' @param kernel A character string specifying the kernel function used to compute spatial weights.
#' Options are \code{"global"}, \code{"gaussian"}, \code{"exponential"}, \code{"bisquare"}, \code{"tricube"}, and \code{"boxcar"}.
#' @param adaptive Logical. If \code{TRUE} (default), uses adaptive bandwidth. If \code{FALSE}, uses fixed bandwidth.
#' @param data A data frame containing all variables used in the model.
#' @param n_permutations Number of Monte Carlo permutations. Default is 99.
#' @param parallel Logical. If \code{TRUE} (default), enables parallel computation to speed up local estimation. Default is \code{TRUE}.
#' @param core Integer. The number of CPU cores to use when \code{parallel = TRUE}.
#' If \code{NULL}, the number of cores is set to the total number of detected cores minus 2.
#'
#' @return A list with one element:
#' \describe{
#'    \item{\code{p_values}}{A named vector of p-values for each coefficient, indicating whether its spatial variation is statistically significant.}
#' }
#'
#' @details
#' In addition to standard binary response variables of numeric \code{0/1} type,
#' this function and the associated GWGLM series uniquely support proportion-type responses
#' expressed as success counts and total trials (e.g., \code{y/n}).
#'
#' When success/trial matrices are provided, the function automatically computes the success proportions internally,
#' eliminating the need for manual preprocessing.
#'
#' This feature broadens the flexibility of modeling binary and proportion data,
#' addressing a common limitation in existing geographically weighted modeling packages.
#'
#'
#' @examples
#' \dontrun{
#' data(tokyo)
#'
#' tokyo$lnoff = log(tokyo$Exp_2564)
#'
#' cordxy <- cbind(tokyo$X, tokyo$Y)
#' formula <- Mort2564~Professl+OwnHome+Elderly+Unemply+offset(lnoff)
#' family <- "poisson"
#'
#' bw_cv <- bw_gwglm_cv(
#'   formula = formula,
#'   family = family,
#'   cordxy = cordxy,
#'   kernel = "gaussian",
#'   adaptive = FALSE,
#'   data = tokyo,
#'   parallel = TRUE
#' )
#'
#' montecarlo_result <- montecarlo_gwglm(
#'   formula = formula,
#'   family = family,
#'   cordxy = cordxy,
#'   h = bw_cv,
#'   kernel = "gaussian",
#'   adaptive= FALSE,
#'   data = tokyo,
#'   n_permutations = 99,
#'   parallel = T,
#'   core = 4
#' )
#' montecarlo_result
#' }
#' @export
montecarlo_gwglm <- function(formula, family, cordxy, distmethod="euclidean",
                             h, kernel, adaptive, data,
                             n_permutations = 99, parallel = TRUE, core = NULL) {
  start_time <- proc.time()

  if (!is.matrix(cordxy)) {
    stop("cordxy must be a matrix.")
  }

  cvar <- 0

  # ------------ Step1: GWGLM variability -------------
  gwglm_result <- gwglm(formula = formula, family = family,
                        cordxy = cordxy, distmethod = distmethod,
                        h = h, kernel = kernel, adaptive = adaptive,
                        data = data, cvar = cvar,
                        parallel = parallel, core = core)
  beta_j_obs <- gwglm_result$beta
  nu_j_obs <- apply(beta_j_obs, 2, var)

  var_names <- colnames(beta_j_obs)

  nu_j_permuted <- matrix(NA, n_permutations, length(nu_j_obs))

  # ------------ Step2~4: Monte Carlo Test ------------
  pb <- progress::progress_bar$new(total = n_permutations)

  for (i in 1:n_permutations) {
    idx <- sample(1:nrow(data), replace = FALSE)
    permuted_coords <- cordxy[idx, ]

    gwglm_permuted <- gwglm(formula = formula, family = family,
                          cordxy = permuted_coords, distmethod = distmethod,
                          h = h, kernel = kernel, adaptive = adaptive,
                          data = data, cvar = cvar,
                          parallel = parallel, core = core)

    beta_j_perm <- gwglm_permuted$beta
    nu_j_permuted[i, ] <- apply(beta_j_perm, 2, var)

    pb$tick()
  }

  # Step 5 : Compute p-value
  p_values <- setNames(sapply(1:length(nu_j_obs), function(j) {
    (sum(nu_j_permuted[ ,j] >= nu_j_obs[j])) / (n_permutations + 1)
  }), var_names)

  monte_time <- proc.time()-start_time
  print(monte_time)

  out <- list(
    p_values = p_values
  )

  return(out)
}
