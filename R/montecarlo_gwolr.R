#' Monte Carlo Test for GWOLR
#'
#' This function performs a Monte Carlo significance test for spatial variability of local coefficients in the Geographically Weighted Ordinal Logistic Regression (GWOLR) Model,
#' based on the method proposed by Dong, Nakaya, and Brunsdon (2020). This implementation also supports parallel computing to enhance computational efficiency.
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
#' @param link A character string specifying the link function for the ordinal regression. Options are \code{"probit"} (default) and \code{"logit"}.
#' @param fixed_vars Optional. A character vector of predictor variable names to be treated as fixed (i.e., non-spatially varying) in the GWOLR model.
#' Default is \code{NULL}, in which case all predictors are treated as spatially varying.
#' @param n_permutations Integer. Number of Monte Carlo permutations to perform. Default is 99.
#' @param parallel Logical. If \code{TRUE} (default), enables parallel computation during the permutation procedure to improve speed.
#' @param core Integer. The number of CPU cores to use when \code{parallel = TRUE}.
#' If \code{NULL}, the number of cores is set to the total number of detected cores minus 2.
#'
#' @return A data frame with p-values for each local coefficient (excluding fixed variables, if specified), representing the statistical significance of spatial variability.
#'
#' @details
#' This function follows the Monte Carlo testing procedure proposed by Dong, Nakaya, and Brunsdon (2020) for assessing spatial variability in Geographically Weighted Ordinal Logistic Regression (GWOLR) Model.
#'
#' In particular, it adopts their approach for handling fixed (global) variables
#' — any variable specified in \code{fixed_vars} is excluded from the Monte Carlo test and treated as spatially stationary.
#'
#' Additionally, this implementation introduces parallel computation via the \code{parallel} and \code{core} arguments
#' to significantly reduce computation time when conducting a large number of permutations.
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
#'                    beta = 0, model = vgm(psill = 1, model = "Exp", range = 50), nmax = 20)
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
#'    monte <- montecarlo_gwolr(formula = formula, cordxy = cordxy, h = bw,
#'                               kernel = "bisquare", adaptive = FALSE,
#'                               data = data, link = "logit", fixed_vars = NULL,
#'                               parallel = TRUE, core = 4)
#' }
#' }
#' @export
montecarlo_gwolr <- function(formula, cordxy, distmethod = "euclidean", h,
                             kernel, adaptive = TRUE, data, link = "probit",
                             fixed_vars = NULL, n_permutations = 99,
                             parallel = TRUE, core = NULL) {
  start <- proc.time()

  if (!is.matrix(cordxy)) {
    stop("cordxy must be a matrix.")
  }

  # ---- Prepare data ----
  x <- as.matrix(model.matrix(formula, data = data))
  if (attr(terms(formula), "intercept")) x <- x[, -1] else x <- x

  # ----- fixed regression coefficients ----
  if (is.null(fixed_vars)) {
    index_fixed <- NA
  } else {
    index_fixed <- match(fixed_vars, colnames(x))
  }

  # ------------ Step1: GWOLR variability -------------
  gwolr_result <- gwolr(
    formula = formula, cordxy = cordxy,
    distmethod = distmethod, h = h, kernel = kernel,
    adaptive = adaptive, data = data,
    link = link, fixed_vars = fixed_vars,
    parallel = parallel, core = core
  )

  if (all(is.na(index_fixed))) {
    beta_j_obs <- gwolr_result$beta
  } else {
    beta_j_obs <- gwolr_result$beta[, -index_fixed]
  }
  nu_j_obs <- apply(beta_j_obs, 2, var)

  var_names <- colnames(beta_j_obs)

  # ------------ Step2~4: Monte Carlo Test ------------
  pb <- progress::progress_bar$new(total = n_permutations)
  nu_j_permuted <- matrix(NA, n_permutations, length(nu_j_obs))

  for (i in 1:n_permutations) {
    idx <- sample(1:nrow(data), replace = FALSE)
    permuted_coords <- cordxy[idx, ]

    gwolr_permuted <- gwolr(
      formula = formula, cordxy = permuted_coords,
      distmethod = distmethod, h = h, kernel = kernel,
      adaptive = adaptive, data = data,
      link = link, fixed_vars = fixed_vars,
      parallel = parallel, core = core
    )

    if (all(is.na(index_fixed))) {
      beta_j_perm <- gwolr_permuted$beta
    } else {
      beta_j_perm <- gwolr_permuted$beta[, -index_fixed]
    }
    nu_j_permuted[i, ] <- apply(beta_j_perm, 2, var)
    #cat("run", i, "now", "\n")
    pb$tick()
  }

  # Step 5 : Compute p-value
  p_values <- sapply(seq_along(nu_j_obs), function(j) {
    sum(nu_j_permuted[, j] >= nu_j_obs[j]) / (n_permutations + 1)
  })

  p_value_df <- as.data.frame(p_values)
  rownames(p_value_df) <- var_names
  colnames(p_value_df) <- "p-value"

  monte_time <- proc.time() - start
  print(monte_time)

  return(p_value_df)
}
