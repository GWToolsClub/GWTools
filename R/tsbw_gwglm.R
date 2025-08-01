#' Two-Stage Bandwidth Selection for Geographically Weighted Generalized Linear Model (GWGLM)
#'
#' This function implements a two-stage bandwidth selection procedure for GWGLM, where the bandwidth at each stage is selected
#' based on either cross-validation (CV) or corrected Akaike Information Criterion (AICc).
#' The function first estimates a global model and selects the optimal bandwidth using either cross-validation (CV) or AICc.
#' The resulting coefficients are used to construct an offset term for the second stage, where the local bandwidth is chosen.
#'
#' @param formula A formula object specifying the model structure (response ~ predictors).
#' @param family A character string specifying the distribution family. Supported options are \code{"gaussian"}, \code{"poisson"}, and \code{"binomial"}.
#' @param cordxy A matrix of spatial coordinates, with columns representing X and Y.
#' @param distmethod A character string indicating the distance calculation method to use.
#' Supported options are \code{"euclidean"} (default), \code{"greatcircle"} (great-circle distance), and \code{"manhattan"} (city-block distance).
#' @param kernel A character string specifying the kernel function used to compute spatial weights.
#' Options are \code{"global"}, \code{"gaussian"}, \code{"exponential"}, \code{"bisquare"}, \code{"tricube"}, and \code{"boxcar"}.
#' @param adaptive Logical. If \code{TRUE} (default), uses adaptive bandwidth. If \code{FALSE}, uses fixed bandwidth.
#' @param criterion A character string specifying the bandwidth selection criterion. Options are \code{"CV"} (cross-validation) and \code{"AIC"} (AICc).
#' @param offset A character string giving the name of a column in the data to be used as offset (optional). If not provided, only the internal offset from fixed coefficients is used.
#' @param local_var A character vector of variable names to be treated as spatially varying (local) coefficients.
#' @param global_var A character vector of variable names to be treated as spatially invariant (global) coefficients.
#' @param data A data frame containing all variables used in the model.
#' @param bw_min Minimum bandwidth value for both stages.
#' If \code{NA} (default), it is automatically set to the 10th percentile of non-zero distances (for fixed bandwidth) or 10% of the number of observations (for adaptive bandwidth).
#' @param bw_max Maximum bandwidth value for both stages.
#' If \code{NA} (default), it is automatically set to the 90th percentile of non-zero distances (for fixed bandwidth) or 90% of the number of observations (for adaptive bandwidth).
#' @param parallel Logical. If \code{TRUE} (default), enables parallel computation during bandwidth search.
#' @param core Integer. The number of CPU cores to use when \code{parallel = TRUE}.
#' If \code{NULL}, the number of cores is set to the total number of detected cores minus 2.
#'
#' @return
#' A list with two elements:
#' \itemize{
#'   \item \code{h1}: The selected bandwidth for the first stage.
#'   \item \code{h2}: The selected bandwidth for the second stage.
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
#' @references
#' Li, D. and Mei, C. (2018). A two-stage estimation method with bootstrap inference for semi-parametric geographically weighted generalized linear models.
#' International Journal of Geographical Information Science, 32(9):1860–1883.
#'
#' Zheng, G.-T. (2024). \emph{An application of semi-parametric geographically weighted logistic regression in real estate transaction data analysis}.
#' Master's thesis, Department of Statistics, National Chengchi University, Taiwan.
#'
#' @importFrom stats as.formula
#' @importFrom stringr str_count
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
#' offset <- "lnoff"
#' local_var <- c("Professl","Unemply")
#' global_var <- c("(Intercept)","OwnHome","Elderly")
#'
#' ts_bw_cv <- tsbw_gwglm(
#'   formula = formula,
#'   family = family,
#'   cordxy = cordxy,
#'   kernel = "gaussian",
#'   adaptive = FALSE,
#'   criterion = "CV",
#'   offset = offset,
#'   local_var = local_var,
#'   global_var = global_var,
#'   data = tokyo,
#'   parallel = TRUE
#' )
#' ts_bw_cv
#'
#' ts_bw_aic <- tsbw_gwglm(
#'   formula = formula,
#'   family = family,
#'   cordxy = cordxy,
#'   kernel = "gaussian",
#'   adaptive = FALSE,
#'   criterion = "AIC",
#'   offset = offset,
#'   local_var = local_var,
#'   global_var = global_var,
#'   data = tokyo,
#'   parallel = TRUE
#' )
#' ts_bw_aic
#' }
#' @export
tsbw_gwglm <- function(formula, family, cordxy, distmethod = "euclidean",
                       kernel, adaptive = TRUE,
                       criterion, offset = NULL,
                       local_var, global_var, data,
                       bw_min = NA, bw_max = NA,
                       parallel = TRUE, core = NULL) {
  start <- proc.time()

  if (!is.matrix(cordxy)) {
    stop("cordxy must be a matrix.")
  }

  # first stage
  # first bandwidth
  cvar <- 0
  print("============ First stage ============")
  if (criterion == "CV") {
    bw_global <- bw_gwglm_cv(
      formula = formula, family = family,
      cordxy = cordxy, distmethod = distmethod,
      kernel = kernel, adaptive = adaptive, data = data,
      bw_min = bw_min, bw_max = bw_max,
      parallel = parallel, core = core
    )
  } else if (criterion == "AIC") {
    bw_global <- bw_gwglm_aic(
      formula = formula, family = family,
      cordxy = cordxy, distmethod = distmethod,
      kernel = kernel, adaptive = adaptive,
      cvar = cvar, data = data,
      bw_min = bw_min, bw_max = bw_max,
      parallel = parallel, core = core
    )
  } else {
    stop("Unsupported criterion. Please use 'CV' or 'AIC'.")
  }
  h1 <- bw_global

  # constant coefficients
  fit_global <- gwglm(
    formula = formula, family = family, cordxy = cordxy,
    distmethod = distmethod, h = h1,
    kernel = kernel, adaptive = adaptive, data = data,
    cvar = cvar, parallel = FALSE, core = NULL
  )
  constant_beta <- apply(as.data.frame(fit_global$beta[, global_var]), 2, mean)
  names(constant_beta) <- c(global_var)

  # second stage
  # second bandwidth
  if (length(global_var) == 0) {
    stop("Please enter the global variables")
  } else if (length(local_var) == 0) {
    stop("Please enter the local variables")
  } else {
    new_global_var <- c()
    for (i in seq_along(global_var)) {
      if (global_var[i] == "(Intercept)") {
        new_global_var[i] <- "Intercept"
      } else {
        new_global_var[i] <- paste0(global_var[i], "_C")
      }
    }
    newvar_df <- as.data.frame(matrix(
      nrow = nrow(data),
      ncol = length(global_var)
    ))
    colnames(newvar_df) <- new_global_var
    for (i in seq_along(global_var)) {
      if (new_global_var[i] == "Intercept") {
        newvar_df[, i] <- rep(
          constant_beta[which(names(constant_beta) == global_var[i])],
          nrow(data)
        )
      } else if (length(global_var) == 1) {
        newvar_df[, i] <- data[, global_var[i]] * constant_beta
      } else {
        newvar_df[, i] <- data[, global_var[i]] * constant_beta[which(names(constant_beta) == global_var[i])]
      }
    }
    data2 <- cbind(data, newvar_df) # 計算XC並合併至data
    if (sum(str_count(new_global_var, pattern = "Intercept")) != 0) {
      if (length(offset) == 0) {
        formula2 <- as.formula(paste(
          formula[2], formula[1], "-1+",
          paste(local_var, collapse = "+"),
          "+offset(",
          paste(new_global_var, collapse = "+"), ")"
        ))
      } else {
        formula2 <- as.formula(paste(
          formula[2], formula[1], "-1+",
          paste(local_var, collapse = "+"),
          "+offset(", offset, "+",
          paste(new_global_var, collapse = "+"), ")"
        ))
      }
    } else {
      if (length(offset) == 0) {
        formula2 <- as.formula(paste(
          formula[2], formula[1],
          paste(local_var[which(local_var != "(Intercept)")],
            collapse = "+"
          ), "+offset(",
          paste(new_global_var, collapse = "+"), ")"
        ))
      } else {
        formula2 <- as.formula(paste(
          formula[2], formula[1],
          paste(local_var[which(local_var != "(Intercept)")],
            collapse = "+"
          ),
          "+offset(", offset, "+",
          paste(new_global_var, collapse = "+"), ")"
        ))
      }
    }
    cvar <- length(global_var)
    print("============ Second stage ============")
    if (criterion == "CV") {
      bw_local <- bw_gwglm_cv(
        formula = formula2, family = family,
        cordxy = cordxy, distmethod = distmethod,
        kernel = kernel, adaptive = adaptive, data = data2,
        bw_min = bw_min, bw_max = bw_max,
        parallel = parallel, core = core
      )
    } else if (criterion == "AIC") {
      bw_local <- bw_gwglm_aic(
        formula = formula2, family = family,
        cordxy = cordxy, distmethod = distmethod,
        kernel = kernel, adaptive = adaptive,
        cvar = cvar, data = data2,
        bw_min = bw_min, bw_max = bw_max,
        parallel = parallel, core = core
      )
    } else {
      stop("Unsupported criterion. Please use 'CV' or 'AIC'.")
    }
    h2 <- bw_local
  }
  bw_time <- proc.time() - start
  #print(bw_time)
  out <- list(h1 = h1, h2 = h2)
  return(out)
}
