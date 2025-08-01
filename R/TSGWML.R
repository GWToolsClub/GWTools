#' Two-Stage Geographically Weighted Maximum Likelihood (TSGWML)
#'
#' Fits a Two-Stage Geographically Weighted Maximum Likelihood (TSGWML) model that separates
#' global (fixed) and local (spatially varying) variables into two sequential estimation stages.
#'
#' @param formula A formula object specifying the model structure (response ~ predictors).
#' @param family A character string specifying the distribution family. Supported options are \code{"gaussian"}, \code{"poisson"}, and \code{"binomial"}.
#' @param cordxy A matrix of spatial coordinates, with columns representing X and Y.
#' @param distmethod A character string indicating the distance calculation method to use.
#' Supported options are \code{"euclidean"} (default), \code{"greatcircle"} (great-circle distance), and \code{"manhattan"} (city-block distance).
#' @param kernel A character string specifying the kernel function used to compute spatial weights.
#' Options are \code{"global"}, \code{"gaussian"}, \code{"exponential"}, \code{"bisquare"}, \code{"tricube"}, and \code{"boxcar"}.
#' @param adaptive Logical. If \code{TRUE} (default), uses adaptive bandwidth. If \code{FALSE}, uses fixed bandwidth.
#' @param h1 Numeric. Bandwidth used in the first-stage global fitting.
#' @param h2 Numeric. Bandwidth used in the second-stage local fitting.
#' @param offset A character string giving the name of a column in the data to be used as offset (optional).
#' If not provided, only the internal offset from fixed coefficients is used.
#' @param local_var A character vector of variable names to be treated as spatially varying (local) coefficients.
#' @param global_var A character vector of variable names to be treated as spatially invariant (global) coefficients.
#' @param data A data frame containing all variables used in the model.
#' @param parallel Logical. If \code{TRUE} (default), enables parallel computing for faster estimation. Default is \code{TRUE}.
#' @param core Integer. The number of CPU cores to use when \code{parallel = TRUE}.
#' If \code{NULL}, the number of cores is set to the total number of detected cores minus 2.
#'
#' @return A list of class \code{"TSGWML_model"} containing:
#' \describe{
#'    \item{\code{formula}}{The formula used in the first stage.}
#'    \item{\code{muhat}}{Fitted values from the second stage model.}
#'    \item{\code{global_beta}}{A matrix of global (fixed) coefficients.}
#'    \item{\code{local_beta}}{A matrix of local (spatially varying) coefficients.}
#'    \item{\code{LRT_statistic}}{The log-likelihood ratio (LRT) statistic computed as the difference between the first-stage (global-only) and second-stage (local+global) models.}
#'    \item{\code{new_data}}{A data frame combining the original input data with estimated global coefficients, local coefficients, and corresponding standard errors, t-values, and p-values.}
#'    \item{\code{call}}{The matched function call.}
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
#' @references
#' Li, D. and Mei, C. (2018). A two-stage estimation method with bootstrap inference for semi-parametric geographically weighted generalized linear models.
#' International Journal of Geographical Information Science, 32(9):1860–1883.
#'
#' Zheng, G.-T. (2024). \emph{An application of semi-parametric geographically weighted logistic regression in real estate transaction data analysis}.
#' Master's thesis, Department of Statistics, National Chengchi University, Taiwan.
#'
#' @importFrom stats as.formula
#'
#' @examples
#' \dontrun{
#' data(tokyo)
#'
#' tokyo[,6:9] <- scale(tokyo[,6:9])
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
#' TSGWML_result <- TSGWML(
#'   formula = formula,
#'   family = family,
#'   cordxy = cordxy,
#'   kernel = "gaussian",
#'   adaptive = FALSE,
#'   h1 = ts_bw_cv$h1,
#'   h2 = ts_bw_cv$h2,
#'   offset = offset,
#'   local_var = local_var,
#'   global_var = global_var,
#'   data = tokyo,
#'   parallel = TRUE
#' )
#' TSGWML_result
#' }
#' @export
TSGWML <- function(formula, family, cordxy, distmethod = "euclidean",
                   kernel, adaptive = TRUE,
                   h1, h2, offset = NULL, local_var, global_var, data,
                   parallel = TRUE, core = NULL) {
  start <- proc.time()

  if (!is.matrix(cordxy)) {
    stop("cordxy must be a matrix.")
  }

  # first stage
  # constant coefficients
  cvar <- 0
  fit_global <- gwglm(
    formula = formula, family = family, cordxy = cordxy,
    distmethod = distmethod, h = h1,
    kernel = kernel, adaptive = adaptive,
    data = data, cvar = cvar, parallel = parallel, core = core
  )
  constant_beta <- apply(as.data.frame(fit_global$beta[, global_var]), 2, mean)
  names(constant_beta) <- c(global_var)
  # constant.std <- apply(as.data.frame(fit_global$std[,global_var]),2,mean)
  # constant.tvalue <- apply(as.data.frame(fit_global$tvalue[,global_var]),2,mean)
  # constant.pvalue <- apply(as.data.frame(fit_global$pvalue[,global_var]),2,mean)

  # second stage
  # second bandwidth
  if (length(global_var) == 0) {
    print("please enter the global variables")
    local_beta <- fit_global$beta
    local_beta_name <- c()
    for (i in 1:length(local_var)) {
      if (local_var[i] == "(Intercept)") {
        local_beta_name[i] <- "Intercept_beta"
      } else {
        local_beta_name[i] <- paste0(local_var[i], "_beta")
      }
    }
    colnames(local_beta) <- local_beta_name
    rownames(local_beta) <- 1:nrow(local_beta)
    output_data <- cbind(data, local_beta)

    TS_time <- proc.time() - start
    #print(TS_time)

    out <- list(
      formula = formula,
      muhat = fit_global$muhat,
      local_beta = local_beta,
      new_data = output_data
    )
    return(out)
  } else if (length(local_var) == 0) {
    print("please enter the local variables")
    global_beta <- as.data.frame(matrix(
      rep(constant_beta, each = nrow(data)),
      nrow(data), length(global_var)
    ))
    global_beta_name <- c()
    for (i in 1:length(global_var)) {
      if (global_var[i] == "(Intercept)") {
        global_beta_name[i] <- "Intercept_beta"
      } else {
        global_beta_name[i] <- paste0(global_var[i], "_beta")
      }
    }
    colnames(global_beta) <- global_beta_name
    output_data <- cbind(data, global_beta)

    TS_time <- proc.time() - start
    #print(TS_time)

    out <- list(
      formula = formula,
      muhat = fit_global$muhat,
      global_beta = constant_beta,
      new_data = output_data
    )
    return(out)
  } else {
    new_global_var <- c()
    for (i in 1:length(global_var)) {
      if (global_var[i] == "(Intercept)") {
        new_global_var[i] <- "Intercept"
      } else {
        new_global_var[i] <- paste0(global_var[i], "_C")
      }
    }
    newvar_df <- as.data.frame(matrix(nrow = nrow(data), ncol = length(global_var)))
    colnames(newvar_df) <- new_global_var
    for (i in 1:length(global_var)) {
      if (new_global_var[i] == "Intercept") {
        newvar_df[, i] <- rep(constant_beta[which(names(constant_beta) == global_var[i])], nrow(data))
      } else if (length(global_var) == 1) {
        newvar_df[, i] <- data[, global_var[i]] * constant_beta
      } else {
        newvar_df[, i] <- data[, global_var[i]] * constant_beta[which(names(constant_beta) == global_var[i])]
      }
    }
    data2 <- cbind(data, newvar_df) # 計算XC並合併至data
    if (sum(stringr::str_count(new_global_var, pattern = "Intercept")) != 0) {
      if (is.null(offset)) {
        formula2 <- as.formula(paste(
          formula[2], formula[1], "-1+", paste(local_var, collapse = "+"),
          "+offset(", paste(new_global_var, collapse = "+"), ")"
        ))
      } else {
        formula2 <- as.formula(paste(
          formula[2], formula[1], "-1+", paste(local_var, collapse = "+"),
          "+offset(", offset, "+", paste(new_global_var, collapse = "+"), ")"
        ))
      }
    } else {
      if (is.null(offset)) {
        formula2 <- as.formula(paste(
          formula[2], formula[1],
          paste(local_var[which(local_var != "(Intercept)")], collapse = "+"),
          "+offset(", paste(new_global_var, collapse = "+"), ")"
        ))
      } else {
        formula2 <- as.formula(paste(
          formula[2], formula[1],
          paste(local_var[which(local_var != "(Intercept)")], collapse = "+"),
          "+offset(", offset, "+", paste(new_global_var, collapse = "+"), ")"
        ))
      }
    }
    cvar <- length(global_var)

    # varing coefficients
    fit_local <- gwglm(
      formula = formula2, family = family, cordxy = cordxy,
      distmethod = distmethod, h = h2,
      kernel = kernel, adaptive = adaptive,
      data = data2, cvar = cvar, parallel = parallel, core = core
    )
    varing_beta <- fit_local$beta
    local_beta <- fit_local$beta
    local_std <- fit_local$std
    local_tvalue <- fit_local$tvalue
    local_pvalue <- fit_local$pvalue

    local_beta_name <- c()
    local_std_name <- c()
    local_tvalue_name <- c()
    local_pvalue_name <- c()
    for (i in 1:length(local_var)) {
      if (local_var[i] == "(Intercept)") {
        local_beta_name[i] <- "Intercept_beta"
        local_std_name[i] <- "Intercept_std"
        local_tvalue_name[i] <- "Intercept_tvalue"
        local_pvalue_name[i] <- "Intercept_pvalue"
      } else {
        local_beta_name[i] <- paste0(local_var[i], "_beta")
        local_std_name[i] <- paste0(local_var[i], "_std")
        local_tvalue_name[i] <- paste0(local_var[i], "_tvalue")
        local_pvalue_name[i] <- paste0(local_var[i], "_pvalue")
      }
    }
    colnames(local_beta) <- local_beta_name
    colnames(local_std) <- local_std_name
    colnames(local_tvalue) <- local_tvalue_name
    colnames(local_pvalue) <- local_pvalue_name
    rownames(local_beta) <- 1:nrow(local_beta)
    rownames(local_std) <- 1:nrow(local_beta)
    rownames(local_tvalue) <- 1:nrow(local_beta)
    rownames(local_pvalue) <- 1:nrow(local_beta)
    global_beta <- as.data.frame(matrix(
      rep(constant_beta, each = nrow(data)),
      nrow(data), length(global_var)
    ))
    global_beta_name <- c()
    for (i in 1:length(global_var)) {
      if (global_var[i] == "(Intercept)") {
        global_beta_name[i] <- "Intercept_beta"
      } else {
        global_beta_name[i] <- paste0(global_var[i], "_beta")
      }
    }
    colnames(global_beta) <- global_beta_name

    muhat <- fit_local$muhat
    llk_H1 <- fit_global$logLik
    llk_H0 <- fit_local$logLik
    LRT_statistic <- llk_H1 - llk_H0

    output_data <- cbind(
      data, global_beta, local_beta, local_std,
      local_tvalue, local_pvalue
    )

    TS_time <- proc.time() - start
    #print(TS_time)

    out <- list(
      formula = formula, muhat = muhat, global_beta = constant_beta,
      local_beta = varing_beta, LRT_statistic = LRT_statistic,
      new_data = output_data,
      call = match.call()
    )
    class(out) <- "TSGWML_model"
    return(out)
  }
}
