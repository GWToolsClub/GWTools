#' Bootstrap Significance Test for TSGWML Models
#'
#' This function evaluates whether the spatially varying (local) coefficients in a two-stage
#' geographically weighted model significantly improve model fit.
#' It performs a parametric bootstrap-based likelihood ratio test (LRT),
#' by regenerating the response variable based on predicted values from the fitted model,
#' to assess the significance of spatially varying coefficients.
#' Additionally, the function estimates bootstrap-based standard errors for all coefficients.
#'
#' @param formula A formula object specifying the model structure (response ~ predictors).
#' @param family A character string specifying the distribution family. Supported options are \code{"gaussian"}, \code{"poisson"}, and \code{"binomial"}.
#' @param cordxy A matrix of spatial coordinates, with columns representing X and Y.
#' @param distmethod A character string indicating the distance calculation method to use.
#' Supported options are \code{"euclidean"} (default), \code{"greatcircle"} (great-circle distance), and \code{"manhattan"} (city-block distance).
#' @param h1 Bandwidth used for the first-stage (fully local) estimation.
#' @param h2 Bandwidth used for the second-stage (local + global) estimation.
#' @param kernel A character string specifying the kernel function used to compute spatial weights.
#' Options are \code{"global"}, \code{"gaussian"}, \code{"exponential"}, \code{"bisquare"}, \code{"tricube"}, and \code{"boxcar"}.
#' @param adaptive Logical. If \code{TRUE} (default), uses adaptive bandwidth. If \code{FALSE}, uses fixed bandwidth.
#' @param offset A character string giving the name of a column in the data to be used as offset (optional).
#' If not provided, only the internal offset from fixed coefficients is used.
#' @param local_var A character vector of variable names to be treated as spatially varying (local) coefficients.
#' @param global_var A character vector of variable names to be treated as spatially invariant (global) coefficients.
#' @param data A data frame containing all variables used in the model.
#' @param R Integer. Number of bootstrap replications. Default is 99.
#' @param parallel Logical. If \code{TRUE} (default), enables parallel computation. Default is \code{TRUE}.
#' @param core Integer. The number of CPU cores to use when \code{parallel = TRUE}.
#' If \code{NULL}, the number of cores is set to the total number of detected cores minus 2.
#'
#' @return A list of class \code{"bootstrap_TSGWML"} with the following components:
#' \itemize{
#'   \item \code{formula}: A formula object.
#'   \item \code{local_var}: A character vector of spatially varying variable names.
#'   \item \code{global_var}: A character vector of spatially invariant variable names.
#'   \item \code{pvalue}: A numeric scalar. The bootstrap p-value from the likelihood ratio test.
#'   \item \code{c_tvalue}: A numeric scalar. The observed LRT statistic from the original data.
#'   \item \code{constant_se}: A named numeric vector of bootstrap standard errors for global coefficients.
#'   \item \code{varying_se}: A data frame of bootstrap standard errors for local coefficients.
#'   \item \code{call}: The matched function call that produced this object.
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
#' @importFrom stats model.frame model.response terms sd rpois rbinom rnorm
#' @importFrom stringr str_c
#' @importFrom progress progress_bar
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
#' result <- bootstrap_TSGWML(
#'   formula = formula,
#'   family = family,
#'   cordxy = cordxy,
#'   h1 = ts_bw_cv$h1,
#'   h2 = ts_bw_cv$h2,
#'   kernel = "gaussian",
#'   adaptive = FALSE,
#'   offset = offset,
#'   local_var = local_var,
#'   global_var = global_var,
#'   data = tokyo,
#'   parallel = TRUE
#' )
#' result
#' }
#' @export
bootstrap_TSGWML <- function(formula, family, cordxy, distmethod = "euclidean",
                             h1, h2, kernel, adaptive, offset = NULL,
                             local_var, global_var, data, R = 99,
                             parallel = TRUE, core = NULL){
  #start_time <- proc.time()

  if (!is.matrix(cordxy)) {
    stop("cordxy must be a matrix.")
  }

  # ------------ Step1: TSGWML variability -------------
  model_obs <- TSGWML(formula = formula, family = family,
                      cordxy = cordxy, distmethod = distmethod,
                      kernel = kernel, adaptive = adaptive,
                      h1 = h1, h2 = h2, offset = offset,
                      local_var = local_var, global_var = global_var,
                      data = data, parallel = parallel, core = core)
  mu <- model_obs$muhat
  c_tvalue <- model_obs$LRT_statistic
  constant_beta <- model_obs$global_beta

  # ------------ Step2~4: Bootstrap Test ------------
  T_value <- c()
  constant_list <- matrix(NA, nrow = R, ncol = length(global_var))
  colnames(constant_list) <- global_var
  varying_list <- vector("list", length(local_var))
  names(varying_list) <- local_var

  pb <- progress_bar$new(total = R)
  y <- model.response(model.frame(formula, data = data), "numeric")
  yname <- all.vars(formula)[attr(terms(formula), "response")]
  sigma_hat <- sd(y - mu)

  for (i in 1:R) {
    data_boot <- data

    for (j in 1:nrow(data_boot)){
      if(family == "poisson"){
        data_boot[, yname][j] <- rpois(1, mu[j])
      }else if(family == "binomial" && ncol(as.data.frame(y))==1){
        data_boot[, yname][j] <- rbinom(1, 1, mu[j])
      }else if(family == "binomial" && ncol(as.data.frame(y))==2){
        data_boot[, yname][j] <- rbinom(1, apply(y,1,sum)[j], mu[j])
      }else if(family == "gaussian"){
        data_boot[, yname][j] <- rnorm(1, mean = mu[j], sd = sigma_hat)
      }
    }

    model_boot <- TSGWML(formula = formula, family = family,
                        cordxy = cordxy, distmethod = distmethod,
                        kernel = kernel, adaptive = adaptive,
                        h1 = h1, h2 = h2, offset = offset,
                        local_var = local_var, global_var = global_var,
                        data = data_boot, parallel = parallel, core = core)

    T_value <- c(T_value, model_boot$LRT_statistic)
    constant_list[i, ] <- model_boot$global_beta

    for (m in seq_along(local_var)) {
      varying_list[[m]] <- cbind(varying_list[[m]], model_boot$local_beta[, m])
    }

    pb$tick()
  }

  # Step 5 : Compute p-value and standard-error
  pvalue <- sum(T_value >= c_tvalue) / R

  constant_se <- apply(as.data.frame(constant_list), 2, sd)
  varying_se <- rapply(varying_list, function(l) as.data.frame(apply(l, 1, sd)), how = "list")
  for(i in 1:length(local_var)) colnames(varying_se[[i]]) <- stringr::str_c(local_var,"_se")[i]#
  varying_se_bind <- do.call(cbind, varying_se)

  #bootstrap_time <- proc.time()-start_time
  #print(bootstrap_time)

  out <- list(
    formula = formula, local_var = local_var, global_var = global_var,
    pvalue = pvalue, c_tvalue = c_tvalue,
    constant_se = constant_se, varying_se = varying_se_bind,
    call = match.call()
  )
  class(out) <- "bootstrap_TSGWML"
  return(out)
}
