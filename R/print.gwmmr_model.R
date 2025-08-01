#' Print Method for GWMMR Model Objects
#'
#' Prints a brief summary of a \code{gwmmr_model} object,
#' including the model formula, estimated local coefficients, standard errors,
#' Pillai Trace Statistics, local p-values,
#' and model fit statistics such as AIC, AICc, and effective degree of freedom.
#'
#' @param x An object of class \code{"gwmmr_model"} returned by \code{\link{gwmmr}}.
#' @param ... Further arguments passed to or from other methods. Currently ignored.
#'
#' @return
#' Invisibly returns \code{NULL}. Called for its side effect of printing model output to the console.
#'
#' @seealso \link{gwmmr}
#'
#' @importFrom utils head
#'
#' @method print gwmmr_model
#' @export
print.gwmmr_model <- function(x, ...) {
  cat("========== GWMMR Model Results: ==========\n")

  cat("\nFormula:\n")
  print(x$formula)

  cat("\nCoefficients (First few observations):\n")
  print(head(x$test[, grep("^estx", names(x$test))], 3))

  cat("\nStandard Errors (First few observations):\n")
  print(head(x$test[, grep("^se_estx", names(x$test))], 3))

  cat("\nPillai Trace Statistics (First few observations):\n")
  print(head(x$test[, grep("^Pillai", names(x$test))], 3))

  cat("\nLocal p-values (first few observations):\n")
  print(head(x$test[, grep("^P_value", names(x$test))], 3))

  cat("\nModel Fit Statistics:\n")
  cat(sprintf("   - AIC  : %.4f\n", x$AIC))
  cat(sprintf("   - AICc : %.4f\n", x$AICc))
  cat(sprintf("   - enp  : %.4f\n", x$enp))
  cat(sprintf("   - enp2 : %.4f\n", x$enp2))
  cat(sprintf("   - df2  : %.4f\n", x$df2))

  cat("\nIf you want to see more details, use $ to extract them.\n")
  cat("===========================================\n")
}
