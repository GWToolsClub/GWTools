#' Print Method for GWGLM Model Objects
#'
#' Prints a brief summary of a \code{gwglm_model} object, including the function call, model formula,
#' a summary of the local coefficients, and diagnostic statistics.
#'
#' @param x An object of class \code{"gwglm_model"} returned by \code{\link{gwglm}}.
#' @param ... Further arguments passed to or from other methods. Currently ignored.
#'
#' @details
#' This print method displays a summary of a geographically weighted generalized linear model,
#' including its call, formula, coefficient summary, and AIC-type diagnostics.
#'
#' @return
#' Invisibly returns \code{NULL}. Called for its side effect of printing model output to the console.
#'
#' @seealso \link{gwglm}
#'
#' @importFrom utils head
#'
#' @method print gwglm_model
#' @export
print.gwglm_model <- function(x, ...) {
  cat("========== GWGLM Model Results: ==========\n")

  cat("Call:\n")
  print(x$call)

  cat("\nFormula:\n")
  print(x$formula)

  cat("\nCoefficients (Summary):\n")
  print(summary(x$beta))

  cat("\nDiagnostic Statistics:\n")
  cat(sprintf("   - Deviance : %.4f\n", x$deviance))
  cat(sprintf("   - LogLik   : %.4f\n", x$logLik))
  cat(sprintf("   - AIC      : %.4f\n", x$AIC))
  cat(sprintf("   - AICc     : %.4f\n", x$AICc))
  cat(sprintf("   - BIC      : %.4f\n", x$BIC))
  cat(sprintf("   - EDF (k)  : %.4f\n", x$k))

  cat("\nIf you want to see more details, use $ to extract them.\n")
  cat("===========================================================\n")
}
