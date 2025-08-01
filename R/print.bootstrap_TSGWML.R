#' Print Method for TSGWML Bootstrap Test Results
#'
#' Prints a brief summary of a \code{bootstrap_TSGWML} object,
#' including the function call, model formula, variable types, the observed likelihood ratio test (LRT) statistic,
#' the associated bootstrap p-value, and the bootstrap-based standard errors for both global and local coefficients.
#'
#' @param x An object of class \code{"bootstrap_TSGWML"} returned by \code{\link{bootstrap_TSGWML}}.
#' @param ... Further arguments passed to or from other methods. Currently ignored.
#'
#' @return
#' Invisibly returns \code{NULL}. Called for its side effect of printing model output to the console.
#'
#' @seealso \code{\link{bootstrap_TSGWML}}
#'
#' @importFrom utils head
#'
#' @method print bootstrap_TSGWML
#' @export
print.bootstrap_TSGWML <- function(x, ...) {
  cat("========== Bootstrap Test Result for TSGWML ==========\n")

  cat("Call:\n")
  print(x$call)

  cat("\nFormula:\n")
  print(x$formula)

  cat("\nLocal Variables:\n")
  print(x$local_var)

  cat("\nGlobal Variables:\n")
  print(x$global_var)

  cat("\nLRT Statistic (Observed):", round(x$c_tvalue, 4), "\n")
  cat("Bootstrap p-value:", round(x$pvalue, 5), "\n")

  cat("\nStandard Errors for Constant (Global) Coefficients:\n")
  print(round(x$constant_se, 4))

  cat("\nFirst 5 rows of Standard Errors for Varying (Local) Coefficients:\n")
  print(round(head(x$varying_se, 5), 4))

  cat("\nIf you want to see more details, use $ to extract them.\n")
  cat("===========================================================\n")
}
