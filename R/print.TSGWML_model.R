#' Print Method for TSGWML Model Objects
#'
#' Prints a brief summary of a \code{TSGWML} object, including the function call, model formula,
#' global and local variables, global coefficient estimates, summary statistics
#' of local coefficients, and the likelihood ratio test (LRT) result.
#'
#' @param x An object of class \code{"gwglm_model"} returned by \code{\link{gwglm}}.
#' @param ... Further arguments passed to or from other methods. Currently ignored.
#'
#' @return
#' Invisibly returns \code{NULL}. Called for its side effect of printing model output to the console.
#'
#' @seealso \link{TSGWML}
#'
#' @method print TSGWML_model
#' @export
print.TSGWML_model <- function(x, ...) {
  cat("========== TSGWML Model Results: ==========\n")

  cat("Call:\n")
  print(x$call)

  cat("\nFormula:\n")
  print(x$formula)

  cat("\nGlobal Coefficient Estimates:\n")
  print(round(x$global_beta, 4))

  cat("\nLocal Coefficient Estimates (Summary):\n")
  print(summary(x$local_beta, 4))

  cat("\nLikelihood Ratio Test:", x$LRT_statistic, "\n")


  cat("\nIf you want to see more details, use $ to extract them.\n")
  cat("===========================================================\n")
}
