#' Print Method for GWQR Model Objects
#'
#' Prints a brief summary of a \code{gwqr_model} object, including the function call, model formula,
#' a summary of the local coefficients, and the check function loss.
#'
#' @param x An object of class \code{"gwqr_model"} returned by \code{\link{gwqr}}.
#' @param ... Further arguments passed to or from other methods. Currently ignored.
#'
#' @details
#' This print method displays a summary of a geographically weighted quantile regression (GWQR) model,
#' including its call, formula, a five-number summary of local coefficient estimates,
#' and the total check function loss used for model fitting.
#'
#' @return
#' Invisibly returns \code{NULL}. Called for its side effect of printing model output to the console.
#'
#' @seealso \link{gwqr}
#'
#' @method print gwqr_model
#' @export
print.gwqr_model <- function(x, ...) {
  cat("========== GWQR Model Results: ==========\n")
  cat("Call:\n")
  print(x$call)

  cat("\nFormula: ")
  print(x$formula)

  cat("\nCoefficients (Summary):\n")
  print(summary(x$beta))

  cat("\nCheck Function:", x$wsum, "\n")

  cat("\nIf you want to see more details, use $ to extract them.\n")
  cat("===========================================================\n")
}
