#' Print Method for GWOLR Model Objects
#'
#' Prints a brief summary of a \code{gwolr_model} object, including the function call, model formula, and a summary of the local coefficients.
#'
#' @param x An object of class \code{"gwolr_model"} returned by \code{\link{gwolr}}.
#' @param ... Further arguments passed to or from other methods. Currently ignored.
#'
#' @return
#' Invisibly returns \code{NULL}. Called for its side effect of printing model output to the console.
#'
#' @seealso \link{gwolr}
#'
#' @method print gwolr_model
#' @export
print.gwolr_model <- function(x, ...) {
  cat("========== GWOLR Model Results: ==========\n")

  cat("Call:\n")
  print(x$call)

  cat("\nFormula:\n")
  print(x$formula)

  cat("\nCoefficients (Summary):\n")
  print(summary(x$beta))

  cat("\nIf you want to see more details, use $ to extract them.\n")
  cat("===========================================================\n")
}
