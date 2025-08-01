#' Print Method for GWMMR MANOVA Objects
#'
#' Prints a brief summary of a \code{gwmmr_manova} object,
#' including the degrees of freedom used in the Pillai trace F-approximation
#' and the first few rows of the local multivariate test results.
#'
#' @param x An object of class \code{"gwmmr_manova"} returned by \code{\link{manova_gwmmr}}.
#' @param ... Further arguments passed to or from other methods. Currently ignored.
#'
#' @return
#' Invisibly returns \code{NULL}. Called for its side effect of printing model output to the console.
#'
#' @seealso \link{manova_gwmmr}
#'
#' @importFrom utils head
#'
#' @method print gwmmr_manova
#' @export
print.gwmmr_manova <- function(x, ...) {
  cat("========== GWMMR MANOVA Results: ==========\n")
  cat("\nDegrees of Freedom:\n")
  cat("  df1 (numerator): ", x$df1, "\n")
  cat("  df2 (denominator): ", round(x$df2), "\n")

  cat("\nLocal Multivariate Test Summary (first 5 rows):\n")
  print(head(x$test, 5))


  cat("\nIf you want to see more details, use $ to extract them.\n")
  cat("===========================================================\n")
}
