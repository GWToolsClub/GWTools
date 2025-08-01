#' Infant Mortality Data - Georgia Counties
#'
#' This dataset contains county-level infant mortality data for the 159 counties
#' in Georgia, USA, covering the period 2002–2004, as described in Chen & Yang (2011).
#'
#' @docType data
#' @name infant
#' @format A data frame with 159 observations and 13 variables:
#' \describe{
#'    \item{COFIPS}{Integer. County FIPS code.}
#'    \item{AVEDEATH03}{Integer. Three-year average number of infant deaths (under age 1) as the "success event".}
#'    \item{AVEBIR03}{Integer. Three-year average number of live births as the "total trials'.}
#'    \item{AVELBW03}{Integer. Three‑year average number of low‑birth‑weight (< 2500g) infants.}
#'    \item{BLACK}{Numeric (%). Percentage of county population identifying as Black.}
#'    \item{HISPANIC}{Numeric (%). Percentage identifying as Hispanic.}
#'    \item{OTHERS}{Numeric (%). Percentage identifying as races other than Black or Hispanic.}
#'    \item{GINI}{Numeric. Gini coefficient of income inequality.}
#'    \item{STABILITY}{Numeric. Residential‑stability index}
#'    \item{Xcoord}{Numeric. X-coordinate (Easting) of the county centroid.}
#'    \item{Ycoord}{Numeric. Y-coordinate (Northing) of the county centroid.}
#' }
#'
#' @details
#' \describe{
#'    \item{Response / trials}{When fitting a GW logistic model, \code{AVEDEATH03}
#'   is treated as the number of “successes” (infant deaths) and \code{AVEBIR03}
#'   as the total trials.}
#'    \item{Source}{All variables are derived from the 2008 Area Resource File (ARF),
#'    a county-level database maintained by the Bureau of Health Professions.}
#' }
#' For full variable definitions and modelling context, see reference below.
#'
#' @source Chen, V. Y-J., & Yang, T.-C. (2012).
#' SAS macro programs for geographically weighted generalized linear modeling
#' with spatial point data Applications to health research.
#' \emph{Computer Methods and Programs in Biomedicine}, 107 (2), 262-273.
"infant"
