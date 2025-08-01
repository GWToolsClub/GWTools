#' Tokyo Mortality Dataset
#'
#' A dataset containing mortality and socio-economic statistics for 262 municipalities
#' in the Tokyo metropolitan area, as described in Nakaya et al. (2005).
#'
#' @docType data
#' @name tokyo
#' @format A data frame with 262 observations and 8 variables:
#' \describe{
#'    \item{ID}{Integer. Municipality ID (0 to 261).}
#'    \item{X}{Numeric. X-coordinate (Easting) of the municipality centroid.}
#'    \item{Y}{Numeric. Y-coordinate (Northings) of the municipality centroid.}
#'    \item{Mort2564}{Integer. Observed number of deaths for ages 25-64 in 1990.}
#'    \item{Exp_2564}{Numeric. Expected number of deaths for ages 25-64, calculated using national age- and sex-specific mortality rates.}
#'    \item{Professl}{Numeric. Proportion of residents employed in professional and technical occupations.}
#'    \item{Elderly}{Numeric. Proportion of population aged 65 or above (OLD).}
#'    \item{OwnHome}{Numeric. Proportion of households that own their home (OWNH)}
#'    \item{Unemply}{Numeric. Unemployment rate (UNEMP)}
#' }
#'
#' @details
#' The dataset includes all 262 municipalities within an approximately 70 km radius
#' from Chiyoda ward (the localtion of the Imperial Palace), as defined in Naakaya et al. (2005)
#' These municipalities constitute the functional urban region of Tokyo. Each explanatory variable
#' was standardized (mean = 0, SD = 1) prior to model fitting.
#'
#' @source Nakaya, T., Fotheringham, A. S., Brunsdon, C., & Charlton, M. (2005).
#' Geographically weighted Poisson regression for disease association mapping.
#' \emph{Statistics in Medicine}, 24(17), 2695–2717.
#'
#' @examples
#' data(tokyo)
#' str(tokyo)
#' head(tokyo)
"tokyo"
