#' Compute Distance Matrix
#'
#' Computes either a full pairwise distance matrix or the distances from a single target point to all others
#' based on specified coordinate pairs and a selected distance metric.
#'
#' @param cordxy A matrix of spatial coordinates, with columns representing X and Y.
#' @param method A character string indicating the distance calculation method to use.
#' Supported options are \code{"euclidean"} (default), \code{"greatcircle"} (great-circle distance), and \code{"manhattan"} (city-block distance).
#' @param target_idx Optional integer. If specified, the function returns a vector of distances from the \code{target_idx}-th point to all others.
#' If \code{NULL} (default), a full pairwise distance matrix is returned.
#'
#' @return
#' Returns a full symmetric distance matrix of size \code{n x n} if \code{target_idx} is not specified.
#' Returns a numeric vector of distances if \code{target_idx} is specified.
#'
#' @details
#' The great-circle distance is calculated assuming a spherical Earth with a radius of 6371 kilometers.
#' It is suitable for coordinates expressed in latitude and longitude.
#'
#' @examples
#' coords <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
#' gwdist(coords)  # Euclidean distance matrix
#' gwdist(coords, method = "manhattan")  # Manhattan distance matrix
#' gwdist(coords, target_idx = 2)  # Distances from second point to all others
#'
#' @export
gwdist <- function(cordxy, method = "euclidean", target_idx = NULL) {
  # ---- check input ----
  if (ncol(cordxy) != 2) {
    stop("cordxy must be a matrix or data frame
         with two columns for coordinates.")
  }

  n <- nrow(cordxy)
  if (!is.null(target_idx)) {
    if (target_idx < 1 || target_idx > n) {
      stop("target_idx is out of range.")
    }
  }

  if (!method %in% c("euclidean", "greatcircle", "manhattan")) {
    stop("Unsupported method. Please choose
         'euclidean', 'greatcircle', or 'manhattan'.")
  }

  # ---- distance common compute ----
  compute_distance <- function(x1, y1, x2, y2, method) {
    # Euclidean Distance
    if (method == "euclidean") {
      return(sqrt((x1 - x2)^2 + (y1 - y2)^2))
      # Manhattan Distance
    } else if (method == "manhattan") {
      return(abs(x1 - x2) + abs(y1 - y2))
      # great circle Distance
    } else if (method == "greatcircle") {
      cos_angle <- sin(y1) * sin(y2) + cos(y1) * cos(y2) * cos(x1 - x2)
      cos_angle <- pmin(pmax(cos_angle, -1), 1)
      return(6371 * acos(cos_angle))
    } else {
      stop("Unsupported method.")
    }
  }

  # ---- Coordinates ----
  cordx <- cordxy[, 1]
  cordy <- cordxy[, 2]

  if (method == "greatcircle") {
    cordx <- cordx * pi / 180
    cordy <- cordy * pi / 180
  }

  # ---- Single point to all points ----
  if (!is.null(target_idx)) {
    xstart <- cordx[target_idx]
    ystart <- cordy[target_idx]
    d <- sapply(1:n, function(i) {
      compute_distance(xstart, ystart, cordx[i], cordy[i], method)
    })
    return(d)
  }

  # ---- Global pairwise distance matrix ----
  dist_mat <- matrix(0, n, n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      dist <- compute_distance(cordx[i], cordy[i], cordx[j], cordy[j], method)
      dist_mat[i, j] <- dist
      dist_mat[j, i] <- dist
    }
  }

  return(dist_mat)
}
