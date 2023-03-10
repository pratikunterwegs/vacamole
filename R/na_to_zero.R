#' @title Convert NAs to zeros.
#' @noRd
na_to_zero <- function(x) {
  ifelse(is.na(x), 0, x)
}
