# ------------------------------------------------------------------
# Helper data wrangling functions ----------------------------------
# ------------------------------------------------------------------
#' Convert data from wide to long
#' @param x data frame to convert
#' @param times vector of time points. should be the same length as
#' nrow(x)
#' @return data frame in long format
#' @keywords vacamole
#' @export
wide_to_long <- function(x, times) {
  x %>%
    dplyr::mutate(time = times) %>%
    tidyr::pivot_longer(
      cols = !.data$time,
      names_to = c("state", "age_group"),
      names_sep = -1,
      values_to = "value"
    )
}

#' Post-process raw results from model run
#' @param x output from model run
#' @return data frame in long format
#' @importFrom stats quantile
#' @keywords vacamole
#' @export
wrangle_results <- function(x) {
  rtn_mle <- x %>%
    dplyr::filter(.data$sim == 0)

  rtn_bounds <- x %>%
    dplyr::filter(.data$sim != 0) %>%
    dplyr::group_by(.data$time, .data$state, .data$age_group) %>%
    dplyr::summarise(
      lower = stats::quantile(.data$value, probs = 0.025),
      upper = stats::quantile(.data$value, probs = 0.975)
    )

  rtn <- dplyr::left_join(
    rtn_mle, rtn_bounds,
    by = c("time", "state", "age_group")
  ) %>%
    dplyr::rename(mean = .data$value) %>%
    dplyr::select(-.data$sim)

  return(rtn)
}
