#' Make a lagged copy of a vector
#' @param x Input vector
#' @return y A lagged copy of the input vector x, with the first entry repeated
#' @details Not normally called directly, but via the cp_wrapper function instead.

lag_v <- function(x) c(x[1], x[1:(length(x) - 1)])

#' Make a lead copy of a vector
#' @param x Input vector of cumulative response counts
#' @return y A lead copy of the input vector x, with the last entry repeated
#' @details Not normally called directly, but via the cp_wrapper function instead.
#'
lead_v <- function(x) c(x[2:(length(x))], x[length(x)])
