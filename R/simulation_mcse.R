#' @title Monte Carlo Standard Error of a Mean
#' @description Estimate the Monte Carlo standard error of a performance measure
#'   that is reported as a mean across simulation repetitions, following the
#'   formulas in Morris, White & Crowther (2019), Statistics in Medicine 38:
#'   2074-2102 (\doi{10.1002/sim.8086}). The Monte Carlo SE of a mean is
#'   `sd(x) / sqrt(n)`.
#' @param x Numeric vector of per-repetition estimates.
#' @param na.rm Logical; drop missing values before computing (default `TRUE`).
#' @return A single numeric value, or `NA_real_` when fewer than two
#'   non-missing values are available.
#' @export
mcse_mean <- function(x, na.rm = TRUE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  n <- length(x)
  if (n < 2L) {
    return(NA_real_)
  }
  stats::sd(x) / sqrt(n)
}

#' @title Monte Carlo Standard Error of a Proportion
#' @description Monte Carlo standard error of a proportion (a binary performance
#'   measure such as a pass rate, power, or type I error rate), computed as
#'   `sqrt(p * (1 - p) / n)` with `p = mean(x)` as in Morris, White & Crowther
#'   (2019).
#' @param x Numeric vector of 0/1 per-repetition indicators, or a logical vector.
#' @param na.rm Logical; drop missing values before computing (default `TRUE`).
#' @return A single numeric value, or `NA_real_` when no non-missing values are
#'   available.
#' @export
mcse_prop <- function(x, na.rm = TRUE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  n <- length(x)
  if (n < 1L) {
    return(NA_real_)
  }
  p <- mean(as.numeric(x))
  sqrt(p * (1 - p) / n)
}

#' @title Monte Carlo Standard Error of a Bias Estimate
#' @description Monte Carlo standard error of a bias estimate. Because
#'   `bias = mean(x) - truth` differs from `mean(x)` by a constant, its Monte
#'   Carlo SE is the same as the Monte Carlo SE of the mean, `sd(x) / sqrt(n)`
#'   (Morris, White & Crowther 2019).
#' @param estimate Numeric vector of per-repetition estimates.
#' @param truth Numeric scalar: the true value of the estimand.
#' @param na.rm Logical; drop missing values before computing (default `TRUE`).
#' @return A single numeric value, or `NA_real_` when fewer than two
#'   non-missing values are available.
#' @export
mcse_bias <- function(estimate, truth, na.rm = TRUE) {
  mcse_mean(estimate, na.rm = na.rm)
}

#' @title Simulation Repetitions Needed for a Target Monte Carlo SE
#' @description Compute the number of simulation repetitions `n_sim` required to
#'   estimate a proportion to a target Monte Carlo standard error, using
#'   Equation (1) of Morris, White & Crowther (2019):
#'   `n_sim = p * (1 - p) / mcse_target^2`. The worst case is `p = 0.5`; for a
#'   target Monte Carlo SE of 0.5% this gives 10,000 repetitions, and for a 95%
#'   coverage target it gives 1,900.
#' @param p Anticipated value of the proportion (between 0 and 1). Use `0.5`
#'   for the worst case when the proportion is unknown.
#' @param mcse_target Target Monte Carlo standard error (on the proportion
#'   scale; e.g. `0.005` for 0.5%).
#' @return A single integer: the ceiling of the required number of repetitions.
#' @export
n_sim_for_mcse <- function(p, mcse_target = 0.005) {
  if (!is.numeric(p) || any(p < 0 | p > 1)) {
    stop("`p` must be between 0 and 1")
  }
  if (!is.numeric(mcse_target) || mcse_target <= 0) {
    stop("`mcse_target` must be a positive number")
  }
  as.integer(ceiling(p * (1 - p) / mcse_target^2 - 1e-9))
}

#' @title Summarise Simulation Performance with Monte Carlo Error
#' @description Group-wise summary of a performance measure across simulation
#'   repetitions, reporting the estimate, its Monte Carlo standard error, and a
#'   Monte Carlo confidence interval (Morris, White & Crowther 2019). This is
#'   the ADEMP "performance measures" dataset builder.
#' @param data A data frame of per-repetition results.
#' @param value Name (string) of the column holding the per-repetition
#'   performance measure.
#' @param group_cols Character vector of column names to group by (e.g.
#'   `"version"`). Use `character(0)` for an overall summary.
#' @param type Type of performance measure: `"mean"` for a continuous measure,
#'   `"prop"` for a proportion (binary measure), or `"bias"` for a bias
#'   estimate (requires `truth`).
#' @param truth Numeric scalar: true value of the estimand; required when
#'   `type = "bias"`.
#' @param level Confidence level for the Monte Carlo interval (default 0.95).
#' @return A data frame with the grouping columns plus `n`, `estimate`, `mcse`,
#'   `conf_low`, `conf_high`, `type`, and `level`.
#' @export
ademp_summarise <- function(data, value, group_cols = character(0),
                            type = c("mean", "prop", "bias"),
                            truth = NULL, level = 0.95) {
  type <- match.arg(type)
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame")
  }
  if (!is.character(value) || length(value) != 1L || !value %in% names(data)) {
    stop("`value` must name a single column of `data`")
  }
  if (!all(group_cols %in% names(data))) {
    stop("all `group_cols` must be columns of `data`")
  }
  if (type == "bias" && is.null(truth)) {
    stop("`truth` must be supplied when `type = \"bias\"`")
  }
  z <- stats::qnorm(1 - (1 - level) / 2)

  compute <- function(x) {
    x <- x[!is.na(x)]
    n <- length(x)
    if (n == 0L) {
      return(data.frame(
        n = 0L, estimate = NA_real_, mcse = NA_real_,
        conf_low = NA_real_, conf_high = NA_real_
      ))
    }
    est <- if (type == "bias") mean(x) - truth else mean(x)
    se <- if (type == "prop") {
      p <- mean(as.numeric(x))
      sqrt(p * (1 - p) / n)
    } else if (n < 2L) {
      NA_real_
    } else {
      stats::sd(x) / sqrt(n)
    }
    data.frame(
      n = n, estimate = est, mcse = se,
      conf_low = est - z * se, conf_high = est + z * se
    )
  }

  if (length(group_cols) == 0L) {
    out <- compute(data[[value]])
    out$type <- type
    out$level <- level
    return(out)
  }

  keys <- data[group_cols]
  groups <- split(seq_len(nrow(data)), keys, drop = TRUE)
  rows <- lapply(groups, function(idx) {
    cbind(keys[idx[1L], group_cols, drop = FALSE], compute(data[[value]][idx]))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$type <- type
  out$level <- level
  out
}
