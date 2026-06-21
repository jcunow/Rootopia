



#' Fit a sine curve to data with optional fixed period
#'
#' Fits a sine curve model: 
#' \eqn{ y = amp * sin(2 * pi / period * (x + phase)) + offset }
#' to data points (x, y) using nonlinear least squares.
#'
#' @param x Numeric vector of independent variable values (e.g., time).
#' @param y Numeric vector of dependent variable values.
#' @param par_start Named list of starting values for parameters:
#'   - amp: Amplitude (default 3)
#'   - phase: Phase shift (default 0)
#'   - offset: Vertical offset (default 0)
#'   - period: Period (default 24)
#' @param fix_period Numeric or NULL. If numeric, period is fixed to this value 
#'   and not estimated. If NULL (default), period is estimated.
#'
#' @return A list containing:
#'   - amp: Estimated amplitude
#'   - phase: Estimated phase (modulo period)
#'   - offset: Estimated offset
#'   - period: Estimated or fixed period
#'   - tss: Total sum of squares
#'   - rss: Residual sum of squares
#'   - R2: Coefficient of determination
#'   - residual_se: Residual standard error
#'   - df_residual: Residual degrees of freedom
#'   - predicted: Fitted values
#'   - residuals: Residuals (y - predicted)
#'
#' @importFrom minpack.lm nls.lm
#' @export
#' @examples
#' set.seed(1)
#' x <- seq(0, 48, length.out = 100)
#' period = 30.7
#' y <- 2 * sin(2 * pi / period * (x + 3)) + 5 + rnorm(100, 0, 0.8)
#' fit_3parameters <- fit_sine_curve(x, y, fix_period = 24)
#' fit_3parameters$R2
#' fit_3parameters$period
#' # estimate period
#' fit_4parameters <- fit_sine_curve(x, y, fix_period = NULL)
#' fit_4parameters$R2 
#' fit_4parameters$period
#' 
#'  plot(x,y)
#'  lines(x,fit_3parameters$predicted, col = "steelblue")
#'  lines(x,fit_4parameters$predicted, col = "coral")
fit_sine_curve <- function(x, y, par_start = list(amp=3, phase=0, offset=0, period=24), fix_period = NULL) {
  # Check input length
  if(length(x) != length(y)) stop("x and y must be of equal length")
  
  if (!is.null(fix_period)) {
    # Fix period to this value, estimate amp, phase, offset
    period <- fix_period
    
    residFun <- function(p, y, x) {
      amp <- p[1]
      phase <- p[2]
      offset <- p[3]
      y - (amp * sin(2 * pi / period * (x + phase)) + offset)
    }
    
    par_startVec <- c(par_start$amp, par_start$phase, par_start$offset)
    
    nls.out <- minpack.lm::nls.lm(par=par_startVec, fn=residFun, y=y, x=x)
    
    amp <- nls.out$par[1]
    phase <- nls.out$par[2] %% period
    offset <- nls.out$par[3]
    
  } else {
    # Estimate all four parameters: amp, phase, offset, period
    residFun <- function(p, y, x) {
      amp <- p["amp"]
      phase <- p["phase"]
      offset <- p["offset"]
      period <- p["period"]
      y - (amp * sin(2 * pi / period * (x + phase)) + offset)
    }
    
    nls.out <- minpack.lm::nls.lm(par=unlist(par_start), fn=residFun, y=y, x=x)
    
    amp <- nls.out$par["amp"]
    phase <- nls.out$par["phase"] %% nls.out$par["period"]
    offset <- nls.out$par["offset"]
    period <- nls.out$par["period"]
  }
  
  n <- length(y)
  k <- ifelse(is.null(fix_period), 4, 3)  # number of estimated parameters
  
  predicted <- amp * sin(2 * pi / period * (x + phase)) + offset
  residuals <- y - predicted
  
  rss <- sum(residuals^2)
  tss <- sum((y - mean(y))^2)
  R2 <- 1 - rss / tss
  df_residual <- n - k
  residual_se <- sqrt(rss / df_residual)
  
  list(
    amp = amp,
    phase = phase,
    offset = offset,
    period = period,
    tss = tss,
    rss = rss,
    R2 = R2,
    residual_se = residual_se,
    df_residual = df_residual,
    predicted = predicted,
    residuals = residuals
  )
}





















#' Assess rhythmicity via sine curve fitting and model comparison
#'
#' Fits a sine curve to time-series data and computes rhythmicity statistics.
#' The amplitude significance is tested via various statistical methods.
#'
#' @param x Numeric vector of time or phase values.
#' @param y Numeric vector of observed values (same length as \code{x}).
#' @param fix_period Either numeric (e.g., 24) to fix the period, or \code{NULL} to estimate it.
#' @param method Type of test to use: \code{"F"} for F-test (default), \code{"FLR"} for finite sample 
#'   likelihood ratio test, or \code{"LR"} for large sample likelihood ratio test.
#' @param par_start Named list with starting values for parameters: \code{amp}, \code{phase}, \code{offset}, 
#'   and \code{period}. Only used when \code{fix_period} is \code{NULL}.
#'
#' @return A named list with the following elements:
#' \itemize{
#'   \item \code{amplitude}: Fitted amplitude
#'   \item \code{phase}: Fitted phase
#'   \item \code{offset}: Fitted offset
#'   \item \code{period}: Fitted or fixed period
#'   \item \code{peakTime}: Peak time (phase + period/4)
#'   \item \code{R2}: Coefficient of determination
#'   \item \code{residual_se}: Residual standard error (F-test only)
#'   \item \code{df_residual}: Residual degrees of freedom (F-test only)
#'   \item \code{sigma02}: Variance under null (LR methods only)
#'   \item \code{sigmaA2}: Variance under alternative (LR methods only)
#'   \item \code{l0}: Log-likelihood under null (LR methods only)
#'   \item \code{l1}: Log-likelihood under alternative (LR methods only)
#'   \item \code{df}: Degrees of freedom for the test
#'   \item \code{stat}: Test statistic
#'   \item \code{p_value}: P-value for the test
#' }
#'
#' @examples
#' set.seed(123)
#' x <- seq(0, 48, length.out = 100)
#' y <- 1.5 * sin(2 * pi / 24 * (x + 4)) + 4 + rnorm(100, 0, 1.5)
#'
#' # Fixed period with F-test
#' rhythmicity(x, y, fix_period = 24, method = "F")
#'
#' # Fixed period with finite sample LR test
#' rhythmicity(x, y, fix_period = 24, method = "FLR")
#'
#' # Estimate period with F-test
#' rhythmicity(x, y, fix_period = NULL, method = "F")
#'
#' @export
#' 
rhythmicity <- function(x, y, fix_period = 24, method = "F", 
                        par_start = list(amp = 3, phase = 0, offset = 0, period = 12)) {
  
  if (!is.numeric(fix_period) && !is.null(fix_period)) {
    stop("fix_period must be numeric (e.g. 24) or NULL")
  }
  
  if (!method %in% c("F", "FLR", "LR")) {
    stop("method must be one of 'F', 'FLR', or 'LR'")
  }
  
  # Fit the sine model
  if (is.numeric(fix_period)) {
    fit <- fit_sine_curve(x, y, fix_period = fix_period)
  } else {
    fit <- fit_sine_curve(x, y, par_start = par_start)
  }
  
  n <- length(y)
  k <- if (is.numeric(fix_period)) 3 else 4
  
  rss_model <- sum((y - fit$predicted)^2)
  rss_null <- sum((y - mean(y))^2)
  R2 <- 1 - rss_model / rss_null
  
  peakTime <- (fit$period / 4 - fit$phase) %% fit$period
  
  result <- list(
    amplitude = fit$amp,
    phase = fit$phase,
    offset = fit$offset,
    period = fit$period,
    peakTime = peakTime,
    R2 = R2
  )
  
  if (method == "F") {
    df_residual <- n - k
    residual_se <- sqrt(rss_model / df_residual)
    
    F_value <- ((rss_null - rss_model) / (k - 1)) / (rss_model / df_residual)
    p_value <- stats::pf(F_value, df1 = k - 1, df2 = df_residual, lower.tail = FALSE)
    
    result$residual_se <- residual_se
    result$df_residual <- df_residual
    result$df <- data.frame(df1 = k - 1, df2 = df_residual)
    result$stat <- data.frame(F_stat = F_value)
    result$p_value <- p_value
    
  } else {
    sigma02 <- sum((y - mean(y))^2) / n
    sigmaA2 <- sum((y - fit$amp * sin(2 * pi / fit$period * (x + fit$phase)) - fit$offset)^2) / n
    
    if (sigma02 <= 0 || sigmaA2 <= 0) {
      warning("Variance estimate is zero or negative, which might indicate issues with the fit.")
    }
    
    l0 <- -n/2 * log(2 * pi * sigma02) - 1/(2 * sigma02) * sum((y - mean(y))^2)
    l1 <- -n/2 * log(2 * pi * sigmaA2) - 1/(2 * sigmaA2) * sum((y - fit$amp * sin(2 * pi / fit$period * (x + fit$phase)) - fit$offset)^2)
    
    dfdiff <- (n - 1) - (n - k)
    LR_stat <- -2 * (l0 - l1)
    
    result$sigma02 <- sigma02
    result$sigmaA2 <- sigmaA2
    result$l0 <- l0
    result$l1 <- l1
    
    if (method == "LR") {
      result$df <- data.frame(n = n, dfdiff = dfdiff)
      result$stat <- data.frame(LR_stat = LR_stat)
      result$p_value <- stats::pchisq(LR_stat, dfdiff, lower.tail = FALSE)
      
    } else if (method == "FLR") {
      r <- k - 1
      LR_stat_adj <- (exp(LR_stat / n) - 1) * (n - k) / r
      result$df <- data.frame(df1 = r, df2 = n - k)
      result$stat <- data.frame(FLR_stat = LR_stat_adj)
      result$p_value <- stats::pf(LR_stat_adj, df1 = r, df2 = n - k, lower.tail = FALSE)
    }
  }
  
  return(result)
}
