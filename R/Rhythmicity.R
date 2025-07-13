

#' Fit a sine curve using non-linear least squares
#'
#' This function fits a sine curve of the form  
#' \eqn{y = \text{amp} \cdot \sin\left(\frac{2\pi}{\text{period}} \cdot (x + \text{phase})\right) + \text{offset}}  
#' to time series data using Levenberg-Marquardt optimization from `minpack.lm::nls.lm`.
#'
#' @param x Numeric vector of time values (e.g., hours or days).
#' @param y Numeric vector of response values.
#' @param parStart Named list with starting values for parameters: `amp`, `phase`, `offset`, and `period`.
#'
#' @return A named list with:
#' \describe{
#'   \item{amp}{Amplitude (peak distance from baseline).}
#'   \item{phase}{Phase shift (horizontal time shift/ scanner rotation) in same units as `tt`.}
#'   \item{offset}{Vertical shift (baseline).}
#'   \item{peak}{Time position of the peak / scanner rotation position, relative to phase and period.}
#'   \item{period}{Distance between peaks}
#'   \item{tss}{Total sum of squares.}
#'   \item{rss}{Residual sum of squares.}
#'   \item{R2}{Coefficient of determination (pseudo R²).}
#'   \item{residual_se}{Residual standard error.}
#'   \item{df_residual}{Degrees of freedom of the residual.}
#'   \item{predicted}{Fitted values.}
#'   \item{residuals}{Model residuals.}
#' }
#'
#' @export
fit_sine_curve <- function(x, y,  parStart = list(amp = 3, phase = 0, offset = 0, period = 24)) {

   # df = data.frame(tt = tt, yy = yy)
   # df = dplyr::arrange(df,df$tt)
   # tt <- df$tt
   # yy <- df$yy
  
  
  getPred <- function(parS, x) {
    parS$amp * sin(2 * pi / parS$period * (x + parS$phase)) + parS$offset
  }
  
  
  residFun <- function(p, y, x) y - getPred(p, x)
  
  nls.out <- try(
    minpack.lm::nls.lm(par = parStart, fn = residFun, y = y, x = x),
    silent = TRUE
  )
  
  if (inherits(nls.out, "try-error")) {
    stop("Nonlinear least squares fitting failed. Check the input data and parameters.")
  }
  
  apar <- nls.out$par
  
  # Adjust amplitude and phase to be positive amplitude and correct time shift
  amp0 <- apar$amp
  asign <- sign(amp0)
  amp <- amp0 * asign
  period = apar$period
  
  phase0 <- apar$phase
  phase <- (phase0 + ifelse(asign == 1, 0, period / 2)) %% period
  offset <- apar$offset
  
  # Calculate peak time (when y reaches max)
  peak <- (period / 2 * sign(amp0) - period / 4 - phase) %% period
  if (peak > period * 0.75) peak <- peak - period
  
  
  # Residuals and model fit stats
  predicted <- getPred(apar, x)
  residuals <- y - predicted
  
  rss <- sum(residuals^2)
  tss <- sum((y - mean(y))^2)
  R2 <- 1 - rss / tss
  
  n <- length(y)
  df_resid <- n - length(apar)
  residual_se <- sqrt(rss / df_resid)
  
  list(
    amp = amp,
    phase = phase,
    peak = peak,
    offset = offset,
    period = period,
    tss = tss,
    rss = rss,
    R2 = R2,
    residual_se = residual_se,
    df_residual = df_resid,
    predicted = predicted,
    residuals = residuals
  )
}







# LRTest with validation and edge case handling
#'@description
#' Test the significance of circadian curve fitting using finite sample likelihood ratio test
#' @title Finite sample/ large sample Likelihood ratio test for circadian pattern detection
#' @param tt Time vector / scanner rotation position
#' @param yy Value
#' @param FN Type of test to use, TRUE = "FN" (finite) for 'F-like Test or FALSE = "LS" (large samples) for Likilihood Ratio Test. Default is finite sample.
#' @param parStart Named list with starting values for parameters: `amp`, `phase`, `offset`, and `period`.
#' @return A list of amp, phase, offset, sigma02, sigmaA2, l0, l1, df, stat, and pvalue.
#' model: \eqn{y= A * sin(2*pi*x+B)+C}
#' \item{y}{a 1*n vector of data y}
#' \item{A}{estimated A^hat from fit_sine_curve}
#' \item{B}{estimated B^hat from fit_sine_curve}
#' \item{C}{estimated C^hat from fit_sine_curve}
#' \item{sigma0}{sigma0^hat under H0}
#' \item{sigmaA}{sigmaA^hat under H1}
#' \item{n}{length of data y}
#' \item{df0}{df under H0}
#' \item{df1}{df under H1}
#'
#'  \item{amp}{Amplitude based on formula 1}
#'  \item{phase}{Phase based on formula 1, phase is restricted within (0, period)}
#'  \item{peakTime}{Phase based on formula 1, peakTime is restricted within (0, period). phase + peakTime = period/4}
#'  \item{offset}{Basal level(vertical shift) based on formula 1 or on formula 2}
#'  \item{sigma02}{Variance estimate under the null (intercept only)}
#'  \item{sigmaA2}{Variance estimate under the alternative (since curve fitting)}
#'  \item{l0}{Log likelihood under the null (intercept only)}
#'  \item{l1}{Log likelihood under the alternative (since curve fitting)}
#'  \item{df}{Degree of freedom for the LR test}
#'  \item{stat}{LR statistics}
#'  \item{pvalue}{P-value from the LR test}
#'  \item{R2}{Pseudo R2 defined as (tss - rss)/tss}
#' @details
#' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
#' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
#' @import dplyr
#' @importFrom stats pchisq
#' @importFrom stats pf
#' @keywords internal
#' @author Caleb (Zhiguang Huo)
#' @details
#' sourced from: https://rdrr.io/github/diffCircadian/diffCircadian/ under GPL3 license
#' author: Zhiguang Huo, Haocheng Ding
#' Feb. 24, 2023, 9:07 a.m.; Version 0.0.0
#' journal article: https://doi.org/10.1093/bib/bbab224
#' tutorial (for circadian data): http://htmlpreview.github.io/?https://github.com/diffCircadian/diffCircadian/blob/master/vignettes/diffCircadian_tutorial.html
#' @examples
#' set.seed(32608)
#' n <- 50
#' Period <- 24
#' tt <- runif(n,0,Period)
#' Amp <- 2
#' Phase <- 6
#' Offset <- 3
#' yy <- Amp * sin(2*pi/Period * (tt + Phase)) + Offset + rnorm(n,0,1)
#' rhytmicity_test(tt, yy,  FN =TRUE, parStart = list(amp = 3, phase = 0, offset = 0, period = 20) )
rhytmicity_test <- function(tt, yy, parStart = list(amp = 3, phase = 0, offset = 0, period = 12) , FN = TRUE) {

  # fit sine curve
  fitCurveOut <- fit_sine_curve(x = tt, y = yy, parStart)
  n <- length(yy)
  
  # Check for zero RSS to avoid issues in likelihood calculation
  if (fitCurveOut$rss == 0) {
    warning("RSS is zero, which might indicate issues with the fit.")
  }
  
  rss <- fitCurveOut$rss
  tss <- fitCurveOut$tss
  
  amp <- fitCurveOut$amp
  phase <- fitCurveOut$phase
  offset <- fitCurveOut$offset
  period <- fitCurveOut$period
  
  sigma02 <- 1/(n)*sum((yy-mean(yy))^2)
  sigmaA2 <- 1/(n)*sum((yy-amp*sin(2*pi/period*(tt+phase))-offset)^2)
  
  l0 <- -n/2*log(2*pi*sigma02)-1/(2*sigma02)*sum((yy-mean(yy))^2)
  l1 <- -n/2*log(2*pi*sigmaA2)-1/(2*sigmaA2)*sum((yy-amp*sin(2*pi/period*(tt+phase))-offset)^2)
  


  if(FN==FALSE){
    # Likelihood Ratio Test  
    dfdiff <- (n-1)-(n-3)
    df = data.frame(n)
    LR_stat <- -2*(l0-l1)
    stat = data.frame(LR_stat)
    pval <- stats::pchisq(LR_stat,dfdiff,lower.tail = F)
  }
  else if(FN==TRUE){
    # F-test 
    
    tss <- fitCurveOut$tss         # null model RSS
    rss <- fitCurveOut$rss         # full model RSS
    df1 <- 2                               # parameters added
    df2 <- n - 3                           # residual df full model
    
    F_stat <- ((tss - rss) / df1) / (rss / df2)
    stat = data.frame(F_stat)
    pval <- stats::pf(F_stat, df1, df2, lower.tail = FALSE)
    
    df = data.frame(df1 = df1, df2 = df2)
  }
  
  
  R2 <- 1-rss/tss
  peakTime = (period / 4 - phase) %% period
  res <- list(
    amp = amp,
    phase = phase,
    offset = offset,
    period = period,
    peakTime = peakTime,
    sigma02=sigma02, 
    sigmaA2=sigmaA2,
    l0=l0,
    l1=l1,
    df = df,
    stat= stat,
    pvalue=pval,
    R2=R2)
  
  return(res)
}









