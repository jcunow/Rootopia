

#############################
# Enhanced tail weight function with error handling
#' Calculate weights for probability distribution comparison
#'
#' @param index Numeric vector specifying the positions for weight calculation
#' @param parameter List containing:
#'   - lambda: Shape parameter (0 = constant weighting)
#'   - x0: Curve offset/inflection point
#' @param index.spacing Character, either "equal" or "custom" for index spacing type
#' @param method Character, weighting function type:
#'   "constant", "asymptotic", "linear", "exponential", "sigmoid", "gompertz", "step"
#' @param baseline.weight Numeric between 0-1, minimum weight value
#' @param inverse Logical, if TRUE reverses weight distribution (left vs right tail)
#' @return Normalized weight vector
#'
#' @keywords internal
#' @export
tail_weight_function <- function(index = NULL, parameter = list(lambda = 0.2, x0=5),
                                 index.spacing = "equal", method = "sigmoid",
                                 baseline.weight = 0, inverse = FALSE) {
  tryCatch({
    # Validate method parameter
    valid_methods <- c("constant", "asympotic", "exponential", "sigmoid",
                       "linear", "gompertz", "step")
    if (!method %in% valid_methods) {
      stop("Invalid method. Must be one of: ",
           paste(valid_methods, collapse = ", "))
    }

    # Validate index spacing
    if (!index.spacing %in% c("equal", "custom")) {
      stop("index.spacing must be either 'equal' or 'custom'")
    }

    # Validate baseline weight
    if (!is.numeric(baseline.weight) || baseline.weight < 0 || baseline.weight > 1) {
      stop("baseline.weight must be a numeric value between 0 and 1")
    }

    # Validate parameter list
    if (!is.list(parameter)) {
      stop("parameter must be a list")
    }
    if (is.null(parameter$lambda) || is.null(parameter$x0)) {
      stop("parameter list must contain 'lambda' and 'x0'")
    }
    if (!is.numeric(parameter$lambda) || !is.numeric(parameter$x0)) {
      stop("parameter$lambda and parameter$x0 must be numeric")
    }

    # Process index based on spacing type
    if (index.spacing == "equal") {
      if (is.null(index)) {
        stop("index cannot be NULL when index.spacing is 'equal'")
      }
      n <- length(index)
      index <- 1:n
    } else {
      if (is.null(index)) {
        stop("index cannot be NULL when index.spacing is 'custom'")
      }
      n <- max(index, na.rm = TRUE)
      if (!is.finite(n)) {
        stop("Invalid index values")
      }
    }

    # Calculate weights based on method
    weights <- switch(method,
                      "constant" = rep(1, length(index)),
                      "asympotic" = {
                        exp_vals <- exp(parameter$lambda / n * abs((index) - n))
                        ((-exp_vals + max(exp_vals)) / max(exp_vals)) *
                          (1-baseline.weight) + baseline.weight
                      },
                      "exponential" = {
                        exp(-parameter$lambda / n * abs((index) - n)) *
                          (1-baseline.weight) + baseline.weight
                      },
                      "sigmoid" = {
                        1 / (1 + exp(-parameter$lambda * (index - parameter$x0))) *
                          (1-baseline.weight) + baseline.weight
                      },
                      "linear" = {
                        (index / n) * (1-baseline.weight) + baseline.weight
                      },
                      "gompertz" = {
                        exp(-parameter$x0 * exp(-parameter$lambda * index)) *
                          (1-baseline.weight) + baseline.weight
                      },
                      "step" = {
                        w <- index
                        w[index <= parameter$x0] <- 0
                        w[index > parameter$x0] <- 1
                        w * (1-baseline.weight) + baseline.weight
                      }
    )

    # Validate weights calculation
    if (any(is.na(weights)) || any(!is.finite(weights))) {
      stop("Invalid weights calculated. Check your parameters.")
    }

    if (inverse) {
      weights <- rev(weights)
    }

    # Normalize weights
    sum_weights <- sum(weights)
    if (sum_weights == 0) {
      stop("Sum of weights is zero")
    }
    out <- weights / sum_weights

    return(out)
  }, error = function(e) {
    stop("Error in tail_weight_function: ", e$message)
  })
}

# Enhanced KL divergence function with error handling
#' Calculate tail-weighted KL divergence for discrete distributions
#'
#' @param Q probability vector 1
#' @param P probability vector 2
#' @param index a positive numeric vector containing probability spacing e.g., depth
#' @param index.spacing whether index intervals are equally distant i.e., c(1,2,3,4....n), if "equal" than index is c(1,n)
#' @param inverse changes from right tail to left tail if TRUE
#' @param parameter list with lambda -> shape parameter (0 = constant weighting) & x0 -> curve offset (= inflexion point )
#' @param method weighting function along index. Available options are: c("constant", "asymptotic", "linear, "exponential", "sigmoid", "gompertz","step")
#' @param alignPQ if TRUE, index end values will be cut off in case of unequal length of P & Q so that length of P & Q is equal
#' @param cut if FALSE, 0 will be added to the shorter vector. If TRUE, the longer vector will be shortened at the end.
#'
#' @details
#' Kullback-Leibler Divergence
#'
#'
#' @return KL divergence, not symmetrical - changing the input order will change the result
#' @export
tail_weighted_kl_divergence <- function(P, Q, index = 1:min(c(length(Q),length(P))),
                                        index.spacing = "equal",
                                        parameter = list(lambda = 0.2, x0=30),
                                        cut = FALSE, inverse = FALSE,
                                        method = "constant", alignPQ = TRUE) {
  tryCatch({
    # Validate inputs
    if (!is.numeric(P) || !is.numeric(Q)) {
      stop("P and Q must be numeric vectors")
    }

    if (any(P < 0) || any(Q < 0)) {
      stop("P and Q must contain non-negative values")
    }

    if (sum(P) == 0 || sum(Q) == 0) {
      stop("P and Q must not sum to zero")
    }

    # Normalize P and Q if they don't sum to 1
    P <- P / sum(P)
    Q <- Q / sum(Q)

    # Handle length matching
    if (length(P) != length(Q)) {
      if (!alignPQ) {
        stop("Distributions P and Q must be of the same length when alignPQ is FALSE")
      }

      n <- min(c(length(Q), length(P)))
      m <- max(c(length(Q), length(P)))

      if (cut) {
        P <- P[1:n]
        Q <- Q[1:n]
      } else {
        vl1 <- abs(length(P) - m)
        vl2 <- abs(length(Q) - m)
        P <- c(P[1:(m-vl1)], rep(0, vl1))
        Q <- c(Q[1:(m-vl2)], rep(0, vl2))
      }
    }

    # Validate index
    if (!is.numeric(index) || any(index <= 0)) {
      stop("index must be a positive numeric vector")
    }
    if (length(index) != length(P)) {
      stop("index length must match distribution length")
    }

    # Calculate weights
    weight <- tail_weight_function(
      index = index,
      index.spacing = index.spacing,
      parameter = parameter,
      method = method,
      inverse = inverse
    )

    # Apply weights
    P <- (P * weight) / sum(P * weight)
    Q <- (Q * weight) / sum(Q * weight)

    # Calculate divergence
    kl_divergence <- sum(sapply(seq_along(index), function(i) {
      p_x <- P[i]
      q_x <- Q[i]
      if (p_x <= 0 || q_x <= 0) return(0)
      p_x * log(p_x / q_x)
    }))

    if (!is.finite(kl_divergence)) {
      stop("Invalid KL divergence calculated")
    }

    return(kl_divergence)
  }, error = function(e) {
    stop("Error in tail_weighted_kl_divergence: ", e$message)
  })
}

# Enhanced JS divergence function with error handling
#' Calculate tail-weighted Jensen-Shannon divergence
#'
#' @param Q probability vector 1
#' @param P probability vector 2
#' @param index a positive numeric vector containing probability spacing e.g., depth
#' @param index.spacing whether index intervals are equally distant i.e., c(1,2,3,4....n), if "equal" than index is c(1,n)
#' @param inverse changes from right tail to left tail if TRUE
#' @param parameter list with lambda -> shape parameter (0 = constant weighting) & x0 -> curve offset (= inflexion point )
#' @param method weighting function along index. Available options are: c("constant", "asymptotic", "linear, "exponential", "sigmoid", "gompertz","step")
#' @param alignPQ if TRUE, index end values will be cut off in case of unequal length of P & Q so that length of P & Q is equal
#' @param cut if FALSE, 0 will be added to the shorter vector. If TRUE, the longer vector will be shortened at the end.
#' @return Numeric JS divergence value
#'
#' @export
#'
#' @examples
#' P <- c(0.025,0.05,0.1,0.15, 0.2, 0.3,0.4, 0.5,0.3,0.1)  # Distribution P
#' Q <- c(0.025,0.05,0.1,0.15, 0.2, 0.3,0.4, 0.5,0.3,0.1)**6  # Distribution Q
#'
#' # Ensure the distributions are valid (non-negative and sum to 1)
#' P <- P / sum(P)
#' Q <- Q / sum(Q)
#'
#' tail_weighted_js_divergence(P,Q,parameter = list(lambda = 0.2,x0=30))
tail_weighted_js_divergence <- function(P, Q, parameter = list(lambda = 0.2, x0=30),
                                        method = "constant", inverse = FALSE,
                                        alignPQ = TRUE, cut = FALSE,
                                        index = 1:min(c(length(Q),length(P))),
                                        index.spacing = "equal") {
  tryCatch({
    # Validate inputs
    if (!is.numeric(P) || !is.numeric(Q)) {
      stop("P and Q must be numeric vectors")
    }

    if (any(P < 0) || any(Q < 0)) {
      stop("P and Q must contain non-negative values")
    }

    if (sum(P) == 0 || sum(Q) == 0) {
      stop("P and Q must not sum to zero")
    }

    # Normalize P and Q
    P <- P / sum(P)
    Q <- Q / sum(Q)

    # Handle length matching
    if (length(P) != length(Q)) {
      if (!alignPQ) {
        stop("Distributions P and Q must be of the same length when alignPQ is FALSE")
      }

      n <- min(c(length(Q), length(P)))
      m <- max(c(length(Q), length(P)))

      if (cut) {
        P <- P[1:n]
        Q <- Q[1:n]
      } else {
        vl1 <- abs(length(P) - m)
        vl2 <- abs(length(Q) - m)
        P <- c(P[1:(m-vl1)], rep(0, vl1))
        Q <- c(Q[1:(m-vl2)], rep(0, vl2))
      }
    }

    # Compute average distribution
    M <- (P + Q) / 2

    # Compute KL divergences with error checking
    KL_PM <- tail_weighted_kl_divergence(
      P, M, index = index,
      index.spacing = index.spacing,
      parameter = parameter,
      method = method,
      inverse = inverse,
      alignPQ = alignPQ
    )

    KL_QM <- tail_weighted_kl_divergence(
      Q, M, index = index,
      index.spacing = index.spacing,
      parameter = parameter,
      method = method,
      inverse = inverse,
      alignPQ = alignPQ
    )

    # Calculate final divergence
    js_divergence <- (KL_PM + KL_QM) / 2
    js_divergence <- round(js_divergence, 4)

    if (!is.finite(js_divergence)) {
      stop("Invalid JS divergence calculated")
    }

    return(js_divergence)
  }, error = function(e) {
    stop("Error in tail_weighted_js_divergence: ", e$message)
  })
}

# Enhanced Wasserstein distance function with error handling
#' A tailweighted Version of 1 dimensional Wasserstein distance betwwen two probability vectors
#'
#' @param Q probability vector 1
#' @param P probability vector 2
#' @param index a positive numeric vector containing probability spacing e.g., depth
#' @param index.spacing whether index intervals are equally distant i.e., c(1,2,3,4....n), if "equal" than index is c(1,n)
#' @param inverse changes from right tail to left tail if TRUE
#' @param parameter list with lambda -> shape parameter (0 = constant weighting) & x0 -> curve offset (= inflexion point )
#' @param method weighting function along index. Available options are: c("constant", "asymptotic", "linear, "exponential", "sigmoid", "gompertz","step")
#' @param baseline.weight Numeric between 0-1
#'
#' @return Numeric Wasserstein distance
#' @export
#'
#' @examples
#' P <- c(0.025,0.05,0.1,0.15, 0.2, 0.3,0.4, 0.5,0.3,0.1)  # Distribution P
#' Q <- c(0.025,0.05,0.1,0.15, 0.2, 0.3,0.4, 0.5,0.3,0.1)**6  # Distribution Q
#'
#' # Ensure the distributions are valid (non-negative and sum to 1)
#' P <- P / sum(P)
#' Q <- Q / sum(Q)
#'
#' tail_weighted_wasserstein_distance(P,Q,
#'   inverse=FALSE,method="constant",parameter = list(lambda = 0.2,x0=3))
tail_weighted_wasserstein_distance <- function(Q, P, inverse = FALSE,
                                               parameter = list(lambda = 0.2, x0=10),
                                               method = "step", baseline.weight = 0,
                                               index = 1:min(c(length(Q),length(P))),
                                               index.spacing = "equal") {
  tryCatch({
    # Validate inputs
    if (!is.numeric(P) || !is.numeric(Q)) {
      stop("P and Q must be numeric vectors")
    }

    if (any(P < 0) || any(Q < 0)) {
      stop("P and Q must contain non-negative values")
    }

    if (sum(P) == 0 || sum(Q) == 0) {
      stop("P and Q must not sum to zero")
    }

    # Validate baseline weight
    if (!is.numeric(baseline.weight) || baseline.weight < 0 || baseline.weight > 1) {
      stop("baseline.weight must be between 0 and 1")
    }

    # Calculate weights
    weights <- tail_weight_function(
      index = index,
      parameter = parameter,
      index.spacing = index.spacing,
      method = method,
      baseline.weight = baseline.weight,
      inverse = inverse
    )

    # Apply weights and normalize
    P <- (P * weights) / sum(P * weights)
    Q <- (Q * weights) / sum(Q * weights)

    # Calculate Wasserstein distance
    distance <- tryCatch({
      transport::wasserstein1d(a = P, b = Q)
    }, error = function(e) {
      stop("Error calculating Wasserstein distance: ", e$message)
    })

    if (!is.finite(distance)) {
      stop("Invalid Wasserstein distance calculated")
    }

    return(distance)
  }, error = function(e) {
    stop("Error in tail_weighted_wasserstein_distance: ", e$message)
  })
}

# Enhanced MRD function with error handling
#' Calculate Mean Rooting Depth
#'
#' @param w Numeric vector of weights (typically depths)
#' @param roots Numeric vector of root coverage values
#' @return Numeric MRD value
#' @export
#' @examples
#' w <- seq(5, 25, 5)
#' roots <- c(0, 10, 7, 3, 1)
#' mean_rooting_depth <- MRD(w, roots)
MRD <- function(w, roots) {
  tryCatch({
    # Validate inputs
    if (!is.numeric(w) || !is.numeric(roots)) {
      stop("w and roots must be numeric vectors")
    }

    if (length(w) != length(roots)) {
      stop("w and roots must have the same length")
    }

    if (any(is.na(w)) || any(is.na(roots))) {
      stop("w and roots cannot contain NA values")
    }

    if (any(roots < 0)) {
      stop("roots values must be non-negative")
    }

    if (sum(roots) == 0) {
      stop("sum of roots cannot be zero")
    }

    # Calculate MRD with small epsilon to avoid division by zero
    mrd <- sum((w + .Machine$double.eps) * roots) / sum(roots)

    if (!is.finite(mrd)) {
      stop("Invalid MRD calculated")
    }

    return(mrd)
  }, error = function(e) {
    stop("Error in MRD: ", e$message)
  })
}

# Enhanced RPI function with error handling

#' Calculate Root Penetration Index
#'
#' @param roots Numeric vector of root coverage values
#' @param w Numeric vector of weights (typically depths)
#' @return Numeric RPI value between -1 and 1
#' @export
#' @examples
#' w <- seq(5, 25, 5)
#' roots <- c(0, 10, 7, 3, 1)
#' rpi <- RPI(roots, w)
RPI <- function(roots, w) {
  tryCatch({
    # Validate inputs
    if (!is.numeric(w) || !is.numeric(roots)) {
      stop("w and roots must be numeric vectors")
    }

    if (length(w) != length(roots)) {
      stop("w and roots must have the same length")
    }

    if (any(is.na(w)) || any(is.na(roots))) {
      stop("w and roots cannot contain NA values")
    }

    if (any(roots < 0)) {
      stop("roots values must be non-negative")
    }

    if (sum(roots) == 0) {
      stop("sum of roots cannot be zero")
    }

    if (sum(w) == 0) {
      stop("sum of weights cannot be zero")
    }

    # Calculate RPI
    rpi <- 1 - 2 * sum(roots/sum(roots) * ((w + .Machine$double.eps)/(sum(w))))

    if (!is.finite(rpi)) {
      stop("Invalid RPI calculated")
    }

    # Validate result is in correct range
    if (rpi < -1 || rpi > 1) {
      stop("RPI calculation resulted in value outside valid range [-1, 1]")
    }

    return(rpi)
  }, error = function(e) {
    stop("Error in RPI: ", e$message)
  })
}
