############################################
# Tail weight function
############################################

#' Calculate weights for probability distribution comparison
#'
#' Constructs a flexible weighting vector used to emphasize or de-emphasize
#' portions of a probability distribution along an index (e.g., depth).
#'
#' @param index Numeric vector specifying positions.
#' @param parameter List with:
#'   \itemize{
#'     \item \code{lambda}: shape parameter controlling steepness
#'     \item \code{x0}: inflection or threshold parameter
#'   }
#' @param index.spacing Character. Either \code{"equal"} or \code{"custom"}.
#' @param weighting Character. Weighting function:
#' \code{"constant"}, \code{"asymptotic"}, \code{"linear"},
#' \code{"exponential"}, \code{"sigmoid"}, \code{"gompertz"}, \code{"step"}.
#' @param baseline.weight Numeric in `[0,1]`. Minimum weight applied.
#' @param inverse Logical. If TRUE, reverses weighting direction.
#'
#' @return Numeric vector of normalized weights.
#' @keywords internal
tail_weight_function <- function(
    index = NULL,
    parameter = list(lambda = 0.2, x0 = 5),
    index.spacing = "equal",
    weighting = "sigmoid",
    baseline.weight = 0,
    inverse = FALSE
) {
  tryCatch({
    
    valid_weightings <- c(
      "constant", "asymptotic", "linear",
      "exponential", "sigmoid", "gompertz", "step"
    )
    
    if (!weighting %in% valid_weightings) {
      stop("Invalid weighting function")
    }
    
    if (!index.spacing %in% c("equal", "custom")) {
      stop("index.spacing must be 'equal' or 'custom'")
    }
    
    if (!is.numeric(baseline.weight) || baseline.weight < 0 || baseline.weight > 1) {
      stop("baseline.weight must be in [0,1]")
    }
    
    if (!is.list(parameter) || is.null(parameter$lambda) || is.null(parameter$x0)) {
      stop("parameter must contain lambda and x0")
    }
    
    if (index.spacing == "equal") {
      if (is.null(index)) stop("index required for equal spacing")
      n <- length(index)
      index <- 1:n
    } else {
      if (is.null(index)) stop("index required for custom spacing")
      n <- max(index, na.rm = TRUE)
    }
    
    weights <- switch(
      weighting,
      
      constant = rep(1, length(index)),
      
      asymptotic = {
        exp_vals <- exp(parameter$lambda / n * abs(index - n))
        ((-exp_vals + max(exp_vals)) / max(exp_vals)) *
          (1 - baseline.weight) + baseline.weight
      },
      
      exponential = {
        exp(-parameter$lambda / n * abs(index - n)) *
          (1 - baseline.weight) + baseline.weight
      },
      
      linear = {
        (index / n) * (1 - baseline.weight) + baseline.weight
      },
      
      sigmoid = {
        1 / (1 + exp(-parameter$lambda * (index - parameter$x0))) *
          (1 - baseline.weight) + baseline.weight
      },
      
      gompertz = {
        exp(-parameter$x0 * exp(-parameter$lambda * index)) *
          (1 - baseline.weight) + baseline.weight
      },
      
      step = {
        w <- ifelse(index > parameter$x0, 1, 0)
        w * (1 - baseline.weight) + baseline.weight
      }
    )
    
    if (inverse) weights <- rev(weights)
    
    weights / sum(weights)
    
  }, error = function(e) {
    stop("tail_weight_function error: ", e$message)
  })
}


############################################
# Alignment of P and Q
############################################

#' Bring two mass vectors to a common length
#'
#' With \code{cut = FALSE} the shorter vector is padded with zeros at the end;
#' with \code{cut = TRUE} both are cut to the shorter length.
#'
#' @param P,Q Numeric mass vectors.
#' @param alignPQ Logical. If \code{FALSE}, unequal lengths are an error.
#' @param cut Logical. Truncate instead of padding.
#' @return A list with the aligned \code{P} and \code{Q}.
#' @keywords internal
#' @noRd
.align_pq <- function(P, Q, alignPQ = TRUE, cut = FALSE) {
  if (length(P) == length(Q)) return(list(P = P, Q = Q))
  if (!alignPQ) stop("P and Q differ in length; set alignPQ = TRUE to align them")
  n <- if (cut) min(length(P), length(Q)) else max(length(P), length(Q))
  fit <- function(x) if (length(x) >= n) x[seq_len(n)] else c(x, rep(0, n - length(x)))
  list(P = fit(P), Q = fit(Q))
}

#' Default or check the bin positions after alignment
#' @keywords internal
#' @noRd
.resolve_index <- function(index, n) {
  if (is.null(index)) return(seq_len(n))
  if (length(index) != n)
    stop("`index` has ", length(index), " position(s) but P and Q have ", n,
         " bin(s) after alignment; give one position per bin.", call. = FALSE)
  index
}


############################################
# KL divergence
############################################

#' Tail-weighted Kullback-Leibler divergence
#'
#' Computes KL divergence between two probability distributions with optional
#' depth-dependent weighting along an index (e.g., soil depth).
#'
#' @param P Numeric probability vector.
#' @param Q Numeric probability vector.
#' @param index Numeric positions of the bins (e.g. depths), one per bin
#'   after alignment. \code{NULL} (default) uses \code{1, 2, ..., n}.
#' @param index.spacing Character. \code{"equal"} or \code{"custom"}.
#' @param parameter List with \code{lambda} and \code{x0}.
#' @param weighting Character. Weighting function applied before comparison.
#' @param baseline.weight Numeric in `[0,1]`. Minimum weight.
#' @param inverse Logical. Reverse weighting direction.
#' @param alignPQ Logical. If TRUE, aligns vectors of unequal length;
#'   if FALSE, unequal lengths are an error.
#' @param cut Logical. If TRUE, truncates the longer vector instead of padding
#'   the shorter one with zeros.
#'
#' @details
#' KL divergence is asymmetric:
#' \deqn{D_{KL}(P \parallel Q) = \sum P \log(P / Q)}
#'
#' Weighting modifies contribution of each index position prior to normalization.
#'
#' @return Numeric KL divergence value.
#' @keywords internal
tail_weighted_kl_divergence <- function(
    P, Q,
    index = NULL,
    index.spacing = "equal",
    parameter = list(lambda = 0.2, x0 = 30),
    weighting = "constant",
    baseline.weight = 0,
    inverse = FALSE,
    alignPQ = TRUE,
    cut = FALSE
) {
  tryCatch({
    
    pq <- .align_pq(P, Q, alignPQ = alignPQ, cut = cut)
    P <- pq$P / sum(pq$P)
    Q <- pq$Q / sum(pq$Q)
    index <- .resolve_index(index, length(P))

    w <- tail_weight_function(
      index = index,
      parameter = parameter,
      index.spacing = index.spacing,
      weighting = weighting,
      baseline.weight = baseline.weight,
      inverse = inverse
    )
    
    P <- (P * w) / sum(P * w)
    Q <- (Q * w) / sum(Q * w)
    
    sum(P * log(P / Q), na.rm = TRUE)
    
  }, error = function(e) {
    stop("KL error: ", e$message)
  })
}


############################################
# JS divergence
############################################

#' Tail-weighted Jensen-Shannon divergence
#'
#' Computes a symmetric divergence between two distributions using a
#' weighted KL decomposition.
#'
#' @inheritParams tail_weighted_kl_divergence
#'
#' @details
#' Jensen-Shannon divergence is defined as:
#' \deqn{JS(P,Q) = 1/2 KL(P || M) + 1/2 KL(Q || M)}
#' where \eqn{M = (P + Q)/2}.
#'
#' Weighting is applied prior to normalization of P, Q, and M.
#'
#' @return Numeric JS divergence value.
#' @keywords internal
tail_weighted_js_divergence <- function(
    P, Q,
    index = NULL,
    index.spacing = "equal",
    parameter = list(lambda = 0.2, x0 = 30),
    weighting = "constant",
    baseline.weight = 0,
    inverse = FALSE,
    alignPQ = TRUE,
    cut = FALSE
) {
  
  pq <- .align_pq(P, Q, alignPQ = alignPQ, cut = cut)
  P <- pq$P / sum(pq$P)
  Q <- pq$Q / sum(pq$Q)
  index <- .resolve_index(index, length(P))
  
  M <- (P + Q) / 2
  
  KL_PM <- tail_weighted_kl_divergence(
    P, M,
    index = index,
    index.spacing = index.spacing,
    parameter = parameter,
    weighting = weighting,
    baseline.weight = baseline.weight,
    inverse = inverse
  )
  
  KL_QM <- tail_weighted_kl_divergence(
    Q, M,
    index = index,
    index.spacing = index.spacing,
    parameter = parameter,
    weighting = weighting,
    baseline.weight = baseline.weight,
    inverse = inverse
  )
  
  (KL_PM + KL_QM) / 2
}


############################################
# Wasserstein distance
############################################

#' Tail-weighted Wasserstein distance (1D)
#'
#' Computes the 1D Wasserstein distance between two probability distributions
#' with optional depth-dependent weighting.
#'
#' @inheritParams tail_weighted_kl_divergence
#'
#' @details
#' Wasserstein distance measures the minimal transport cost required to
#' transform one distribution into another along a one-dimensional index: the
#' mass-weighted mean distance the probability has to travel, in the units of
#' \code{index}. For depth profiles indexed in centimetres it is therefore a
#' depth, and a profile shifted \code{d} cm returns exactly \code{d}.
#'
#' \code{transport::wasserstein1d()} does the work. Its \code{a} and \code{b}
#' are positions and its \code{wa} and \code{wb} the mass at those positions,
#' so \code{index} is passed as the positions and \code{P} and \code{Q} as
#' the weights. Passing the mass vectors as positions instead would sort them
#' and discard \code{index}, which measures how differently concentrated the
#' two profiles are rather than how far apart they sit.
#'
#' @return Numeric Wasserstein distance, in the units of \code{index}.
#' @keywords internal
tail_weighted_wasserstein_distance <- function(
    P, Q,
    index = NULL,
    index.spacing = "equal",
    parameter = list(lambda = 0.2, x0 = 10),
    weighting = "constant",
    baseline.weight = 0,
    inverse = FALSE,
    alignPQ = TRUE,
    cut = FALSE
) {
  
  # `index` carries the positions the mass sits at and P/Q the mass at each,
  # so the cost of moving mass between two bins is the distance between their
  # index values. Zero mass on a padded bin moves nothing.
  pq <- .align_pq(P, Q, alignPQ = alignPQ, cut = cut)
  P <- pq$P; Q <- pq$Q
  index <- .resolve_index(index, length(P))

  w <- tail_weight_function(
    index = index,
    parameter = parameter,
    index.spacing = index.spacing,
    weighting = weighting,
    baseline.weight = baseline.weight,
    inverse = inverse
  )
  
  P <- (P * w) / sum(P * w)
  Q <- (Q * w) / sum(Q * w)

  tryCatch({
    transport::wasserstein1d(a = index, b = index, wa = P, wb = Q)
  }, error = function(e) {
    stop("Wasserstein error: ", e$message)
  })
}


############################################
# Wrapper
############################################

#' Compare depth distributions using multiple metrics
#'
#' Unified interface for comparing two probability distributions using
#' Wasserstein distance, Jensen-Shannon divergence, or KL divergence,
#' with optional depth-dependent weighting.
#'
#' @param P Numeric probability vector.
#' @param Q Numeric probability vector.
#' @param metric Character. One of \code{"wasserstein"} (default),
#' \code{"js"}, or \code{"kl"}.
#' @param tail_weight Logical. If TRUE, enables depth-dependent weighting.
#' @param weighting Character. Weighting function used when
#' \code{tail_weight = TRUE}.
#' @param parameter List with weighting parameters \code{lambda} and \code{x0}.
#' @param inverse Logical. Reverse weighting direction.
#' @param index Numeric positions of the bins (e.g. depths), one per bin
#'   after alignment. \code{NULL} (default) uses \code{1, 2, ..., n}.
#' @param index.spacing Character. \code{"equal"} or \code{"custom"}.
#' @param alignPQ Logical. Aligns unequal-length vectors; if FALSE, unequal
#'   lengths are an error.
#' @param cut Logical. If TRUE, truncates the longer vector instead of padding
#'   the shorter one with zeros.
#' @param baseline.weight Numeric in `[0,1]`. Minimum weight applied.
#'
#' @section Data requirements:
#'
#' \code{P} and \code{Q} are non-negative vectors representing discrete
#' probability mass functions over a shared ordered index (e.g. depth bins).
#'
#' They are internally normalized to sum to 1 if necessary.
#'
#' Typical construction:
#'
#' \preformatted{
#' counts_P <- c(5, 10, 20, 15)
#' counts_Q <- c(2, 8, 25, 20)
#'
#' P <- counts_P / sum(counts_P)
#' Q <- counts_Q / sum(counts_Q)
#' }
#'
#' Continuous data must be discretised (e.g. binning) before use.
#'
#' @section Alignment:
#'
#' If \code{P} and \code{Q} differ in length, they are made comparable as follows:
#'
#' \strong{cut = TRUE:} both vectors are truncated to the shared minimum length.
#'
#' \strong{cut = FALSE} (default): the shorter vector is padded with zeros.
#'
#' The same rule applies to all three metrics. No interpolation or re-binning
#' is performed. Alignment only enforces equal vector length; consistent bin
#' definitions are assumed. Padding gives \code{"kl"} a value of \code{Inf}
#' whenever \code{P} has mass in a bin that was padded onto \code{Q}, since
#' \eqn{P \log(P / 0)} is infinite; \code{"js"} and \code{"wasserstein"} stay
#' finite.
#'
#' @section Interpretation:
#'
#' Elements of \code{P} and \code{Q} correspond to probability mass at
#' positions given by \code{index}. The same index must apply to both
#' distributions for valid comparison.
#'
#' @details
#' This function wraps:
#' \itemize{
#'   \item \code{tail_weighted_wasserstein_distance()}
#'   \item \code{tail_weighted_js_divergence()}
#'   \item \code{tail_weighted_kl_divergence()}
#' }
#'
#' If \code{tail_weight = FALSE}, uniform weighting is used.
#'
#' @return Numeric distance or divergence value.
#'
#' @examples
#' counts_P <- c(5, 10, 20, 15)
#' counts_Q <- c(2, 8, 25, 20)
#'
#' P <- counts_P / sum(counts_P)
#' Q <- counts_Q / sum(counts_Q) 
#' compare_depth_distribution(P, Q)
#' compare_depth_distribution(P, Q, metric = "js")
#' compare_depth_distribution(P, Q, tail_weight = TRUE, weighting = "sigmoid")
#'
#' @export
compare_depth_distribution <- function(
    P, Q,
    metric = "wasserstein",
    tail_weight = FALSE,
    weighting = "sigmoid",
    parameter = list(lambda = 0.2, x0 = 30),
    inverse = FALSE,
    index = NULL,
    index.spacing = "equal",
    alignPQ = TRUE,
    cut = FALSE,
    baseline.weight = 0
) {
  
  metric <- match.arg(metric, c("wasserstein", "js", "kl"))
  
  if (!tail_weight) weighting <- "constant"
  
  switch(
    metric,
    
    wasserstein = tail_weighted_wasserstein_distance(
      P, Q,
      index = index,
      index.spacing = index.spacing,
      parameter = parameter,
      weighting = weighting,
      baseline.weight = baseline.weight,
      inverse = inverse,
      alignPQ = alignPQ,
      cut = cut
    ),
    
    js = tail_weighted_js_divergence(
      P, Q,
      index = index,
      index.spacing = index.spacing,
      parameter = parameter,
      weighting = weighting,
      baseline.weight = baseline.weight,
      inverse = inverse,
      alignPQ = alignPQ,
      cut = cut
    ),
    
    kl = tail_weighted_kl_divergence(
      P, Q,
      index = index,
      index.spacing = index.spacing,
      parameter = parameter,
      weighting = weighting,
      baseline.weight = baseline.weight,
      inverse = inverse,
      alignPQ = alignPQ,
      cut = cut
    )
  )
}








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
