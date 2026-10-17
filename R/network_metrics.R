#' Connected components of an edge list
#'
#' Union-find with path compression. Iterative, so a long chain cannot exhaust
#' the evaluation stack the way a recursive \code{find()} would.
#'
#' @param fi,ti Integer node indices, one pair per edge.
#' @param n Number of nodes.
#' @return Integer vector of component labels, one per node.
#' @keywords internal
#' @noRd
.net_components <- function(fi, ti, n) {
  uf <- seq_len(n)
  find <- function(x) {
    r <- x
    while (uf[r] != r) r <- uf[r]
    while (uf[x] != r) { nx <- uf[x]; uf[x] <<- r; x <- nx }
    r
  }
  for (e in seq_along(fi)) {
    a <- find(fi[e]); b <- find(ti[e])
    if (a != b) uf[a] <- b
  }
  as.integer(factor(vapply(seq_len(n), find, integer(1))))
}


#' Topology metrics for a spatial network
#'
#' Summarizes an edge table as a graph: node and edge counts, components, tips,
#' branch points, and the connectivity coefficients used in network ecology.
#' Works on any network whose segments are listed as node pairs, not only root
#' systems.
#'
#' @section Why the topology has to be declared:
#' The same quantity means three different things depending on whether loops in
#' the data are real.
#'
#' \describe{
#'   \item{\code{"reticulate"}}{Loops are biology: fungal mycelium with
#'     anastomoses, leaf venation, the endoplasmic reticulum, vasculature. The
#'     cyclomatic number and the alpha coefficient measure redundancy, so they
#'     are returned as traits.}
#'   \item{\code{"tree"}}{Loops are impossible: a root system, a bronchial tree,
#'     a river network. Then \eqn{E = V - G} holds by definition, so the
#'     cyclomatic number is identically zero and the alpha coefficient carries
#'     no information. They are returned under \code{$qc}, where a value above
#'     zero means the input is not actually a forest.}
#'   \item{\code{"projected_tree"}}{A tree photographed in two dimensions, which
#'     is what a root scan is. Every loop is two roots crossing in projection,
#'     so the cyclomatic number counts unresolved overlaps and belongs under
#'     \code{$qc}. \code{artefact_rate} normalizes that count by network length
#'     so it is comparable between images.}
#' }
#'
#' Reporting a meshedness of 0.002 as a biological result, when it is a count of
#' root crossings the segmentation failed to separate, is the mistake this
#' separation is meant to prevent.
#'
#' @param edges An edge table with character columns \code{from} and \code{to},
#'   as returned by \code{\link{root_graph_pipeline}}, or the object returned by
#'   \code{\link{branch_order_map}}, whose \code{$edges} element is used.
#'   \code{length} is used for \code{total_length} and for
#'   \code{artefact_rate} when present.
#' @param topology One of \code{"auto"}, \code{"projected_tree"},
#'   \code{"tree"}, \code{"reticulate"}. \code{"auto"} calls a network with no
#'   cycles a tree, and anything else a projected tree, because an edge table
#'   alone cannot distinguish a real loop from a crossing. Declare the topology
#'   whenever you know it.
#' @return A list with \code{topology}, \code{traits}, \code{qc} and
#'   \code{note}. An empty or zero-row input returns the same structure with
#'   empty \code{traits} rather than an error.
#' @export
#' @seealso \code{\link{root_graph_pipeline}}, \code{\link{branch_order_map}}
#' @examples
#' # A three-segment tree: one basal edge splitting into two laterals
#' et <- data.frame(from = c("A", "B", "B"), to = c("B", "C", "D"),
#'                  length = c(10, 6, 7), stringsAsFactors = FALSE)
#' network_metrics(et, topology = "tree")$traits
#'
#' # The same edges plus one closing a loop: a crossing on a root scan
#' et2 <- rbind(et, data.frame(from = "C", to = "D", length = 4))
#' network_metrics(et2, topology = "projected_tree")$qc
network_metrics <- function(edges,
                            topology = c("auto", "projected_tree",
                                         "tree", "reticulate")) {
  topology <- match.arg(topology)

  if (inherits(edges, "branchOrderMap")) edges <- edges$edges
  empty <- list(topology = topology, traits = list(), qc = list(),
                note = "empty network: no edges to summarize")
  if (is.null(edges) || !is.data.frame(edges) || nrow(edges) == 0L) return(empty)
  if (!all(c("from", "to") %in% names(edges)))
    stop("`edges` needs `from` and `to` columns naming the nodes of each segment.",
         call. = FALSE)

  from <- as.character(edges$from); to <- as.character(edges$to)
  if (anyNA(from) || anyNA(to))
    stop("`from` and `to` must not contain NA.", call. = FALSE)

  nodes <- unique(c(from, to))
  idx <- stats::setNames(seq_along(nodes), nodes)
  fi <- as.integer(idx[from]); ti <- as.integer(idx[to])
  V <- length(nodes); E <- nrow(edges)
  G <- length(unique(.net_components(fi, ti, V)))

  # A self-loop adds two to its node's degree, which is the usual convention
  # and keeps the handshake identity sum(degree) = 2E true.
  deg <- tabulate(c(fi, ti), nbins = V)
  n_tips <- sum(deg == 1L); n_branch <- sum(deg >= 3L)

  has_len <- !is.null(edges$length) && is.numeric(edges$length)
  total_len <- if (has_len) sum(edges$length, na.rm = TRUE) else NA_real_

  mu <- E - V + G                                   # cyclomatic number
  # Both coefficients are ratios to the maximum a planar graph can reach. Those
  # maxima are only positive once the graph is big enough to have them, so the
  # denominators are guarded rather than allowed to go negative.
  alpha <- if (2 * V - 5 > 0) mu / (2 * V - 5) else NA_real_
  gamma <- if (3 * V - 6 > 0) E / (3 * V - 6) else NA_real_

  if (topology == "auto")
    topology <- if (mu == 0L) "tree" else "projected_tree"

  base <- list(
    n_nodes = V, n_edges = E, n_components = G,
    n_tips = n_tips, n_branch_points = n_branch,
    tip_ratio = n_tips / V,
    beta = E / V,
    total_length = total_len,
    mean_edge_length = if (has_len) mean(edges$length, na.rm = TRUE) else NA_real_
  )
  loops <- list(cyclomatic = mu, alpha_meshedness = alpha, gamma = gamma)

  # Two inputs that inflate the cycle count without being overlaps, reported so
  # a surprising cyclomatic number can be traced rather than guessed at.
  dup <- sum(duplicated(t(apply(cbind(pmin(fi, ti), pmax(fi, ti)), 1L, sort))))
  qc_extra <- list(n_self_loops = sum(fi == ti), n_duplicate_edges = dup)

  if (topology == "reticulate")
    return(list(topology = topology, traits = c(base, loops), qc = qc_extra,
                note = "loops counted as structure: cyclomatic number and alpha are traits"))

  if (topology == "tree")
    return(list(topology = topology, traits = base,
                qc = c(loops, qc_extra, list(is_forest = mu == 0L)),
                note = paste("on a tree E = V - G, so the cyclomatic number and alpha are",
                             "identities rather than measurements; a value above zero means",
                             "the input is not a forest")))

  list(topology = topology, traits = base,
       qc = c(loops, qc_extra, list(
         is_forest = mu == 0L,
         artefact_cycles = mu,
         artefact_rate = if (is.finite(total_len) && total_len > 0)
           mu / total_len else NA_real_)),
       note = paste("loops in a projected tree are crossings the segmentation did not",
                    "separate; artefact_rate divides that count by network length, in",
                    "whatever unit the edge table carries, so images can be compared"))
}
