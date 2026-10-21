#' Root each component and orient every edge parent -> child
#'
#' A skeleton is not rooted: the base of a root system is a free end exactly
#' like a tip. Each component is therefore rooted at its thickest degree-1
#' node, which on a scan is the proximal end of the main axis. This is the same
#' convention \code{compute_strahler_order()} uses, so orders and transport
#' quantities refer to the same base.
#'
#' @param et Edge table with \code{from}/\code{to}.
#' @return list with per-node parent edge/node, BFS order, per-edge child and
#'   parent node, component labels, degrees, and any edges that close a cycle.
#' @keywords internal
#' @noRd
.root_forest <- function(et) {
  nodes <- unique(c(et$from, et$to))
  ni <- stats::setNames(seq_along(nodes), nodes)
  fi <- as.integer(ni[et$from]); ti <- as.integer(ni[et$to])
  V <- length(nodes); E <- nrow(et)

  # adjacency, built with split() rather than by growing a list per edge
  node  <- c(fi, ti); other <- c(ti, fi); edge <- rep(seq_len(E), 2L)
  o <- order(node)
  lv <- factor(node[o], levels = seq_len(V))
  nb <- split(other[o], lv); eg <- split(edge[o], lv)

  comp <- .net_components(fi, ti, V)
  deg  <- tabulate(c(fi, ti), nbins = V)
  w <- if (!is.null(et$mean_diameter)) et$mean_diameter else
       if (!is.null(et$length)) et$length else rep(1, E)
  w[!is.finite(w)] <- 0

  parent_edge <- integer(V); parent_node <- integer(V); depth <- integer(V)
  root_of <- integer(V); seen <- logical(V); used <- logical(E)
  bfs <- integer(V); nb_seen <- 0L

  for (cp in sort(unique(comp))) {
    members <- which(comp == cp)
    free <- members[deg[members] == 1L]
    s <- if (length(free)) {
      free[which.max(vapply(free, function(v) max(w[eg[[v]]]), numeric(1)))]
    } else members[1L]                       # a pure cycle has no free end
    seen[s] <- TRUE; root_of[s] <- s
    q <- integer(length(members)); q[1L] <- s; qh <- 1L; qt <- 1L
    while (qh <= qt) {
      u <- q[qh]; qh <- qh + 1L
      nb_seen <- nb_seen + 1L; bfs[nb_seen] <- u
      nbs <- nb[[u]]; egs <- eg[[u]]
      for (k in seq_along(nbs)) {
        e <- egs[k]; v <- nbs[k]
        if (used[e] || seen[v]) next
        seen[v] <- TRUE; used[e] <- TRUE
        parent_edge[v] <- e; parent_node[v] <- u
        depth[v] <- depth[u] + 1L; root_of[v] <- s
        qt <- qt + 1L; q[qt] <- v
      }
    }
  }

  child_of_edge <- integer(E); parent_of_edge <- integer(E)
  for (v in seq_len(V)) {
    e <- parent_edge[v]
    if (e != 0L) { child_of_edge[e] <- v; parent_of_edge[e] <- parent_node[v] }
  }
  list(nodes = nodes, fi = fi, ti = ti, V = V, E = E, comp = comp, deg = deg,
       parent_edge = parent_edge, parent_node = parent_node, depth = depth,
       root_of = root_of, bfs = bfs[seq_len(nb_seen)],
       child_of_edge = child_of_edge, parent_of_edge = parent_of_edge,
       cycle_edges = which(!used))
}


#' Per-segment transport quantities
#'
#' Adds, for every segment, how far and through how much resistance water would
#' travel from the base of its root system to that segment, and how much root
#' the segment serves. The result is a per-segment vector, so it can be painted
#' onto the image with \code{\link{render_order_overlay}}.
#'
#' @section What the columns mean:
#' \describe{
#'   \item{\code{resistance}}{Hagen-Poiseuille resistance of the segment alone,
#'     \eqn{8 \eta L / (\pi r^4)}. A geometric proxy: it treats the root as an
#'     open pipe, where real resistance is set by xylem anatomy. The viscosity
#'     cancels in every comparison, so the default \code{eta = 1} leaves the
#'     values relative.}
#'   \item{\code{path_length}, \code{path_resistance}}{Cumulative from the base
#'     of the component to the far end of this segment.}
#'   \item{\code{downstream_tips}, \code{downstream_length},
#'     \code{downstream_surface}}{How many tips, how much root length, and how
#'     much lateral surface area lie beyond this segment -- what it serves.
#'     These peak at the base.}
#' }
#'
#' Lengths and diameters are used in whatever unit the edge table carries, so
#' run \code{\link{convert_root_units}} first if you want centimeters.
#'
#' @param et Edge table from \code{\link{root_graph_pipeline}} or
#'   \code{\link{branch_order_map}}, needing \code{from}, \code{to},
#'   \code{length} and \code{mean_diameter}.
#' @param eta Dynamic viscosity. Cancels in any ratio; change it only if you
#'   want resistance in physical units.
#' @param min_radius Radius floor, applied before the fourth power so a
#'   zero-radius segment cannot give infinite resistance. \code{NULL}
#'   (default) uses the smallest positive radius in the table, so the floor
#'   carries whatever unit the table is in; a segment with no usable diameter
#'   is then treated as the thinnest root present. Give a number to set the
#'   floor yourself, in the table's length unit.
#' @return \code{et} with the six columns above added. Segments on a cycle get
#'   \code{NA} for the downstream columns, since they have no subtree, and a
#'   warning names how many there were. A component with no free end at all --
#'   a pure ring, which a root system does not produce -- has no proximal end
#'   to root at, so an arbitrary node is used and its \code{path_*} values are
#'   measured from a meaningless start.
#' @export
#' @seealso \code{\link{render_order_overlay}} to paint any of these columns.
#' @examples
#' \dontrun{
#' res <- branch_order_map(skel, mask, unit = "cm", dpi = 300)
#' et  <- root_transport(res$edges)
#'
#' # a transport map: color every segment by its path resistance
#' render_order_overlay(attr(et, "segments"), et$path_resistance,
#'                      attr(et, "dims"), "transport.png",
#'                      value_label = "path resistance")
#' }
root_transport <- function(et, eta = 1, min_radius = NULL) {
  if (is.null(et) || !is.data.frame(et) || nrow(et) == 0L) return(et)
  need <- c("from", "to", "length", "mean_diameter")
  miss <- setdiff(need, names(et))
  if (length(miss))
    stop(sprintf("`et` is missing: %s", paste(miss, collapse = ", ")), call. = FALSE)

  if (anyNA(et$from) || anyNA(et$to))
    stop("`from` and `to` must not contain NA.", call. = FALSE)

  E <- nrow(et)
  r <- et$mean_diameter / 2
  bad <- !is.finite(r) | r <= 0
  if (is.null(min_radius)) {
    # The table may be in pixels, centimeters or inches, so a fixed floor would
    # mean something different in each. The thinnest root actually measured is
    # in the right unit by construction.
    min_radius <- if (all(bad)) NA_real_ else min(r[!bad])
  }
  r[bad | r < min_radius] <- min_radius
  et$resistance <- 8 * eta * et$length / (pi * r^4)

  rf <- .root_forest(et)
  V <- rf$V
  len <- et$length; res <- et$resistance
  surf <- pi * et$mean_diameter * et$length

  # downward pass: accumulate from each component's base
  pl <- rep(NA_real_, V); pr <- rep(NA_real_, V)
  for (u in rf$bfs) {
    pe <- rf$parent_edge[u]
    if (pe == 0L) { pl[u] <- 0; pr[u] <- 0 }
    else {
      pn <- rf$parent_node[u]
      pl[u] <- pl[pn] + len[pe]; pr[u] <- pr[pn] + res[pe]
    }
  }

  # upward pass: a node's subtree totals are finished only after its children,
  # which reversing the BFS order guarantees
  kids <- vector("list", V)
  for (e in seq_len(E)) {
    cv <- rf$child_of_edge[e]
    if (cv != 0L) kids[[rf$parent_of_edge[e]]] <- c(kids[[rf$parent_of_edge[e]]], e)
  }
  sub_tips <- integer(V); sub_len <- numeric(V); sub_surf <- numeric(V)
  for (u in rev(rf$bfs)) {
    ke <- kids[[u]]
    if (!length(ke)) { sub_tips[u] <- 1L; next }
    for (e in ke) {
      cv <- rf$child_of_edge[e]
      sub_tips[u] <- sub_tips[u] + sub_tips[cv]
      sub_len[u]  <- sub_len[u]  + len[e]  + sub_len[cv]
      sub_surf[u] <- sub_surf[u] + surf[e] + sub_surf[cv]
    }
  }

  cv <- rf$child_of_edge
  has <- cv != 0L
  et$path_length        <- ifelse(has, pl[pmax(cv, 1L)], NA_real_)
  et$path_resistance    <- ifelse(has, pr[pmax(cv, 1L)], NA_real_)
  et$downstream_tips    <- ifelse(has, sub_tips[pmax(cv, 1L)], NA_integer_)
  et$downstream_length  <- ifelse(has, len + sub_len[pmax(cv, 1L)], NA_real_)
  et$downstream_surface <- ifelse(has, surf + sub_surf[pmax(cv, 1L)], NA_real_)

  if (length(rf$cycle_edges))
    warning(sprintf(paste("%d segment(s) lie on a cycle and have no subtree;",
                          "their downstream_* and path_* values are NA. On a root",
                          "scan a cycle is usually two roots crossing -- see",
                          "network_metrics(..., 'projected_tree')$qc."),
                    length(rf$cycle_edges)), call. = FALSE)
  et
}
