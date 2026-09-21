# ============================================================
# ROOT GRAPH PIPELINE  (SKELETON -> SEGMENT GRAPH -> TIP-ORDER METRICS)
# + junction-cluster contraction + order-colored validation PNG
# ------------------------------------------------------------
# For FRAGMENTED ("broken") roots: no base/origin assumed; each
# component handled independently; tip_order via leaf-peeling.
#
# Contraction: connected clumps of junction pixels (degree >= 3)
# are merged into a single graph node BEFORE tracing. This
# dissolves the 2-px "ladder" segments and spurious micro-cycles
# that an imperfectly-thinned skeleton produces (which otherwise
# leave tip_order = NA).
#
# Deps: base R + grDevices/graphics. terra only for SpatRaster
#       input. imager used for the DT if installed (else pure-R fallback).
# ============================================================


#' Coerce an image to a 0/1 integer matrix
#'
#' @param x A single-layer \code{SpatRaster}, matrix, or data.frame.
#' @return A 0/1 integer matrix (NA treated as 0).
#' @keywords internal
#' @noRd
.to_binary_matrix <- function(x) {
  if (inherits(x, "SpatRaster")) {
    if (terra::nlyr(x) > 1L) x <- x[[1L]]
    m <- terra::as.matrix(x, wide = TRUE)
  } else if (is.matrix(x)) {
    m <- x
  } else if (is.data.frame(x)) {
    m <- as.matrix(x)
  } else stop("Provide a matrix or single-layer SpatRaster.")
  m[is.na(m)] <- 0
  (m > 0) * 1L
}


#' Flip a matrix for image-style (x right, y up) plotting with graphics::image
#'
#' @param m A matrix.
#' @return \code{m} transposed and with rows reversed.
#' @keywords internal
#' @noRd
.flip_matrix <- function(m) t(apply(m, 2L, rev))


#' Euclidean distance transform (distance to nearest background)
#'
#' Dispatches to \pkg{imager} when available, falling back to an exact pure-R
#' separable transform.
#'
#' @param mask 0/1 matrix; distance is measured from each pixel to the nearest 0.
#' @param backend One of \code{"auto"}, \code{"imager"}, \code{"baseR"}.
#' @return Numeric matrix of Euclidean distances (inscribed radius at root pixels).
#' @keywords internal
#' @noRd
.distance_transform_edt <- function(mask, backend = "auto") {
  nr <- nrow(mask); nc <- ncol(mask)
  if (all(mask > 0)) warning("DT: no background pixels; radii unbounded.")
  res <- NULL
  if (backend %in% c("auto", "imager") && requireNamespace("imager", quietly = TRUE)) {
    res <- tryCatch({
      d <- imager::distance_transform(imager::as.cimg(mask * 1.0), value = 0, metric = 2L)
      matrix(as.numeric(d), nr, nc)
    }, error = function(e) NULL)
  }
  if (is.null(res)) res <- .distance_transform_edt_baseR(mask)
  res
}

#' One-dimensional squared Euclidean distance transform
#'
#' Felzenszwalb & Huttenlocher lower-envelope, O(n); used by the base-R EDT.
#' @keywords internal
#' @noRd
.edt_sq_1d <- function(f) {
  n <- length(f); if (n <= 1L) return(f)
  d <- numeric(n); v <- integer(n); v[1L] <- 1L
  z <- numeric(n + 1L); z[1L] <- -Inf; z[2L] <- Inf; k <- 1L
  for (q in 2:n) {
    repeat { s <- ((f[q] + q*q) - (f[v[k]] + v[k]*v[k])) / (2*q - 2*v[k]); if (s > z[k]) break; k <- k - 1L }
    k <- k + 1L; v[k] <- q; z[k] <- s; z[k + 1L] <- Inf
  }
  k <- 1L
  for (q in 1:n) { while (z[k + 1L] < q) k <- k + 1L; dq <- q - v[k]; d[q] <- dq*dq + f[v[k]] }
  d
}

#' Exact pure-R Euclidean distance transform (separable, O(N))
#'
#' Fallback used when no compiled backend is installed.
#' @keywords internal
#' @noRd
.distance_transform_edt_baseR <- function(mask) {
  nr <- nrow(mask); nc <- ncol(mask); BIG <- (nr + nc)^2 + 1
  f <- matrix(ifelse(mask > 0, BIG, 0), nr, nc)
  for (cc in seq_len(nc)) f[, cc] <- .edt_sq_1d(f[, cc])
  for (rr in seq_len(nr)) f[rr, ] <- .edt_sq_1d(f[rr, ])
  sqrt(f)
}


#' Warn when the traced graph covers little of the skeleton
#'
#' Junction contraction only dissolves a handful of pixels per branch point, so
#' a real 1-px skeleton traces to >95% coverage. A low value almost always means
#' the "skeleton" is actually a filled mask (or a badly thinned one): nearly
#' every pixel then has degree >= 3, is contracted away, and the pipeline
#' silently returns a near-empty graph.
#'
#' @param segs Traced segment list.
#' @param skel The 0/1 skeleton matrix the segments were traced from.
#' @param min_frac Coverage below which to warn.
#' @return Invisibly, the covered fraction.
#' @keywords internal
#' @noRd
.check_skeleton_coverage <- function(segs, skel, min_frac = 0.75) {
  npx <- sum(skel > 0)
  if (npx == 0L) return(invisible(NA_real_))
  hit <- matrix(FALSE, nrow(skel), ncol(skel))
  for (s in segs) hit[s$coords] <- TRUE
  frac <- sum(hit) / npx
  if (frac < min_frac)
    warning(sprintf(paste0("Only %.1f%% of skeleton pixels were traced into segments. ",
                           "Is 'skel' a 1-pixel-wide skeleton? Passing a filled mask ",
                           "as 'skel' yields a near-empty graph."), 100 * frac),
            call. = FALSE)
  invisible(frac)
}


#' Trace skeleton segments with junction-cluster contraction
#'
#' Walks degree-2 pixel chains between nodes (tips and junctions). Connected
#' clumps of junction pixels (degree >= 3) are contracted to a single graph node
#' so that imperfect 1-px thinning does not fragment the graph.
#'
#' @param skel Binary skeleton (matrix or single-layer \code{SpatRaster}).
#' @return A list of segments; each is \code{list(coords, from, to)} where
#'   \code{coords} is an Nx2 (row, col) matrix and \code{from}/\code{to} are node
#'   labels (\code{"J<id>"} for junction clusters, \code{"<row>_<col>"} for tips).
#'   \code{attr(., "dims")} holds the skeleton dimensions and
#'   \code{attr(., "node_xy")} the (row, col) centroid of each junction cluster.
#' @keywords internal
#' @noRd
trace_segments <- function(skel) {
  skel <- (.to_binary_matrix(skel) > 0) * 1L
  nr <- nrow(skel); nc <- ncol(skel)
  lin <- which(skel > 0); np <- length(lin)
  if (np == 0L) {
    out <- list(); attr(out, "dims") <- c(nr, nc)
    attr(out, "node_xy") <- matrix(numeric(0), 0L, 2L,
                                   dimnames = list(NULL, c("row", "col")))
    return(out)
  }
  
  id_mat <- integer(nr * nc); id_mat[lin] <- seq_len(np); dim(id_mat) <- c(nr, nc)
  rr0 <- ((lin - 1L) %% nr) + 1L
  cc0 <- ((lin - 1L) %/% nr) + 1L
  
  doff  <- c(-1L,-1L, 0L, 1L, 1L, 1L, 0L,-1L)
  coff  <- c( 0L, 1L, 1L, 1L, 0L,-1L,-1L,-1L)
  rev_d <- c(  5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L)
  
  NB <- matrix(0L, np, 8L)
  for (d in 1:8) {
    r2 <- rr0 + doff[d]; c2 <- cc0 + coff[d]
    ok <- r2 >= 1L & r2 <= nr & c2 >= 1L & c2 <= nc
    nb <- integer(np); nb[ok] <- id_mat[(c2[ok] - 1L) * nr + r2[ok]]
    NB[, d] <- nb
  }
  deg <- rowSums(NB > 0L)
  
  # union-find over adjacent junction (deg>=3) pixels.
  # NOTE: 'uf' (not 'parent') to avoid clashing with any package's parent().
  # Inside find(), <<- reaches this frame's 'uf'; at body level we use plain <-.
  uf <- seq_len(np)
  find <- function(x) {
    r <- x; while (uf[r] != r) r <- uf[r]
    while (uf[x] != r) { nx <- uf[x]; uf[x] <<- r; x <- nx }; r
  }
  is_j <- deg >= 3L
  for (p in which(is_j)) for (d in 1:8) {
    q <- NB[p, d]
    if (q > 0L && is_j[q]) { a <- find(p); b <- find(q); if (a != b) uf[a] <- b }
  }
  
  # Absorb degree-2 pixels that bridge two pixels of the SAME junction cluster.
  # A 90-degree bend in an 8-connected skeleton leaves such a pixel wedged
  # between two contracted junction pixels; tracing it produced a 3-px J->J
  # self-loop (a spurious micro-cycle: NA tip_order + a phantom branch point).
  # Folding it into the cluster removes the artifact at its source.
  d2 <- which(deg == 2L)
  if (length(d2)) {
    pos <- NB[d2, , drop = FALSE] > 0L
    i1  <- max.col(pos, ties.method = "first")
    pos[cbind(seq_along(d2), i1)] <- FALSE
    i2  <- max.col(pos, ties.method = "first")
    n1  <- NB[d2, , drop = FALSE][cbind(seq_along(d2), i1)]
    n2  <- NB[d2, , drop = FALSE][cbind(seq_along(d2), i2)]
    repeat {
      cand <- which(!is_j[d2] & is_j[n1] & is_j[n2])
      if (!length(cand)) break
      r1 <- vapply(n1[cand], find, integer(1L))
      r2 <- vapply(n2[cand], find, integer(1L))
      take <- cand[r1 == r2]
      if (!length(take)) break
      for (i in take) { p <- d2[i]; is_j[p] <- TRUE; uf[p] <- find(n1[i]) }
    }
  }
  
  node_label <- function(p) if (is_j[p]) paste0("J", find(p)) else paste0(rr0[p], "_", cc0[p])
  
  consumed <- matrix(FALSE, np, 8L)
  segs <- list()
  # start at every node pixel: tips/isolated (deg != 2) and every junction pixel
  # (which after absorption may itself have degree 2).
  for (p in which(deg != 2L | is_j)) {
    for (d in 1:8) {
      q <- NB[p, d]
      if (q == 0L || consumed[p, d]) next
      if (is_j[q]) { consumed[p, d] <- TRUE; next }   # junction<->junction: contracted away
      consumed[p, d] <- TRUE
      
      buf <- integer(64L); buf[1L] <- p; buf[2L] <- q; len <- 2L
      prev <- p; cur <- q; din <- rev_d[d]
      while (!is_j[cur] && deg[cur] == 2L) {
        dirs <- which(NB[cur, ] > 0L)
        out_dir <- dirs[dirs != din][1L]
        nxt <- NB[cur, out_dir]
        if (len == length(buf)) buf <- c(buf, integer(len))
        len <- len + 1L; buf[len] <- nxt
        prev <- cur; cur <- nxt; din <- rev_d[out_dir]
      }
      consumed[cur, din] <- TRUE
      ids <- buf[seq_len(len)]
      segs[[length(segs) + 1L]] <- list(
        coords = cbind(rr0[ids], cc0[ids]),
        from = node_label(p), to = node_label(cur)
      )
    }
  }
  
  # Junction-cluster centroids. Contraction dissolves the cluster interior, so
  # the traced arms stop on the cluster rim; build_edge_table() re-attaches the
  # rim-to-centroid stub, which is otherwise lost from every root length.
  node_xy <- matrix(numeric(0), 0L, 2L, dimnames = list(NULL, c("row", "col")))
  jp <- which(is_j)
  if (length(jp)) {
    jr  <- vapply(jp, find, integer(1L))
    cr  <- tapply(rr0[jp], jr, mean)
    cc_ <- tapply(cc0[jp], jr, mean)
    node_xy <- cbind(row = as.numeric(cr), col = as.numeric(cc_))
    rownames(node_xy) <- paste0("J", names(cr))
  }
  
  attr(segs, "dims")    <- c(nr, nc)
  attr(segs, "node_xy") <- node_xy
  segs
}


#' Distance from a segment's terminal pixel to its junction-cluster center
#'
#' Junction contraction stops the traced arms on the cluster rim, so this stub is
#' the piece of root length the contraction removed. Tip nodes are not in
#' \code{node_xy} (the pixel already is the node) and contribute nothing.
#'
#' @param lab Node label of that segment end.
#' @param px Length-2 (row, col) of the segment's terminal pixel.
#' @param node_xy Junction-cluster centroids, as attached by \code{trace_segments}.
#' @return A single non-negative length in pixels.
#' @keywords internal
#' @noRd
.stub_length <- function(lab, px, node_xy) {
  if (is.null(node_xy) || nrow(node_xy) == 0L) return(0)
  k <- match(lab, rownames(node_xy))
  if (is.na(k)) return(0)                         # tip node: the pixel is the node
  sqrt((node_xy[k, 1L] - px[1L])^2 + (node_xy[k, 2L] - px[2L])^2)
}


#' Polyline length of one segment, junction stubs included
#'
#' The same definition \code{build_edge_table()} reports as \code{length}, so
#' anything that thresholds on length (pruning) and the table the user reads
#' cannot drift apart.
#'
#' @param s One segment, \code{list(coords, from, to)}.
#' @param node_xy Junction-cluster centroids.
#' @return Length in pixels.
#' @keywords internal
#' @noRd
.segment_length <- function(s, node_xy) {
  p <- s$coords
  sum(sqrt(rowSums(diff(p)^2))) +
    .stub_length(s$from, p[1L, ], node_xy) +
    .stub_length(s$to,   p[nrow(p), ], node_xy)
}


#' Map each node label to the (segment, side) pairs incident on it
#'
#' Three passes of the pipeline need this table. Building it by \code{rbind()}
#' into a growing list slot re-allocated on every segment end; \code{split()}
#' builds it in one pass.
#'
#' Label order is the order of first occurrence, not sorted: the union-find in
#' \code{assign_root_order()} merges roots in the order it walks the nodes, so a
#' different order would renumber \code{root_id}.
#'
#' @param segs Segment list.
#' @return Named list; each element an Nx2 matrix of \code{(segment, side)} rows,
#'   side 1 = \code{from}, side 2 = \code{to}.
#' @keywords internal
#' @noRd
.node_incidence <- function(segs) {
  ns <- length(segs)
  if (ns == 0L) return(list())
  ends <- character(2L * ns)
  ends[seq(1L, 2L * ns, by = 2L)] <- vapply(segs, function(s) s$from, character(1))
  ends[seq(2L, 2L * ns, by = 2L)] <- vapply(segs, function(s) s$to,   character(1))
  seg  <- rep(seq_len(ns), each = 2L)
  side <- rep(c(1L, 2L), times = ns)
  idx  <- split(seq_len(2L * ns), factor(ends, levels = unique(ends)))
  lapply(idx, function(k) cbind(seg[k], side[k]))
}


#' Per-segment measurements from traced segments
#'
#' @param segs Output of \code{trace_segments} / \code{\link{resolve_crossings}}.
#' @param DT Distance-transform matrix on the (cropped) mask.
#' @param node_xy Junction-cluster centroids, as attached by
#'   \code{trace_segments()}. Each segment end that lands on a contracted
#'   junction is extended to its cluster centroid, which restores the
#'   cluster-interior length that contraction removes. Defaults to the
#'   \code{"node_xy"} attribute of \code{segs}; pass \code{NULL} for the raw
#'   rim-to-rim length.
#' @return A data.frame with \code{edge_id}, \code{from}, \code{to}, \code{n_px},
#'   \code{n_orth}, \code{n_diag}, \code{length} (sqrt(2) polyline, px),
#'   \code{length_kimura} (px), and \code{mean}/\code{median}/\code{min_diameter} (px).
#' @keywords internal
#' @noRd
build_edge_table <- function(segs, DT, node_xy = attr(segs, "node_xy")) {
  n <- length(segs)
  if (n == 0L) return(NULL)

  # Filling plain vectors and calling data.frame() once at the end, rather than
  # building one data.frame per segment and rbind()-ing them, is what keeps this
  # from dominating the pipeline: on a root system of a few thousand segments the
  # per-segment version cost more than the tracing it summarizes.
  from <- character(n); to <- character(n)
  n_px <- integer(n); n_orth <- integer(n); n_diag <- integer(n)
  len_poly <- numeric(n); len_kimura <- numeric(n)
  d_mean <- numeric(n); d_median <- numeric(n); d_min <- numeric(n)

  for (i in seq_len(n)) {
    s <- segs[[i]]; p <- s$coords
    dr <- diff(p[, 1]); dc <- diff(p[, 2]); seglen <- sqrt(dr*dr + dc*dc)
    adr <- abs(dr); adc <- abs(dc)
    orth <- (adr + adc) == 1L            # horizontal/vertical unit step
    diag <- (adr == 1L & adc == 1L)      # diagonal unit step
    no <- sum(orth); nd <- sum(diag)
    gap <- sum(seglen[!(orth | diag)])   # jumps across contracted-junction gaps
    ends <- .stub_length(s$from, p[1L, ], node_xy) +
            .stub_length(s$to,   p[nrow(p), ], node_xy)
    rad <- DT[p]

    from[i] <- s$from; to[i] <- s$to
    n_px[i] <- nrow(p); n_orth[i] <- no; n_diag[i] <- nd
    len_poly[i]   <- sum(seglen) + ends                              # sqrt(2) chain code
    len_kimura[i] <- sqrt(nd^2 + (nd + no/2)^2) + no/2 + gap + ends  # Kimura per segment
    d_mean[i]   <- mean(2*rad)
    d_median[i] <- stats::median(2*rad)
    d_min[i]    <- min(2*rad)
  }

  data.frame(edge_id = seq_len(n), from = from, to = to, n_px = n_px,
             n_orth = n_orth, n_diag = n_diag,
             length = len_poly, length_kimura = len_kimura,
             mean_diameter = d_mean, median_diameter = d_median,
             min_diameter = d_min, stringsAsFactors = FALSE)
}


#' Topological tip order by leaf-peeling
#'
#' Assigns \code{tip_order} per segment: terminal segments are 1; peeling
#' terminals away round by round gives \code{1 + max(child orders)} inward.
#' Cyclic components without endpoints receive \code{NA}.
#'
#' @param edge_tbl Edge table from \code{\link{build_edge_table}}.
#' @return \code{edge_tbl} with an added integer \code{tip_order} column.
#' @keywords internal
#' @noRd
#' @name compute_tip_order
compute_tip_order <- function(edge_tbl) {
  if (is.null(edge_tbl) || nrow(edge_tbl) == 0L) return(edge_tbl)
  nodes <- unique(c(edge_tbl$from, edge_tbl$to))
  ni <- stats::setNames(seq_along(nodes), nodes)
  fi <- as.integer(ni[edge_tbl$from]); ti <- as.integer(ni[edge_tbl$to])
  M <- length(nodes); E <- nrow(edge_tbl)
  active <- rep(TRUE, E); ord <- rep(NA_integer_, E); peel <- 0L
  repeat {
    inc  <- tabulate(c(fi[active], ti[active]), nbins = M)
    leaf <- active & (inc[fi] == 1L | inc[ti] == 1L)
    if (!any(leaf)) break
    peel <- peel + 1L; ord[leaf] <- peel; active[leaf] <- FALSE
  }
  edge_tbl$tip_order <- ord
  if (any(is.na(ord)))
    warning(sprintf("%d segment(s) NA tip_order (genuine cycle, no endpoints).", sum(is.na(ord))))
  edge_tbl
}


#' Strahler order, rooted at the base of each root system
#'
#' Assigns \code{strahler_order} per segment: every terminal segment is 1, and a
#' segment's order rises only where two of its children share the highest order
#' -- otherwise it inherits that highest order unchanged. This is the convention
#' the fine-root literature means by root order (Pregitzer et al. 2002; Fitter's
#' topological ordering), where first-order roots are the distal, unbranched
#' absorptive roots.
#'
#' It shares the peeling loop with \code{compute_tip_order()} but not its rule:
#' there, order climbs at every junction, so an axis carrying ten laterals
#' reaches order 11; here it stays 1 until two equal-order roots actually meet.
#'
#' \strong{Rooting.} Strahler ordering is defined on a \emph{rooted} tree, and a
#' skeleton is not rooted: the base of a root system is a free end like any tip,
#' so peeling from every free end at once scores the basal axis 1, alongside the
#' finest laterals. Each connected component is therefore rooted at its thickest
#' free end -- the degree-1 node whose segment has the largest mean diameter,
#' which on a root scan is the proximal end of the main axis. That node is held
#' back from the peeling, so its segment is scored last, from its children, as a
#' rooted traversal would. Components with no free end at all (a pure cycle)
#' receive \code{NA}, as for \code{tip_order}.
#'
#' @param edge_tbl Edge table from \code{\link{build_edge_table}}. Uses
#'   \code{mean_diameter} to root each component, falling back to \code{length}
#'   and then to the first free end found.
#' @return \code{edge_tbl} with an added integer \code{strahler_order} column.
#' @keywords internal
#' @noRd
compute_strahler_order <- function(edge_tbl) {
  if (is.null(edge_tbl) || nrow(edge_tbl) == 0L) return(edge_tbl)
  nodes <- unique(c(edge_tbl$from, edge_tbl$to))
  ni <- stats::setNames(seq_along(nodes), nodes)
  fi <- as.integer(ni[edge_tbl$from]); ti <- as.integer(ni[edge_tbl$to])
  M <- length(nodes); E <- nrow(edge_tbl)

  # ---- one root node per connected component --------------------------------
  ufr  <- seq_len(M)
  find <- function(x) { r <- x; while (ufr[r] != r) r <- ufr[r]
    while (ufr[x] != r) { nx <- ufr[x]; ufr[x] <<- r; x <- nx }; r }
  for (e in seq_len(E)) {
    a <- find(fi[e]); b <- find(ti[e])
    if (a != b) ufr[a] <- b
  }
  comp <- vapply(seq_len(M), find, integer(1))

  deg0 <- tabulate(c(fi, ti), nbins = M)
  w <- if (!is.null(edge_tbl$mean_diameter)) edge_tbl$mean_diameter
       else if (!is.null(edge_tbl$length))   edge_tbl$length
       else rep(1, E)
  w[!is.finite(w)] <- 0

  is_root  <- rep(FALSE, M)
  end_node <- c(fi, ti)                       # every segment end, both columns
  end_w    <- c(w,  w)
  free     <- deg0[end_node] == 1L
  if (any(free)) {
    en <- end_node[free]; o <- order(end_w[free], decreasing = TRUE)
    en <- en[o]
    is_root[en[!duplicated(comp[en])]] <- TRUE   # thickest free end per component
  }

  # ---- peel toward the root -------------------------------------------------
  active <- rep(TRUE, E)
  ord    <- rep(NA_integer_, E)
  # Orders of the segments already peeled away at each node. A segment's
  # children are whatever was removed at the end it became a leaf on; a genuine
  # tip has none, which is what makes it order 1.
  kids <- vector("list", M)

  strahler <- function(k) {
    if (length(k) == 0L) return(1L)
    m <- max(k)
    if (sum(k == m) >= 2L) m + 1L else m
  }

  repeat {
    inc <- tabulate(c(fi[active], ti[active]), nbins = M)
    # A root node is never a leaf, so its segment waits until everything above
    # it has been scored.
    leaf_at_f <- inc[fi] == 1L & !is_root[fi]
    leaf_at_t <- inc[ti] == 1L & !is_root[ti]
    leaf <- which(active & (leaf_at_f | leaf_at_t))
    if (length(leaf) == 0L) break

    # Score first, then remove and push upward, so segments peeled in the same
    # round cannot see each other's orders.
    leaf_end <- ifelse(leaf_at_f[leaf], fi[leaf], ti[leaf])
    far_end  <- ifelse(leaf_end == fi[leaf], ti[leaf], fi[leaf])
    new_ord  <- vapply(leaf_end, function(u) strahler(kids[[u]]), integer(1))

    ord[leaf]    <- new_ord
    active[leaf] <- FALSE
    for (k in seq_along(leaf))
      kids[[far_end[k]]] <- c(kids[[far_end[k]]], new_ord[k])
  }

  edge_tbl$strahler_order <- ord
  edge_tbl
}


#' Write an order-colored validation overlay (PNG)
#'
#' Paints each segment by its order value over the mask; unordered (cyclic)
#' segments are drawn in red so failures are visible.
#'
#' @param segs Segment list (kept via \code{keep_segments = TRUE}).
#' @param tip_order Numeric vector of the order value per segment to color by.
#' @param dims Integer length-2 vector, the (cropped) raster dimensions.
#' @param file Output PNG path.
#' @param mask Optional 0/1 mask drawn as a gray background.
#' @param max_side Long-side cap (px) for the output image.
#' @return Invisibly, \code{file}.
#' @export
render_order_overlay <- function(segs, tip_order, dims, file,
                                 mask = NULL, max_side = 2000) {
  nr <- dims[1]; nc <- dims[2]
  ord_img <- matrix(NA_real_, nr, nc); na_img <- matrix(NA_real_, nr, nc)
  for (i in seq_along(segs)) {
    o <- tip_order[i]; p <- segs[[i]]$coords
    if (is.na(o)) na_img[p] <- 1 else ord_img[p] <- o
  }
  maxo <- suppressWarnings(max(tip_order, na.rm = TRUE)); if (!is.finite(maxo)) maxo <- 1L
  ord_f <- .flip_matrix(ord_img); na_f <- .flip_matrix(na_img)
  mk_f  <- if (!is.null(mask)) .flip_matrix((.to_binary_matrix(mask) > 0) * 1) else .flip_matrix(matrix(0, nr, nc))
  xs <- seq(0, 1, length.out = nrow(ord_f)); ys <- seq(0, 1, length.out = ncol(ord_f))
  
  scale <- min(1, max_side / max(nr, nc))
  W <- max(1L, round(nc * scale)); H <- max(1L, round(nr * scale))
  grDevices::png(file, width = W, height = H, units = "px", res = 96)
  on.exit(grDevices::dev.off())
  graphics::par(mar = c(0, 0, 0, 0))
  
  graphics::image(xs, ys, mk_f, col = c("white", "gray88"), zlim = c(0, 1),
                  axes = FALSE, useRaster = TRUE)
  pal <- grDevices::hcl.colors(max(maxo, 1L), "viridis")
  # image() derives breaks from the data, so an all-NA layer must be skipped
  if (any(!is.na(ord_f)))
    graphics::image(xs, ys, ord_f, col = pal, breaks = seq(0.5, maxo + 0.5, by = 1),
                    add = TRUE, useRaster = TRUE)
  if (any(!is.na(na_f)))
    graphics::image(xs, ys, na_f, col = "red", zlim = c(0, 1),
                    add = TRUE, useRaster = TRUE)
  
  for (k in seq_len(maxo)) {
    yt <- 0.99 - (k - 1) * 0.035
    graphics::rect(0.01, yt - 0.03, 0.035, yt, col = pal[k], border = NA)
    graphics::text(0.04, yt - 0.015, sprintf("order %d", k), adj = 0, cex = 0.7)
  }
  if (any(is.na(tip_order))) {
    yt <- 0.99 - maxo * 0.035
    graphics::rect(0.01, yt - 0.03, 0.035, yt, col = "red", border = NA)
    graphics::text(0.04, yt - 0.015, "unordered", adj = 0, cex = 0.7)
  }
  invisible(file)
}


#' Skeleton-to-graph root ordering pipeline (pixel units)
#'
#' Core engine: crops to the foreground, computes the distance transform, traces
#' segments with junction contraction, resolves crossings, optionally prunes weak
#' tips, and assigns all three order schemes. Returns a per-segment edge table in
#' \emph{pixels}; use \code{\link{branch_order_map}} for a unit-aware wrapper.
#'
#' @details
#' Three order columns are produced (see \code{\link{branch_order_map}} for the
#' full rules): \code{tip_order} (per-segment leaf-peeling), \code{root_order}
#' (per-root max tip_order), \code{branch_order} (per-root centrifugal generation
#' from the thickest root). \code{color_by} selects which one the overlay uses.
#'
#' @param skel Binary skeleton: single-layer \code{SpatRaster} (preferred) or 0/1
#'   matrix. Do not pre-convert a raster with \code{as.matrix()}. If
#'   \code{NULL}, it is computed from \code{mask} via
#'   \code{skeletonize_image()}.
#' @param mask Filled root mask on the same grid for the distance transform.
#'   Required if \code{skel} is \code{NULL}.
#' @param verbose Print progress.
#' @param dt_backend Distance-transform backend: \code{"auto"}, \code{"imager"}, or \code{"baseR"}.
#' @param crop Crop to the foreground bounding box (with a 1-px pad) first.
#' @param overlay_png Optional path for the validation PNG.
#' @param max_side Long-side cap (px) for the overlay.
#' @param keep_segments Attach the traced segments as \code{attr(., "segments")}
#'   (needed for re-plotting and classification maps).
#' @param resolve_overlaps Resolve degree-4 crossings by continuity.
#' @param splice_passthrough Dissolve contracted junctions that carry only two
#'   segment ends (bends and thinning artifacts, not branch points) by splicing
#'   the two arms into one segment. Leave \code{TRUE} unless you specifically
#'   want one edge per traced skeleton chain.
#' @param crossing_straight Straightness threshold for crossing resolution.
#' @param crossing_diam_ratio Opt-in thickness test at a degree-4 node, off by
#'   default (\code{0} = resolve on geometry alone).
#'   When set (0.5 is a reasonable value), a node is treated as a crossing only
#'   if the two candidate through-roots are within this thickness ratio of each
#'   other; otherwise it is read as a bilateral branch and left intact. This is
#'   a genuine trade-off, not a free win: a thick axis crossed by a thin root and
#'   a thick axis with two thin laterals are the same shape in outline, so
#'   raising the threshold fixes the second case and breaks the first. Set it
#'   only if your material has more bilateral branching than fine-over-coarse
#'   crossing.
#' @param color_by Which order column the overlay colors by.
#' @param diam_weight Diameter-vs-angle weight for the root-continuation choice.
#' @param prune_min_length,prune_min_diameter,prune_iter Optional terminal-segment
#'   pruning (off when \code{prune_iter = 0}); see \code{\link{prune_terminal_segments}}.
#' @return A per-segment edge table (data.frame) with order and diameter columns
#'   in pixels, carrying \code{attr}s \code{crop_offset}, \code{dims}, and
#'   (if requested) \code{segments}.
#' @seealso \code{\link{branch_order_map}}
#' @export
root_graph_pipeline <- function(skel = NULL, mask = NULL, verbose = TRUE,
                                dt_backend = "auto", crop = TRUE,
                                overlay_png = NULL, max_side = 2000,
                                keep_segments = FALSE,
                                resolve_overlaps = TRUE, splice_passthrough = TRUE,
                                crossing_straight = -0.5, crossing_diam_ratio = 0,
                                color_by = c("strahler_order", "branch_order", "root_order", "tip_order"),
                                diam_weight = 0.5,
                                prune_min_length = 0, prune_min_diameter = 0, prune_iter = 0L) {
  color_by <- match.arg(color_by)
  t0 <- Sys.time()
  if (is.null(skel)) {
    if (is.null(mask)) stop("Either 'skel' or 'mask' must be supplied.")
    if (verbose) cat("No 'skel' supplied; computing skeleton via skeletonize_image()...\n")
    skel <- skeletonize_image(mask, verbose = FALSE)
  }
  skel <- .to_binary_matrix(skel)
  if (is.null(mask)) {
    warning("No 'mask'; using skeleton for the DT (diameters ~1 px).")
    mask <- skel
  } else {
    mask <- .to_binary_matrix(mask)
    if (!all(dim(mask) == dim(skel))) stop("'skel' and 'mask' dims must match.")
  }
  
  ro <- 0L; co <- 0L
  if (crop) {
    fg <- (skel > 0) | (mask > 0)
    if (!any(fg)) { warning("Empty image."); return(NULL) }
    rs <- which(rowSums(fg) > 0); cs <- which(colSums(fg) > 0)
    r1 <- min(rs); r2 <- max(rs); c1 <- min(cs); c2 <- max(cs)
    # 1-px background pad so foreground at the bbox edge keeps a nearby background
    r1 <- max(1L, r1 - 1L); r2 <- min(nrow(skel), r2 + 1L)
    c1 <- max(1L, c1 - 1L); c2 <- min(ncol(skel), c2 + 1L)
    skel <- skel[r1:r2, c1:c2, drop = FALSE]; mask <- mask[r1:r2, c1:c2, drop = FALSE]
    ro <- r1 - 1L; co <- c1 - 1L
    if (verbose) cat(sprintf("Cropped to %dx%d (offset %d,%d)\n", nrow(skel), ncol(skel), ro, co))
  }
  dims <- c(nrow(skel), ncol(skel))
  
  if (verbose) cat("Distance transform...\n")
  DT <- .distance_transform_edt(mask, backend = dt_backend)
  if (verbose) cat("Tracing + contracting junctions...\n")
  segs <- trace_segments(skel)
  if (verbose) cat(sprintf("  %d segment(s)\n", length(segs)))
  .check_skeleton_coverage(segs, skel)
  
  if (splice_passthrough) {
    n0 <- length(segs)
    segs <- .splice_passthrough(segs)
    if (verbose) cat(sprintf("  spliced pass-through nodes: %d -> %d segment(s)\n",
                             n0, length(segs)))
  }
  
  if (resolve_overlaps) {
    n0 <- length(segs)
    segs <- resolve_crossings(segs, DT = DT, straight_dot = crossing_straight,
                              diam_ratio = crossing_diam_ratio, verbose = verbose)
    if (verbose) cat(sprintf("  resolved crossings: %d -> %d segment(s)\n", n0, length(segs)))
  }
  
  if (prune_iter > 0L && (prune_min_length > 0 || prune_min_diameter > 0)) {
    n0 <- length(segs)
    segs <- prune_terminal_segments(segs, DT, min_length = prune_min_length,
                                    min_diameter = prune_min_diameter, iter = prune_iter)
    if (verbose) cat(sprintf("  pruned weak tips: %d -> %d segment(s)\n", n0, length(segs)))
  }
  
  edge_tbl <- build_edge_table(segs, DT)
  if (is.null(edge_tbl)) { warning("No segments."); return(edge_tbl) }
  if (verbose) cat("Tip order...\n")
  edge_tbl <- compute_tip_order(edge_tbl)
  edge_tbl <- compute_strahler_order(edge_tbl)
  edge_tbl <- assign_root_order(segs, edge_tbl, diam_weight = diam_weight)
  
  attr(edge_tbl, "crop_offset") <- c(row = ro, col = co)
  attr(edge_tbl, "dims") <- dims
  if (keep_segments) attr(edge_tbl, "segments") <- segs
  
  if (!is.null(overlay_png)) {
    if (verbose) cat(sprintf("Overlay (%s) -> %s\n", color_by, overlay_png))
    render_order_overlay(segs, edge_tbl[[color_by]], dims, overlay_png,
                         mask = mask, max_side = max_side)
  }
  if (verbose) {
    ok <- !is.na(edge_tbl$tip_order)
    cat(sprintf("Done in %.2fs. %d segments, max order %d, total length %.1f px\n",
                as.numeric(difftime(Sys.time(), t0, units = "secs")), nrow(edge_tbl),
                if (any(ok)) max(edge_tbl$tip_order[ok]) else 0L, sum(edge_tbl$length)))
  }
  edge_tbl
}


# ---- example ---------------------------------------------------------------
# skel_rast <- skeletonize_image(img)
# mask_rast <- load_flexible_image(img, output_format = "spatrast", scale = "binary")
# edge_tbl  <- root_graph_pipeline(skel_rast, mask_rast,
#                                  dt_backend  = "imager",
#                                  overlay_png = "order_validation.png")
# # Re-plot later (needs keep_segments = TRUE on the call above):
# # render_order_overlay(attr(edge_tbl,"segments"), edge_tbl$tip_order,
# #                      attr(edge_tbl,"dims"), "order_validation.png", mask = mask_rast)


#' Native-resolution validation of a sub-window
#'
#' Renders a small window at \code{scale}x magnification (no downsampling) so the
#' order-colored graph can be checked against the skeleton. Skeleton pixels are
#' gray; traced/ordered pixels are colored on top, so bare gray reveals skeleton
#' the graph missed. Coordinates are in the original image frame.
#'
#' @param et An \code{edges} table carrying \code{attr(., "segments")} and
#'   \code{attr(., "crop_offset")} (run with \code{keep_segments = TRUE}).
#' @param skel The original skeleton (\code{SpatRaster} or matrix) for the gray background.
#' @param r_range,c_range Integer length-2 row/column ranges (original coordinates).
#' @param scale Magnification factor (device px per image px).
#' @param file Output PNG path.
#' @param order_col Which order column to color by (default \code{"root_order"}).
#' @return Invisibly, \code{file}.
#' @examples
#' \dontrun{
#' res <- branch_order_map(skel, mask, order = "root_order", unit = "px")
#' plot_order_window(res$edges, skel, r_range = c(1, 500), c_range = c(1, 500),
#'                   file = tempfile(fileext = ".png"))
#' }
#' @export
plot_order_window <- function(et, skel, r_range, c_range, scale = 3, file = "window.png",
                              order_col = "root_order") {
  segs <- attr(et, "segments")
  if (is.null(segs)) stop("Re-run root_graph_pipeline(..., keep_segments = TRUE).")
  if (!order_col %in% names(et)) order_col <- "tip_order"
  ord_vec <- et[[order_col]]
  off <- attr(et, "crop_offset"); ro <- as.integer(off["row"]); co <- as.integer(off["col"])
  smat <- .to_binary_matrix(skel)
  r0 <- r_range[1]; r1 <- r_range[2]; c0 <- c_range[1]; c1 <- c_range[2]
  hh <- r1 - r0 + 1L; ww <- c1 - c0 + 1L
  skel_w <- smat[r0:r1, c0:c1]
  ord <- matrix(NA_real_, hh, ww)
  for (i in seq_along(segs)) {
    p <- segs[[i]]$coords
    R <- p[, 1] + ro; C <- p[, 2] + co                 # cropped frame -> original frame
    keep <- R >= r0 & R <= r1 & C >= c0 & C <= c1
    if (!any(keep)) next
    ord[cbind(R[keep] - r0 + 1L, C[keep] - c0 + 1L)] <- ord_vec[i]
  }
  maxo <- suppressWarnings(max(ord_vec, na.rm = TRUE)); if (!is.finite(maxo)) maxo <- 1L
  pal <- grDevices::hcl.colors(max(maxo, 1L), "viridis")
  grDevices::png(file, width = ww * scale, height = hh * scale); on.exit(grDevices::dev.off())
  graphics::par(mar = c(0, 0, 0, 0))
  skel_f <- .flip_matrix(skel_w); ord_f <- .flip_matrix(ord)
  xs <- seq(0, 1, length.out = nrow(skel_f)); ys <- seq(0, 1, length.out = ncol(skel_f))
  graphics::image(xs, ys, skel_f, col = c("white", "gray75"), axes = FALSE, useRaster = TRUE)
  graphics::image(xs, ys, ord_f, col = pal, breaks = seq(0.5, maxo + 0.5, by = 1),
                  add = TRUE, useRaster = TRUE)
  invisible(file)
}


#' Unit tangent of a segment at one endpoint
#'
#' Direction pointing away from the node into the segment, averaged over the
#' first/last \code{look} pixels.
#' @keywords internal
#' @noRd
.endpoint_tangent <- function(coords, side, look = 5L) {
  n <- nrow(coords); k <- min(look, n)
  v <- if (side == 1L) coords[k, ] - coords[1, ] else coords[n - k + 1L, ] - coords[n, ]
  nrm <- sqrt(sum(v^2)); if (nrm == 0) return(c(0, 0)); v / nrm
}

#' Resolve root crossings by continuity
#'
#' At clean degree-4 nodes (an X overlap of two roots in projection) the two
#' arms whose tangents are most collinear are spliced into one through-root,
#' dissolving the crossing-induced cycle. Genuine branches (degree 3) and
#' ambiguous nodes are left untouched.
#'
#' @param segs Segment list from \code{trace_segments}.
#' @param DT Distance-transform matrix. Optional but recommended: it enables the
#'   thickness test that separates a crossing from a bilateral branch (see
#'   \code{diam_ratio}). Without it, geometry alone decides.
#' @param straight_dot A crossing is resolved only if both candidate pairs have
#'   tangent dot product below this (more negative = straighter; -0.5 ~ >120deg).
#' @param look Pixels used to estimate endpoint tangents and local radii.
#' @param diam_ratio Minimum thickness ratio between the two through-roots, or
#'   \code{0} (default) to disable the test. Two thin laterals leaving a thick
#'   axis in opposite directions look exactly like an X in outline, and so does
#'   a thin root crossing a thick one; no threshold gets both right, which is
#'   why this is off unless the caller asks for it.
#' @param verbose Report nodes with more than four arms, which this function
#'   cannot resolve and \code{assign_root_order()} will read as one continuation
#'   plus laterals.
#' @return A segment list with crossings spliced; \code{attr(., "dims")} preserved.
#' @keywords internal
#' @noRd
resolve_crossings <- function(segs, DT = NULL, straight_dot = -0.5, look = 5L,
                              diam_ratio = 0, verbose = FALSE) {
  ns <- length(segs)
  if (ns < 2L) return(segs)
  
  tang <- lapply(seq_len(ns), function(i) {
    p <- segs[[i]]$coords
    list(.endpoint_tangent(p, 1L, look), .endpoint_tangent(p, 2L, look))
  })
  side_tan <- function(i, s) tang[[i]][[s]]
  # Representative radius per segment, for the crossing-vs-bilateral-branch test.
  # Measured over the whole segment, not at the node: every arm is swollen where
  # it enters a junction, so an end-window makes a thin lateral look thick.
  seg_rad <- if (is.null(DT) || diam_ratio <= 0) NULL else
    vapply(seq_len(ns), function(i) stats::median(DT[segs[[i]]$coords]), numeric(1))
  
  node_inc <- .node_incidence(segs)

  partner_seg  <- matrix(0L, ns, 2)
  partner_side <- matrix(0L, ns, 2)
  
  n_highdeg <- 0L
  for (lab in names(node_inc)) {
    inc <- node_inc[[lab]]
    if (nrow(inc) > 4L) n_highdeg <- n_highdeg + 1L
    if (nrow(inc) != 4L) next                       # only clean X crossings
    tg  <- lapply(seq_len(4), function(m) side_tan(inc[m, 1], inc[m, 2]))
    dot <- function(a, b) sum(tg[[a]] * tg[[b]])
    pairings <- list(c(1, 2, 3, 4), c(1, 3, 2, 4), c(1, 4, 2, 3))
    best <- NULL; best_tot <- Inf
    for (pr in pairings) {
      d <- c(dot(pr[1], pr[2]), dot(pr[3], pr[4]))
      if (sum(d) < best_tot) { best_tot <- sum(d); best <- list(pr = pr, d = d) }
    }
    if (max(best$d) >= straight_dot) next            # both pairs must be straight
    if (!is.null(seg_rad)) {
      # a crossing joins two roots of comparable thickness; a bilateral branch
      # pairs a thick axis against two thin laterals, and must be left alone
      pr0 <- best$pr
      rad <- seg_rad[inc[, 1]]
      pd  <- c(mean(rad[pr0[1:2]]), mean(rad[pr0[3:4]]))
      if (max(pd) > 0 && min(pd) / max(pd) < diam_ratio) next
    }
    pr <- best$pr
    for (q in c(1L, 3L)) {                            # link (pr[q], pr[q+1])
      a <- inc[pr[q], ]; b <- inc[pr[q + 1L], ]
      partner_seg[a[1], a[2]]  <- b[1]; partner_side[a[1], a[2]]  <- b[2]
      partner_seg[b[1], b[2]]  <- a[1]; partner_side[b[1], b[2]]  <- a[2]
    }
  }
  
  if (verbose && n_highdeg > 0L)
    cat(sprintf(paste0("  %d node(s) with >4 arms left unresolved; each is read as ",
                       "one continuation plus laterals, not as a crossing\n"), n_highdeg))

  if (!any(partner_seg > 0L)) return(segs)
  .stitch_chains(segs, partner_seg, partner_side)
}


#' Follow partner links and concatenate segments into single polylines
#'
#' Shared back end for \code{\link{resolve_crossings}} and
#' \code{\link{.splice_passthrough}}. Chains are always entered from a free end
#' (a side with no partner); only after every open chain is consumed are the
#' remaining closed rings walked. Starting mid-chain would emit the two halves
#' as separate segments and leave the node unresolved.
#'
#' @param segs Segment list.
#' @param partner_seg,partner_side \code{length(segs) x 2} integer matrices; the
#'   segment/side that each (segment, side) is spliced to, or 0 for a free end.
#' @return A segment list with chains concatenated; \code{"dims"} and
#'   \code{"node_xy"} attributes are preserved.
#' @keywords internal
#' @noRd
.stitch_chains <- function(segs, partner_seg, partner_side) {
  ns <- length(segs)
  side_lab <- function(i, s) if (s == 1L) segs[[i]]$from else segs[[i]]$to
  # open chains first, closed rings last
  order_starts <- c(which(partner_seg[, 1] == 0L | partner_seg[, 2] == 0L),
                    which(partner_seg[, 1] != 0L & partner_seg[, 2] != 0L))
  visited <- logical(ns); out <- list()
  for (s0 in order_starts) {
    if (visited[s0]) next
    start_side <- if (partner_seg[s0, 1] == 0L) 1L
    else if (partner_seg[s0, 2] == 0L) 2L else 1L     # last case: closed loop fallback
    coords_list <- list()
    cur <- s0; enter_side <- start_side
    start_lab <- side_lab(s0, start_side)
    last_seg <- s0; last_exit <- if (start_side == 1L) 2L else 1L
    guard <- 0L
    repeat {
      visited[cur] <- TRUE
      p <- segs[[cur]]$coords
      if (enter_side == 2L) p <- p[nrow(p):1, , drop = FALSE]
      coords_list[[length(coords_list) + 1L]] <- p
      exit_side <- if (enter_side == 1L) 2L else 1L
      last_seg <- cur; last_exit <- exit_side
      nxt <- partner_seg[cur, exit_side]
      if (nxt == 0L || visited[nxt]) break
      enter_side <- partner_side[cur, exit_side]
      cur <- nxt
      guard <- guard + 1L; if (guard > ns) break
    }
    full <- do.call(rbind, coords_list)
    out[[length(out) + 1L]] <- list(coords = full,
                                    from = start_lab,
                                    to   = side_lab(last_seg, last_exit))
  }
  attr(out, "dims")    <- attr(segs, "dims")
  attr(out, "node_xy") <- attr(segs, "node_xy")
  out
}


#' Splice segments across pass-through nodes
#'
#' A contracted junction cluster carrying only \emph{two} segment ends is not a
#' branch point -- it is a bend or a thinning artifact that fragmented one root
#' into two edges. Splicing the pair back together restores contiguous segments,
#' so \code{n_segments} and \code{mean_segment_length} describe roots rather than
#' skeleton staircases, and the cluster the two arms straddle is re-counted as
#' polyline length instead of being dropped.
#'
#' @param segs Segment list from \code{\link{trace_segments}}.
#' @return A segment list with pass-through nodes dissolved; attributes preserved.
#' @keywords internal
#' @noRd
.splice_passthrough <- function(segs) {
  ns <- length(segs)
  if (ns < 2L) return(segs)
  node_inc <- .node_incidence(segs)
  partner_seg  <- matrix(0L, ns, 2)
  partner_side <- matrix(0L, ns, 2)
  for (lab in names(node_inc)) {
    inc <- node_inc[[lab]]
    if (nrow(inc) != 2L) next
    if (inc[1, 1] == inc[2, 1]) next                  # self-loop: leave it alone
    a <- inc[1, ]; b <- inc[2, ]
    partner_seg[a[1], a[2]] <- b[1]; partner_side[a[1], a[2]] <- b[2]
    partner_seg[b[1], b[2]] <- a[1]; partner_side[b[1], b[2]] <- a[2]
  }
  if (!any(partner_seg > 0L)) return(segs)
  .stitch_chains(segs, partner_seg, partner_side)
}


#' Group segments into roots and assign root_order and branch_order
#'
#' Segments are grouped into continuous roots: at each junction the continuation
#' pair is the two arms maximizing \code{straightness + diam_weight * diameter
#' consistency}, where consistency combines how alike the two arms are with how
#' thick they are relative to the node; the odd arm out is a lateral. \code{root_order} is the maximum
#' \code{tip_order} along each root. \code{branch_order} is a centrifugal
#' generation: the thickest root (length-weighted mean diameter) in each
#' connected component is 1, counting branching hops outward.
#'
#' @param segs Segment list (post crossing-resolution).
#' @param edge_tbl Edge table with \code{tip_order} already assigned.
#' @param diam_weight Diameter-vs-angle weight for the continuation choice (>= 0).
#' @param look Pixels used to estimate endpoint tangents.
#' @return \code{edge_tbl} with added \code{root_order}, \code{branch_order},
#'   and \code{root_id} columns.
#' @keywords internal
#' @noRd
assign_root_order <- function(segs, edge_tbl, diam_weight = 0.5, look = 5L) {
  ns <- length(segs)
  if (ns == 0L) return(edge_tbl)
  diam <- edge_tbl$mean_diameter
  tang <- lapply(seq_len(ns), function(i) {
    p <- segs[[i]]$coords
    list(.endpoint_tangent(p, 1L, look), .endpoint_tangent(p, 2L, look))
  })
  
  node_inc <- .node_incidence(segs)

  ufr <- seq_len(ns)
  find <- function(x) { r <- x; while (ufr[r] != r) r <- ufr[r]
  while (ufr[x] != r) { nx <- ufr[x]; ufr[x] <<- r; x <- nx }; r }
  unite <- function(i, j) { ri <- find(i); rj <- find(j); if (ri != rj) ufr[ri] <<- rj }
  
  node_cont <- list(); node_lats <- list()
  for (lab in names(node_inc)) {
    inc <- node_inc[[lab]]; k <- nrow(inc)
    if (k < 2L) next
    if (k == 2L) {                                      # pure pass-through, no lateral
      unite(inc[1, 1], inc[2, 1])
      node_cont[[lab]] <- c(inc[1, 1], inc[2, 1]); node_lats[[lab]] <- integer(0)
      next
    }
    # k >= 3: pick the continuation pair (straightest + most similar diameter)
    # Diameter similarity alone would pick the wrong pair: two identical thin
    # laterals score a perfect 1 while two near-identical thick arms score
    # slightly less, so at a 4-way node the laterals outscore the parent axis.
    # The diameter term is therefore scaled by how thick the pair is relative
    # to the node (`dmag`), which breaks that tie in favour of the thick axis
    # among equally straight, equally consistent candidates.
    dmax <- max(diam[inc[, 1]], na.rm = TRUE)
    bi <- NA_integer_; bj <- NA_integer_; best <- -Inf
    for (m in 1:(k - 1L)) for (n in (m + 1L):k) {
      ti <- tang[[inc[m, 1]]][[inc[m, 2]]]; tj <- tang[[inc[n, 1]]][[inc[n, 2]]]
      di <- diam[inc[m, 1]]; dj <- diam[inc[n, 1]]
      straight <- -sum(ti * tj)
      dsim <- if (max(di, dj) > 0) min(di, dj) / max(di, dj) else 1
      dmag <- if (is.finite(dmax) && dmax > 0) min(di, dj) / dmax else 1
      sc <- straight + diam_weight * dsim * dmag
      if (sc > best) { best <- sc; bi <- m; bj <- n }
    }
    unite(inc[bi, 1], inc[bj, 1])
    node_cont[[lab]] <- c(inc[bi, 1], inc[bj, 1])
    node_lats[[lab]] <- inc[-c(bi, bj), 1]              # remaining arms = laterals
  }
  
  root <- vapply(seq_len(ns), find, integer(1))
  roots_u <- sort(unique(root))
  rid <- match(root, roots_u)                          # 1..n_roots per segment
  n_roots <- length(roots_u)
  
  # (2) root_order: max tip_order along each root (thick main axis stays high)
  grp_max <- tapply(edge_tbl$tip_order, rid,
                    function(v) if (all(is.na(v))) NA_integer_ else max(v, na.rm = TRUE))
  edge_tbl$root_order <- as.integer(grp_max[as.character(rid)])
  
  # (3) branch_order: generation, THICKEST root in each connected component = 1,
  # counting up by branching hops outward. parent = continuation root at a
  # junction; child = each lateral root. Component connectivity is undirected.
  ep <- integer(0); ec <- integer(0)
  for (lab in names(node_cont)) {
    cont <- node_cont[[lab]]; lats <- node_lats[[lab]]
    if (length(cont) == 0L || length(lats) == 0L) next
    Rp <- rid[cont[1]]
    for (L in lats) { ep <- c(ep, Rp); ec <- c(ec, rid[L]) }
  }
  
  # per-root diameter (length-weighted) to choose the thickest = seed
  root_dia <- numeric(n_roots)
  for (g in seq_len(n_roots)) {
    ix <- which(rid == g); w <- edge_tbl$length[ix]; d <- edge_tbl$mean_diameter[ix]
    root_dia[g] <- if (sum(w) > 0) sum(d * w) / sum(w) else mean(d)
  }
  
  # undirected root adjacency
  adj <- vector("list", n_roots)
  for (e in seq_along(ep)) {
    adj[[ep[e]]] <- c(adj[[ep[e]]], ec[e]); adj[[ec[e]]] <- c(adj[[ec[e]]], ep[e])
  }
  
  gen <- rep(NA_integer_, n_roots); seen <- logical(n_roots)
  for (start in seq_len(n_roots)) {
    if (seen[start]) next
    comp <- integer(0); q <- start; seen[start] <- TRUE; qi <- 1L   # gather component
    while (qi <= length(q)) {
      v <- q[qi]; qi <- qi + 1L; comp <- c(comp, v)
      for (w in adj[[v]]) if (!seen[w]) { seen[w] <- TRUE; q <- c(q, w) }
    }
    seed <- comp[which.max(root_dia[comp])]                          # thickest root
    d <- rep(NA_integer_, n_roots); d[seed] <- 0L; q2 <- seed; qj <- 1L
    while (qj <= length(q2)) {                                       # hops from seed
      v <- q2[qj]; qj <- qj + 1L
      for (w in adj[[v]]) if (is.na(d[w])) { d[w] <- d[v] + 1L; q2 <- c(q2, w) }
    }
    gen[comp] <- d[comp] + 1L
  }
  edge_tbl$branch_order <- as.integer(gen[rid])
  
  # (4) architecture counts: tips (degree-1 endpoints) and branch points.
  # A branch point is a branching junction on the parent (continuation) root;
  # it is attributed to the parent segment so summing by order aggregates
  # correctly. For an ordinary 3-way fork this equals the number of laterals.
  node_deg <- vapply(node_inc, nrow, integer(1))            # incidence per node label
  edge_tbl$n_tips <- as.integer((node_deg[edge_tbl$from] == 1L) +
                                  (node_deg[edge_tbl$to]   == 1L))
  n_bp <- integer(ns)
  for (lab in names(node_cont)) {
    lats <- node_lats[[lab]]; cont <- node_cont[[lab]]
    if (length(lats) > 0L && length(cont) > 0L) n_bp[cont[1]] <- n_bp[cont[1]] + 1L
  }
  edge_tbl$n_branch_points <- as.integer(n_bp)
  
  edge_tbl$root_id <- rid
  edge_tbl
}


# ============================================================
# PACKAGE-FACING WRAPPERS
# ============================================================

#' Prune short or thin terminal segments
#'
#' Iteratively removes terminal (degree-1) segments below a length or diameter
#' threshold. Segments are only deleted, never re-routed, so the ordering stays
#' valid. Skeleton-level pruning before the pipeline is an alternative, fully
#' modular approach.
#'
#' @details
#' Between passes the junctions left behind are re-checked: once a lateral is
#' gone its branch point carries only two arms, so the parent's two pieces are
#' spliced back into one segment before the next pass measures them. Without
#' that step a multi-pass prune eats healthy axes, because each parent piece
#' between two former branch points is short enough on its own to look like a
#' spur. Set \code{splice = FALSE} to delete laterals without splicing.
#'
#' @param segs Segment list from \code{trace_segments}.
#' @param DT Distance-transform matrix (for the diameter test).
#' @param min_length Minimum segment length (px) to keep a terminal segment.
#' @param min_diameter Minimum segment diameter (px) to keep a terminal segment.
#' @param iter Number of pruning passes.
#' @param splice Re-join the parent across a junction that a removed lateral has
#'   left with only two arms, before the next pass.
#' @return The pruned segment list.
#' @export
prune_terminal_segments <- function(segs, DT, min_length = 0, min_diameter = 0,
                                    iter = 1L, splice = TRUE) {
  for (it in seq_len(iter)) {
    if (length(segs) == 0L) break
    if (splice && it > 1L) segs <- .splice_passthrough(segs)
    labs <- unlist(lapply(segs, function(s) c(s$from, s$to)))
    deg <- table(labs)
    node_xy <- attr(segs, "node_xy")
    drop <- logical(length(segs))
    for (i in seq_along(segs)) {
      s <- segs[[i]]
      if (deg[s$from] != 1L && deg[s$to] != 1L) next        # only terminal segments
      # Measured the same way build_edge_table() reports it, junction stub
      # included -- otherwise a segment can be pruned at a min_length it clears
      # in the table the user reads.
      L  <- .segment_length(s, node_xy)
      dm <- 2 * min(DT[s$coords])
      if (L < min_length || dm < min_diameter) drop[i] <- TRUE
    }
    if (!any(drop)) break
    dims <- attr(segs, "dims"); node_xy <- attr(segs, "node_xy")
    segs <- segs[!drop]                                  # `[` drops attributes
    attr(segs, "dims") <- dims; attr(segs, "node_xy") <- node_xy
  }
  segs
}

#' Prune short or thin spurs from a skeleton or segmentation image
#'
#' Standalone clean-up step that removes short/thin terminal branches ("spurs")
#' that survive segmentation. Spurs are defined on the skeleton graph (that is
#' where connectivity lives), but the result can be returned either as the
#' cleaned skeleton (\code{output = "skeleton"}) or propagated back into
#' segmentation space (\code{output = "mask"}).
#'
#' @details
#' The skeleton is traced into segments and \code{\link{prune_terminal_segments}}
#' removes terminal segments below \code{min_length} / \code{min_diameter},
#' iterated \code{iter} times. For \code{output = "mask"}, each surviving
#' skeleton pixel is regrown by its local distance-transform radius (the medial
#' axis reconstruction) and intersected with the original \code{mask}, so the
#' spur's body is removed from the segmentation while the real roots are
#' restored to their original thickness. A \code{mask} is therefore strongly
#' recommended (and required for trustworthy \code{output = "mask"} and the
#' \code{min_diameter} test); without it diameters collapse to ~1 px.
#'
#' @param skel Binary skeleton: single-layer \code{SpatRaster} or 0/1 matrix.
#' @param mask Filled root mask on the same grid, for the distance transform and
#'   (when \code{output = "mask"}) reconstruction. Optional but recommended.
#' @param min_length Minimum terminal-segment length (px) to keep.
#' @param min_diameter Minimum terminal-segment diameter (px) to keep.
#' @param iter Number of pruning passes.
#' @param output \code{"skeleton"} (default) returns the pruned skeleton;
#'   \code{"mask"} returns the pruned binary segmentation.
#' @param verbose Print progress.
#' @return A cleaned image matching the class of \code{skel} (\code{SpatRaster}
#'   in, \code{SpatRaster} out; matrix in, matrix out).
#' @seealso \code{\link{prune_terminal_segments}}, \code{\link{clean_image}}
#' @export
#' @examples
#' \dontrun{
#' skl <- skeletonize_image(mask)
#' # remove spurs shorter than 15 px, return a cleaned skeleton
#' skl2 <- prune_skeleton(skl, mask, min_length = 15, iter = 2)
#' # ... or clean the segmentation itself
#' seg2 <- prune_skeleton(skl, mask, min_length = 15, output = "mask")
#' }
prune_skeleton <- function(skel, mask = NULL,
                           min_length = 0, min_diameter = 0, iter = 1L,
                           output = c("skeleton", "mask"), verbose = FALSE) {
  output <- match.arg(output)
  is_rast  <- inherits(skel, "SpatRaster")
  template <- if (is_rast) skel[[1]] else NULL

  skel_m <- .to_binary_matrix(skel)
  if (is.null(mask)) {
    if (output == "mask")
      warning("No 'mask' supplied; reconstructed mask will be ~1 px wide. Pass 'mask' for seg-space pruning.")
    mask_m <- skel_m
  } else {
    mask_m <- .to_binary_matrix(mask)
    if (!all(dim(mask_m) == dim(skel_m))) stop("'skel' and 'mask' dims must match.")
  }

  DT   <- .distance_transform_edt(mask_m)
  segs <- trace_segments(skel_m)
  n0   <- length(segs)
  segs <- prune_terminal_segments(segs, DT, min_length = min_length,
                                  min_diameter = min_diameter, iter = iter)
  if (verbose) cat(sprintf("prune_skeleton: %d -> %d segment(s)\n", n0, length(segs)))

  nr <- nrow(skel_m); nc <- ncol(skel_m)
  kept <- matrix(0L, nr, nc)
  for (s in segs) kept[s$coords] <- 1L          # coords are (row, col)

  if (output == "skeleton") {
    out_m <- kept
  } else {
    # Medial-axis reconstruction: stamp a DT-radius disk at each kept pixel,
    # then clip to the original mask so we never grow beyond the real root.
    out_m <- matrix(0L, nr, nc)
    idx   <- which(kept > 0L, arr.ind = TRUE)
    if (nrow(idx) > 0L) {
      rad <- DT[kept > 0L]
      disk_cache <- list()
      for (k in seq_len(nrow(idx))) {
        rr <- max(0L, as.integer(floor(rad[k])))
        key <- as.character(rr)
        off <- disk_cache[[key]]
        if (is.null(off)) {                      # offsets for a disk of radius rr
          g <- expand.grid(dr = -rr:rr, dc = -rr:rr)
          off <- as.matrix(g[g$dr^2 + g$dc^2 <= rr^2, , drop = FALSE])
          disk_cache[[key]] <- off
        }
        R <- idx[k, 1] + off[, 1]; C <- idx[k, 2] + off[, 2]
        ok <- R >= 1L & R <= nr & C >= 1L & C <= nc
        out_m[cbind(R[ok], C[ok])] <- 1L
      }
    }
    out_m <- out_m * mask_m                      # clip to original segmentation
  }

  if (is_rast) {
    out <- terra::rast(template)
    terra::values(out) <- as.vector(t(out_m))
    names(out) <- if (output == "mask") "pruned_mask" else "pruned_skeleton"
    out
  } else out_m
}

#' Rasterise a per-segment value onto the image grid
#'
#' Writes any per-segment column (default the order class) back onto the full
#' image grid, aligned to \code{template}, for masking and zonal statistics.
#' Background pixels are \code{NA}.
#'
#' @details
#' Only pixels that belong to a segment are painted. Junction contraction
#' dissolves the interior of each branch-point cluster, so the map carries one
#' unpainted pixel per branch point and a pixel count taken off the map runs
#' that much below the skeleton. Read lengths from \code{et$length}, which
#' accounts for the contracted pixels; the map is for masking and zonal
#' statistics, not for measuring.
#'
#' @param et An \code{edges} table carrying \code{attr(., "segments")} and
#'   \code{attr(., "crop_offset")} (run with \code{keep_segments = TRUE}).
#' @param template \code{SpatRaster} (or matrix) defining the output grid/extent.
#' @param value Column of \code{et} to rasterise (e.g. \code{"branch_order"},
#'   \code{"mean_diameter"}).
#' @return A \code{SpatRaster} (or matrix) of \code{value} per root pixel, \code{NA}
#'   elsewhere, aligned to \code{template}.
#' @examples
#' \dontrun{
#' res <- branch_order_map(skel, mask, order = "branch_order", unit = "px")
#' cmap <- order_classification_map(res$edges, template = skel,
#'                                  value = "branch_order")
#' }
#' @export
order_classification_map <- function(et, template, value = "branch_order") {
  segs <- attr(et, "segments")
  if (is.null(segs)) stop("Run with keep_segments = TRUE (branch_order_map does this).")
  if (!value %in% names(et)) stop(sprintf("'%s' is not a column of et.", value))
  off <- attr(et, "crop_offset"); ro <- as.integer(off["row"]); co <- as.integer(off["col"])
  is_rast <- inherits(template, "SpatRaster")
  nr <- if (is_rast) terra::nrow(template) else nrow(template)
  nc <- if (is_rast) terra::ncol(template) else ncol(template)
  m <- matrix(NA_real_, nr, nc)
  vv <- et[[value]]
  for (i in seq_along(segs)) {
    p <- segs[[i]]$coords
    R <- p[, 1] + ro; C <- p[, 2] + co
    ok <- R >= 1L & R <= nr & C >= 1L & C <= nc
    if (any(ok)) m[cbind(R[ok], C[ok])] <- vv[i]
  }
  if (is_rast) {
    out <- terra::rast(template[[1]])
    terra::values(out) <- as.vector(t(m))               # row-major to match terra
    names(out) <- value
    out
  } else m
}

#' Per-order summary of length and diameter
#'
#' Aggregates the per-segment edge table into one row per order class. All
#' lengths and diameters are in the unit of \code{et} (see
#' \code{\link{convert_root_units}}); counts are integers.
#'
#' @param et An \code{edges} table (ideally already in real units via
#'   \code{\link{convert_root_units}} or \code{\link{branch_order_map}}).
#' @param order_col Which order column to group by.
#' @return A data.frame with one row per order:
#'   \describe{
#'     \item{\code{order}}{The order class value.}
#'     \item{\code{n_segments}}{Number of graph segments (edges between nodes).}
#'     \item{\code{n_tips}}{Number of root tips/apices (degree-1 endpoints).}
#'     \item{\code{n_branch_points}}{Number of branching junctions on roots of this
#'       order (= number of laterals departing, for ordinary 3-way forks).}
#'     \item{\code{total_length}}{Summed root length.}
#'     \item{\code{mean_segment_length}}{Mean length per segment.}
#'     \item{\code{branching_frequency}}{Branch points per unit length
#'       (\code{n_branch_points / total_length}).}
#'     \item{\code{mean_diameter}}{Length-weighted mean diameter.}
#'     \item{\code{median_diameter}}{Median of segment median diameters.}
#'   }
#'   Unordered (NA) segments are excluded and counted in \code{attr(., "n_unordered")}.
#' @export
summarize_orders <- function(et, order_col = "branch_order") {
  # Thin wrapper over order_metrics() (focal = NULL) so the per-order
  # aggregation lives in one place. The columns are reshaped to this
  # function's output format: a leading integer `order` column and no
  # `length_fraction` (which order_metrics adds).
  if (is.null(et)) return(data.frame())
  out <- order_metrics(et, order_col = order_col, focal = NULL)
  if (nrow(out) == 0L) return(out)

  n_unordered <- attr(out, "n_unordered")
  out$length_fraction <- NULL
  names(out)[1] <- "order"
  out$order           <- as.integer(out$order)
  out$n_segments      <- as.integer(out$n_segments)
  out$n_tips          <- as.integer(out$n_tips)
  out$n_branch_points <- as.integer(out$n_branch_points)
  out <- out[, c("order", "n_segments", "n_tips", "n_branch_points",
                 "total_length", "mean_segment_length", "branching_frequency",
                 "mean_diameter", "median_diameter")]
  rownames(out) <- NULL
  attr(out, "order_col")   <- order_col
  attr(out, "n_unordered") <- n_unordered
  out
}

#' Convert edge-table lengths and diameters to real units
#'
#' Pixels convert as \code{inch = px / dpi} and \code{cm = px * 2.54 / dpi}.
#' \code{length} is set to the chosen method; \code{length_poly} and
#' \code{length_kimura} are both retained in the same unit.
#'
#' @param et Edge table in pixels (from \code{\link{root_graph_pipeline}}).
#' @param unit One of \code{"cm"}, \code{"inch"}, \code{"px"}.
#' @param dpi Scan resolution (dots per inch).
#' @param length_method \code{"polyline"} (sqrt(2) chain code) or \code{"kimura"}.
#' @return \code{et} with length/diameter columns in \code{unit}; records
#'   \code{unit}, \code{dpi}, \code{length_method} as attributes.
#' @examples
#' \dontrun{
#' et <- root_graph_pipeline(skel, mask)
#' et_cm <- convert_root_units(et, unit = "cm", dpi = 300)
#' }
#' @export
convert_root_units <- function(et, unit = c("cm", "inch", "px"), dpi = 300,
                               length_method = c("polyline", "kimura")) {
  if (is.null(et) || nrow(et) == 0L) return(et)
  unit <- match.arg(unit); length_method <- match.arg(length_method)
  px_factor <- function(u, d) switch(u, px = 1, inch = 1 / d, cm = 2.54 / d)
  # Convert from whatever unit the table is already in, not blindly from pixels,
  # so that applying this twice is idempotent rather than scaling twice. The
  # source unit is read back from the `unit`/`dpi` attributes the last call set.
  from_unit <- attr(et, "unit"); from_dpi <- attr(et, "dpi")
  if (is.null(from_unit)) { from_unit <- "px"; from_dpi <- dpi }
  if (is.null(from_dpi)) from_dpi <- dpi
  f <- px_factor(unit, dpi) / px_factor(from_unit, from_dpi)
  poly_px <- if (is.null(et$length_poly)) et$length else et$length_poly
  et$length_poly     <- poly_px * f
  et$length_kimura   <- et$length_kimura * f
  et$length          <- if (length_method == "kimura") et$length_kimura else et$length_poly
  et$mean_diameter   <- et$mean_diameter * f
  et$median_diameter <- et$median_diameter * f
  et$min_diameter    <- et$min_diameter * f
  attr(et, "unit") <- unit; attr(et, "dpi") <- dpi
  attr(et, "length_method") <- length_method
  et
}

#' Branch-order classification of a root skeleton
#'
#' Main entry point. Converts a binary root skeleton into a per-segment graph,
#' assigns each segment a branching order, and returns the per-segment table, a
#' per-order summary, and a classification raster aligned to the input. Lengths
#' and diameters are reported in real units.
#'
#' @details
#' All three order schemes are always computed and stored on \code{$edges};
#' \code{order} only selects which labels \code{$class_map} and \code{$summary}.
#' \describe{
#'   \item{\code{tip_order} (per segment)}{Topological leaf-peeling depth. Every
#'     terminal segment is order 1; peeling terminals away round by round, a
#'     segment's order is \code{1 + max(child orders)}. Order rises toward the
#'     interior, so the distal end of even a thick root is 1. This is centrifugal
#'     depth, not strict Strahler ordering: Strahler only increments where two
#'     children share the highest order, while this increments at every junction,
#'     so a long chain of single laterals keeps climbing.}
#'   \item{\code{root_order} (per root)}{Segments are grouped into continuous roots
#'     (continuation rule below); each root takes the \emph{maximum} \code{tip_order}
#'     along it, so a thick main axis keeps its high order out to its tip. Sensitive
#'     to how deep the deepest subtree runs.}
#'   \item{\code{branch_order} (per root)}{Centrifugal generation. The thickest root
#'     (length-weighted mean diameter) in each connected component is order 1;
#'     counting branching hops outward, its laterals are 2, theirs 3, and so on.
#'     Independent of subtree depth, so one heavily branched root does not inflate
#'     the rest.}
#' }
#' \strong{Continuation rule} (groups segments into roots for \code{root_order} and
#' \code{branch_order}): at each junction the two arms forming the same root are
#' chosen by \code{straightness + diam_weight * diameter_similarity}.
#' \code{diam_weight = 0} uses angle only; larger values let thickness dominate.
#' The odd arm out is a lateral. \code{tip_order} does not use this rule.
#'
#' @param skel Binary skeleton: single-layer \code{SpatRaster} (preferred) or 0/1
#'   matrix. Pass the raster directly; do not pre-convert with \code{as.matrix()}.
#'   If \code{NULL}, it is computed from \code{mask} via
#'   \code{skeletonize_image()}.
#' @param mask Filled (un-thinned) root mask on the same grid, used for the
#'   distance transform / diameters. If \code{NULL}, the skeleton is used and
#'   diameters collapse to ~1 px. Required if \code{skel} is \code{NULL}.
#' @param order Which scheme labels \code{$class_map}/\code{$summary}:
#'   \code{"branch_order"}, \code{"root_order"}, or \code{"tip_order"}.
#' @param unit Reporting unit: \code{"cm"}, \code{"inch"}, or \code{"px"}.
#' @param dpi Scan resolution (dots per inch); required for cm/inch.
#' @param length_method \code{"polyline"} (sqrt(2) chain code, follows curves) or
#'   \code{"kimura"} (per-segment Kimura correction, better for straight segments).

#' @param template \code{SpatRaster} defining the \code{$class_map} grid; defaults
#'   to \code{skel} when it is a raster.
#' @param overlay_png Optional path; writes the order-colored validation PNG.
#' @param return_map Build \code{$class_map}.
#' @param ... Passed to \code{\link{root_graph_pipeline}} (e.g. \code{dt_backend},
#'   \code{crossing_straight}, \code{prune_iter}, and \code{diam_weight} —
#'   the diameter-vs-angle weight for the continuation rule, >= 0).
#' @return An object of class \code{"branchOrderMap"}: a list with \code{$edges},
#'   \code{$summary}, \code{$class_map} (\code{SpatRaster}; chosen order per pixel,
#'   \code{NA} off-root), and \code{$order}, \code{$unit}, \code{$dpi},
#'   \code{$length_method}.
#' @seealso \code{\link{root_graph_pipeline}}, \code{\link{order_classification_map}},
#'   \code{\link{summarize_orders}}, \code{\link{convert_root_units}}
#' @examples
#' \dontrun{
#' skel <- skeletonize_image(mask)
#' res  <- branch_order_map(skel, mask, order = "branch_order",
#'                          unit = "cm", dpi = 300, dt_backend = "imager")
#' res$summary
#' terra::plot(res$class_map, maxcell = Inf)
#' }
#' @export
branch_order_map <- function(skel = NULL, mask = NULL,
                             order = c("strahler_order", "branch_order", "root_order", "tip_order"),
                             unit = "cm", dpi = 300, length_method = "polyline",
                             template = NULL, overlay_png = NULL, return_map = TRUE, ...) {
  order <- match.arg(order)
  et <- root_graph_pipeline(skel, mask, color_by = order, keep_segments = TRUE,
                            overlay_png = overlay_png, ...)
  et <- convert_root_units(et, unit = unit, dpi = dpi, length_method = length_method)
  res <- list(edges = et, summary = summarize_orders(et, order),
              order = order, unit = unit, dpi = dpi, length_method = length_method)
  if (return_map && (is.null(et) || nrow(et) == 0L)) {
    warning("No root segments found; returning an empty result without a class_map.")
    return(structure(res, class = "branchOrderMap"))
  }
  if (return_map) {
    tmpl <- if (!is.null(template)) template
    else if (inherits(skel, "SpatRaster")) skel
    else if (inherits(mask, "SpatRaster")) mask else NULL
    if (is.null(tmpl)) {
      warning("No SpatRaster template available; skipping class_map (pass template=).")
    } else {
      res$class_map <- order_classification_map(et, tmpl, value = order)
    }
  }
  structure(res, class = "branchOrderMap")
}


#' Aggregate root architecture by order, or split focal-vs-rest
#'
#' Flexible query over a \code{branch_order_map} result (or its \code{$edges}
#' table): length, diameter, and branching architecture, either per order class
#' or split into a focal group versus the rest. Diameter is length-weighted so it
#' reflects how much root length sits at each width. Lengths/diameters are in the
#' unit of \code{x}; counts are integers; \code{branching_frequency} is per unit
#' length.
#'
#' @param x A \code{"branchOrderMap"} object or an \code{edges} data.frame.
#' @param order_col Order column to aggregate by; defaults to the object's
#'   \code{$order}, else \code{"branch_order"}.
#' @param focal Controls the split:
#'   \itemize{
#'     \item \code{NULL} (default): one row per order class.
#'     \item \code{"thinnest"} / \code{"thickest"} (aliases \code{"finest"} /
#'       \code{"coarsest"}): the order class with the smallest / largest
#'       length-weighted diameter, versus all others. Selected by \emph{diameter},
#'       so it is independent of the order numbering; the \code{orders} column
#'       reports which order number(s) fell in each group.
#'     \item a numeric order value or vector (e.g. \code{1} or \code{c(4, 5)}):
#'       those orders as the focal group versus the rest.
#'   }
#' @return A data.frame with one row per group:
#'   \describe{
#'     \item{\code{group}}{Order value (no split) or \code{"focal"}/\code{"rest"}.}
#'     \item{\code{orders}}{Comma-separated order numbers contained in the group
#'       (so a thinnest/thickest split shows which orders it picked).}
#'     \item{\code{n_segments}}{Graph segments in the group.}
#'     \item{\code{n_tips}}{Root tips/apices (degree-1 endpoints).}
#'     \item{\code{n_branch_points}}{Branching junctions on roots of the group.}
#'     \item{\code{total_length}}{Summed root length.}
#'     \item{\code{length_fraction}}{Share of total ordered length.}
#'     \item{\code{mean_segment_length}}{Mean length per segment.}
#'     \item{\code{branching_frequency}}{Branch points per unit length.}
#'     \item{\code{mean_diameter}}{Length-weighted mean diameter.}
#'     \item{\code{median_diameter}}{Median of segment median diameters.}
#'   }
#'   \code{attr(., "n_unordered")} counts NA-order segments; for a split,
#'   \code{attr(., "focal_orders")} lists the focal order number(s).
#' @examples
#' \dontrun{
#' res <- branch_order_map(skel, mask, order = "branch_order", unit = "cm", dpi = 300)
#' order_metrics(res)                       # per-order table
#' order_metrics(res, focal = "thinnest")   # thinnest order class vs all others
#' order_metrics(res, focal = c(1))         # order 1 vs the rest
#' }
#' @export
order_metrics <- function(x, order_col = NULL, focal = NULL) {
  et <- if (inherits(x, "branchOrderMap")) x$edges else x
  # An empty pipeline (e.g. an image with no roots) yields NULL / 0-row edges;
  # return an empty summary rather than erroring on the missing order column.
  if (is.null(et) || nrow(et) == 0L) return(data.frame())
  if (is.null(order_col))
    order_col <- if (inherits(x, "branchOrderMap")) x$order else "branch_order"
  if (!order_col %in% names(et)) stop(sprintf("'%s' is not a column of the edge table.", order_col))
  
  o <- et[[order_col]]; keep <- !is.na(o)
  n_unordered <- sum(!keep)
  o <- o[keep]; len <- et$length[keep]; dia <- et$mean_diameter[keep]; mdn <- et$median_diameter[keep]
  ntip <- if (!is.null(et$n_tips)) et$n_tips[keep] else rep(NA_integer_, length(o))
  nbp  <- if (!is.null(et$n_branch_points)) et$n_branch_points[keep] else rep(NA_integer_, length(o))
  if (length(o) == 0L) return(data.frame())
  total_len <- sum(len)
  
  agg <- function(g) {
    gl <- if (is.numeric(g)) sort(unique(g)) else unique(g)
    do.call(rbind, lapply(gl, function(k) {
      ix <- g == k; w <- len[ix]; L <- sum(w)
      data.frame(
        group               = k,
        orders              = paste(sort(unique(o[ix])), collapse = ","),
        n_segments          = sum(ix),
        n_tips              = sum(ntip[ix]),
        n_branch_points     = sum(nbp[ix]),
        total_length        = L,
        length_fraction     = L / total_len,
        mean_segment_length = mean(w),
        branching_frequency = if (L > 0) sum(nbp[ix]) / L else NA_real_,
        mean_diameter       = if (L > 0) sum(dia[ix] * w) / L else mean(dia[ix]),
        median_diameter     = stats::median(mdn[ix]),
        stringsAsFactors = FALSE
      )
    }))
  }
  
  if (is.null(focal)) {
    res <- agg(o); res$orders <- NULL; names(res)[1] <- order_col   # group == order here
  } else {
    if (is.character(focal) && tolower(focal[1]) %in% c("thinnest","thickest","finest","coarsest")) {
      want_thin <- tolower(focal[1]) %in% c("thinnest", "finest")
      lv <- sort(unique(o))
      wdia <- vapply(lv, function(k) { ix <- o == k; w <- len[ix]
      if (sum(w) > 0) sum(dia[ix] * w) / sum(w) else mean(dia[ix]) }, numeric(1))
      focal_set <- lv[if (want_thin) which.min(wdia) else which.max(wdia)]
    } else {
      focal_set <- as.numeric(focal)
    }
    g <- ifelse(o %in% focal_set, "focal", "rest")
    res <- agg(g)
    attr(res, "focal_orders") <- focal_set
  }
  attr(res, "order_col")   <- order_col
  attr(res, "n_unordered") <- n_unordered
  res
}
