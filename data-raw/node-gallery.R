# ---------------------------------------------------------------------------
# What Rootopia does at a node
#
# Ten synthetic nodes, three to six arms, each run through
# root_graph_pipeline() at its default settings. Nothing here needs a scan.
#
#   source("data-raw/node-gallery.R")
#
# That writes both figures the overview vignette uses and prints what came
# back for each node:
#
#   vignettes/figures/node_gallery.png       the ten nodes
#   vignettes/figures/node_alternatives.png  the runner-up reading of four
#
# Paths are relative to the working directory, so run it from the package
# root. gallery(which = 8:10) redraws a subset larger, for a closer look.
#
# THE TEN NODES
#
#   arms  node
#   3     1  a wide Y: a thin lateral leaves a thick axis at 60 degrees
#   3     2  a narrow fork: two thick arms leaving 20 degrees apart
#   3     3  a thin arm continuing straight, a thick arm bending 50 degrees
#   4     4  an axis with two thin laterals both sweeping forward
#   4     5  two roots crossing square
#   4     6  two roots crossing at 20 degrees, nearly parallel
#   4     7  the same 20 degree crossing, one root 21 px wide and the other 9
#   5     8  a crossing with a lateral at 70 degrees, nearer one of the two
#   5     9  a crossing with a lateral at 45 degrees, between the two roots
#   6    10  three roots through one point, each arm carrying a stub
#
# Diameters are the drawn widths in pixels.
#
# A ROOT HERE MEANS AN AXIS: a path the tracer read as going on through its
# nodes. Its laterals are branches of it, not roots in their own right, so the
# counts under each panel give the two separately, along with the junction
# points: node 1 is one root with one branch, node 5 is two roots with one
# junction, node 10 is three roots with six branches and three junctions.
#
# HOW THE ARMS ARE READ
#
# At a node of four arms or more, resolve_crossings() pairs arms wherever a
# pair runs through the node. Each arm's direction comes from its last look = 5
# pixels pointing away; a pair scores the dot product of the two, -1 running
# straight through and 0 turning a right angle, and a pair may be joined when
# it scores below crossing_straight (default -0.5, so a root may bend by at
# most 60 degrees). Of the ways to pair the arms, the one pairing the MOST arms
# wins, ties going to the straightest total. An arm left over is a lateral,
# and the pair it leaves is held back from splicing so that the lateral stays
# attached and the node is read as an ordinary junction. A node of three arms
# is a junction to begin with, and assign_root_order() picks which two arms are
# one root going on, by straightness and diameter together through diam_weight.
#
# WHAT EACH NODE IS FOR
#
#   1-3  the junction rule, from easy to contested. In 1 the straight arm is
#        also the thick one. In 2 both arms are thick and nearly aligned, so
#        the two candidate pairs score 1.39 and 1.23. In 3 they disagree: the
#        arm continuing straight is 5 px and the arm bending 50 degrees is 19,
#        as thick as the stem, and diam_weight decides. The thick arm wins
#        above about 0.38, so at the default 0.5 the root is read as turning 50
#        degrees. That is the one node where the default is a judgement call.
#   4    four arms that are not a crossing: both laterals sweep forward, 90
#        degrees apart, so the only pair that runs through is the axis, and
#        holding it back for the laterals leaves the node a junction.
#   5-7  crossings. The angle BETWEEN the two roots is not what is read: 5 and
#        6 differ only in that angle and come out alike. 7 adds a width
#        difference, and comes out alike again, because the pairing reads
#        directions only. Its runner-up pairing, which keeps each root at one
#        width and kinks both by 20 degrees, is 0.12 behind.
#   8    an arm left over. The two crossing roots pair at -1 each, the lateral
#        at 70 degrees pairs with nothing straight enough, so it is held as a
#        lateral: 0.94 to the vertical arm against 0.34 to the horizontal, and
#        the angle settles which root it leaves on its own.
#   9    the same node with the lateral at 45 degrees, where the angle ties at
#        0.71 against both crossing roots and says nothing. The thinner pair
#        breaks the tie, so the lateral leaves the thin root.
#   10   three roots crossing at one pixel, each with stubs: three pairs, all
#        straight, nothing left over.
#
# HOW SHALLOW A CROSSING SURVIVES
#
# Two one-pixel centerlines crossing at t degrees share pixels for about
# 1/tan(t) of their length either side of the crossing: 1.4 px at 20 degrees,
# 1.9 at 15, 2.8 at 10. Below roughly 20 degrees the shared run is long enough
# that the crossing is no longer one point in the skeleton, and the tracer
# returns two three-arm junctions a few pixels apart joined by a short segment.
# Such a node never reaches resolve_crossings() at all. The scoring margin
# narrows with the angle too -- 0.12 at 20 degrees, 0.07 at 15, 0.03 at 10 --
# but the tracing gives way first.
#
# WHAT THE PANELS SHOW
#
# Each hue is one root. A lateral is drawn as a darker tone of the root it
# leaves and narrower than it, so an axis and its branches are told apart by
# both. The gray underneath is the segmentation, so stroke thickness stays
# visible, and the dashed lines are the input centerlines.
# ---------------------------------------------------------------------------

if (!requireNamespace("Rootopia", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(Rootopia)
stopifnot(requireNamespace("terra", quietly = TRUE))
if (!interactive()) grDevices::graphics.off()   # leave a screen device alone

# --- drawing ---------------------------------------------------------------
# A stroke is a straight piece in image coordinates: x is the column, y the row
# counted from the top. It carries the diameter it is drawn with and the name
# of the root it belongs to, which is the truth we compare against.

line_px <- function(p0, p1) {                    # integer points along a line
  n <- max(abs(p1 - p0)) + 1
  cbind(x = round(seq(p0[1], p1[1], length.out = n)),
        y = round(seq(p0[2], p1[2], length.out = n)))
}

disc_px <- function(m, pts, r) {                 # disc of radius r at each point
  k <- ceiling(r)
  off <- as.matrix(expand.grid(dx = -k:k, dy = -k:k))
  off <- off[off[, 1]^2 + off[, 2]^2 <= r^2 + 1e-9, , drop = FALSE]
  for (i in seq_len(nrow(pts))) {
    xx <- pts[i, 1] + off[, 1]; yy <- pts[i, 2] + off[, 2]
    ok <- xx >= 1 & xx <= ncol(m) & yy >= 1 & yy <= nrow(m)
    m[cbind(yy[ok], xx[ok])] <- 1
  }
  m
}

# The skeleton and the segmentation come from the same centrelines, so the
# skeleton is one pixel wide by construction and sits inside the segmentation.
paint_node <- function(strokes, size = 260L) {
  skel <- mask <- matrix(0, size, size)
  for (s in strokes) {
    pts  <- line_px(s$from, s$to)
    skel <- disc_px(skel, pts, 0)
    mask <- disc_px(mask, pts, (s$d - 1) / 2)
  }
  list(skel = skel, mask = mask, strokes = strokes)
}

# arms radiating from the canvas centre, given as angle (degrees, y up),
# diameter and truth label
star <- function(arms, L = 108, size = 260L) {
  c0 <- size / 2
  lapply(arms, function(a) list(
    from = c(c0, c0),
    to   = c(c0 + L * cos(a$ang * pi / 180), c0 - L * sin(a$ang * pi / 180)),
    d = a$d, root = a$root, ang = a$ang))
}

arm <- function(ang, d, root) list(ang = ang, d = d, root = root)

# a short side branch leaving one arm partway out, for a shape that needs
# branches of its own rather than only a central node
stub <- function(ang, at = 0.62, len = 36, d = 5, turn = 45, L = 108,
                 size = 260L, root = "stub") {
  c0 <- size / 2
  b  <- c(c0 + at * L * cos(ang * pi / 180), c0 - at * L * sin(ang * pi / 180))
  a2 <- ang + turn
  list(from = b, to = c(b[1] + len * cos(a2 * pi / 180),
                        b[2] - len * sin(a2 * pi / 180)),
       d = d, root = root)
}

# --- the eight nodes -------------------------------------------------------

shapes <- list(
  list(label = "1  wide Y, 3 arms",
       note = "thin lateral at 60 deg; straightest arm is also thickest",
       arms = list(arm(270, 17, "axis"), arm( 90, 15, "axis"),
                   arm( 30,  5, "lat"))),

  list(label = "2  narrow fork, 3 arms",
       note = "two thick arms only 20 deg apart",
       arms = list(arm(270, 17, "axis"), arm( 90, 15, "axis"),
                   arm( 70, 13, "lat"))),

  list(label = "3  thin straight vs thick bent, 3 arms",
       note = "straight arm 5 px wide, bent arm 19 px, same as the stem",
       arms = list(arm(270, 19, "stem"), arm( 85,  5, "thin"),
                   arm(140, 19, "thick"))),

  list(label = "4  axis with two forward laterals, 4 arms",
       note = "best pairing turns 90 deg, so it is refused",
       arms = list(arm( 90, 19, "axis"), arm(270, 19, "axis"),
                   arm(225,  5, "lat1"), arm(315,  5, "lat2"))),

  list(label = "5  square crossing, 4 arms",
       note = "both roots run straight across the node",
       arms = list(arm(  0, 17, "A"), arm(180, 17, "A"),
                   arm( 90, 17, "B"), arm(270, 17, "B"))),

  list(label = "6  20 deg crossing, 4 arms",
       note = "the same case as shape 5, drawn as a narrow X",
       arms = list(arm(  0, 15, "A"), arm(180, 15, "A"),
                   arm( 20, 15, "B"), arm(200, 15, "B"))),

  list(label = "7  20 deg crossing, thick by thin, 4 arms",
       note = "21 px and 9 px roots; straightness and thickness disagree",
       arms = list(arm(  0, 21, "thick"), arm(200, 21, "thick"),
                   arm(180,  9, "thin"),  arm( 20,  9, "thin"))),

  list(label = "8  70 deg lateral on a crossing, 5 arms",
       note = "lateral nearer the vertical root; angle alone decides",
       arms = list(arm( 90, 19, "axis"),  arm(270, 19, "axis"),
                   arm(  0, 15, "cross"), arm(180, 15, "cross"),
                   arm( 70,  5, "lat"))),

  list(label = "9  45 deg lateral on a crossing, 5 arms",
       note = "lateral exactly between the two roots; thickness decides",
       arms = list(arm( 90, 21, "thick"), arm(270, 21, "thick"),
                   arm(  0,  9, "thin"),  arm(180,  9, "thin"),
                   arm( 45,  5, "lat"))),

  list(label = "10  three roots through one point, 6 arms",
       note = "each arm carries a stub, so axes and branches both show",
       arms = list(arm(  0, 15, "A"), arm(180, 15, "A"),
                   arm( 60, 15, "B"), arm(240, 15, "B"),
                   arm(120, 15, "C"), arm(300, 15, "C")),
       extra = lapply(c(0, 60, 120, 180, 240, 300), stub))
)

for (i in seq_along(shapes)) {
  extra <- if (is.null(shapes[[i]]$extra)) list() else shapes[[i]]$extra
  shapes[[i]]$img  <- paint_node(c(star(shapes[[i]]$arms), extra))
  shapes[[i]]$narm <- length(shapes[[i]]$arms)   # arms at the central node
}

# --- run -------------------------------------------------------------------

run <- function(img) {
  Rootopia::root_graph_pipeline(
    skel = terra::rast(img$skel), mask = terra::rast(img$mask),
    crop = FALSE, keep_segments = TRUE, verbose = FALSE,
    color_by = "order_from_base")
}

# --- what the package returns ----------------------------------------------

orders_of <- function(et) {
  ob <- et[["order_from_base"]]
  if (is.null(ob)) rep(1L, length(attr(et, "segments"))) else as.integer(ob)
}

# root_id groups the segments the ordering rule united into one continuous root
roots_of <- function(et) {
  r <- et[["root_id"]]
  if (is.null(r)) seq_along(attr(et, "segments")) else as.integer(r)
}

# A junction point is a place where two roots overlap. The package does not
# report them -- a resolved crossing leaves no node behind, the two segments
# simply share pixels -- so they are counted here from the segment list: a
# pixel carried by m segments is choose(m, 2) overlaps. On these synthetic
# nodes two centerlines cross in exactly one pixel, so the count is the number
# of crossings; on a real scan a crossing can share a short run of pixels and
# this would count it more than once.
junction_points <- function(et) {
  sg  <- attr(et, "segments")
  key <- unlist(lapply(sg, function(s)
    unique(paste(s$coords[, 1], s$coords[, 2]))))
  m   <- table(key)
  sum(choose(as.integer(m[m > 1L]), 2))
}

# What the pipeline came back with: which segments belong to one root, how deep
# each is, and the two counts.
grouping <- function(et) {
  list(rid      = roots_of(et),
       ord      = orders_of(et),
       n_branch = sum(et$n_branch_points, na.rm = TRUE),
       n_junc   = junction_points(et))
}

# Segment lines are one pixel wide, so they have to be grown to be seen at
# panel size. An axis is drawn wide and a lateral narrow, which separates a
# main root from its branches by width as well as by shade.
stamp_val <- function(m, pts, r, v) {
  k <- ceiling(r)
  off <- as.matrix(expand.grid(dx = -k:k, dy = -k:k))
  off <- off[off[, 1]^2 + off[, 2]^2 <= r^2 + 1e-9, , drop = FALSE]
  for (i in seq_len(nrow(pts))) {
    xx <- pts[i, 2] + off[, 1]; yy <- pts[i, 1] + off[, 2]
    ok <- xx >= 1 & xx <= ncol(m) & yy >= 1 & yy <= nrow(m)
    m[cbind(yy[ok], xx[ok])] <- v
  }
  m
}

overlay <- function(et, dims, value, ord) {
  sg <- attr(et, "segments")
  m  <- matrix(NA_real_, dims[1], dims[2])
  for (i in order(ord))                # axes first, laterals painted over them
    m <- stamp_val(m, sg[[i]]$coords, if (ord[i] <= 1L) 3.2 else 1.4, value[i])
  m
}

# The Okabe-Ito qualitative set, ordered strongest first. It avoids pairing red
# with green, which is the one combination a red-green color vision deficiency
# cannot separate, and the hues stay apart in grayscale print.
root_hues <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
               "#56B4E9", "#D55E00", "#F0E442")

# A lateral is drawn as a darker tone of the root it leaves: scaling the
# channels keeps the hue, and darkening rather than lightening keeps it
# readable against the pale gray of the segmentation.
shade <- function(col, k) {
  v <- grDevices::col2rgb(col) / 255
  f <- max(0.42, 1 - 0.32 * k)          # 1 for an axis, 0.68 for its laterals
  grDevices::rgb(v[1] * f, v[2] * f, v[3] * f)
}

# One hue per root system, one shade per branch order inside it. A root's
# parent is the root it shares a node with that sits one order shallower, and
# the hue comes from the order-1 ancestor at the top of that chain.
segment_colors <- function(et, grp) {
  rid <- grp$rid; ob <- grp$ord
  roots <- sort(unique(rid))
  r_ord <- vapply(roots, function(r) min(ob[rid == r]), numeric(1))

  # A lateral's parent can be handed in, which the alternative-reading figure
  # does; otherwise it is the root it shares a node with that sits one order
  # shallower. Since a resolved crossing leaves no node behind, the only roots
  # meeting at a node are a continuation and its laterals, so that is
  # unambiguous.
  if (!is.null(grp$parent)) {
    parent <- grp$parent
  } else {
    parent <- rep(NA_integer_, length(roots))
    nodes  <- c(et$from, et$to); nrid <- rep(rid, 2)
    for (nd in unique(nodes)) {
      here <- unique(nrid[nodes == nd])
      if (length(here) < 2L) next
      o <- r_ord[match(here, roots)]
      pa <- here[which.min(o)]
      for (q in here[o > min(o)]) {
        k <- match(q, roots)
        if (is.na(parent[k])) parent[k] <- pa
      }
    }
  }

  # walk up to the shallowest root of each family; the hop limit is there
  # because an unresolved tangle can leave a loop in the parent links
  anc <- vapply(seq_along(roots), function(k) {
    for (hop in 1:10) {
      if (is.na(parent[k])) break
      k <- match(parent[k], roots)
    }
    roots[k]
  }, numeric(1))

  fam  <- unique(anc)
  base <- root_hues[(seq_along(fam) - 1L) %% length(root_hues) + 1L]
  cols <- vapply(seq_along(roots), function(k)
    shade(base[match(anc[k], fam)], r_ord[k] - 1), character(1))
  cols[match(rid, roots)]
}

# Every panel carries the same two numbers, computed here from the drawn
# geometry so they cannot drift from the figure: the score of the arrangement
# the rule picks, and the score of the runner-up. What is being scored depends
# only on how many arms meet, which is what separates the node classes.
#
#   three arms      the continuation score, straight + diam_weight * dsim *
#                   dmag, higher being better: which two arms are one root
#                   going on through the node.
#   four or more    the total dot product of a pairing of the arms, lower being
#                   better, taken over every way of pairing them up whether or
#                   not a pair is straight enough to be allowed. Acceptance is
#                   a separate test, on the WORST single pair of the winning
#                   pairing, which is printed as well: above -0.5 and the
#                   pairing is refused however good its total.
#   an arm left over  where the winning pairing leaves an arm unpaired, the
#                   parent score of the best and second-best candidate root,
#                   higher being better. Shapes whose arms all pair up have no
#                   such choice to make, so the line is absent.
#
# Diameters here are the drawn widths in pixels. The pipeline measures them
# from the distance transform instead, so its numbers can differ slightly.
spec_scores <- function(sh, diam_weight = 0.5, parent_thin_weight = 0.5,
                        pair_diam_weight = 0) {
  ang <- vapply(sh$arms, function(a) a$ang, numeric(1))
  dia <- vapply(sh$arms, function(a) a$d,   numeric(1))
  k   <- length(ang)
  u   <- cbind(cos(ang * pi / 180), sin(ang * pi / 180))
  dm  <- u %*% t(u)

  if (k == 3L) {
    sc <- c()
    for (a in 1:2) for (b in (a + 1L):3L) {
      di <- dia[a]; dj <- dia[b]
      sc <- c(sc, -dm[a, b] + diam_weight * (min(di, dj) / max(di, dj)) *
                (min(di, dj) / max(dia)))
    }
    sc <- sort(sc, decreasing = TRUE)
    return(sprintf("continuation %.2f, next %.2f", sc[1], sc[2]))
  }

  mats <- list()
  rec <- function(free, pairs, tot) {
    if (length(free) < 2L) {
      mats[[length(mats) + 1L]] <<- list(n = length(pairs), tot = tot,
                                         pairs = pairs)
      return(invisible(NULL))
    }
    a <- free[1]
    for (b in free[-1])
      rec(setdiff(free, c(a, b)), c(pairs, list(c(a, b))),
          tot + dm[a, b] - pair_diam_weight *
            min(dia[a], dia[b]) / max(dia[a], dia[b]))
    rec(free[-1], pairs, tot)                            # leave this arm over
  }
  rec(seq_len(k), list(), 0)

  keep <- Filter(function(m) m$n == max(vapply(mats, `[[`, numeric(1), "n")),
                 mats)
  tot  <- vapply(keep, `[[`, numeric(1), "tot")
  o    <- order(tot)
  best <- keep[[o[1]]]
  line <- sprintf("pairing %.2f, next %.2f; worst pair %.2f", tot[o[1]],
                  if (length(o) > 1L) tot[o[2]] else NA_real_,
                  max(vapply(best$pairs, function(q) dm[q[1], q[2]],
                             numeric(1))))

  over <- setdiff(seq_len(k), unlist(best$pairs))
  if (length(over) == 1L) {
    cand <- unlist(best$pairs)
    ps <- sort(dm[over, cand] + parent_thin_weight * (1 - dia[cand] / max(dia)),
               decreasing = TRUE)
    line <- c(line, sprintf("parent %.2f, next %.2f", ps[1], ps[2]))
  }
  line
}

panel <- function(sh, et, grp, lines = NULL, head = sh$label,
                  note = sh$note) {
  img  <- sh$img
  dims <- dim(img$mask)
  bg <- img$mask; bg[bg == 0] <- NA
  terra::plot(terra::rast(bg), col = "grey84", breaks = c(0.5, 1.5),
              legend = FALSE, axes = FALSE,
              mar = c(0.3, 0.3, if (is.null(lines)) 3.0 else 5.4, 0.3))

  rid  <- grp$rid
  ob   <- grp$ord
  cols <- segment_colors(et, grp)
  uc   <- unique(cols)
  sm   <- overlay(et, dims, match(cols, uc), ob)
  terra::plot(terra::rast(sm), add = TRUE, legend = FALSE,
              breaks = c(0.5, seq_along(uc) + 0.5), col = uc)

  for (s in img$strokes)                         # truth, in image coordinates
    graphics::lines(c(s$from[1], s$to[1]) - 0.5,
                    dims[1] - c(s$from[2], s$to[2]) + 0.5,
                    lty = 3, lwd = 0.7, col = "black")

  # A root is an axis: a path the tracer read as going on through its nodes.
  # Its laterals are branches of it, not roots of their own, so they are
  # counted separately.
  n_root <- length(unique(rid[ob <= 1L]))
  n_br   <- if (is.null(grp$n_branch)) length(unique(rid[ob > 1L])) else grp$n_branch
  txt <- sprintf("%d root%s", n_root, if (n_root == 1L) "" else "s")
  if (n_br > 0L)
    txt <- paste0(txt, sprintf(", %d branch point%s", n_br,
                               if (n_br == 1L) "" else "s"))
  if (!is.null(grp$n_junc) && grp$n_junc > 0L)
    txt <- paste0(txt, sprintf(", %d junction%s", grp$n_junc,
                               if (grp$n_junc == 1L) "" else "s"))
  graphics::text(dims[2] / 2, 6, adj = c(0.5, 0), cex = 0.7, col = "grey20",
                 labels = txt)
  top <- if (is.null(lines)) 1.5 else 4.0
  graphics::mtext(head, side = 3, line = top, cex = 0.72, font = 2)
  graphics::mtext(note, side = 3, line = top - 0.9, cex = 0.58, col = "grey30")
  for (j in seq_along(lines))
    graphics::mtext(lines[j], side = 3, line = 2.0 - (j - 1L), cex = 0.56,
                    col = "grey15")
}

# ---------------------------------------------------------------------------
# THE RUNNER-UP READING
#
# Where the winning arrangement is only a little better than the next one, the
# node is close to ambiguous and one panel tells half the story. ambiguity()
# draws the alternative for those shapes: the second-best pairing where two
# pairings are close, or, where the pairing is settled but the lateral could as
# well have gone to the other root, the second-best parent. Each panel gives
# both scores and the gap between them, so how firm the reading is can be seen
# rather than assumed.
# ---------------------------------------------------------------------------

open_png <- function(file, width, height) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(outdir))
    stop("cannot create '", outdir, "' from '", getwd(),
         "'. Run this from the package root, or set outdir to an ",
         "absolute path.", call. = FALSE)
  grDevices::png(file.path(outdir, file), width = width, height = height,
                 res = 150)
}

# which arm each segment is, from the direction of its far end
arm_index <- function(sh, et, size = 260L) {
  c0  <- size / 2
  ang <- vapply(sh$arms, function(a) a$ang, numeric(1))
  sg  <- attr(et, "segments")
  vapply(seq_along(sg), function(i) {
    cc  <- sg[[i]]$coords
    far <- which.max((cc[, 1] - c0)^2 + (cc[, 2] - c0)^2)
    th  <- atan2(c0 - cc[far, 1], cc[far, 2] - c0) * 180 / pi
    which.min(abs(((ang - th + 180) %% 360) - 180))
  }, integer(1))
}

# every way of pairing the arms, scored as P2 scores them, best first
arm_pairings <- function(sh, pair_diam_weight = 0) {
  ang <- vapply(sh$arms, function(a) a$ang, numeric(1))
  dia <- vapply(sh$arms, function(a) a$d,   numeric(1))
  u   <- cbind(cos(ang * pi / 180), sin(ang * pi / 180))
  dm  <- u %*% t(u)
  out <- list()
  rec <- function(free, pairs, tot) {
    if (length(free) < 2L) {
      out[[length(out) + 1L]] <<- list(n = length(pairs), tot = tot,
                                       pairs = pairs, over = free)
      return(invisible(NULL))
    }
    a <- free[1]
    for (b in free[-1])
      rec(setdiff(free, c(a, b)), c(pairs, list(c(a, b))),
          tot + dm[a, b] -
            pair_diam_weight * min(dia[a], dia[b]) / max(dia[a], dia[b]))
    rec(free[-1], pairs, tot)
  }
  rec(seq_along(ang), list(), 0)
  keep <- Filter(function(m) m$n == max(vapply(out, `[[`, numeric(1), "n")), out)
  keep[order(vapply(keep, `[[`, numeric(1), "tot"))]
}

# the candidate parents of a pairing's leftover arm, best first
parent_rank <- function(sh, m, parent_thin_weight = 0.5) {
  if (length(m$over) != 1L) return(NULL)
  ang  <- vapply(sh$arms, function(a) a$ang, numeric(1))
  dia  <- vapply(sh$arms, function(a) a$d,   numeric(1))
  cand <- unlist(m$pairs)
  sc   <- cos((ang[m$over] - ang[cand]) * pi / 180) +
    parent_thin_weight * (1 - dia[cand] / max(dia))
  o <- order(sc, decreasing = TRUE)
  list(arm = cand[o], score = sc[o])
}

# the grouping one pairing implies, as root ids over the arms
grouping_from_pairs <- function(sh, m, parent_arm = NA_integer_) {
  rid <- integer(length(sh$arms))
  for (i in seq_along(m$pairs)) rid[m$pairs[[i]]] <- i
  nxt <- length(m$pairs)
  for (a in m$over) { nxt <- nxt + 1L; rid[a] <- nxt }
  ord <- rep(1L, nxt); parent <- rep(NA_integer_, nxt)
  if (length(m$over) > 0L) ord[rid[m$over]] <- 2L
  if (length(m$over) == 1L && !is.na(parent_arm))
    parent[rid[m$over]] <- rid[parent_arm]
  list(rid_arm = rid, ord_root = ord, parent = parent,
       n_junc = as.integer(choose(length(m$pairs), 2)),
       n_branch = length(m$over))
}

ambiguity <- function(which = 6:9) {
  use <- shapes[which]
  open_png("node_alternatives.png", 520 * length(use), 600)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, length(use)))

  tab <- list()
  for (sh in use) {
    et <- run(sh$img)
    ai <- arm_index(sh, et)
    pr  <- arm_pairings(sh)
    gap <- pr[[2]]$tot - pr[[1]]$tot
    pa  <- parent_rank(sh, pr[[1]])

    if (!is.null(pa)) {        # an arm was left over, so a parent was chosen
      g    <- grouping_from_pairs(sh, pr[[1]], pa$arm[2])
      what <- "parent"
      lines <- c(sprintf("this parent %.2f, the chosen one %.2f",
                         pa$score[2], pa$score[1]),
                 sprintf("parent gap %.2f, pairing gap %.2f",
                         pa$score[1] - pa$score[2], gap))
    } else {                   # every arm paired, so the pairing is the doubt
      alt  <- pr[[2]]
      g    <- grouping_from_pairs(sh, alt)
      what <- "pairing"
      lines <- c(sprintf("this pairing %.2f, the chosen one %.2f",
                         alt$tot, pr[[1]]$tot),
                 sprintf("pairing gap %.2f", gap))
    }

    grp <- list(rid = g$rid_arm[ai], ord = g$ord_root[g$rid_arm[ai]],
                parent = g$parent, n_junc = g$n_junc, n_branch = g$n_branch)
    panel(sh, et, grp, lines = lines,
          head = paste0(sub("^([0-9]+).*", "\\1", sh$label),
                        "  the other ", what),
          note = sh$note)
    tab[[length(tab) + 1L]] <- data.frame(shape = sh$label, doubt = what,
                                          gap = round(gap, 3))
  }
  cat("\nthe runner-up reading -- wrote ",
      file.path(outdir, "node_alternatives.png"), "\n", sep = "")
  print(do.call(rbind, tab), row.names = FALSE)
  invisible(tab)
}

outdir <- "vignettes/figures"      # where the vignette looks for them

gallery <- function(which = seq_along(shapes), ncol = 5L) {
  use  <- shapes[which]
  fits <- lapply(use, function(sh) run(sh$img))
  grps <- lapply(fits, grouping)

  nr <- ceiling(length(use) / ncol)
  open_png("node_gallery.png", 520 * min(ncol, length(use)), 540 * nr)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(nr, min(ncol, length(use))))
  for (i in seq_along(use)) panel(use[[i]], fits[[i]], grps[[i]])

  tab <- do.call(rbind, lapply(seq_along(use), function(i) {
    ob <- grps[[i]]$ord; rid <- grps[[i]]$rid
    data.frame(shape = use[[i]]$label, arms = use[[i]]$narm,
               roots = length(unique(rid[ob <= 1L])),
               branch_pts = if (is.null(grps[[i]]$n_branch)) NA_integer_ else
                 grps[[i]]$n_branch,
               junctions = if (is.null(grps[[i]]$n_junc)) NA_integer_ else
                 grps[[i]]$n_junc,
               segments = length(attr(fits[[i]], "segments")),
               orders = paste(sort(unique(ob)), collapse = "-"))
  }))
  cat("\nwrote ", file.path(outdir, "node_gallery.png"), "\n", sep = "")
  print(tab, row.names = FALSE)
  cat("\nscores from the drawn geometry\n")
  for (sh in use)
    cat(sprintf("  %-40s %s\n", sh$label,
                paste(spec_scores(sh), collapse = "  |  ")))
  invisible(tab)
}

gallery()
ambiguity()
