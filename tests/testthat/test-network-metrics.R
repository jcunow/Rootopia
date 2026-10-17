test_that("counts and identities are right on a known tree", {
  et <- data.frame(from = c("A", "B", "B"), to = c("B", "C", "D"),
                   length = c(10, 6, 7), stringsAsFactors = FALSE)
  r <- network_metrics(et, "tree")
  expect_equal(r$traits$n_nodes, 4L)
  expect_equal(r$traits$n_edges, 3L)
  expect_equal(r$traits$n_components, 1L)
  expect_equal(r$traits$n_tips, 3L)          # A, C, D
  expect_equal(r$traits$n_branch_points, 1L) # B
  expect_equal(r$traits$total_length, 23)
  expect_equal(r$traits$beta, 3 / 4)
  expect_equal(r$qc$cyclomatic, 0L)
  expect_true(r$qc$is_forest)
})

test_that("E = V - G holds on a forest of separate edges", {
  et <- data.frame(from = c("A", "C", "E"), to = c("B", "D", "F"),
                   length = rep(1, 3), stringsAsFactors = FALSE)
  r <- network_metrics(et, "tree")
  expect_equal(r$traits$n_edges,
               r$traits$n_nodes - r$traits$n_components)
  expect_equal(r$qc$cyclomatic, 0L)
})

test_that("loop metrics are traits only when the network is reticulate", {
  et <- data.frame(from = c("A", "B", "B", "C"), to = c("B", "C", "D", "D"),
                   length = c(10, 6, 7, 4), stringsAsFactors = FALSE)
  expect_false(is.null(network_metrics(et, "reticulate")$traits$alpha_meshedness))
  expect_null(network_metrics(et, "tree")$traits$alpha_meshedness)
  expect_null(network_metrics(et, "projected_tree")$traits$alpha_meshedness)
  expect_false(is.null(network_metrics(et, "tree")$qc$alpha_meshedness))
  expect_false(is.null(network_metrics(et, "projected_tree")$qc$artefact_cycles))
})

test_that("auto calls an acyclic graph a tree and anything else projected", {
  acyc <- data.frame(from = c("A", "B"), to = c("B", "C"),
                     length = c(1, 1), stringsAsFactors = FALSE)
  cyc  <- rbind(acyc, data.frame(from = "C", to = "A", length = 1))
  expect_equal(network_metrics(acyc, "auto")$topology, "tree")
  expect_equal(network_metrics(cyc,  "auto")$topology, "projected_tree")
})

test_that("artefact_rate normalizes the cycle count by network length", {
  et <- data.frame(from = c("A", "B", "C"), to = c("B", "C", "A"),
                   length = c(10, 10, 5), stringsAsFactors = FALSE)
  r <- network_metrics(et, "projected_tree")
  expect_equal(r$qc$artefact_cycles, 1L)
  expect_equal(r$qc$artefact_rate, 1 / 25)
  no_len <- et; no_len$length <- NULL
  expect_true(is.na(network_metrics(no_len, "projected_tree")$qc$artefact_rate))
})

test_that("degenerate input is handled without erroring", {
  et <- data.frame(from = "A", to = "B", length = 3, stringsAsFactors = FALSE)
  expect_type(network_metrics(NULL), "list")
  expect_length(network_metrics(NULL)$traits, 0L)
  expect_length(network_metrics(et[0, ])$traits, 0L)
  expect_length(network_metrics(list(a = 1))$traits, 0L)
  # two nodes: both planar maxima are non-positive, so the ratios are undefined
  r <- network_metrics(et, "reticulate")
  expect_true(is.na(r$traits$alpha_meshedness))
  expect_true(is.na(r$traits$gamma))
})

test_that("malformed input fails with a message that names the problem", {
  expect_error(network_metrics(data.frame(a = 1, b = 2)), "from.*to")
  expect_error(network_metrics(data.frame(from = NA, to = "B")), "NA")
})

test_that("self-loops and parallel edges are reported, not silently counted", {
  sl <- data.frame(from = c("A", "B"), to = c("B", "B"),
                   length = c(1, 1), stringsAsFactors = FALSE)
  expect_equal(network_metrics(sl, "tree")$qc$n_self_loops, 1L)
  dp <- data.frame(from = c("A", "A"), to = c("B", "B"),
                   length = c(1, 1), stringsAsFactors = FALSE)
  expect_equal(network_metrics(dp, "tree")$qc$n_duplicate_edges, 1L)
  # both raise the cycle count, which is why they are worth surfacing
  expect_equal(network_metrics(sl, "tree")$qc$cyclomatic, 1L)
  expect_equal(network_metrics(dp, "tree")$qc$cyclomatic, 1L)
})
