# root_transport on a hand-computed tree.
#
#   A -e1- B -e2- C -e3- D        e1 is the thickest segment, so A is the base
#          |             |        of the component and every quantity below is
#          e4            e5       measured outward from A.
#          E             F
#
# All five segments are 10 long; e1 has diameter 4, the rest 2. The expected
# values are written out rather than recomputed from the function.
tree5 <- function() {
  data.frame(from = c("A", "B", "C", "B", "C"),
             to   = c("B", "C", "D", "E", "F"),
             length = rep(10, 5),
             mean_diameter = c(4, 2, 2, 2, 2),
             stringsAsFactors = FALSE)
}

test_that("downstream totals are the subtree beyond each segment", {
  out <- root_transport(tree5())
  expect_equal(out$downstream_tips, c(3L, 2L, 1L, 1L, 1L))
  expect_equal(out$downstream_length, c(50, 30, 10, 10, 10))
  # the basal segment serves the whole component
  expect_equal(out$downstream_length[1], sum(tree5()$length))
  expect_equal(out$downstream_tips[1], 3L)
  expect_equal(out$downstream_surface[1], pi * 4 * 10 + 4 * (pi * 2 * 10))
})

test_that("path quantities accumulate from the base outward", {
  out <- root_transport(tree5())
  expect_equal(out$path_length, c(10, 20, 30, 20, 30))

  r1 <- 8 * 10 / (pi * 2^4)   # e1: radius 2
  r2 <- 8 * 10 / (pi * 1^4)   # the rest: radius 1
  expect_equal(out$resistance, c(r1, r2, r2, r2, r2))
  expect_equal(out$path_resistance,
               c(r1, r1 + r2, r1 + 2 * r2, r1 + r2, r1 + 2 * r2))
  # strictly increasing along the main axis A -> D
  expect_true(all(diff(out$path_resistance[1:3]) > 0))
})

test_that("eta scales resistance but not length", {
  a <- root_transport(tree5(), eta = 1)
  b <- root_transport(tree5(), eta = 3)
  expect_equal(b$resistance, 3 * a$resistance)
  expect_equal(b$path_length, a$path_length)
})

test_that("the thickest free end is the base", {
  et <- tree5()
  et$mean_diameter <- c(2, 2, 4, 2, 2)   # now e3 is thickest, so D is the base
  out <- root_transport(et)
  expect_equal(out$downstream_tips[3], 3L)
  expect_equal(out$path_length[3], 10)
})

test_that("components are rooted independently", {
  et <- rbind(tree5(),
              data.frame(from = "X", to = "Y", length = 5, mean_diameter = 1))
  out <- root_transport(et)
  expect_equal(out$downstream_tips[6], 1L)
  expect_equal(out$downstream_length[6], 5)
  expect_equal(out$downstream_tips[1:5], root_transport(tree5())$downstream_tips)
})

test_that("a segment on a cycle gets NA and a warning", {
  et <- rbind(tree5(),
              data.frame(from = "E", to = "F", length = 7, mean_diameter = 2))
  expect_warning(out <- root_transport(et), "cycle")
  expect_true(is.na(out$downstream_tips[6]))
  expect_true(is.na(out$path_resistance[6]))
  # resistance is per-segment geometry, so it is still defined
  expect_true(is.finite(out$resistance[6]))
})

test_that("the radius floor keeps resistance finite and is unit-free", {
  et <- tree5()
  et$mean_diameter[2] <- 0
  out <- root_transport(et)
  expect_true(all(is.finite(out$resistance)))
  # the default floor is the thinnest radius present, here 1
  expect_equal(out$resistance[2], 8 * 10 / (pi * 1^4))

  # an explicit floor is taken in the table's own unit
  expect_equal(root_transport(et, min_radius = 0.5)$resistance[2],
               8 * 10 / (pi * 0.5^4))

  # converting the whole table to another unit rescales resistance by a power
  # of the factor instead of hitting a fixed floor
  cm <- tree5(); f <- 2.54 / 300
  cm$length <- cm$length * f; cm$mean_diameter <- cm$mean_diameter * f
  expect_equal(root_transport(cm)$resistance,
               root_transport(tree5())$resistance / f^3)
})

test_that("a table with no usable diameter gives NA resistance, not Inf", {
  et <- tree5()
  et$mean_diameter <- 0
  out <- root_transport(et)
  expect_true(all(is.na(out$resistance)))
  expect_true(all(is.na(out$path_resistance)))
  # the length-based columns are unaffected
  expect_equal(out$downstream_tips, c(3L, 2L, 1L, 1L, 1L))
  expect_equal(out$path_length, c(10, 20, 30, 20, 30))
})

test_that("degenerate and malformed input", {
  et <- tree5()
  expect_equal(root_transport(et[0, ]), et[0, ])
  expect_null(root_transport(NULL))
  expect_error(root_transport(et[, c("from", "to")]), "mean_diameter")
})

test_that("a single segment is its own base and tip", {
  et <- data.frame(from = "A", to = "B", length = 4, mean_diameter = 2,
                   stringsAsFactors = FALSE)
  out <- root_transport(et)
  expect_equal(out$downstream_tips, 1L)
  expect_equal(out$downstream_length, 4)
  expect_equal(out$path_length, 4)
})
