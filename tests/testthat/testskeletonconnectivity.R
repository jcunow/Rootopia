test_that("skeletonization preserves connectivity (no fragmentation, no vanishing)", {
  # connected-component count of a 0/1 matrix.
  # `conn = 8` (default) for the foreground skeleton; connectivity duality means
  # the background must be counted with `conn = 4`, else inside/outside of a
  # closed 8-connected loop leak through the corners.
  n_components <- function(m, conn = 8) {
    nr <- nrow(m); nc <- ncol(m)
    offs <- if (conn == 4) list(c(-1,0), c(1,0), c(0,-1), c(0,1))
            else { o <- list(); for (dy in -1:1) for (dx in -1:1) if (dy||dx) o[[length(o)+1]] <- c(dy,dx); o }
    seen <- matrix(FALSE, nr, nc); n <- 0L
    for (i in seq_len(nr)) for (j in seq_len(nc)) {
      if (m[i, j] == 1 && !seen[i, j]) {
        n <- n + 1L; stack <- list(c(i, j)); seen[i, j] <- TRUE
        while (length(stack)) {
          p <- stack[[length(stack)]]; stack[[length(stack)]] <- NULL
          for (d in offs) {
            y <- p[1] + d[1]; x <- p[2] + d[2]
            if (y >= 1 && y <= nr && x >= 1 && x <= nc &&
                m[y, x] == 1 && !seen[y, x]) { seen[y, x] <- TRUE; stack[[length(stack) + 1]] <- c(y, x) }
          }
        }
      }
    }
    n
  }

  skel_mat <- function(m) {
    r <- terra::rast(m)
    s <- skeletonize_image(r, verbose = FALSE)
    matrix(terra::values(s), nrow = terra::nrow(s), byrow = TRUE)
  }

  # 2px-wide line: must survive as one connected component (regression: it vanished)
  a <- matrix(0L, 9, 20); a[4:5, 3:18] <- 1L
  s <- skel_mat(a)
  expect_gt(sum(s), 0)                     # not erased
  expect_equal(n_components(s), 1L)        # still connected

  # thick ring: one component, and the central hole is preserved
  a <- matrix(0L, 21, 21); a[4:18, 4:18] <- 1L; a[8:14, 8:14] <- 0L
  s <- skel_mat(a)
  expect_gt(sum(s), 0)
  expect_equal(n_components(s), 1L)
  # the loop encloses a hole, so the background (counted 4-connected, per
  # connectivity duality) splits into outside + inside
  expect_equal(n_components(1L - s, conn = 4), 2L)

  # cross stays one component
  a <- matrix(0L, 21, 21); a[10:12, 3:19] <- 1L; a[3:19, 10:12] <- 1L
  s <- skel_mat(a)
  expect_equal(n_components(s), 1L)
})
