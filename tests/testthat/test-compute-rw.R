# ---- scgas_compute random walk ---------------------------------------------

test_that(".batched_rw on a column-stochastic matrix keeps unit column mass", {
  set.seed(4)
  n <- 60
  A <- Matrix::rsparsematrix(n, n, 0.15)
  A <- (abs(A) + Matrix::t(abs(A))) > 0
  A <- A * 1
  diag(A) <- 0
  A <- Matrix::drop0(A)
  M <- A %*% Matrix::Diagonal(x = 1 / pmax(Matrix::colSums(A), 1))
  # every node has at least one edge in this fixture
  skip_if(any(Matrix::colSums(A) == 0))
  P0 <- matrix(0, n, 2); P0[1:6, 1] <- 1 / 6; P0[7:12, 2] <- 1 / 6
  P  <- scGAS:::.batched_rw(M, P0, gamma = 0.05)
  expect_equal(unname(colSums(P)), c(1, 1), tolerance = 1e-6)
  expect_true(all(P >= 0))
  # different seeds must give different stationary profiles
  expect_lt(cor(P[, 1], P[, 2]), 0.99)
})
