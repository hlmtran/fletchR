context("Testing pruneCol")

col_test <- matrix(
  c(
    1, 0, 0,
    1, 1, 0,
    1, 1, 1
  ),
  nrow = 3
)
col_test

testthat::expect_equal(
  pruneCols(col_test, prune = 2),
  c(TRUE, FALSE, FALSE)
)