context("Test bad input")



test_that("Mirror throws an error if it is not converging to a solution", {
  A = matrix(rep(1, 3), nrow = 1)
  x0 = c(-1, -1, -1)
  
  expect_error(mirror(A, x0 = x0, n = 10))
  
})