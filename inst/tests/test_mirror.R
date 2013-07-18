context("Test mirror")



test_that("Single sample", {
  A <- matrix(rep(1, 3), nrow = 1)
  x0 <- c(1, 1, 1)
  
  set.seed(10)
  expect_that(round(mirror(A, x0, 1), 2), equals(c(1.87, 0.36, 0.77)))
   
})
