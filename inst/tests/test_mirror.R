context("Test mirror")

test_that("Single sample", {
  A <- matrix(rep(1, 3), nrow = 1)
  x0 <- c(1, 1, 1)
  
  set.seed(10)
  expect_that(round(mirror(A, x0, 1), 2), equals(c(1.87, 0.36, 0.77)))
   
})

test_that("make sure 100 samples meet constraint", {
  A <- matrix(rep(1, 3), nrow = 1)
  x0 <- c(2, 2, 2)
  
  ## All solutions, each of which is a 3-tuple, should have components that add
  ## to 6, just as the initial point does.
  
  set.seed(5)
  expect_that(apply(mirror(A, x0, 100), 2, sum), equals(rep(6, 100)))
  
})



