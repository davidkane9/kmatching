context("Hit and Run with NA")

test_that("Hit and run fails for missing values", {
  
  Amat = matrix(c(1,1,NA), nrow = 1, ncol = 3)
  b = 1
  
  set.seed(40)
  expect_error(hitandrun(A = Amat, b = b, n = 1000))
  
})

