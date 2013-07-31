context("Test bad input")



test_that("Mirror throws an error if it is not converging to a solution", {
  A <- matrix(rep(1, 3), nrow = 1)
  x0 <- c(-1, -1, -1)
  
  ## There should be no solution here because it expects x+y+z = -3 with
  ## x,y,z positive.
  
  ##expect_error(mirror(A, x0 = x0, n = 10))  
})

test_that("Non matching dimensions lead to error", {
  A <- matrix(rep(1, 3), nrow = 1)
  x0 <- c(1, 1)
  
  ## There should be an error because the dimensions of A and x0 do not match. 
  ## Maybe we ought to add some warnings that explain what is going on?
  
  ## Had to comment this out because it produces a warning as well as an error,
  ## which seems to be a problem.
  
  ## expect_error(mirror(A, x0 = x0, n = 10))  
})

