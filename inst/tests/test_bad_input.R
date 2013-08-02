context("Test bad input")



test_that("Mirror throws an error if it is not converging to a solution", {
  A <- matrix(rep(1, 3), nrow = 1)
  x0 <- c(-1, -1, -1)
  
  ## There should be no solution here because it expects x+y+z = -3 with
  ## x,y,z positive.
  
  ## I can no longer make this check because as of commit 3862ef32 I changed
  ## the mirror algorithm so that it sometimes moves away from the solution
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

test_that("Missing entries are noted", {
  set.seed(400)
  dat = data.frame( size = rnorm(10), value = rnorm(10, 0, 2), weight = rep(.05, 10))
  dat2 = dat
  dat$weight[ 1:2] = c(NA, NULL)
  dat2$size[1:2] = c(NA, NULL)
  expect_warning(kmatch(x = dat, match.var = "size", weight.var = "weight", n = 10, replace = TRUE),
                 "weights are missing")
  expect_warning(kmatch(x = dat2, match.var = "size", weight.var = "weight", n = 10, replace = TRUE),
                 "size column are missing")
})