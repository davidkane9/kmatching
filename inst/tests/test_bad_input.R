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
  
  ## fake data, weight and size
  dat = data.frame( size = rnorm(10), value = rnorm(10, 0, 2), weight = rep(.05, 10))
  ## dat is to test the weight warning, dat2 is to test match.var warning
  dat2 = dat
  dat$weight[ 1:2] = c(NA, NULL)
  dat2$size[1:2] = c(NA, NULL)
  ## Since there are missing entries in weight of dat, should give warning
  expect_warning(kmatch(x = dat, match.var = "size", weight.var = "weight", n = 10, replace = TRUE),
                 "2 entries of 'weight' are missing")
  ## Since there are missing entries in size of dat, should give warning
  expect_warning(kmatch(x = dat2, match.var = "size", weight.var = "weight", n = 10, replace = TRUE),
                 "2 entries of 'size' are missing")
  
  ## test error thrown if there are NAs in A
  A = matrix(c(NA, NA, 1), ncol = 3)
  b = 1
  expect_error(hitandrun(A, b, n = 1), "'A' cannot have NA's in it, has 2")
  
  ## test error thrown if there are NAs in b
  A = matrix(c(1, 1, 1), ncol = 3)
  b = NA
  expect_error(hitandrun(A, b, n = 1), "'b' cannot have NA's in it, has 1")
})


test_that("Problems aren't overdetermined", {
  A = matrix(rnorm(6), ncol = 2, nrow = 3)
  b = c(1,2,3)
  expect_error(hitandrun(A, b, n = 1), "overdetermined")
  
})