context("Test bad input with mirror")

test_that("Missing entries are noted in mirror", {
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
  Amat = matrix(c(NA, NA, 1), ncol = 3)
  x0 <- c(.3, .3, .4)
  expect_error(mirror(Amat, x0, 1), "'Amat' cannot have NA's in it. It currently has 2.")
  
  ## test error thrown if there are NAs in b
  Amat = matrix(c(1, 1, 1), ncol = 3)
  x0 <-c(.3, .7, NA)
  expect_error(mirror(Amat, x0, 1), "'x0' cannot have NA's in it. It currently has 1.")
})


test_that("Problems aren't overdetermined in mirror", {
  Amat = matrix(rnorm(6), ncol = 2, nrow = 3)
  x0 = c(.1,.2,.3)
  expect_error(mirror(Amat, x0, 1), "overdetermined")
  
})