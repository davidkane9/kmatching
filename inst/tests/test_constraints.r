context("Testing constraints")

set.seed(44)

orig <- rexp(100)
orig <- orig/sum(orig)
data <- data.frame(value = rnorm(100), growth = rnorm(100), orig = orig, country = sample(c("Cuba", "USA", "Mexico", "Canada"), 100, replace = T))

k <- kmatch(data = data, match.var = c("value", "growth", "country"), weight.var = "orig", n = 100, replace = TRUE)
Amat <- matrix(c(data$value, data$growth), ncol = nrow(data), byrow = TRUE)
Amat <- rbind(Amat, .dummy(data$country))

#'Constraints:
#'-must add up to same as weight var
#'-must have same exposures

#'Errors:
#'-Must throw error if no solution
  
test_that("Generated Weights match the constraints", {
    ## match Ax = b, within tolerance (division error approximately 1e-12)
    expect_true(all(apply(k, 2, function(x) abs(Amat %*% x - Amat %*% data$orig) < 1e-12)))
    ## match x > 0
    expect_that(all(k > 0), is_true())
    ## expect that sums match
    expect_true(all(apply(k, 2, function(x) abs(sum(x) - sum(data$orig)) < 1e-12)))
    ## expect that there are the same number of row in output
    expect_that(nrow(k), equals(nrow(data)))
    })

orig2 = rexp(60)
orig2 = orig2/sum(orig2)

data$orig = c(orig2, rep(0, 40))
k <- kmatch(data = data, match.var = c("value", "growth", "country"), weight.var = "orig", n = 100)

test_that("replace = FALSE works correctly", {
  ## match Ax = b, within tolerance (division error approximately 1e-12)
  expect_true(all(apply(k, 2, function(x) abs(Amat %*% x - Amat %*% data$orig) < 1e-12)))
  ## match x > 0
  expect_that(all(k >= 0), is_true())
  ## expect that sums match
  expect_true(all(apply(k, 2, function(x) abs(sum(x) - sum(data$orig)) < 1e-12)))
  ## expect that there are the same number of row in output
  expect_that(nrow(k), equals(nrow(data)))
  ## expect that none of the weighted rows in the original are weighted in the output
  expect_true(all(apply(k[which(data$orig > 0),], 1, function(x) all(x==0))))
})

