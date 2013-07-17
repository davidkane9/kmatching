library(kmatching)
library(testthat)
context("Testing constraints")

orig <- rexp(1000)
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
    ## match Ax = b
    expect_true(all(apply(k, 2, function(x) abs(Amat %*% x - Amat %*% data$orig) < 1e-12)))
    ## match x > 0
    expect_that(all(k > 0), is_true())
    ## expect that sums match
    expect_true(all(apply(k, 2, function(x) abs(sum(x) - sum(data$orig)) < 1e-12)))
    ## expect that there are the same number of row in output
    expect_that(nrow(k), equals(nrow(data)))
    })
