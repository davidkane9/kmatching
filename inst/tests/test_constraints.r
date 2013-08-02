context("Testing constraints")
library(stats)
library(utils)

set.seed(44)

orig <- rexp(100)
orig <- orig/sum(orig)
data <- data.frame(value = rnorm(100), growth = rnorm(100), orig = orig, country = sample(c("Cuba", "USA", "Mexico", "Canada"), 100, replace = T))

k <- kmatch(data, match.var = c("value", "growth", "country"), weight.var = "orig", n = 100, replace = TRUE)[[1]]
Amat <- matrix(c(data$value, data$growth), ncol = nrow(data), byrow = TRUE)
Amat <- rbind(Amat, dummy(data$country))

#'Constraints:
#' Test constraints for all chains
#' 
#' 
  
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
k <- kmatch(x = data, match.var = c("value", "growth", "country"), weight.var = "orig", n = 100)[[1]]

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

test_that("samples are generated uniformly", {

  ## 2 variables:
  ## A sample of 1000 coordinates for 2 variables should follow binomial
  ## distribution with respect to being above/below .5
  A = matrix(1,ncol = 2)
  b = 1
  set.seed(1)
  k = hitandrun(A, b, n = 1000)[[1]]
  expect_true(all(apply(k, 1, function(x) length(which(x> .5)) < qbinom(.99, 1000, .5))))

  ## 3 variables: 
  A = matrix(1, ncol = 3)
  b = 1
  set.seed(2)
  k = hitandrun(A,b, n = 1000)[[1]]
  bins = apply(k, 2, function(x) {
    ## cool 'hack'
    bin = (x[1] >= x[3]) + 2*(x[1] >= x[2]) + 4*(x[2] >= x[3])
    return(bin)
  })
  ChiSquared = chisq.test(table(bins))
  expect_true(ChiSquared$p.value > .1) 
  
  
  ## n variables, these bins will grow like 2^(n^2), unfortunately, suffering
  ## from an extreme curse of dimensionality. High dimensions will make most bins=0.
  n = 5
  A = matrix(1, ncol = n)
  b = 1
  set.seed(11)
  k = hitandrun(A,b, n = 1000)[[1]]
  
  ## this creates a list of combinations of the 10 variables we could have
  ## we will compare each combination, and bin them based on whether one is
  ## greater or less than the other, there are roughly n^2 combinations
  combin = combn(1:n, 2)
  bins = apply(k, 2, function(x) {
    bin = "b"
    ## digit i is 1 if ith combination returns true else 0, there are roughly
    ## 2^n ways each weighted rows 
    for(i in 1:ncol(combin))
      bin = paste(bin, as.numeric(x[combin[1,i]] >= x[combin[2,i]]), sep = "")
    return(bin)    
  })
  ChiSquared = chisq.test(table(bins))
  expect_true(ChiSquared$p.value > .1)
})

test_that("hitandrun knows when there is only 1 solution", {
  A = matrix(c(1,1,1,1,0,-1,1,-1,0), ncol = 3, byrow = T)
  b = c(1,0,0)
  expect_warning((h = hitandrun(A, b, n = 10)[[1]]), "unique")
  expect_that(nrow(h), equals(3))
  expect_that(ncol(h), equals(1))
  expect_that(h[,1], equals(c(solve(A) %*% b)))
})

test_that("kmatch knows when there is only 1 solutions", {
  A = matrix(c(1,1,1,1,0,-1,1,-1,0), ncol = 3, byrow = T)
  set.seed(87)
  dat = as.data.frame(t(A[2:3,]))
  dat$x0 = c(1/3,1/3,1/3)
  expect_warning((
    k = kmatch(dat, weight.var= "x0",match.var = c("V1", "V2"), n = 10, replace= TRUE)[[1]]
  ), "unique")
  expect_that(k, equals(solve(A) %*% A %*% dat$x0))
  
  set.seed(20)
  dat2 = data.frame(a = rnorm(3), weight = c(1,0,0))
  dat2[1,1] = 0
  expect_warning((
    k = kmatch(dat2, weight.var = "weight", match.var = "a", n = 10)
  ), "unique")
  expect_that(round(k[[1]][,1],2), equals(c(0, .75, .25 ) ))
})
