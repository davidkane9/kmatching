context("Test kmatch")

## Ought to consider pulling out Lalonde data and adding it to the kmatching
## package, but not today.


test_that("Using Lalonde produces results that match perfectly", {
  library(MatchIt)
  data(lalonde)
  
  ## Single 0/1 variable.
  
  z <- kmatch(lalonde, match.var = c("hispan"), weight.var = "treat", n = 10)
    
  for(i in 1:10){
  expect_that(sum(lalonde$hispan[lalonde$treat == 1]), 
              equals(apply(lalonde$hispan * z, 2, sum)[i]))
  }  
  
  ## Double check that it works for a factor variable as well.
  
  lalonde2 <- lalonde
  lalonde2$hispan <- factor(lalonde2$hispan, labels = c("Non-hispanic", "Hispanic"))
  
  z <- kmatch(lalonde2, match.var = c("hispan"), weight.var = "treat", n = 10)
  for(i in 1:10){
    expect_that(sum(lalonde$hispan[lalonde$treat == 1]), 
                equals(apply(lalonde$hispan * z, 2, sum)[i]))
  } 
  
  ## Two 0/1 variables.
  
  z <- kmatch(lalonde, match.var = c("hispan", "black"), weight.var = "treat", n = 10)
  
  for(i in 1:10){
    expect_that(sum(lalonde$hispan[lalonde$treat == 1]), 
                equals(apply(lalonde$hispan * z, 2, sum)[i]))
    expect_that(sum(lalonde$black[lalonde$treat == 1]), 
                equals(apply(lalonde$black * z, 2, sum)[i]))
  }
   
  ## Two 0/1 variables and continuous variable.
  
  z <- kmatch(lalonde, match.var = c("hispan", "black", "educ"), weight.var = "treat", n = 10)
  
  ## Be careful in the use of sum versus mean. You want 'mean' to work as well,
  ## but it doesn't because the number of observation is different when
  ## comparing the subset which are Hispanic to a matching vector which includes
  ## the entire sample.
  
  for(i in 1:10){
    expect_that(sum(lalonde$hispan[lalonde$treat == 1]), 
                equals(apply(lalonde$hispan * z, 2, sum)[i]))
    expect_that(sum(lalonde$black[lalonde$treat == 1]), 
                equals(apply(lalonde$black * z, 2, sum)[i]))
    expect_that(sum(lalonde$educ[lalonde$treat == 1]), 
                equals(apply(lalonde$educ * z, 2, sum)[i]))
  }
  
  
  
})

test_that("make sure 100 samples meet constraint", {
  A <- matrix(rep(1, 3), nrow = 1)
  x0 <- c(2, 2, 2)
  
  ## All solutions, each of which is a 3-tuple, should have components that add
  ## to 6, just as the initial point does.
  
  set.seed(5)
  expect_that(apply(mirror(A, x0, 100), 2, sum), equals(rep(6, 100)))
  
})



